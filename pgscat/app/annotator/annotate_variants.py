#!/usr/bin/env python3
"""
Variant Annotator — PGS/PRS Platform  v1.3
============================================
Annotates each variant in a betamap file using:
  - GENCODE GFF3 annotation (required)
  - Reference FASTA (optional — enables coding-consequence + splice GT/AG logic)
  - dbNSFP5 tabix index (optional — CADD, REVEL, SIFT, PolyPhen2, ClinVar)
  - Regulatory BED files (optional — ENCODE, Ensembl Regulatory, custom)
  - dbSNP population frequency VCF (optional — rsid, global AF, population AFs)

Usage (inside Apptainer):
    python annotate_variants.py \\
        --betamap  /data/PGS000001/PGS000001_hmPOS_GRCh38.betamap.tsv.gz \\
        --gff3     /ref/gencode.v49.basic.annotation.gff3 \\
        --outdir   /annotations/PGS000001 \\
        --pgs-id   PGS000001 \\
        [--fasta   /fasta/hg38.fa] \\
        [--dbnsfp  /dbnsfp/dbNSFP5.0a_grch38.gz] \\
        [--dbsnp   /dbsnp_freq/freq.vcf.gz] \\
        [--regulatory-bed /regulatory/encode.bed.gz \\
                          /regulatory/ensembl_reg.bed.gz] \\
        [--no-parquet] [--verbose]

Outputs written to --outdir:
    {PGS_ID}_variants_annotated.tsv.gz          full annotation table
    {PGS_ID}_variants_annotated.parquet          same, columnar format
    {PGS_ID}_annotation_summary.json            per-region-class counts + metadata

Input (betamap) columns (must be present, order-independent):
    PRS_ID  CHROM  POS  ID  EFFECT_ALLELE  OTHER_ALLELE  BETA  IS_FLIP

Annotation columns added (v1.0):
    gene_name  gene_id  gene_type  transcript_id  transcript_type
    feature_type  region_class  consequence
    is_coding  is_regulatory  is_intergenic
    n_overlapping_genes  strand  distance_nearest_gene

Annotation columns added (v1.1):
    distance_to_splice_site  consequence_priority
    all_overlapping_genes  all_overlapping_transcripts  all_region_classes
    regulatory_source

Annotation columns added (v1.2, FASTA/dbNSFP-backed):
    codon_ref  codon_alt  aa_ref  aa_alt  aa_ref_3  aa_alt_3
    splice_type  cadd_phred  revel_score  sift_pred  polyphen2_pred
    clinvar_clnsig  is_missense  is_synonymous  is_lof

Annotation columns added (v1.3, dbSNP population frequency):
    rsid  af_global  af_max_population  af_population_summary  rarity_class
"""
from __future__ import annotations

import argparse
import gzip
import json
import logging
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import List, Optional

import pandas as pd

# ── Local imports (same package) ─────────────────────────────────────────────
# When running inside Apptainer the working directory is /app/annotator
# so relative imports work.
try:
    from gff3_parser import GFF3Index
    from classify import classify, summarise, AnnotationResult
    from regulatory import RegulatoryBEDIndex
    from fasta_engine import FASTAEngine
    from dbnsfp import DbNSFP5Fetcher
    from dbsnp_freq import DbSNPFreqFetcher
except ImportError:
    import sys as _sys
    _sys.path.insert(0, str(Path(__file__).parent))
    from gff3_parser import GFF3Index
    from classify import classify, summarise, AnnotationResult
    from regulatory import RegulatoryBEDIndex
    from fasta_engine import FASTAEngine
    from dbnsfp import DbNSFP5Fetcher
    from dbsnp_freq import DbSNPFreqFetcher

# ── Logging ───────────────────────────────────────────────────────────────────
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    datefmt="%Y-%m-%dT%H:%M:%S",
)
logger = logging.getLogger("variant_annotator")

# ── Constants ─────────────────────────────────────────────────────────────────
# Single source of truth for version — propagated to summary JSON,
# log messages, and the annotation_tool / schema_version fields.
TOOL_VERSION    = "variant_annotator/1.6"
SCHEMA_VERSION  = TOOL_VERSION.split("/", 1)[-1]   # "1.6"

BETAMAP_COLUMNS = ["PRS_ID", "CHROM", "POS", "ID",
                   "EFFECT_ALLELE", "OTHER_ALLELE", "BETA", "IS_FLIP"]

# v1.0 annotation columns (backward-compatible)
_ANNOTATION_COLUMNS_V1 = [
    "gene_name", "gene_id", "gene_type",
    "transcript_id", "transcript_type",
    "feature_type", "region_class", "consequence",
    "is_coding", "is_regulatory", "is_intergenic",
    "n_overlapping_genes", "strand", "distance_nearest_gene",
]

# v1.1 new annotation columns
_ANNOTATION_COLUMNS_V1_1 = [
    "distance_to_splice_site",
    "consequence_priority",
    "all_overlapping_genes",
    "all_overlapping_transcripts",
    "all_region_classes",
    "regulatory_source",
]

# v1.2 new annotation columns (FASTA + dbNSFP5 backed)
_ANNOTATION_COLUMNS_V1_2 = [
    "codon_ref",
    "codon_alt",
    "aa_ref",
    "aa_alt",
    "aa_ref_3",
    "aa_alt_3",
    "splice_type",
    "cadd_phred",
    "revel_score",
    "sift_pred",
    "polyphen2_pred",
    "clinvar_clnsig",
    "is_missense",
    "is_synonymous",
    "is_lof",
]

# v1.3 new annotation columns (dbSNP population frequency)
_ANNOTATION_COLUMNS_V1_3 = [
    "rsid",
    "af_global",
    "af_max_population",
    "af_population_summary",
    "rarity_class",
]

# v1.5 columns: regulatory element detail (real BED overlap or nearest-element distance)
_ANNOTATION_COLUMNS_V1_5 = [
    "regulatory_element_type",   # promoter | enhancer | insulator | TFBS | open_chromatin | …
    "regulatory_element_id",     # TF name, OREG ID, or element name
    "distance_to_regulatory",    # 0 = overlapping; >0 = bp to nearest element; None = no BED
]

# v1.4 columns: gnomAD 4.1.1 + authoritative ClinVar (added at annotation time
# when --gnomad-dir and/or --clinvar flags are provided)
_ANNOTATION_COLUMNS_V1_4 = [
    "af_gnomad",
    "af_gnomad_afr",
    "af_gnomad_amr",
    "af_gnomad_eas",
    "af_gnomad_nfe",
    "af_gnomad_sas",
    "af_gnomad_mid",
    "an_gnomad",
    "nhomalt_gnomad",
    "rarity_class_gnomad",
    "clinvar_clndn",
    "clinvar_clnrevstat",
    "clinvar_alleleid",
    "clinvar_clnvc",
]

# v1.6 columns: MCPS population AFs + effective AF (best available across sources)
_ANNOTATION_COLUMNS_V1_6 = [
    "af_mcps",           # MCPS global AF (AF_RAW)
    "af_mcps_mex",       # MCPS Mexican population AF (primary for MCPS)
    "af_mcps_eur",       # MCPS European AF
    "af_mcps_afr",       # MCPS African AF
    "an_mcps",           # MCPS allele number
    "ac_mcps",           # MCPS alt allele count
    "rarity_class_mcps", # rarity based on MCPS AF
    # Effective AF: highest-priority available AF across all sources
    "af_effective",          # MCPS MEX › MCPS RAW › gnomAD › dbSNP
    "af_source",             # "mcps_mex" | "mcps_raw" | "gnomad" | "dbsnp" | None
    "rarity_class_effective",# rarity classification of af_effective
]

ANNOTATION_COLUMNS = (
    _ANNOTATION_COLUMNS_V1
    + _ANNOTATION_COLUMNS_V1_1
    + _ANNOTATION_COLUMNS_V1_2
    + _ANNOTATION_COLUMNS_V1_3
    + _ANNOTATION_COLUMNS_V1_4
    + _ANNOTATION_COLUMNS_V1_5
    + _ANNOTATION_COLUMNS_V1_6
)
OUTPUT_COLUMNS     = BETAMAP_COLUMNS + ANNOTATION_COLUMNS


# ── Betamap loader ────────────────────────────────────────────────────────────

def load_betamap(path: Path) -> pd.DataFrame:
    """
    Load betamap TSV.GZ.
    Validates required columns; normalises CHROM to string.
    """
    logger.info("Loading betamap: %s", path)
    df = pd.read_csv(
        path,
        sep="\t",
        dtype={
            "PRS_ID": str,
            "CHROM": str,
            "POS": "int64",
            "ID": str,
            "EFFECT_ALLELE": str,
            "OTHER_ALLELE": str,
            "BETA": "float64",
            "IS_FLIP": "int8",
        },
        compression="infer",
    )

    missing = [c for c in BETAMAP_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(
            f"Betamap is missing required columns: {missing}\n"
            f"Found: {list(df.columns)}"
        )

    # Normalise CHROM: remove 'chr' prefix, map M→MT
    df["CHROM"] = (
        df["CHROM"]
        .astype(str)
        .str.lstrip("chr")
        .replace({"M": "MT"})
    )

    logger.info(
        "Betamap loaded: %d variants across %d chromosomes",
        len(df),
        df["CHROM"].nunique(),
    )
    return df


# ── Main annotation loop ──────────────────────────────────────────────────────

def annotate(
    df: pd.DataFrame,
    idx: GFF3Index,
    reg_idx: Optional[RegulatoryBEDIndex] = None,
    fasta_engine: Optional["FASTAEngine"] = None,
    dbnsfp_fetcher: Optional["DbNSFP5Fetcher"] = None,
    dbsnp_fetcher: Optional["DbSNPFreqFetcher"] = None,
) -> pd.DataFrame:
    """
    For each variant row, query the GFF3 index (and optional regulatory index,
    FASTA engine, dbNSFP5 fetcher, and dbSNP frequency fetcher) and classify.

    Returns the input DataFrame with ANNOTATION_COLUMNS appended.

    v1.2 additions:
      - fasta_engine: if provided, calls get_codon_consequence() for CDS variants
        and get_splice_type() for near-splice intronic variants.
      - dbnsfp_fetcher: if provided, queries CADD/REVEL/SIFT/PolyPhen2/ClinVar
        for all SNP variants.

    v1.3 additions:
      - dbsnp_fetcher: if provided, queries population allele frequencies from
        dbSNP freq.vcf.gz (rsid, af_global, af_max_population, rarity_class).
        dbNSFP5 is RETAINED; dbSNP is complementary.
    """
    n = len(df)
    logger.info("Annotating %d variants …", n)
    if reg_idx and not reg_idx.is_empty:
        logger.info("  Regulatory BED sources: %s", ", ".join(reg_idx.sources))
    if fasta_engine:
        logger.info("  FASTA engine active — coding consequences and splice GT/AG enabled")
    if dbnsfp_fetcher:
        logger.info("  dbNSFP5 fetcher active — CADD/REVEL/SIFT/PolyPhen2/ClinVar enabled")
    if dbsnp_fetcher:
        logger.info("  dbSNP freq fetcher active — rsid/AF/rarity_class enabled")

    rows: list[dict] = []
    t0 = time.monotonic()
    log_every = max(1, n // 20)  # log ~20 progress updates

    for i, (_, row) in enumerate(df.iterrows()):
        chrom         = str(row["CHROM"])
        pos           = int(row["POS"])
        effect_allele = str(row["EFFECT_ALLELE"])
        other_allele  = str(row["OTHER_ALLELE"])
        is_flip       = int(row.get("IS_FLIP", 0))

        # ── GFF3 overlap ─────────────────────────────────────────────────────
        features = idx.query(chrom, pos)

        nearest_dist: Optional[int] = None
        if not features:
            nearest_dist = idx.nearest_gene_distance(chrom, pos)

        # ── Splice site distance (all variants) ───────────────────────────────
        splice_dist = idx.nearest_exon_boundary_distance(chrom, pos)

        # ── Regulatory BED overlap + nearest element distance ─────────────────
        reg_hits: list[dict] = []
        reg_nearest_dist: Optional[int] = None
        if reg_idx and not reg_idx.is_empty:
            reg_hits = reg_idx.query(chrom, pos)
            if not reg_hits:
                # No overlap — find nearest element (up to 100 kb)
                near_feat, near_dist = reg_idx.nearest(chrom, pos, max_distance=100_000)
                if near_feat is not None and near_dist is not None:
                    # Annotate as nearest-element hint (distance > 0 → not a hit)
                    reg_nearest_dist = near_dist
                    # Attach distance to hits so classify() can pick it up
                    near_feat = dict(near_feat)
                    near_feat["distance"] = near_dist
                    reg_hits = [near_feat] if near_dist == 0 else []
                    # Only flag as regulatory hit if actually overlapping (distance=0)

        # ── FASTA: coding consequence (CDS variants only) ─────────────────────
        codon_result: Optional[dict] = None
        if fasta_engine:
            # Determine whether this variant is in a CDS before the expensive call.
            # We check quickly: any feature in the hit list is a CDS.
            is_cds = any(f.feature == "CDS" for f in features)
            if is_cds:
                # Get the primary transcript_id from the best CDS feature
                cds_feat = next((f for f in features if f.feature == "CDS"), None)
                if cds_feat and cds_feat.transcript_id:
                    try:
                        codon_result = fasta_engine.get_codon_consequence(
                            chrom=chrom,
                            pos=pos,
                            effect_allele=effect_allele,
                            other_allele=other_allele,
                            transcript_id=cds_feat.transcript_id,
                            gff3_index=idx,
                        )
                    except Exception as exc:
                        logger.debug(
                            "codon_consequence failed %s:%d: %s", chrom, pos, exc
                        )

        # ── FASTA: splice site direction (near-splice intronic variants) ───────
        fasta_splice_type: Optional[str] = None
        if fasta_engine and splice_dist is not None and splice_dist <= 10:
            # Only useful for intronic variants near exon boundaries
            is_intronic_near_splice = (
                not any(f.feature in {"CDS", "exon", "five_prime_UTR", "three_prime_UTR"}
                        for f in features)
                and len(features) > 0
            )
            if is_intronic_near_splice:
                try:
                    fasta_splice_type = fasta_engine.get_splice_type(
                        chrom=chrom,
                        pos=pos,
                        gff3_index=idx,
                    )
                except Exception as exc:
                    logger.debug(
                        "get_splice_type failed %s:%d: %s", chrom, pos, exc
                    )

        # ── Determine REF/ALT for tabix queries (shared by dbNSFP5 + dbSNP) ─────
        # Betamap convention:
        #   IS_FLIP=0 → EFFECT_ALLELE is ALT, OTHER_ALLELE is REF (standard)
        #   IS_FLIP=1 → EFFECT_ALLELE is REF (dosage = 2 − DS), OTHER_ALLELE is ALT
        resolved_ref: Optional[str] = None
        resolved_alt: Optional[str] = None
        if dbnsfp_fetcher or dbsnp_fetcher:
            if fasta_engine and codon_result and codon_result.get("ref_allele"):
                # FASTA-confirmed alleles — most accurate
                resolved_ref = codon_result["ref_allele"]
                resolved_alt = codon_result["alt_allele"]
            elif fasta_engine:
                resolved_ref, resolved_alt = fasta_engine.determine_ref_alt(
                    chrom, pos, effect_allele, other_allele
                )
                if not resolved_ref:
                    # FASTA could not determine — fall back to betamap convention
                    if is_flip:
                        resolved_ref, resolved_alt = effect_allele, other_allele
                    else:
                        resolved_ref, resolved_alt = other_allele, effect_allele
            else:
                # No FASTA — use betamap convention directly
                if is_flip:
                    resolved_ref, resolved_alt = effect_allele, other_allele
                else:
                    resolved_ref, resolved_alt = other_allele, effect_allele

        # ── dbNSFP5 scores (SNPs only — indels not scored) ────────────────────
        dbnsfp_scores: Optional[dict] = None
        if dbnsfp_fetcher and resolved_ref and resolved_alt:
            # Only query for SNPs (ref and alt both single nucleotides)
            if len(resolved_ref) == 1 and len(resolved_alt) == 1:
                try:
                    dbnsfp_scores = dbnsfp_fetcher.fetch(
                        chrom, pos, resolved_ref, resolved_alt
                    )
                except Exception as exc:
                    logger.debug(
                        "dbnsfp_fetch failed %s:%d: %s", chrom, pos, exc
                    )

        # ── dbSNP population frequencies (v1.3) ───────────────────────────────
        dbsnp_freq: Optional[dict] = None
        if dbsnp_fetcher and resolved_ref and resolved_alt:
            try:
                dbsnp_freq = dbsnp_fetcher.fetch(
                    chrom, pos, resolved_ref, resolved_alt
                )
            except Exception as exc:
                logger.debug(
                    "dbsnp_freq fetch failed %s:%d: %s", chrom, pos, exc
                )

        # ── Classify ──────────────────────────────────────────────────────────
        ann = classify(
            features,
            nearest_gene_dist=nearest_dist,
            nearest_exon_dist=splice_dist,
            regulatory_hits=reg_hits,
            codon_result=codon_result,
            splice_type=fasta_splice_type,
            dbnsfp_scores=dbnsfp_scores,
        )

        # Build row dict — start from classification result, then add dbSNP fields
        row_dict = ann.to_dict()

        # v1.5: distance_to_regulatory — override with nearest-element distance
        # if the variant didn't directly overlap (classify() sets it to None in that case)
        if reg_idx and not reg_idx.is_empty:
            if row_dict.get("distance_to_regulatory") is None and reg_nearest_dist is not None:
                row_dict["distance_to_regulatory"] = reg_nearest_dist

        # v1.3: merge dbSNP population frequency fields
        if dbsnp_freq:
            row_dict["rsid"]                  = dbsnp_freq.get("rsid")
            row_dict["af_global"]             = dbsnp_freq.get("af_global")
            row_dict["af_max_population"]     = dbsnp_freq.get("af_max_population")
            row_dict["af_population_summary"] = dbsnp_freq.get("af_population_summary")
            row_dict["rarity_class"]          = dbsnp_freq.get("rarity_class", "novel")
        else:
            row_dict["rsid"]                  = None
            row_dict["af_global"]             = None
            row_dict["af_max_population"]     = None
            row_dict["af_population_summary"] = None
            row_dict["rarity_class"]          = "novel"

        rows.append(row_dict)

        if (i + 1) % log_every == 0:
            elapsed = time.monotonic() - t0
            rate = (i + 1) / elapsed if elapsed > 0 else 0
            eta  = (n - i - 1) / rate if rate > 0 else 0
            logger.info(
                "  %d / %d (%.0f%%)  %.0f var/s  ETA %.0f s",
                i + 1, n, 100 * (i + 1) / n, rate, eta,
            )

    ann_df = pd.DataFrame(rows, index=df.index)
    result = pd.concat([df, ann_df], axis=1)

    elapsed = time.monotonic() - t0
    logger.info(
        "Annotation complete: %d variants in %.1f s (%.0f var/s)",
        n, elapsed, n / elapsed if elapsed > 0 else 0,
    )
    return result


# ── Writers ───────────────────────────────────────────────────────────────────

def write_tsv_gz(df: pd.DataFrame, path: Path) -> None:
    logger.info("Writing TSV.GZ: %s", path)
    df.to_csv(path, sep="\t", index=False, compression="gzip")
    logger.info("  → %s (%.1f MB)", path.name, path.stat().st_size / 1e6)


def write_parquet(df: pd.DataFrame, path: Path) -> None:
    try:
        import pyarrow as pa
        import pyarrow.parquet as pq
    except ImportError:
        logger.warning("pyarrow not available — skipping parquet output.")
        return

    logger.info("Writing Parquet: %s", path)
    # Cast boolean columns explicitly for clean parquet schema
    for col in ["is_coding", "is_regulatory", "is_intergenic"]:
        if col in df.columns:
            df[col] = df[col].astype(bool)

    table = pa.Table.from_pandas(df, preserve_index=False)
    pq.write_table(table, str(path), compression="snappy")
    logger.info("  → %s (%.1f MB)", path.name, path.stat().st_size / 1e6)


def write_summary(
    df: pd.DataFrame,
    pgs_id: str,
    betamap_path: Path,
    gff3_path: Path,
    outdir: Path,
    elapsed_s: float,
    fasta_path: Optional[Path] = None,
    dbnsfp_path: Optional[Path] = None,
    dbsnp_path: Optional[Path] = None,
    regulatory_bed_paths: Optional[List[Path]] = None,
    reg_idx: Optional[RegulatoryBEDIndex] = None,
) -> dict:
    """Build and write annotation_summary.json."""
    # Use only columns that are present in the dataframe (graceful for missing v1.2 cols)
    avail_ann_cols = [c for c in ANNOTATION_COLUMNS if c in df.columns]
    annotations = [
        AnnotationResult(
            gene_name              = row.get("gene_name", ""),
            gene_id                = row.get("gene_id", ""),
            gene_type              = row.get("gene_type", ""),
            transcript_id          = row.get("transcript_id", ""),
            transcript_type        = row.get("transcript_type", ""),
            feature_type           = row.get("feature_type", ""),
            region_class           = row.get("region_class", "intergenic"),
            consequence            = row.get("consequence", "intergenic_variant"),
            is_coding              = bool(row.get("is_coding", False)),
            is_regulatory          = bool(row.get("is_regulatory", False)),
            is_intergenic          = bool(row.get("is_intergenic", True)),
            n_overlapping_genes    = int(row.get("n_overlapping_genes", 0)),
            strand                 = row.get("strand", "."),
            distance_nearest_gene  = row.get("distance_nearest_gene"),
            distance_to_splice_site = row.get("distance_to_splice_site"),
            consequence_priority   = int(row.get("consequence_priority", 900)),
            all_overlapping_genes  = str(row.get("all_overlapping_genes", "") or ""),
            all_overlapping_transcripts = str(row.get("all_overlapping_transcripts", "") or ""),
            all_region_classes     = str(row.get("all_region_classes", "") or ""),
            regulatory_source         = str(row.get("regulatory_source", "") or ""),
            # v1.5 fields
            regulatory_element_type   = str(row.get("regulatory_element_type", "") or ""),
            regulatory_element_id     = str(row.get("regulatory_element_id", "") or ""),
            distance_to_regulatory    = row.get("distance_to_regulatory"),
            # v1.2 fields
            codon_ref       = str(row.get("codon_ref", "") or ""),
            codon_alt       = str(row.get("codon_alt", "") or ""),
            aa_ref          = str(row.get("aa_ref", "") or ""),
            aa_alt          = str(row.get("aa_alt", "") or ""),
            splice_type     = str(row.get("splice_type", "") or ""),
            cadd_phred      = row.get("cadd_phred"),
            revel_score     = row.get("revel_score"),
            sift_pred       = str(row.get("sift_pred", "") or ""),
            polyphen2_pred  = str(row.get("polyphen2_pred", "") or ""),
            clinvar_clnsig  = str(row.get("clinvar_clnsig", "") or ""),
            is_missense     = bool(row.get("is_missense", False)),
            is_synonymous   = bool(row.get("is_synonymous", False)),
            is_lof          = bool(row.get("is_lof", False)),
        )
        for row in df[avail_ann_cols].to_dict(orient="records")
    ]

    stats = summarise(annotations)

    # Per-chromosome breakdown
    chrom_counts: dict = {}
    for chrom, grp in df.groupby("CHROM"):
        chrom_counts[str(chrom)] = {
            "n_variants":      len(grp),
            "n_coding":        int(grp["is_coding"].sum()),
            "n_intergenic":    int(grp["is_intergenic"].sum()),
            "n_splice_site":   int((grp["consequence"] == "splice_site_variant").sum()),
            "n_splice_region": int((grp["consequence"] == "splice_region_variant").sum()),
        }

    # Determine what's resolved vs pending
    resolved_limitations = [
        "Splice site distance computed (distance_to_splice_site column).",
        "Splice site and splice region variants reclassified from intron_variant.",
        "All overlapping genes/transcripts now reported.",
        "Consequence priority numeric ranking added.",
        "Regulatory BED overlap supported via --regulatory-bed (if provided).",
    ]
    if fasta_path:
        resolved_limitations += [
            "Coding consequences computed from FASTA: missense, synonymous, stop_gained,"
            " stop_lost, start_lost, frameshift, inframe_insertion, inframe_deletion.",
            "Splice donor vs acceptor classified by GT/AG FASTA inspection.",
            "REF allele verified against GRCh38 reference FASTA.",
        ]
    if dbnsfp_path:
        resolved_limitations.append(
            "Functional scores fetched from dbNSFP5: CADD, REVEL, SIFT, PolyPhen2, ClinVar."
        )
    if dbsnp_path:
        resolved_limitations.append(
            "Population allele frequencies fetched from dbSNP: rsid, af_global, "
            "af_max_population, af_population_summary, rarity_class."
        )

    remaining_limitations = []
    if not fasta_path:
        remaining_limitations.append(
            "FASTA reference not provided — coding consequences (missense, synonymous,"
            " stop-gained, frameshift) not computed. Provide --fasta."
        )
    if not dbnsfp_path:
        remaining_limitations.append(
            "dbNSFP5 not provided — CADD/REVEL/SIFT/PolyPhen2/ClinVar scores absent."
            " Provide --dbnsfp."
        )
    if not dbsnp_path:
        remaining_limitations.append(
            "dbSNP population frequency VCF not provided — rsid, AF, and rarity_class absent."
            " Provide --dbsnp /dbsnp_freq/freq.vcf.gz."
        )
    if not regulatory_bed_paths:
        remaining_limitations.append(
            "Regulatory annotation is heuristic (±2000 bp from gene boundaries)."
            " Provide --regulatory-bed or set REGULATORY_BED env var with"
            " ENCODE TFBS/DHS or ORegAnno BED files from"
            " /mnt/cephfs/hot/pgscat/work/regulatory_derived/."
        )
    else:
        bed_names = [Path(p).name for p in regulatory_bed_paths]
        resolved_limitations.append(
            "Regulatory BED annotation active: element-level overlap for"
            f" {len(regulatory_bed_paths)} source(s): {', '.join(bed_names)}."
            " Fields: regulatory_element_type, regulatory_element_id,"
            " distance_to_regulatory."
        )

    summary = {
        "pgs_id":           pgs_id,
        "annotated_at":     datetime.now(timezone.utc).isoformat(),
        "annotation_tool":  TOOL_VERSION,
        "schema_version":   SCHEMA_VERSION,
        "gff3_reference":   str(gff3_path),
        "fasta_reference":  str(fasta_path) if fasta_path else None,
        "dbnsfp_reference": str(dbnsfp_path) if dbnsfp_path else None,
        "dbsnp_reference":  str(dbsnp_path) if dbsnp_path else None,
        "regulatory_beds":  [str(p) for p in (regulatory_bed_paths or [])],
        "betamap_input":    str(betamap_path),
        "elapsed_seconds":  round(elapsed_s, 1),
        "stats":            stats,
        "per_chromosome":   chrom_counts,
        "resolved_limitations":   resolved_limitations,
        "remaining_limitations":  remaining_limitations,
    }

    out_path = outdir / f"{pgs_id}_annotation_summary.json"
    logger.info("Writing summary: %s", out_path)
    out_path.write_text(json.dumps(summary, indent=2, default=str), encoding="utf-8")
    return summary


# ── CLI ───────────────────────────────────────────────────────────────────────

def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=f"Annotate PGS betamap variants with GENCODE GFF3. ({TOOL_VERSION})",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--version", action="version", version=TOOL_VERSION)
    p.add_argument("--betamap", required=True, type=Path,
                   help="Path to <PGS_ID>_hmPOS_GRCh38.betamap.tsv.gz")
    p.add_argument("--gff3",    required=True, type=Path,
                   help="Path to GENCODE GFF3 (plain or .gz)")
    p.add_argument("--outdir",  required=True, type=Path,
                   help="Output directory (created if absent)")
    p.add_argument("--pgs-id",  required=True, dest="pgs_id",
                   help="PGS identifier (e.g. PGS000001)")

    # Optional inputs
    p.add_argument(
        "--fasta", type=Path, default=None,
        metavar="FASTA",
        help=(
            "Path to GRCh38 reference FASTA (.fa or .fa.gz). "
            "Enables coding-consequence annotation (missense, synonymous, stop-gained, "
            "frameshift) and splice donor/acceptor classification via GT/AG inspection."
        ),
    )
    p.add_argument(
        "--dbnsfp", type=Path, default=None,
        metavar="DBNSFP",
        help=(
            "Path to dbNSFP5 tabix-indexed file (dbNSFP5.0a_grch38.gz). "
            "Enables CADD, REVEL, SIFT, PolyPhen2, and ClinVar score annotation."
        ),
    )
    p.add_argument(
        "--dbsnp", type=Path, default=None,
        metavar="DBSNP",
        help=(
            "Path to dbSNP population frequency VCF (freq.vcf.gz, tabix-indexed). "
            "Enables rsid, af_global, af_max_population, af_population_summary, and "
            "rarity_class annotation. Complementary to --dbnsfp; both can be provided."
        ),
    )
    p.add_argument(
        "--regulatory-bed", nargs="*", type=Path, default=None,
        dest="regulatory_beds",
        metavar="BED",
        help=(
            "One or more regulatory element BED files (plain or .gz). "
            "Examples: ENCODE TFBS clusters, DNase DHS clusters, ORegAnno. "
            "Can also be set via env var REGULATORY_BED (comma- or space-separated). "
            "Variants overlapping these elements get is_regulatory=True, "
            "regulatory_element_type, regulatory_element_id filled, "
            "and distance_to_regulatory=0. Non-overlapping variants within 100 kb "
            "also get distance_to_regulatory filled."
        ),
    )
    p.add_argument(
        "--gnomad-dir", type=Path, default=None,
        dest="gnomad_dir",
        metavar="GNOMAD_DIR",
        help=(
            "Path to gnomAD 4.1.1 variants directory containing exome/ and genome/ "
            "subdirectories with per-chromosome vcf.bgz + .tbi files. "
            "Enables af_gnomad, af_gnomad_afr, af_gnomad_amr, af_gnomad_eas, "
            "af_gnomad_nfe, af_gnomad_sas, af_gnomad_mid, rarity_class_gnomad. "
            "gnomAD AFs are higher quality than dbSNP for most variants."
        ),
    )
    p.add_argument(
        "--gnomad-parquet-dir", type=Path, default=None,
        dest="gnomad_parquet_dir",
        metavar="GNOMAD_PARQUET_DIR",
        help=(
            "Path to pre-built gnomAD parquet index (from gnomad_index.py). "
            "When present, enables fast DuckDB batch queries (one query/chrom) "
            "instead of per-variant tabix seeks over CephFS. "
            "Build once with: python gnomad_index.py --gnomad-dir ... --out-dir ..."
        ),
    )
    p.add_argument(
        "--mcps-dir", type=Path, default=None,
        dest="mcps_dir",
        metavar="MCPS_DIR",
        help=(
            "Path to MCPS per-chromosome TSV.GZ source directory "
            "(contains chrN.freq.tsv.gz files). "
            "Enables af_mcps, af_mcps_mex, af_mcps_eur, af_mcps_afr, "
            "an_mcps, ac_mcps, rarity_class_mcps. "
            "A per-run parquet index is built in outdir/mcps_idx/ on first use. "
            "Used together with --gnomad-dir to compute af_effective and af_source."
        ),
    )
    p.add_argument(
        "--clinvar", type=Path, default=None,
        dest="clinvar_vcf",
        metavar="CLINVAR_VCF",
        help=(
            "Path to ClinVar GRCh38 VCF (clinvar.vcf.gz, BGZF). "
            "Loads the full ClinVar into memory (~240MB RAM) and annotates "
            "clinvar_clnsig (authoritative), clinvar_clndn, clinvar_clnrevstat, "
            "clinvar_alleleid. Replaces/enriches the dbNSFP5-derived clinvar_clnsig."
        ),
    )

    p.add_argument("--no-parquet", action="store_true",
                   help="Skip Parquet output (write TSV.GZ only)")
    p.add_argument("--verbose",    action="store_true",
                   help="Set log level to DEBUG")
    return p


def main(argv=None) -> int:
    parser = _build_parser()
    args   = parser.parse_args(argv)

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    # ── Validate required inputs ──────────────────────────────────────────────
    if not args.betamap.exists():
        logger.error("Betamap not found: %s", args.betamap)
        return 1
    if not args.gff3.exists():
        logger.error("GFF3 not found: %s", args.gff3)
        return 1

    # ── Validate optional FASTA ───────────────────────────────────────────────
    if args.fasta:
        if not args.fasta.exists():
            logger.error("FASTA not found: %s", args.fasta)
            return 1
    else:
        logger.info("FASTA reference not provided — coding consequences will be"
                    " limited to region-level (coding_sequence_variant).")

    # ── Validate optional dbNSFP5 ─────────────────────────────────────────────
    if args.dbnsfp:
        if not args.dbnsfp.exists():
            logger.error("dbNSFP5 file not found: %s", args.dbnsfp)
            return 1
        tbi = Path(str(args.dbnsfp) + ".tbi")
        if not tbi.exists():
            logger.error("dbNSFP5 tabix index not found: %s", tbi)
            return 1
    else:
        logger.info("dbNSFP5 not provided — CADD/REVEL/SIFT/PolyPhen2/ClinVar will be absent.")

    # ── Validate optional dbSNP frequency VCF (v1.3) ─────────────────────────
    if args.dbsnp:
        if not args.dbsnp.exists():
            logger.error("dbSNP freq VCF not found: %s", args.dbsnp)
            return 1
        tbi_snp = Path(str(args.dbsnp) + ".tbi")
        csi_snp = Path(str(args.dbsnp) + ".csi")
        if not tbi_snp.exists() and not csi_snp.exists():
            logger.error("dbSNP freq VCF tabix index not found: %s(.tbi|.csi)", args.dbsnp)
            return 1
    else:
        logger.info(
            "dbSNP freq VCF not provided — rsid/AF/rarity_class will be absent. "
            "Provide --dbsnp /dbsnp_freq/freq.vcf.gz."
        )

    # ── Validate optional regulatory BEDs ────────────────────────────────────
    # CLI --regulatory-bed takes precedence; fall back to REGULATORY_BED env var.
    regulatory_beds: List[Path] = []
    if args.regulatory_beds:
        for bed in args.regulatory_beds:
            if not bed.exists():
                logger.error("Regulatory BED not found: %s", bed)
                return 1
            regulatory_beds.append(bed)
    elif not regulatory_beds:
        import os as _os
        reg_env = _os.environ.get("REGULATORY_BED", "").strip()
        if reg_env:
            logger.info("REGULATORY_BED env var: %s", reg_env)
            for p_str in reg_env.replace(",", " ").split():
                p = Path(p_str.strip())
                if p.exists():
                    regulatory_beds.append(p)
                else:
                    logger.warning("REGULATORY_BED: file not found — skipping: %s", p)

    args.outdir.mkdir(parents=True, exist_ok=True)

    pgs_id = args.pgs_id.upper()
    logger.info("=== Variant Annotator v1.2  PGS_ID=%s ===", pgs_id)

    t_start = time.monotonic()

    # ── Load betamap ──────────────────────────────────────────────────────────
    try:
        betamap = load_betamap(args.betamap)
    except Exception as exc:
        logger.error("Failed to load betamap: %s", exc)
        return 1

    # ── Build GFF3 index ──────────────────────────────────────────────────────
    idx = GFF3Index()
    try:
        idx.build_from_gff3(args.gff3)
    except Exception as exc:
        logger.error("Failed to build GFF3 index: %s", exc)
        return 1

    if idx.n_features == 0:
        logger.error("GFF3 index is empty — check file format and path.")
        return 1

    # ── Load regulatory BEDs (optional) ──────────────────────────────────────
    reg_idx = RegulatoryBEDIndex()
    for bed_path in regulatory_beds:
        try:
            n = reg_idx.load_bed(bed_path)
            logger.info("Loaded %d regulatory elements from %s", n, bed_path.name)
        except Exception as exc:
            logger.error("Failed to load regulatory BED %s: %s", bed_path, exc)
            return 1

    if reg_idx.is_empty:
        logger.info("No regulatory BED files loaded — using heuristic ±2 kb upstream window.")
    else:
        logger.info(
            "Regulatory BED index: %d elements from %d source(s): %s",
            reg_idx.n_features, len(reg_idx.sources), ", ".join(reg_idx.sources),
        )

    # ── Initialise FASTA engine (optional) ───────────────────────────────────
    fasta_engine: Optional[FASTAEngine] = None
    if args.fasta:
        try:
            fasta_engine = FASTAEngine(args.fasta)
        except Exception as exc:
            logger.error("Failed to open FASTA %s: %s", args.fasta, exc)
            return 1

    # ── Initialise dbNSFP5 fetcher (optional) ─────────────────────────────────
    dbnsfp_fetcher: Optional[DbNSFP5Fetcher] = None
    if args.dbnsfp:
        try:
            dbnsfp_fetcher = DbNSFP5Fetcher(args.dbnsfp)
        except Exception as exc:
            logger.error("Failed to open dbNSFP5 %s: %s", args.dbnsfp, exc)
            return 1

    # ── Initialise dbSNP frequency fetcher (optional, v1.3) ───────────────────
    dbsnp_fetcher: Optional[DbSNPFreqFetcher] = None
    if args.dbsnp:
        try:
            dbsnp_fetcher = DbSNPFreqFetcher(args.dbsnp)
        except Exception as exc:
            logger.error("Failed to open dbSNP freq VCF %s: %s", args.dbsnp, exc)
            return 1

    # ── Initialise gnomAD 4.1.1 annotator (optional, v1.4) ───────────────────
    gnomad_annotator = None
    if getattr(args, "gnomad_dir", None):
        try:
            from gnomad_freq import GnomadAnnotator
            _gnomad_parquet_dir = getattr(args, "gnomad_parquet_dir", None)
            gnomad_annotator = GnomadAnnotator(
                gnomad_dir=str(args.gnomad_dir),
                parquet_dir=str(_gnomad_parquet_dir) if _gnomad_parquet_dir else None,
            )
            if _gnomad_parquet_dir:
                logger.info(
                    "gnomAD 4.1.1 annotator: parquet fast path at %s", _gnomad_parquet_dir
                )
            else:
                logger.info(
                    "gnomAD 4.1.1 annotator: tabix fallback (no parquet index) — "
                    "run gnomad_index.py and set GNOMAD_PARQUET_DIR for faster annotation"
                )
        except Exception as exc:
            logger.warning("gnomAD annotator init failed: %s", exc)

    # ── Initialise MCPS annotator (optional, v1.6) ────────────────────────────
    mcps_annotator = None
    if getattr(args, "mcps_dir", None):
        try:
            from mcps_freq import McpsAnnotator
            mcps_index_dir = args.outdir / "mcps_idx"
            mcps_annotator = McpsAnnotator(
                source_dir=str(args.mcps_dir),
                index_dir=str(mcps_index_dir),
            )
            # Pre-build parquet index for chromosomes present in the betamap
            betamap_chroms = list(betamap["CHROM"].unique())
            logger.info(
                "MCPS: building index for %d chromosome(s): %s",
                len(betamap_chroms), ", ".join(sorted(betamap_chroms)),
            )
            mcps_annotator.build_for_chroms(betamap_chroms)
            logger.info("MCPS annotator active: %s", args.mcps_dir)
        except Exception as exc:
            logger.warning("MCPS annotator init failed: %s", exc)
            mcps_annotator = None

    # ── Initialise ClinVar annotator (optional, v1.4) ────────────────────────
    clinvar_annotator = None
    if getattr(args, "clinvar_vcf", None):
        try:
            from clinvar_anno import ClinVarAnnotator
            clinvar_annotator = ClinVarAnnotator(str(args.clinvar_vcf))
            if clinvar_annotator.is_loaded():
                logger.info(
                    "ClinVar annotator active: %d variants loaded", len(clinvar_annotator)
                )
            else:
                logger.warning("ClinVar annotator failed to load — skipping")
                clinvar_annotator = None
        except Exception as exc:
            logger.warning("ClinVar annotator init failed: %s", exc)

    # ── Annotate ──────────────────────────────────────────────────────────────
    try:
        annotated = annotate(
            betamap,
            idx,
            reg_idx=reg_idx if not reg_idx.is_empty else None,
            fasta_engine=fasta_engine,
            dbnsfp_fetcher=dbnsfp_fetcher,
            dbsnp_fetcher=dbsnp_fetcher,
        )
    except Exception as exc:
        logger.exception("Annotation failed: %s", exc)
        return 1

    # ── Log dbSNP summary ─────────────────────────────────────────────────────
    if dbsnp_fetcher:
        dbsnp_fetcher.log_summary()

    # ── gnomAD 4.1.1 enrichment ───────────────────────────────────────────────
    # Fast path  : parquet index (build once with gnomad_index.py) → DuckDB,
    #              one query per chromosome — eliminates CephFS seek latency.
    # Slow path  : pysam tabix range queries sorted by position (fallback).
    # htslib "index file is older than data file" warnings are non-fatal.
    if gnomad_annotator is not None:
        t_gnomad = time.monotonic()
        gnomad_cols = [
            "af_gnomad", "af_gnomad_afr", "af_gnomad_amr", "af_gnomad_eas",
            "af_gnomad_nfe", "af_gnomad_sas", "af_gnomad_mid",
            "an_gnomad", "nhomalt_gnomad", "rarity_class_gnomad",
        ]
        _g_col_data: dict = {c: [None] * len(annotated) for c in gnomad_cols}

        # Vectorized extraction (avoids slow iterrows per variant)
        _g_rst   = annotated.reset_index(drop=True)
        _g_chr   = _g_rst["CHROM"].astype(str).to_numpy()
        _g_pos   = _g_rst["POS"].astype(int).to_numpy()
        _g_flip  = _g_rst["IS_FLIP"].fillna(0).astype(bool).to_numpy()
        _g_ea    = _g_rst["EFFECT_ALLELE"].fillna("").astype(str).to_numpy()
        _g_oa    = _g_rst["OTHER_ALLELE"].fillna("").astype(str).to_numpy()

        from collections import defaultdict as _ddict
        _g_by_chrom: dict = _ddict(list)
        for _gi in range(len(_g_rst)):
            _g_by_chrom[_g_chr[_gi]].append(_gi)

        _g_n_matched = 0

        for _g_chrom, _g_ilocs in sorted(_g_by_chrom.items()):
            # Build (pos, ref_u, alt_u, iloc) list for this chromosome
            _g_v4c = []
            for _gi in _g_ilocs:
                _gf  = _g_flip[_gi]
                _gea, _goa = _g_ea[_gi], _g_oa[_gi]
                _gr  = (_goa if not _gf else _gea).upper()
                _ga  = (_gea if not _gf else _goa).upper()
                _g_v4c.append((int(_g_pos[_gi]), _gr, _ga, _gi))

            # ── Fast path: DuckDB parquet ────────────────────────────────────
            if gnomad_annotator.has_parquet(_g_chrom):
                logger.info(
                    "  gnomAD chr%s: %d variants → DuckDB parquet (fast path)",
                    _g_chrom, len(_g_v4c),
                )
                _g_lookup = list({(p, r, a) for p, r, a, _ in _g_v4c})
                _g_buf    = gnomad_annotator.fetch_batch_parquet(_g_chrom, _g_lookup)

            # ── Slow path: tabix range queries (sorted by position) ──────────
            else:
                _g_v4c.sort(key=lambda x: x[0])
                _g_MERGE = 100_000
                _g_buf   = {}
                _g_ranges: list = []
                _grs = _gre = _g_v4c[0][0]; _gri = [_g_v4c[0]]
                for _gv in _g_v4c[1:]:
                    if _gv[0] - _gre <= _g_MERGE:
                        _gre = _gv[0]; _gri.append(_gv)
                    else:
                        _g_ranges.append((_grs, _gre, _gri))
                        _grs = _gre = _gv[0]; _gri = [_gv]
                _g_ranges.append((_grs, _gre, _gri))
                logger.info(
                    "  gnomAD chr%s: %d variants → %d tabix range quer%s (slow path — "
                    "run gnomad_index.py to enable fast path)",
                    _g_chrom, len(_g_v4c), len(_g_ranges),
                    "y" if len(_g_ranges) == 1 else "ies",
                )
                for _grstart, _grend, _grvars in _g_ranges:
                    _g_buf.update(gnomad_annotator.fetch_range(_g_chrom, _grstart, _grend))

            # ── Map results → col_data ───────────────────────────────────────
            for _gpos, _gref, _galt, _gi in _g_v4c:
                _gres = _g_buf.get((_gpos, _gref, _galt)) or _g_buf.get((_gpos, _galt, _gref))
                if _gres:
                    _g_n_matched += 1
                    for _gc in gnomad_cols:
                        _g_col_data[_gc][_gi] = _gres.get(_gc)

        for c in gnomad_cols:
            annotated[c] = _g_col_data[c]
        gnomad_annotator.n_matched = _g_n_matched
        gnomad_annotator.n_novel   = len(annotated) - _g_n_matched
        gnomad_annotator.log_summary()
        logger.info(
            "gnomAD enrichment done in %.1f s  (matched %d / %d)",
            time.monotonic() - t_gnomad, _g_n_matched, len(annotated),
        )

    # ── MCPS enrichment (optional, v1.6) ─────────────────────────────────────
    if mcps_annotator is not None:
        logger.info("Enriching with MCPS population frequencies ...")
        mcps_cols = [
            "af_mcps", "af_mcps_mex", "af_mcps_eur", "af_mcps_afr",
            "an_mcps", "ac_mcps", "rarity_class_mcps",
        ]
        mcps_col_data: dict = {c: [None] * len(annotated) for c in mcps_cols}
        for i, row in annotated.iterrows():
            chrom = str(row.get("CHROM", ""))
            pos   = int(row.get("POS", 0))
            flip  = bool(row.get("IS_FLIP", False))
            ea    = str(row.get("EFFECT_ALLELE", "") or "")
            oa    = str(row.get("OTHER_ALLELE", "") or "")
            ref, alt = (oa, ea) if flip else (ea, oa)
            res = mcps_annotator.fetch(chrom, pos, ref, alt)
            if res:
                idx_loc = annotated.index.get_loc(i)
                for c in mcps_cols:
                    mcps_col_data[c][idx_loc] = res.get(c)
        for c in mcps_cols:
            annotated[c] = mcps_col_data[c]
        mcps_annotator.log_summary()

    # ── Compute af_effective / af_source / rarity_class_effective (v1.6) ─────
    # Priority: MCPS MEX > MCPS RAW > gnomAD > dbSNP (af_global)
    # MCPS reflects the Mexican population better than gnomAD for MCPS cohorts.
    # Only computed when at least one frequency source was run.
    if gnomad_annotator is not None or mcps_annotator is not None or dbsnp_fetcher is not None:
        logger.info("Computing af_effective (priority: MCPS MEX > MCPS RAW > gnomAD > dbSNP) …")
        from mcps_freq import classify_rarity_mcps

        af_eff_list:  list = [None] * len(annotated)
        af_src_list:  list = [None] * len(annotated)
        rar_eff_list: list = [None] * len(annotated)

        for iloc, (_, row) in enumerate(annotated.iterrows()):
            af_mex = row.get("af_mcps_mex")
            af_raw = row.get("af_mcps")
            af_g   = row.get("af_gnomad")
            af_db  = row.get("af_global")     # dbSNP column from v1.3

            if af_mex is not None and af_mex >= 0:
                af_eff_list[iloc] = af_mex
                af_src_list[iloc] = "mcps_mex"
            elif af_raw is not None and af_raw >= 0:
                af_eff_list[iloc] = af_raw
                af_src_list[iloc] = "mcps_raw"
            elif af_g is not None and af_g >= 0:
                af_eff_list[iloc] = af_g
                af_src_list[iloc] = "gnomad"
            elif af_db is not None and af_db >= 0:
                af_eff_list[iloc] = af_db
                af_src_list[iloc] = "dbsnp"

            rar_eff_list[iloc] = classify_rarity_mcps(af_eff_list[iloc])

        annotated["af_effective"]           = af_eff_list
        annotated["af_source"]              = af_src_list
        annotated["rarity_class_effective"] = rar_eff_list

    # ── ClinVar enrichment (optional, v1.4) ───────────────────────────────────
    if clinvar_annotator is not None:
        logger.info("Enriching with ClinVar GRCh38 ...")
        cv_col_map = {
            "clinvar_clnsig":     "clinvar_clnsig",
            "clinvar_clndn":      "clinvar_clndn",
            "clinvar_clnrevstat": "clinvar_clnrevstat",
            "clinvar_alleleid":   "clinvar_alleleid",
            "clinvar_clnvc":      "clinvar_clnvc",
        }
        cv_col_data: dict = {c: [None] * len(annotated) for c in cv_col_map}
        n_found = 0
        for i, row in annotated.iterrows():
            chrom = str(row.get("CHROM", ""))
            pos   = int(row.get("POS", 0))
            flip  = bool(row.get("IS_FLIP", False))
            ea    = str(row.get("EFFECT_ALLELE", "") or "")
            oa    = str(row.get("OTHER_ALLELE", "") or "")
            ref, alt = (oa, ea) if flip else (ea, oa)
            res = clinvar_annotator.lookup(chrom, pos, ref, alt)
            if res:
                n_found += 1
                idx_loc = annotated.index.get_loc(i)
                for cv_col, res_key in cv_col_map.items():
                    cv_col_data[cv_col][idx_loc] = res.get(res_key)
        for cv_col in cv_col_map:
            # Overwrite dbNSFP-derived clinvar_clnsig only where ClinVar found a result
            if cv_col == "clinvar_clnsig" and cv_col in annotated.columns:
                # Keep dbNSFP value where ClinVar has nothing
                mask = [v is not None for v in cv_col_data[cv_col]]
                annotated.loc[mask, cv_col] = [v for v in cv_col_data[cv_col] if v is not None]
            else:
                annotated[cv_col] = cv_col_data[cv_col]
        logger.info("ClinVar: %d / %d variants matched", n_found, len(annotated))

    # ── Write outputs ─────────────────────────────────────────────────────────
    # Only include columns that are actually present (new columns may be absent
    # on older DataFrames in future edge cases).
    output_cols = [c for c in OUTPUT_COLUMNS if c in annotated.columns]

    tsv_path = args.outdir / f"{pgs_id}_variants_annotated.tsv.gz"
    write_tsv_gz(annotated[output_cols], tsv_path)

    if not args.no_parquet:
        parquet_path = args.outdir / f"{pgs_id}_variants_annotated.parquet"
        write_parquet(annotated[output_cols], parquet_path)

    elapsed = time.monotonic() - t_start
    summary = write_summary(
        df=annotated,
        pgs_id=pgs_id,
        betamap_path=args.betamap,
        gff3_path=args.gff3,
        outdir=args.outdir,
        elapsed_s=elapsed,
        fasta_path=args.fasta,
        dbnsfp_path=args.dbnsfp,
        dbsnp_path=args.dbsnp,
        regulatory_bed_paths=regulatory_beds if regulatory_beds else None,
        reg_idx=reg_idx,
    )

    # ── Print summary to stdout ───────────────────────────────────────────────
    stats = summary["stats"]
    print("\n── Annotation Summary ──────────────────────────────────────────")
    print(f"  PGS ID:           {pgs_id}")
    print(f"  Total variants:   {stats['total_variants']:>10,}")
    print(f"  Coding:           {stats['n_coding']:>10,}  ({stats['pct_coding']:.1f}%)")
    print(f"  Intergenic:       {stats['n_intergenic']:>10,}  ({stats['pct_intergenic']:.1f}%)")
    print(f"  Regulatory:       {stats['n_regulatory']:>10,}")
    print(f"  Splice site:      {stats['n_splice_site']:>10,}")
    print(f"  Splice donor:     {stats.get('n_splice_donor', 0):>10,}")
    print(f"  Splice acceptor:  {stats.get('n_splice_acceptor', 0):>10,}")
    print(f"  Splice region:    {stats['n_splice_region']:>10,}")
    print(f"  Missense:         {stats.get('n_missense', 0):>10,}  ({stats.get('pct_missense', 0):.1f}%)")
    print(f"  Synonymous:       {stats.get('n_synonymous', 0):>10,}")
    print(f"  LoF:              {stats.get('n_lof', 0):>10,}  ({stats.get('pct_lof', 0):.1f}%)")
    print(f"  Stop gained:      {stats.get('n_stop_gained', 0):>10,}")
    print(f"  Frameshift:       {stats.get('n_frameshift', 0):>10,}")
    print()
    print("  Region class breakdown:")
    for rc, cnt in stats["region_class_counts"].items():
        pct = 100 * cnt / stats["total_variants"] if stats["total_variants"] else 0
        print(f"    {rc:<28s} {cnt:>8,}  ({pct:.1f}%)")
    print(f"\n  Elapsed: {elapsed:.1f} s")
    print(f"  FASTA:   {args.fasta or '(not provided)'}")
    print(f"  dbNSFP5: {args.dbnsfp or '(not provided)'}")
    if args.dbsnp and dbsnp_fetcher:
        print(f"  dbSNP:   {args.dbsnp}")
        print(f"    matched={dbsnp_fetcher.n_matched}  "
              f"af_computed={dbsnp_fetcher.n_af_found}  "
              f"novel={dbsnp_fetcher.n_novel}")
    else:
        print(f"  dbSNP:   (not provided)")
    print(f"  Output:  {args.outdir}")
    print("────────────────────────────────────────────────────────────────\n")

    return 0


if __name__ == "__main__":
    sys.exit(main())
