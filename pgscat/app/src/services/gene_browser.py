"""
GeneBrowserService
==================
Provides gene-centric interpretation data for the gene browser view.

This service reads from existing annotation parquet files (produced by the
Variant Annotator pipeline) and builds structured track data for Plotly-based
genomic visualisation.

Design decision: Plotly native (not LocusZoom.js)
--------------------------------------------------
Rationale for choosing Plotly over LocusZoom.js:
  - Plotly is already loaded in _base.tpl; zero additional dependencies
  - plotly.subplots.make_subplots + shapes covers all required track types:
    scatter for variants, rectangles for exons, lines for introns
  - LocusZoom.js is optimised for GWAS locus-zoom plots and requires
    specific data-source adapters; adapting it to a gene-centric browser
    would require significant JS plumbing with no benefit over Plotly
  - All data is served as JSON from Bottle — Plotly consumes it natively
  - No CDN dependency on external genomics infrastructure

Track architecture (3 main tracks):
  Track 1 — Variants scatter
    x = genomic position (POS)
    y = consequence priority (lower = higher impact) or CADD/REVEL
    color = consequence category
    size = |BETA| or AF-derived

  Track 2 — Gene model
    Derived from annotation data: coding positions → exon blocks
    Intronic positions → connecting lines
    Approximated from parquet (exact model needs GFF3 cache)

  Track 3 — Functional summary
    Per-segment consequence distribution: heatmap / bar chart
    Rarity distribution across gene

NOTE: The track MTR (Missense Tolerance Ratio w31) is NOT implemented.
"""
from __future__ import annotations

import logging
import math
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from services.clinical_genes import get_clinical_genes_service

logger = logging.getLogger(__name__)

# ── Consequence colour map (consistent with the gene browser UI) ─────────────
CONSEQUENCE_COLORS = {
    # LoF
    "stop_gained":           "#d32f2f",
    "frameshift_variant":    "#b71c1c",
    "splice_donor_variant":  "#e64a19",
    "splice_acceptor_variant":"#e64a19",
    "start_lost":            "#c62828",
    "stop_lost":             "#ad1457",
    # Missense
    "missense_variant":      "#f57c00",
    "inframe_insertion":     "#ff8f00",
    "inframe_deletion":      "#ff8f00",
    # Synonymous / coding
    "synonymous_variant":    "#388e3c",
    "coding_sequence_variant":"#43a047",
    # Splice region
    "splice_site_variant":   "#e65100",
    "splice_region_variant": "#fb8c00",
    # Non-coding
    "5_prime_UTR_variant":   "#7b1fa2",
    "3_prime_UTR_variant":   "#8e24aa",
    "intron_variant":        "#546e7a",
    "non_coding_exon_variant":"#78909c",
    "upstream_gene_variant": "#b0bec5",
    "regulatory_region_variant":"#0097a7",
    # Intergenic
    "intergenic_variant":    "#9e9e9e",
}

# Consequence category grouping (for legend)
CONSEQUENCE_CATEGORIES = {
    "pLoF": [
        "stop_gained", "frameshift_variant",
        "splice_donor_variant", "splice_acceptor_variant",
        "start_lost", "stop_lost",
    ],
    "missense": ["missense_variant", "inframe_insertion", "inframe_deletion"],
    "synonymous": ["synonymous_variant", "coding_sequence_variant"],
    "splice_region": ["splice_site_variant", "splice_region_variant"],
    "non_coding": [
        "5_prime_UTR_variant", "3_prime_UTR_variant",
        "non_coding_exon_variant", "intron_variant",
        "upstream_gene_variant", "regulatory_region_variant",
        "intergenic_variant",
    ],
}

# Reverse map: consequence → category
_CSQ_TO_CAT = {}
for _cat, _csqs in CONSEQUENCE_CATEGORIES.items():
    for _csq in _csqs:
        _CSQ_TO_CAT[_csq] = _cat


def _category_color(consequence: str) -> str:
    if consequence in CONSEQUENCE_COLORS:
        return CONSEQUENCE_COLORS[consequence]
    return "#9e9e9e"


def _safe_float(v) -> Optional[float]:
    if v is None:
        return None
    try:
        f = float(v)
        return None if (math.isnan(f) or math.isinf(f)) else f
    except (TypeError, ValueError):
        return None


class GeneBrowserService:
    """
    Gene-centric data service for the genomic browser view.

    Reads from the annotations directory to collect all variants annotated
    for a given gene symbol (across all available PGS parquets). Builds
    structured track data for Plotly subplots rendering.
    """

    def __init__(self, cfg) -> None:
        self.cfg = cfg
        self.annotations_dir: Path = cfg.ANNOTATIONS_DIR

    # ── Public API ──────────────────────────────────────────────────────────────

    def get_gene_info(self, gene_symbol: str) -> Dict[str, Any]:
        """
        Return summary information about a gene: genomic range, variant counts,
        list of PGS IDs that contain variants in this gene.
        """
        df = self._load_gene_variants(gene_symbol)
        if df is None:
            cg_svc = get_clinical_genes_service()
            cg_entry = cg_svc.get_gene_info(gene_symbol) or {}
            return {
                "gene_symbol":       gene_symbol,
                "found":             False,
                "is_clinical_gene":  cg_svc.is_clinical_gene(gene_symbol),
                "clinical_sources":  cg_svc.get_clinical_sources(gene_symbol),
                "clinical_confidence": cg_svc.get_confidence(gene_symbol),
                "clinical_moi":      cg_entry.get("moi"),
                "clinical_evidence": cg_entry.get("evidence"),
                "clinical_disease":  cg_entry.get("disease"),
                "message": (
                    f"No annotated variants found for gene '{gene_symbol}'. "
                    "Run the annotation pipeline and ensure the parquet files are available."
                ),
            }

        import pandas as pd

        chrom = df["CHROM"].mode().iloc[0] if not df.empty else None
        pos_min = int(df["POS"].min()) if not df.empty else None
        pos_max = int(df["POS"].max()) if not df.empty else None
        n_variant_records = len(df)
        pgs_ids = sorted(df["PRS_ID"].unique().tolist()) if "PRS_ID" in df.columns else []

        # Unique genomic positions across all PGS files
        _dedup_cols = [c for c in ("CHROM", "POS", "EFFECT_ALLELE") if c in df.columns]
        n_unique_variants = (
            len(df.drop_duplicates(subset=_dedup_cols))
            if _dedup_cols else n_variant_records
        )

        consequence_counts = {}
        if "consequence" in df.columns:
            consequence_counts = {k: int(v) for k, v in df["consequence"].value_counts().items()}

        rarity_counts = {}
        if "rarity_class" in df.columns:
            rarity_counts = {k: int(v) for k, v in
                             df["rarity_class"].fillna("novel").value_counts().items()}

        # ── Clinical gene annotation ──────────────────────────────────────────
        cg_svc = get_clinical_genes_service()
        is_clinical   = cg_svc.is_clinical_gene(gene_symbol)
        cg_info       = cg_svc.get_gene_info(gene_symbol) or {}
        cg_sources    = cg_svc.get_clinical_sources(gene_symbol)
        cg_confidence = cg_svc.get_confidence(gene_symbol)

        return {
            "gene_symbol":       gene_symbol,
            "found":             True,
            "chrom":             str(chrom) if chrom else None,
            "genomic_start":     pos_min,
            "genomic_end":       pos_max,
            "genomic_span_kb":   round((pos_max - pos_min) / 1000, 1) if pos_min and pos_max else None,
            "n_variants":        n_unique_variants,   # unique positions (CHROM,POS,EFFECT_ALLELE)
            "n_variant_records": n_variant_records,   # raw rows across all PGS files
            "n_pgs_ids":         len(pgs_ids),
            "pgs_ids":           pgs_ids[:20],   # cap for display
            "consequence_counts": consequence_counts,
            "rarity_counts":      rarity_counts,
            # Clinical gene fields
            "is_clinical_gene":   is_clinical,
            "clinical_sources":   cg_sources,
            "clinical_confidence": cg_confidence,
            "clinical_moi":       cg_info.get("moi"),
            "clinical_evidence":  cg_info.get("evidence"),
            "clinical_disease":   cg_info.get("disease"),
        }

    def get_tracks(self, gene_symbol: str) -> Dict[str, Any]:
        """
        Build all track data for the gene browser Plotly figure.
        NOTE: MTR (Missense Tolerance Ratio w31) track is intentionally excluded.
        """
        df = self._load_gene_variants(gene_symbol)
        if df is None:
            return {"gene_symbol": gene_symbol, "found": False, "tracks": {}}

        df = df.sort_values("POS").reset_index(drop=True)

        # Live-enrich with gnomAD / ClinVar / MCPS if not already present
        df = self._enrich_gene_variants(df)

        chrom   = str(df["CHROM"].mode().iloc[0]) if not df.empty else "?"
        pos_min = int(df["POS"].min())
        pos_max = int(df["POS"].max())
        span    = max(pos_max - pos_min, 1_000)
        pad     = max(int(span * 0.04), 500)
        x_min   = pos_min - pad
        x_max   = pos_max + pad

        # Deduplicate by genomic position across PGS: keep highest-CADD row per pos
        df_uniq = self._deduplicate_variants(df)

        tracks = {
            "variants":    self._build_variant_track(df_uniq),
            "gene_model":  self._build_gene_model_track(df_uniq, x_min, x_max),
            "density":     self._build_functional_summary_track(df_uniq, x_min, x_max),
            "af":          self._build_af_track(df_uniq),
        }

        strand = ""
        if "strand" in df_uniq.columns:
            sv = df_uniq["strand"].dropna().unique()
            strand = str(sv[0]) if len(sv) else ""

        return {
            "gene_symbol":    gene_symbol,
            "found":          True,
            "chrom":          chrom,
            "genomic_start":  pos_min,
            "genomic_end":    pos_max,
            "strand":         strand,
            "n_total":        len(df),
            "n_unique":       len(df_uniq),
            "tracks":         tracks,
            "layout_hints": {
                "x_range": [x_min, x_max],
                "x_title": f"chr{chrom} position (GRCh38)",
                # Domains defined in JS — not here — so they can react to window width
                "excluded_tracks": ["Missense Tolerance Ratio (MTR) w31"],
            },
        }

    @staticmethod
    def _deduplicate_variants(df):
        """
        Collapse duplicate genomic positions across PGS IDs.
        For each unique (CHROM, POS, EFFECT_ALLELE) keep the row with the
        highest CADD score (or first row if CADD absent).
        """
        if "cadd_phred" in df.columns:
            df = df.sort_values("cadd_phred", ascending=False, na_position="last")
        return df.drop_duplicates(
            subset=["CHROM", "POS", "EFFECT_ALLELE"], keep="first"
        ).reset_index(drop=True)

    # ── Track builders ──────────────────────────────────────────────────────────

    # Maximum variants to render in the scatter track (priority-sampled)
    _MAX_PLOT_VARIANTS = 3_000

    # Consequence-based CADD fallback values (used when cadd_phred is null)
    # Approximates pathogenicity on the CADD phred scale for grouping
    _CSQ_CADD_FALLBACK = {
        "stop_gained": 40, "frameshift_variant": 38,
        "splice_donor_variant": 35, "splice_acceptor_variant": 35,
        "start_lost": 35, "stop_lost": 32,
        "missense_variant": 22, "inframe_insertion": 18, "inframe_deletion": 18,
        "splice_site_variant": 16, "splice_region_variant": 14,
        "synonymous_variant": 8, "coding_sequence_variant": 8,
        "5_prime_UTR_variant": 5, "3_prime_UTR_variant": 5,
        "intron_variant": 2, "non_coding_exon_variant": 4,
        "upstream_gene_variant": 1, "regulatory_region_variant": 3,
        "intergenic_variant": 0,
    }

    _CAT_COLORS = {
        "pLoF":         "#e53935",
        "missense":     "#fb8c00",
        "synonymous":   "#43a047",
        "splice_region":"#8e24aa",
        "non_coding":   "#78909c",
    }

    def _build_variant_track(self, df) -> Dict[str, Any]:
        """
        Track 1 — Variant scatter (redesigned).

        Y-axis: CADD phred score (primary pathogenicity metric).
          Variants without CADD get a consequence-based fallback value so
          every variant is visible. A dashed reference line at CADD=20 marks
          the "possibly deleterious" threshold.

        Deduplication: caller (_deduplicate_variants) already collapsed
          multi-PGS rows to unique genomic positions.

        Sampling: if more than _MAX_PLOT_VARIANTS unique positions remain,
          we priority-sample: always keep pLoF + missense; sample the rest
          proportionally to stay under the cap.

        One Plotly trace per consequence category → clean shared legend.
        """
        import pandas as pd
        import numpy as np

        if df.empty:
            return {"traces": [], "n_variants": 0, "n_sampled": 0}

        # ── Vectorized column extraction ───────────────────────────────────────
        df = df.copy()

        # Fill consequence
        df["_csq"] = df["consequence"].fillna("intergenic_variant").astype(str)
        df["_cat"] = df["_csq"].map(lambda c: _CSQ_TO_CAT.get(c, "non_coding"))

        # Y = CADD phred when available; consequence-based fallback otherwise
        cadd_col = "cadd_phred" if "cadd_phred" in df.columns else None
        if cadd_col:
            df["_y"] = pd.to_numeric(df[cadd_col], errors="coerce")
        else:
            df["_y"] = pd.Series([None] * len(df), dtype=float)

        df["_y_fallback"] = df["_csq"].map(
            lambda c: self._CSQ_CADD_FALLBACK.get(c, 0)
        ).astype(float)
        df["_y_final"] = df["_y"].fillna(df["_y_fallback"])
        df["_has_cadd"] = df["_y"].notna()

        # AF and rarity — use unified effective AF (MCPS>gnomAD>dbSNP) if available
        if "af_effective" not in df.columns:
            from services.variant_annotation import _compute_effective_af
            df = _compute_effective_af(df)
        df["_af"] = pd.to_numeric(df["af_effective"], errors="coerce") if "af_effective" in df.columns else pd.Series([None]*len(df), dtype=float)
        _rarity_col = "rarity_class_effective" if "rarity_class_effective" in df.columns else "rarity_class"
        _rarity_series = df.get(_rarity_col, pd.Series(["novel"]*len(df), index=df.index))
        df["_rarity"] = _rarity_series.fillna("novel").astype(str).str.replace("-", "_")

        # BETA-based marker size: clamp to [5, 18]
        beta_col = "BETA" if "BETA" in df.columns else None
        if beta_col:
            df["_beta_abs"] = pd.to_numeric(df[beta_col], errors="coerce").abs()
            df["_size"] = (5 + (df["_beta_abs"].fillna(0) * 25)).clip(5, 18).round().astype(int)
        else:
            df["_beta_abs"] = pd.Series([None] * len(df), dtype=float)
            df["_size"] = 8

        # REVEL score (raw, for client-side y-axis switching)
        if "revel_score" in df.columns:
            df["_revel_raw"] = pd.to_numeric(df["revel_score"], errors="coerce")
        else:
            df["_revel_raw"] = pd.Series([None] * len(df), dtype=float)

        # Marker symbol: filled circle for CADD-backed, diamond for fallback
        df["_symbol"] = df["_has_cadd"].map(lambda x: "circle" if x else "diamond")

        # Opacity: higher for impactful categories
        df["_opacity"] = df["_cat"].map({
            "pLoF": 0.95, "missense": 0.90, "splice_region": 0.88,
            "synonymous": 0.75, "non_coding": 0.60,
        }).fillna(0.65)

        # Colour per row (exact consequence colour for richer display)
        df["_color"] = df["_csq"].map(_category_color)

        # ── Tooltip ────────────────────────────────────────────────────────────
        def _make_hover(row) -> str:
            gene  = str(row.get("gene_name", "") or "")
            csq   = row["_csq"]
            pos   = int(row["POS"])
            rsid  = str(row.get("rsid", "") or "")
            beta  = row.get("BETA")
            af    = row["_af"]
            rarity= row["_rarity"]
            cadd  = row.get("cadd_phred")
            revel = row.get("revel_score")
            sift  = str(row.get("sift_pred", "") or "")
            pp2   = str(row.get("polyphen2_pred", "") or "")
            clnv  = str(row.get("clinvar_clnsig", "") or "")
            aa    = ""
            if row.get("aa_ref_3") and row.get("aa_alt_3"):
                aa = f"{row['aa_ref_3']}→{row['aa_alt_3']}"
            elif row.get("aa_ref") and row.get("aa_alt"):
                aa = f"{row['aa_ref']}→{row['aa_alt']}"
            pgs = str(row.get("PRS_ID", "") or "")

            parts = [f"<b>{csq}</b>"]
            if gene:     parts.append(f"Gene: {gene}")
            parts.append(f"chr{row.get('CHROM','')}:{pos:,}  {row.get('EFFECT_ALLELE','')}>{row.get('OTHER_ALLELE','')}")
            if rsid:     parts.append(f"rsID: {rsid}")
            if isinstance(beta, float): parts.append(f"β = {beta:+.4f}")
            if isinstance(af, float):   parts.append(f"AF = {af:.2e}  ({rarity})")
            if isinstance(cadd, float): parts.append(f"CADD phred = {cadd:.1f}")
            if isinstance(revel, float):parts.append(f"REVEL = {revel:.3f}")
            if sift:     parts.append(f"SIFT: {sift}")
            if pp2:      parts.append(f"PolyPhen2: {pp2}")
            if aa:       parts.append(f"AA change: {aa}")
            if clnv and clnv not in (".", ""):
                         parts.append(f"ClinVar: {clnv}")
            # gnomAD population AFs
            gaf_amr = _safe_float(row.get("af_gnomad_amr"))
            gaf_afr = _safe_float(row.get("af_gnomad_afr"))
            gaf_nfe = _safe_float(row.get("af_gnomad_nfe"))
            gaf_eas = _safe_float(row.get("af_gnomad_eas"))
            pop_parts = []
            if gaf_amr is not None: pop_parts.append(f"AMR={gaf_amr:.2e}")
            if gaf_afr is not None: pop_parts.append(f"AFR={gaf_afr:.2e}")
            if gaf_nfe is not None: pop_parts.append(f"NFE={gaf_nfe:.2e}")
            if pop_parts: parts.append("gnomAD: " + "  ".join(pop_parts))
            # MCPS
            mcps_mex = _safe_float(row.get("af_mcps_mex"))
            mcps_af  = _safe_float(row.get("af_mcps"))
            if mcps_mex is not None: parts.append(f"MCPS-MEX={mcps_mex:.2e}")
            elif mcps_af is not None: parts.append(f"MCPS={mcps_af:.2e}")
            if pgs:      parts.append(f"PGS: {pgs}")
            return "<br>".join(parts)

        df["_hover"] = df.apply(_make_hover, axis=1)

        # ── Priority sampling if needed ────────────────────────────────────────
        n_total = len(df)
        if n_total > self._MAX_PLOT_VARIANTS:
            high_impact = df[df["_cat"].isin(["pLoF", "missense", "splice_region"])]
            rest = df[~df["_cat"].isin(["pLoF", "missense", "splice_region"])]
            budget = max(0, self._MAX_PLOT_VARIANTS - len(high_impact))
            if budget > 0 and len(rest) > budget:
                rest = rest.sample(n=budget, random_state=42)
            df = pd.concat([high_impact, rest]).sort_values("POS")
        n_sampled = len(df)

        # ── Build one trace per category ───────────────────────────────────────
        traces = []
        cat_order = ["pLoF", "missense", "splice_region", "synonymous", "non_coding"]
        for cat in cat_order:
            sub = df[df["_cat"] == cat]
            if sub.empty:
                continue
            cat_color = self._CAT_COLORS.get(cat, "#9e9e9e")
            # customdata columns (for client-side y-axis switching):
            #   [0] cadd_final  — CADD with consequence fallback
            #   [1] revel       — raw REVEL score (null if unavailable)
            #   [2] beta_abs    — |BETA| (null if unavailable)
            #   [3] af_global   — allele frequency (null if novel)
            cd = list(zip(
                [_safe_float(v) for v in sub["_y_final"].tolist()],
                [_safe_float(v) for v in sub["_revel_raw"].tolist()],
                [_safe_float(v) for v in sub["_beta_abs"].tolist()],
                [_safe_float(v) for v in sub["_af"].tolist()],
            ))
            traces.append({
                "name":           cat,
                "type":           "scatter",
                "mode":           "markers",
                "x":              sub["POS"].tolist(),
                "y":              sub["_y_final"].tolist(),
                "customdata":     cd,
                "text":           sub["_hover"].tolist(),
                "hovertemplate":  "%{text}<extra></extra>",
                "marker": {
                    "color":    sub["_color"].tolist(),
                    "size":     sub["_size"].tolist(),
                    "symbol":   sub["_symbol"].tolist(),
                    "opacity":  float(sub["_opacity"].mean()),
                    "line":     {"width": 0.8, "color": "rgba(255,255,255,0.7)"},
                },
            })

        # CADD=20 reference line (sent as a shape in layout_hints, applied in JS)
        has_cadd = bool(df["_has_cadd"].any())

        return {
            "traces":      traces,
            "n_variants":  n_total,
            "n_sampled":   n_sampled,
            "has_cadd":    has_cadd,
            "cadd_threshold": 20,  # "possibly deleterious"
            "yaxis": {
                "title": "CADD Phred Score  (● = scored, ◆ = estimated)",
                "rangemode": "tozero",
                "showgrid": True,
                "gridcolor": "#eeeeee",
            },
        }

    def _build_gene_model_track(
        self, df, x_min: int, x_max: int
    ) -> Dict[str, Any]:
        """
        Track 2 — Gene model.

        Exon blocks derived from coding variant positions (clustered within
        500 bp). Blocks are widened to a minimum visible size.
        Strand direction shown as arrow overlay.
        UTR-like regions (5'/3' UTR consequence) shown in lighter colour.
        """
        import pandas as pd
        import numpy as np

        gene_name = str(df["gene_name"].mode().iloc[0]) if not df.empty else "?"
        chrom     = str(df["CHROM"].mode().iloc[0]) if not df.empty else "?"
        strand    = "."
        if "strand" in df.columns:
            sv = df["strand"].dropna().unique()
            if len(sv):
                strand = str(sv[0])

        # Coding positions
        coding_positions = []
        utr_positions    = []
        if "is_coding" in df.columns:
            coding_positions = sorted(
                df.loc[df["is_coding"].fillna(False).astype(bool), "POS"].tolist()
            )
        if "consequence" in df.columns:
            utr_mask = df["consequence"].isin(["5_prime_UTR_variant", "3_prime_UTR_variant"])
            utr_positions = sorted(df.loc[utr_mask, "POS"].tolist())

        def _cluster(positions: List[int], gap: int) -> List[Tuple[int, int]]:
            """Cluster positions within `gap` bp into blocks."""
            blocks = []
            if not positions:
                return blocks
            s, e = positions[0], positions[0]
            for p in positions[1:]:
                if p - e <= gap:
                    e = p
                else:
                    blocks.append((s, e))
                    s = e = p
            blocks.append((s, e))
            return blocks

        exon_blocks = _cluster(coding_positions, gap=500)
        utr_blocks  = _cluster(utr_positions, gap=300)

        # Minimum visible exon width: 0.5% of gene span
        min_width = max(50, (x_max - x_min) // 200)

        shapes = []
        # ── Intron backbone ────────────────────────────────────────────────────
        shapes.append({
            "type": "line",
            "x0": x_min, "x1": x_max,
            "y0": 0.5, "y1": 0.5,
            "line": {"color": "#455a64", "width": 2.5},
        })

        # Direction tick marks (chevrons) along the intron line
        n_ticks = min(12, max(4, (x_max - x_min) // 10_000))
        tick_step = (x_max - x_min) // (n_ticks + 1)
        tick_dx = min_width // 2
        for i in range(1, n_ticks + 1):
            tx = x_min + i * tick_step
            # Short diagonal tick to indicate direction
            tick_x0, tick_x1 = (tx - tick_dx, tx + tick_dx) if strand == "+" else (tx + tick_dx, tx - tick_dx)
            shapes.append({
                "type": "line",
                "x0": tick_x0, "x1": tx,
                "y0": 0.35, "y1": 0.5,
                "line": {"color": "#78909c", "width": 1.2},
            })
            shapes.append({
                "type": "line",
                "x0": tx, "x1": tick_x1,
                "y0": 0.5, "y1": 0.65,
                "line": {"color": "#78909c", "width": 1.2},
            })

        # ── UTR blocks (lighter shade) ─────────────────────────────────────────
        for (us, ue) in utr_blocks:
            w = max(ue - us, min_width)
            shapes.append({
                "type": "rect",
                "x0": us, "x1": us + w,
                "y0": 0.3, "y1": 0.7,
                "fillcolor": "#90caf9",
                "line": {"color": "#64b5f6", "width": 1},
                "opacity": 0.75,
            })

        # ── CDS exon blocks ────────────────────────────────────────────────────
        for (es, ee) in exon_blocks:
            w = max(ee - es, min_width)
            shapes.append({
                "type": "rect",
                "x0": es, "x1": es + w,
                "y0": 0.2, "y1": 0.8,
                "fillcolor": "#1565c0",
                "line": {"color": "#0d47a1", "width": 1},
                "opacity": 0.88,
            })

        # ── Annotations ────────────────────────────────────────────────────────
        strand_sym = "→" if strand == "+" else ("←" if strand == "-" else "")
        span_kb = round((x_max - x_min) / 1000, 1)
        mid_pos = (x_min + x_max) // 2

        annotations_plotly = [{
            "x":         mid_pos,
            "y":         0.95,
            "text":      (
                f"<b>{gene_name}</b> {strand_sym}&nbsp;&nbsp;"
                f"chr{chrom}:{x_min:,}–{x_max:,}&nbsp;"
                f"({span_kb} kb)&nbsp;&nbsp;"
                f"{len(exon_blocks)} exon block(s)"
            ),
            "showarrow": False,
            "font":      {"size": 11, "color": "#1a237e", "family": "monospace"},
            "xanchor":   "center",
            "yanchor":   "top",
        }]

        # Per-exon mini-labels (only if few blocks and gene is wide enough)
        if 1 <= len(exon_blocks) <= 12:
            for idx, (es, ee) in enumerate(exon_blocks, 1):
                mid = (es + ee) // 2
                annotations_plotly.append({
                    "x":         mid,
                    "y":         0.12,
                    "text":      f"E{idx}",
                    "showarrow": False,
                    "font":      {"size": 9, "color": "#546e7a"},
                    "xanchor":   "center",
                })

        return {
            "shapes":        shapes,
            "annotations":   annotations_plotly,
            "exon_blocks":   [{"start": s, "end": e} for s, e in exon_blocks],
            "n_exon_blocks": len(exon_blocks),
            "gene_name":     gene_name,
            "strand":        strand,
        }

    def _build_functional_summary_track(
        self, df, x_min: int, x_max: int, n_bins: int = 25
    ) -> Dict[str, Any]:
        """
        Track 3 — Functional summary (vectorized).

        Uses pd.cut() to bin variants along the genomic axis.
        Produces:
          - Stacked bar: consequence category counts per bin
          - Dotted line: mean CADD per bin (right y-axis)
        """
        import pandas as pd
        import numpy as np

        if df.empty:
            return {"traces": [], "n_bins": 0}

        pos_min = int(df["POS"].min())
        pos_max = int(df["POS"].max())
        if pos_max <= pos_min:
            pos_max = pos_min + 1

        # ── Vectorized binning with pd.cut ─────────────────────────────────────
        bins = pd.cut(df["POS"], bins=n_bins, labels=False, include_lowest=True)
        bin_edges = pd.cut(df["POS"], bins=n_bins, include_lowest=True, retbins=True)[1]
        bin_centers = ((bin_edges[:-1] + bin_edges[1:]) / 2).astype(int).tolist()

        df = df.copy()
        df["_bin"]  = bins
        df["_cat"]  = df.get("consequence", pd.Series(["intergenic_variant"]*len(df))).fillna("intergenic_variant").map(
            lambda c: _CSQ_TO_CAT.get(str(c), "non_coding")
        )

        # Consequence counts per bin per category
        cat_order = ["pLoF", "missense", "splice_region", "synonymous", "non_coding"]
        pivot = (
            df.groupby(["_bin", "_cat"], observed=False)
            .size()
            .unstack(fill_value=0)
        )
        # Reindex to ensure all bins and categories present
        pivot = pivot.reindex(range(n_bins), fill_value=0)
        for cat in cat_order:
            if cat not in pivot.columns:
                pivot[cat] = 0

        traces = []
        for cat in cat_order:
            counts = pivot[cat].tolist()
            if sum(counts) == 0:
                continue
            traces.append({
                "name":          cat,
                "type":          "bar",
                "x":             bin_centers,
                "y":             counts,
                "hovertemplate": f"<b>{cat}</b><br>Pos: %{{x:,.0f}}<br>Count: %{{y}}<extra></extra>",
                "marker":        {"color": self._CAT_COLORS.get(cat, "#9e9e9e"), "opacity": 0.82},
            })

        # CADD mean per bin (secondary y-axis, dotted line)
        if "cadd_phred" in df.columns:
            # Use .astype(float) before .tolist() so pandas.NA → np.nan (a real float)
            # then _safe_float handles np.nan / None / NAType uniformly.
            cadd_by_bin = (
                df.groupby("_bin", observed=False)["cadd_phred"]
                .mean()
                .reindex(range(n_bins))
                .astype(float)   # converts pd.NA → np.nan; NAType → float NaN
                .tolist()
            )
            cadd_vals = [_safe_float(v) for v in cadd_by_bin]
            if any(v is not None for v in cadd_vals):
                traces.append({
                    "name":          "Mean CADD",
                    "type":          "scatter",
                    "mode":          "lines+markers",
                    "x":             bin_centers,
                    "y":             [v if v is not None else 0 for v in cadd_vals],
                    "yaxis":         "y2",
                    "line":          {"color": "#e91e63", "width": 2, "dash": "dot"},
                    "marker":        {"size": 5, "color": "#e91e63"},
                    "hovertemplate": "CADD mean: %{y:.1f}<br>Pos: %{x:,.0f}<extra></extra>",
                })

        return {
            "traces":      traces,
            "bin_centers": bin_centers,
            "n_bins":      n_bins,
            "barmode":     "stack",
        }

    def _build_af_track(self, df) -> Dict[str, Any]:
        """
        Track 4 — Population allele frequency scatter.

        Y-axis: −log₁₀(AF) so rare variants appear higher on the axis.
        Color by rarity class. Skips variants with no AF (novel).

        AF source priority (unified af_effective hierarchy):
          1. af_mcps   (MCPS cohort — best for Mexican/Latin ancestry)
          2. af_gnomad (gnomAD 4.1.1 — highest quality public AF)
          3. af_global (dbSNP GRAF-pop — from annotator pipeline)
        Rarity derived from rarity_class_effective (unified) when available.
        """
        import pandas as pd
        import numpy as np

        # Ensure effective AF columns exist
        if "af_effective" not in df.columns:
            from services.variant_annotation import _compute_effective_af
            df = _compute_effective_af(df)

        if "af_effective" not in df.columns or df["af_effective"].isna().all():
            return {"traces": [], "has_af": False, "af_source": None}

        af_df = df.copy()
        af_df["_af"] = pd.to_numeric(af_df["af_effective"], errors="coerce")
        af_df = af_df[af_df["_af"].notna() & (af_df["_af"] > 0)].copy()

        # Build a per-row source label for the hover (e.g. "MCPS", "gnomAD", "dbSNP")
        # af_source column was set by _compute_effective_af
        if "af_source" in af_df.columns:
            af_source_counts = af_df["af_source"].value_counts().to_dict()
            # Primary label = most frequent source
            af_source_label = max(af_source_counts, key=af_source_counts.get) if af_source_counts else "mixed"
        else:
            af_source_label = "unknown"

        if af_df.empty:
            return {"traces": [], "has_af": False, "af_source": af_source_label}

        af_df["_log_af"] = -np.log10(af_df["_af"].clip(lower=1e-10))
        _rarity_col = "rarity_class_effective" if "rarity_class_effective" in af_df.columns else "rarity_class"
        af_df["_rarity"] = (
            af_df[_rarity_col].fillna("novel").astype(str)
            if _rarity_col in af_df.columns
            else pd.Series(["novel"] * len(af_df), index=af_df.index)
        )
        # Normalise rarity strings (ultra-rare / ultra_rare both accepted)
        af_df["_rarity"] = af_df["_rarity"].str.replace("-", "_")

        rarity_order = ["ultra_rare", "rare", "low_frequency", "common", "novel"]
        rarity_colors = {
            "ultra_rare":   "#f44336",
            "rare":         "#ff9800",
            "low_frequency":"#2196f3",
            "common":       "#4caf50",
            "novel":        "#9e9e9e",
        }
        rarity_labels = {
            "ultra_rare":   "Ultra-rare (<0.01%)",
            "rare":         "Rare (0.01–0.1%)",
            "low_frequency":"Low-freq (0.1–1%)",
            "common":       "Common (≥1%)",
            "novel":        "Novel / Ultra-rare",
        }

        # Per-population AF columns to show in hover (gnomAD)
        _pop_cols = {
            "af_gnomad_amr": "AMR",
            "af_gnomad_afr": "AFR",
            "af_gnomad_eas": "EAS",
            "af_gnomad_nfe": "NFE",
            "af_mcps_mex":   "MCPS-MEX",
            "af_mcps":       "MCPS",
        }

        def _hover_af(row) -> str:
            parts = [f"<b>{rarity_labels.get(row['_rarity'], row['_rarity'])}</b>"]
            row_src = str(row.get("af_source") or af_source_label)
            parts.append(f"AF ({row_src}) = {row['_af']:.2e}")
            parts.append(f"chr{row.get('CHROM', '?')}:{int(row['POS']):,}")
            rsid = str(row.get("rsid", "") or "")
            if rsid and rsid not in ("nan", "None", ""):
                parts.append(f"rsID: {rsid}")
            csq = str(row.get("consequence", "") or "")
            if csq and csq != "nan":
                parts.append(f"Csq: {csq}")
            # Population breakdown
            pop_parts = []
            for col, lbl in _pop_cols.items():
                v = _safe_float(row.get(col))
                if v is not None and v > 0:
                    pop_parts.append(f"{lbl}={v:.2e}")
            if pop_parts:
                parts.append("  ".join(pop_parts[:4]))
            return "<br>".join(parts)

        af_df["_hover"] = af_df.apply(_hover_af, axis=1)

        traces = []
        for rarity in rarity_order:
            sub = af_df[af_df["_rarity"] == rarity]
            if sub.empty:
                continue
            traces.append({
                "name":          rarity_labels.get(rarity, rarity),
                "type":          "scatter",
                "mode":          "markers",
                "x":             sub["POS"].tolist(),
                "y":             [_safe_float(v) or 0.0 for v in sub["_log_af"].tolist()],
                "text":          sub["_hover"].tolist(),
                "hovertemplate": "%{text}<extra></extra>",
                "marker":        {
                    "color":   rarity_colors.get(rarity, "#9e9e9e"),
                    "size":    5,
                    "opacity": 0.75,
                    "line":    {"width": 0.5, "color": "rgba(255,255,255,0.6)"},
                },
            })

        return {
            "traces":    traces,
            "has_af":    True,
            "n_with_af": len(af_df),
            "af_source": af_source_label,
            "af_col":    "af_effective",
        }

    # ── Data loading ─────────────────────────────────────────────────────────────

    @staticmethod
    def _read_file(parquet_path: Path, tsv_path: Path,
                   usecols: Optional[List[str]] = None):
        """
        Load a parquet or TSV annotation file into a pandas DataFrame.

        Priority:
          1. duckdb  — in requirements.txt; reads parquet natively
          2. pyarrow — may be present in some container images
          3. pandas  — TSV.GZ fallback; always available (in requirements via pgscat)

        The TSV fallback is tried whenever the parquet cannot be read,
        regardless of whether the parquet file exists.
        """
        import pandas as pd

        # ── 1. Try duckdb ──────────────────────────────────────────────────────
        if parquet_path.exists():
            try:
                import duckdb
                con = duckdb.connect()
                try:
                    if usecols:
                        # Detect available columns via LIMIT 0 — works across all
                        # DuckDB versions (avoids parquet_schema() API differences).
                        avail = set(
                            con.execute(
                                "SELECT * FROM read_parquet(?) LIMIT 0",
                                [str(parquet_path)]
                            ).df().columns
                        )
                        cols    = [c for c in usecols if c in avail] or list(avail)
                        col_sql = ", ".join(f'"{c}"' for c in cols)
                    else:
                        col_sql = "*"
                    df = con.execute(
                        f"SELECT {col_sql} FROM read_parquet(?)", [str(parquet_path)]
                    ).df()
                finally:
                    con.close()
                # DuckDB may return nullable extension types (Int32, Float64, …).
                # Cast numeric ones to float64 so NaN handling is uniform
                # across duckdb / pyarrow / TSV paths.
                for _c in list(df.columns):
                    if (pd.api.types.is_extension_array_dtype(df[_c].dtype)
                            and pd.api.types.is_numeric_dtype(df[_c].dtype)):
                        df[_c] = df[_c].astype("float64")
                logger.debug("duckdb read %s (%d rows)", parquet_path.name, len(df))
                return df
            except Exception as exc:
                logger.debug("duckdb failed for %s: %s", parquet_path.name, exc)

        # ── 2. Try pyarrow ─────────────────────────────────────────────────────
        if parquet_path.exists():
            try:
                import pyarrow.parquet as pq
                tbl = pq.read_table(str(parquet_path), columns=usecols)
                df = tbl.to_pandas()
                logger.debug("pyarrow read %s (%d rows)", parquet_path.name, len(df))
                return df
            except Exception as exc:
                logger.debug("pyarrow failed for %s: %s", parquet_path.name, exc)

        # ── 3. TSV.GZ fallback ─────────────────────────────────────────────────
        if tsv_path.exists():
            try:
                df = pd.read_csv(
                    tsv_path, sep="\t",
                    usecols=usecols,
                    dtype={"CHROM": str, "POS": "int64"},
                    compression="infer",
                    low_memory=False,
                )
                logger.debug("TSV read %s (%d rows)", tsv_path.name, len(df))
                return df
            except Exception as exc:
                logger.debug("TSV failed for %s: %s", tsv_path.name, exc)

        return None

    @staticmethod
    @staticmethod
    def _is_valid_gene_symbol(token: str) -> bool:
        """Return True if token looks like a gene symbol (not an Ensembl ID or empty)."""
        t = token.strip()
        if not t or t in (".", "NA", "nan", "None"):
            return False
        if t.startswith(("ENSG", "ENST", "ENSP", "ENSE")):
            return False
        return True

    def _load_gene_variants(self, gene_symbol: str):
        """
        Load annotated variants for gene_symbol from all available files.

        For each PGS annotation dir, tries:
          1. duckdb with push-down WHERE filter (fastest for parquets)
          2. pyarrow + pandas filter
          3. pandas TSV.GZ + pandas filter (vectorized — works for large files)

        Filters rows where:
          - gene_name == gene_symbol  (case-insensitive exact match only)

        Note: all_overlapping_genes is intentionally NOT used here so that all
        metrics (header counts, PRS panel, Gene×PGS) share the same definition
        of "variant in this gene".  Overlapping-gene rows inflate counts without
        being directly assigned to the gene and are excluded from primary metrics.
        """
        import pandas as pd

        if not self.annotations_dir.exists():
            logger.warning("Annotations dir not found: %s", self.annotations_dir)
            return None

        frames = []
        gene_upper = gene_symbol.upper()

        for subdir in sorted(self.annotations_dir.iterdir()):
            if not subdir.is_dir():
                continue
            parquet_path = subdir / f"{subdir.name}_variants_annotated.parquet"
            tsv_path     = subdir / f"{subdir.name}_variants_annotated.tsv.gz"

            df_chunk = self._read_filtered(parquet_path, tsv_path, gene_upper)
            if df_chunk is None or df_chunk.empty:
                continue

            logger.debug("%s: %d variants matched for gene %s",
                         subdir.name, len(df_chunk), gene_symbol)
            frames.append(df_chunk)

        if not frames:
            logger.debug("No variants found for gene %s", gene_symbol)
            return None

        try:
            return pd.concat(frames, ignore_index=True)
        except Exception as exc:
            logger.error("Cannot concat gene variant frames: %s", exc)
            return None

    def _enrich_gene_variants(self, df):
        """
        Live-enrich gene variant DataFrame with gnomAD / ClinVar / MCPS.

        Uses batch/region fetch for all sources:
          - gnomAD: fetch_batch() → fetch_region() — ONE tabix seek for the
            entire gene region instead of N individual seeks per variant.
          - MCPS:   fetch_batch() — single DuckDB VALUES JOIN.
          - ClinVar: fetch_batch() — single DuckDB VALUES JOIN.
        Only enriches columns not already present in df.
        """
        import pandas as pd

        required = ["CHROM", "POS", "EFFECT_ALLELE", "OTHER_ALLELE"]
        if not all(c in df.columns for c in required):
            return df
        if df.empty:
            return df

        needs_gnomad  = "af_gnomad" not in df.columns
        needs_mcps    = "af_mcps" not in df.columns
        needs_clinvar = "clinvar_clnsig" not in df.columns

        # Vectorized variant list construction — no iterrows()
        flip_mask = df.get("IS_FLIP", pd.Series(False, index=df.index)).astype(bool)
        ea = df["EFFECT_ALLELE"].fillna("").astype(str)
        oa = df["OTHER_ALLELE"].fillna("").astype(str)
        ref_ser = oa.where(flip_mask, ea)
        alt_ser = ea.where(flip_mask, oa)
        variants = list(zip(
            df["CHROM"].astype(str),
            df["POS"].astype(int),
            ref_ser,
            alt_ser,
        ))

        if needs_gnomad:
            try:
                from services.gnomad_fetcher import get_gnomad_fetcher
                gf = get_gnomad_fetcher()
                if gf.is_available():
                    gnomad_cols = [
                        "af_gnomad", "af_gnomad_afr", "af_gnomad_amr",
                        "af_gnomad_eas", "af_gnomad_nfe", "rarity_class_gnomad",
                    ]
                    # fetch_batch groups by chrom and uses fetch_region for
                    # compact regions → single tabix seek for the whole gene
                    batch = gf.fetch_batch(variants)
                    col_data: dict = {c: [None] * len(df) for c in gnomad_cols}
                    for i, v in enumerate(variants):
                        res = batch.get(v)
                        if res:
                            for c in gnomad_cols:
                                col_data[c][i] = res.get(c)
                    df = df.copy()
                    for c in gnomad_cols:
                        df[c] = col_data[c]
            except Exception as exc:
                logger.warning("GeneBrowser gnomAD enrichment failed: %s", exc)

        if needs_mcps:
            try:
                from services.mcps_fetcher import get_mcps_fetcher
                mf = get_mcps_fetcher()
                if mf.is_available() and mf.index_ready():
                    mcps_cols = ["af_mcps", "af_mcps_afr", "af_mcps_eur", "af_mcps_mex"]
                    # fetch_batch uses a single DuckDB VALUES JOIN
                    batch = mf.fetch_batch(variants)
                    col_data = {c: [None] * len(df) for c in mcps_cols}
                    for i, v in enumerate(variants):
                        res = batch.get(v)
                        if res:
                            for c in mcps_cols:
                                col_data[c][i] = res.get(c)
                    df = df.copy()
                    for c in mcps_cols:
                        df[c] = col_data[c]
            except Exception as exc:
                logger.warning("GeneBrowser MCPS enrichment failed: %s", exc)

        if needs_clinvar:
            try:
                from services.clinvar_fetcher import get_clinvar_fetcher
                cf = get_clinvar_fetcher()
                if cf.is_available():
                    clinvar_map = {
                        "clnsig":       "clinvar_clnsig",
                        "clndn":        "clinvar_clndn",
                        "clnrevstat":   "clinvar_clnrevstat",
                        "alleleid":     "clinvar_alleleid",
                        "clnvc":        "clinvar_clnvc",
                        "clinvar_rsid": "clinvar_rsid",
                    }
                    # fetch_batch uses a single DuckDB VALUES JOIN
                    batch = cf.fetch_batch(variants)
                    col_data = {dest: [None] * len(df) for dest in clinvar_map.values()}
                    for i, v in enumerate(variants):
                        res = batch.get(v)
                        if res:
                            for src, dest in clinvar_map.items():
                                col_data[dest][i] = res.get(src)
                    df = df.copy()
                    for dest in clinvar_map.values():
                        df[dest] = col_data[dest]
            except Exception as exc:
                logger.warning("GeneBrowser ClinVar enrichment failed: %s", exc)

        return df

    def _read_filtered(self, parquet_path: Path, tsv_path: Path,
                       gene_upper: str):
        """
        Read a single annotation file and return only rows matching gene_upper.
        Uses duckdb push-down filter when available; falls back to pandas.
        Filter: gene_name == gene_upper (exact, case-insensitive).
        """
        import pandas as pd

        # ── 1. duckdb with SQL WHERE push-down ───────────────────────────────
        if parquet_path.exists():
            try:
                import duckdb
                con = duckdb.connect()
                # Escape single quotes in gene symbol for SQL safety
                gene_sql = gene_upper.replace("'", "''")
                df = con.execute(f"""
                    SELECT * FROM read_parquet(?)
                    WHERE upper(gene_name) = '{gene_sql}'
                """, [str(parquet_path)]).df()
                con.close()
                # Normalise nullable extension types → float64
                for _c in list(df.columns):
                    if (pd.api.types.is_extension_array_dtype(df[_c].dtype)
                            and pd.api.types.is_numeric_dtype(df[_c].dtype)):
                        df[_c] = pd.to_numeric(df[_c], errors="coerce")
                logger.debug("duckdb filtered %s → %d rows", parquet_path.name, len(df))
                return df
            except Exception as exc:
                logger.debug("duckdb filter failed for %s: %s", parquet_path.name, exc)

        # ── 2. pyarrow full read + pandas filter ─────────────────────────────
        if parquet_path.exists():
            try:
                import pyarrow.parquet as pq
                df = pq.read_table(str(parquet_path)).to_pandas()
                return self._pandas_filter(df, gene_upper)
            except Exception as exc:
                logger.debug("pyarrow failed for %s: %s", parquet_path.name, exc)

        # ── 3. TSV.GZ fallback: chunked read + pandas filter ─────────────────
        if tsv_path.exists():
            try:
                # Use chunksize for large files to avoid OOM; each chunk is filtered
                chunk_frames = []
                for chunk in pd.read_csv(
                    tsv_path, sep="\t",
                    dtype={"CHROM": str, "POS": "int64"},
                    compression="infer",
                    low_memory=False,
                    chunksize=50_000,
                ):
                    filtered = self._pandas_filter(chunk, gene_upper)
                    if not filtered.empty:
                        chunk_frames.append(filtered)
                if chunk_frames:
                    return pd.concat(chunk_frames, ignore_index=True)
                return pd.DataFrame()
            except Exception as exc:
                logger.debug("TSV read failed for %s: %s", tsv_path.name, exc)

        return None

    @staticmethod
    def _pandas_filter(df, gene_upper: str):
        """Vectorized pandas filter: gene_name exact match only.

        all_overlapping_genes is intentionally excluded — see _load_gene_variants.
        """
        import pandas as pd
        if "gene_name" not in df.columns:
            return df.iloc[0:0]  # empty with same columns
        mask = df["gene_name"].fillna("").str.upper() == gene_upper
        return df[mask]

    def list_available_pgs_ids(self) -> List[str]:
        """
        Return sorted list of PGS IDs that have annotation parquet/TSV files.
        Used to populate the PGS selector dropdown in the genes UI.
        """
        if not self.annotations_dir.exists():
            return []
        ids = []
        for subdir in sorted(self.annotations_dir.iterdir()):
            if not subdir.is_dir():
                continue
            parquet_path = subdir / f"{subdir.name}_variants_annotated.parquet"
            tsv_path     = subdir / f"{subdir.name}_variants_annotated.tsv.gz"
            if parquet_path.exists() or tsv_path.exists():
                ids.append(subdir.name)
        return ids

    def list_annotated_genes(
        self,
        min_variants: int = 1,
        pgs_id: Optional[str] = None,
    ) -> List[Dict[str, Any]]:
        """
        Return gene symbols found in annotation files with variant counts.

        Parameters
        ----------
        pgs_id : str | None
            When provided, restrict results to that single PGS annotation file.
            When None, scan all annotation directories (global view).

        - Reads only gene_name + all_overlapping_genes columns (fast even for 1M-row files)
        - Vectorized: uses pandas value_counts() and str.split().explode() — no iterrows()
        - Global results (pgs_id=None) are cached in memory; per-PGS results are NOT cached
          (per-PGS queries are fast single-file reads, caching adds complexity for little gain)
        - Includes genes from both gene_name and all_overlapping_genes (ENSG filtered out)
        """
        import pandas as pd

        # ── Per-PGS fast path: single file, no cache ──────────────────────────
        if pgs_id:
            return self._list_genes_for_pgs(pgs_id, min_variants)

        # ── Global path: check in-memory cache ────────────────────────────────
        # Invalidate if any annotation subdir was modified since last scan
        try:
            newest_mtime = max(
                p.stat().st_mtime
                for p in self.annotations_dir.glob("*/")
                if p.is_dir()
            ) if self.annotations_dir.exists() else 0.0
        except (OSError, ValueError):
            newest_mtime = 0.0

        if (hasattr(self, "_gene_cache")
                and self._gene_cache_mtime == newest_mtime
                and self._gene_cache is not None):
            return self._gene_cache

        # ── Build global gene index ───────────────────────────────────────────
        gene_counts: Counter = Counter()

        if not self.annotations_dir.exists():
            logger.warning("Annotations dir not found: %s", self.annotations_dir)
            return []

        gene_cols = ["gene_name", "all_overlapping_genes"]

        for subdir in sorted(self.annotations_dir.iterdir()):
            if not subdir.is_dir():
                continue
            parquet_path = subdir / f"{subdir.name}_variants_annotated.parquet"
            tsv_path     = subdir / f"{subdir.name}_variants_annotated.tsv.gz"

            df = self._read_file(parquet_path, tsv_path, usecols=gene_cols)
            if df is None or df.empty:
                continue

            self._accumulate_gene_counts(df, gene_counts)
            logger.debug("%s: gene index now has %d symbols", subdir.name, len(gene_counts))

        # ── Annotate with clinical gene info ──────────────────────────────────
        result = self._build_gene_list(gene_counts, min_variants)
        logger.info("list_annotated_genes (global): %d genes indexed", len(result))

        self._gene_cache: List[Dict[str, Any]] = result
        self._gene_cache_mtime: float = newest_mtime
        return result

    def _list_genes_for_pgs(
        self, pgs_id: str, min_variants: int = 1
    ) -> List[Dict[str, Any]]:
        """Build gene list restricted to a single PGS annotation file."""
        subdir       = self.annotations_dir / pgs_id
        parquet_path = subdir / f"{pgs_id}_variants_annotated.parquet"
        tsv_path     = subdir / f"{pgs_id}_variants_annotated.tsv.gz"

        if not subdir.is_dir():
            logger.warning("No annotation dir for %s", pgs_id)
            return []

        gene_cols = ["gene_name", "all_overlapping_genes"]
        df = self._read_file(parquet_path, tsv_path, usecols=gene_cols)
        if df is None or df.empty:
            return []

        gene_counts: Counter = Counter()
        self._accumulate_gene_counts(df, gene_counts)
        result = self._build_gene_list(gene_counts, min_variants)
        logger.info("list_annotated_genes (%s): %d genes", pgs_id, len(result))
        return result

    @staticmethod
    def _accumulate_gene_counts(df, gene_counts: Counter) -> None:
        """Vectorised accumulation of gene_name + all_overlapping_genes into counter."""
        if "gene_name" in df.columns:
            vc = (
                df["gene_name"].dropna().astype(str).str.strip().value_counts()
            )
            for gene, cnt in vc.items():
                if GeneBrowserService._is_valid_gene_symbol(gene):
                    gene_counts[gene] += int(cnt)

        if "all_overlapping_genes" in df.columns:
            exploded = (
                df["all_overlapping_genes"]
                .dropna().astype(str).str.split(",").explode().str.strip()
            )
            for gene, cnt in exploded.value_counts().items():
                if GeneBrowserService._is_valid_gene_symbol(gene):
                    gene_counts[gene] = max(gene_counts[gene], int(cnt))

    @staticmethod
    def _build_gene_list(gene_counts: Counter, min_variants: int) -> List[Dict[str, Any]]:
        """Convert counter to sorted list of gene dicts with clinical annotation."""
        cg_svc = get_clinical_genes_service()
        return [
            {
                "gene_symbol":      gene,
                "n_variants":       count,
                "is_clinical_gene": cg_svc.is_clinical_gene(gene),
            }
            for gene, count in gene_counts.most_common()
            if count >= min_variants
        ]

    # ── Gene summary (heatmap / ranking table) ────────────────────────────────

    # ── Internal: load & aggregate per-gene metrics ───────────────────────────

    def _load_combined_metrics(
        self,
        pgs_id: Optional[str] = None,
    ):
        """
        Read metric columns from one or all annotation parquets and return a
        single concatenated DataFrame. Returns None when no data is found.

        Global results (pgs_id=None) are cached in memory with mtime-based
        invalidation — the same pattern used by list_annotated_genes.
        Per-PGS results are NOT cached (single-file read, fast enough).
        """
        import pandas as pd

        metric_cols = ["gene_name", "BETA", "af_global", "af_gnomad", "af_mcps",
                       "af_effective", "af_source", "consequence", "cadd_phred",
                       "revel_score", "rarity_class", "rarity_class_gnomad",
                       "rarity_class_effective", "is_lof"]

        if not self.annotations_dir.exists():
            return None

        frames: list = []

        if pgs_id:
            # ── Per-PGS fast path: single file, no cache ──────────────────────
            subdir       = self.annotations_dir / pgs_id
            parquet_path = subdir / f"{pgs_id}_variants_annotated.parquet"
            tsv_path     = subdir / f"{pgs_id}_variants_annotated.tsv.gz"
            if not subdir.is_dir():
                return None
            df = self._read_file(parquet_path, tsv_path, usecols=metric_cols)
            if df is not None and not df.empty:
                frames.append(df)
            return pd.concat(frames, ignore_index=True) if frames else None

        # ── Global path: mtime-based in-memory cache ──────────────────────────
        try:
            newest_mtime = max(
                p.stat().st_mtime
                for p in self.annotations_dir.glob("*/")
                if p.is_dir()
            ) if self.annotations_dir.exists() else 0.0
        except (OSError, ValueError):
            newest_mtime = 0.0

        if (hasattr(self, "_metrics_cache")
                and self._metrics_cache_mtime == newest_mtime
                and self._metrics_cache is not None):
            return self._metrics_cache

        for subdir in sorted(self.annotations_dir.iterdir()):
            if not subdir.is_dir():
                continue
            parquet_path = subdir / f"{subdir.name}_variants_annotated.parquet"
            tsv_path     = subdir / f"{subdir.name}_variants_annotated.tsv.gz"
            df = self._read_file(parquet_path, tsv_path, usecols=metric_cols)
            if df is not None and not df.empty:
                frames.append(df)

        result = pd.concat(frames, ignore_index=True) if frames else None
        self._metrics_cache: object = result
        self._metrics_cache_mtime: float = newest_mtime
        return result

    @staticmethod
    def _filter_sort_gene_rows(
        gene_rows: list,
        sort_by: str = "sum_abs_beta",
        sort_dir: str = "desc",
        clinical_only: bool = False,
        clinical_confidence: Optional[str] = None,
        min_variants: int = 1,
        cg_svc=None,
        q: str = "",
    ) -> list:
        """
        Apply client-requested filters and sorting to a list of gene metric dicts.

        sort_by options:
            ranking_score  → ranking_score_mean
            variant_count  → n_variants
            sum_abs_beta   → sum_abs_beta
            mean_cadd      → mean_cadd  (None treated as 0)
            alphabetical   → gene_name
        """
        # Symbol substring filter — applied first so total_filtered reflects it
        if q:
            q_upper = q.strip().upper()
            gene_rows = [g for g in gene_rows if q_upper in g["gene_name"].upper()]

        if min_variants > 1:
            gene_rows = [g for g in gene_rows if g["n_variants"] >= min_variants]

        if clinical_only:
            gene_rows = [g for g in gene_rows if g.get("is_clinical_gene")]

        if clinical_confidence and cg_svc and cg_svc.is_available():
            conf_lower = clinical_confidence.lower()
            allowed = {
                sym for sym in cg_svc.get_all_symbols()
                if str(cg_svc.get_confidence(sym) or "").lower() == conf_lower
            }
            gene_rows = [g for g in gene_rows if g["gene_name"].upper() in allowed]

        # Sorting
        rev = (sort_dir.lower() != "asc")
        key_map = {
            "ranking_score": lambda g: g.get("ranking_score_mean") or 0.0,
            "variant_count":  lambda g: g.get("n_variants") or 0,
            "sum_abs_beta":   lambda g: g.get("sum_abs_beta") or 0.0,
            "mean_cadd":      lambda g: g.get("mean_cadd") or 0.0,
            "alphabetical":   lambda g: g.get("gene_name", ""),
        }
        sort_key = key_map.get(sort_by, key_map["sum_abs_beta"])
        # alphabetical ascending by default
        if sort_by == "alphabetical":
            rev = not rev
        gene_rows.sort(key=sort_key, reverse=rev)
        return gene_rows

    @staticmethod
    def _build_heatmap(gene_rows: list) -> Dict[str, Any]:
        """Build Plotly-ready heatmap matrix from a list of gene metric dicts."""
        import math

        def _safe(v) -> float:
            if v is None:
                return 0.0
            try:
                f = float(v)
                return 0.0 if (math.isnan(f) or math.isinf(f)) else f
            except (TypeError, ValueError):
                return 0.0

        def _norm(vals):
            mx = max(vals) if vals else 1.0
            return [round(v / mx, 4) if mx else 0.0 for v in vals]

        gene_names = [g["gene_name"] for g in gene_rows]
        z_matrix   = [
            _norm([_safe(g["n_variants"])        for g in gene_rows]),
            _norm([_safe(g["sum_abs_beta"])       for g in gene_rows]),
            _norm([_safe(g["mean_cadd"])          for g in gene_rows]),
            [_safe(g["ranking_score_mean"])        for g in gene_rows],
        ]
        return {
            "gene_names": gene_names,
            "metrics":    ["n_variants", "sum_abs_beta", "mean_cadd", "ranking_score"],
            "z":          z_matrix,
        }

    # ── Gene summary (heatmap / ranking table) ────────────────────────────────

    def get_gene_summary(
        self,
        limit: int = 200,
        pgs_id: Optional[str] = None,
        sort_by: str = "sum_abs_beta",
        sort_dir: str = "desc",
        clinical_only: bool = False,
        clinical_confidence: Optional[str] = None,
        min_variants: int = 1,
        page: int = 1,
        page_size: Optional[int] = None,
        q: str = "",
    ) -> Dict[str, Any]:
        """
        Aggregate per-gene metrics across annotation parquets with full
        filtering, sorting, and optional pagination.

        Parameters
        ----------
        pgs_id             : restrict to one PGS (None = global)
        sort_by            : ranking_score | variant_count | sum_abs_beta |
                             mean_cadd | alphabetical
        sort_dir           : asc | desc
        clinical_only      : only genes in clinical dataset
        clinical_confidence: high | medium | low (GenCC tier)
        min_variants       : minimum variant count
        page, page_size    : pagination (page_size=None → return all up to limit)

        Returns
        -------
        {
          pgs_id, total_genes, total_filtered,
          page, page_size, total_pages,
          genes[],        # current page
          heatmap{}       # top-N for heatmap (always limit genes, no pagination)
        }
        """
        from services.variant_ranking import aggregate_gene_metrics

        cg_svc  = get_clinical_genes_service()
        cg_syms = cg_svc.get_all_symbols() if cg_svc.is_available() else set()

        combined = self._load_combined_metrics(pgs_id)
        empty = {"pgs_id": pgs_id, "genes": [], "heatmap": {}, "total_genes": 0,
                 "total_filtered": 0, "page": page, "page_size": page_size,
                 "total_pages": 1}
        if combined is None:
            if pgs_id:
                empty["error"] = f"No annotation directory found for {pgs_id}."
            return empty

        # For global queries, cache the expensive aggregate_gene_metrics result
        # using the same mtime already stored by _load_combined_metrics.
        # Per-PGS queries are fast single-file reads so they are not cached.
        if not pgs_id and (
            hasattr(self, "_all_rows_cache")
            and self._all_rows_cache is not None
            and getattr(self, "_all_rows_cache_mtime", None) == getattr(self, "_metrics_cache_mtime", object())
        ):
            all_rows = self._all_rows_cache
        else:
            all_rows = aggregate_gene_metrics(combined, cg_syms)
            if not pgs_id:
                self._all_rows_cache = all_rows
                self._all_rows_cache_mtime = getattr(self, "_metrics_cache_mtime", 0.0)

        total_genes = len(all_rows)

        # Apply filters + sort (q applied here so total_filtered reflects it)
        filtered = self._filter_sort_gene_rows(
            all_rows,
            sort_by=sort_by,
            sort_dir=sort_dir,
            clinical_only=clinical_only,
            clinical_confidence=clinical_confidence,
            min_variants=min_variants,
            cg_svc=cg_svc,
            q=q,
        )
        total_filtered = len(filtered)

        # Heatmap always shows top `limit` of filtered results
        heatmap_rows = filtered[:limit]
        heatmap      = self._build_heatmap(heatmap_rows)

        # Pagination
        if page_size:
            ps          = max(1, min(page_size, 2000))
            total_pages = max(1, (total_filtered + ps - 1) // ps)
            pg          = max(1, min(page, total_pages))
            start       = (pg - 1) * ps
            page_rows   = filtered[start: start + ps]
        else:
            ps          = total_filtered
            total_pages = 1
            pg          = 1
            page_rows   = filtered[:limit]

        return {
            "pgs_id":         pgs_id,
            "total_genes":    total_genes,
            "total_filtered": total_filtered,
            "page":           pg,
            "page_size":      ps,
            "total_pages":    total_pages,
            "genes":          page_rows,
            "heatmap":        heatmap,
        }

    # ── Gene × PGS matrix ────────────────────────────────────────────────────

    def get_gene_pgs_matrix(
        self,
        gene_symbol: str,
        min_variants: int = 1,
        sort_by: str = "sum_abs_beta",
    ) -> Dict[str, Any]:
        """
        For a given gene, compute per-PGS metrics across all annotation parquets.

        Heatmap 2: Gene × PGS
        ---------------------
        Axis X = PGS IDs, Axis Y = metrics for this gene in that PGS.
        Answers: in how many PGS does this gene appear? In which does it
        contribute most to the score?

        Uses _load_gene_variants (the same data path as get_gene_info) for
        consistency and robustness, then groups by PRS_ID.

        Returns
        -------
        {
          "gene_symbol":  str,
          "total_pgs":    int,
          "pgs_rows": [
            {
              "pgs_id":             str,
              "n_variants":         int,
              "coding_count":       int,
              "splice_count":       int,
              "sum_abs_beta":       float,
              "mean_cadd":          float | None,
              "max_ranking_score":  float,
              "mean_ranking_score": float,
            }, ...
          ],
          "heatmap": {
            "pgs_ids":  [str],    # x-axis
            "metrics":  [str],    # y-axis
            "z":        [[float]],
          }
        }
        """
        import pandas as pd
        import math

        from services.variant_ranking import score_dataframe as _sdf

        _empty = {"gene_symbol": gene_symbol, "total_pgs": 0, "pgs_rows": [], "heatmap": {}}

        # Reuse the same data path as get_gene_info: all variants for this gene
        # across all PGS files, combined into a single DataFrame with PRS_ID column.
        df = self._load_gene_variants(gene_symbol)
        if df is None or df.empty or "PRS_ID" not in df.columns:
            return _empty

        cg_svc  = get_clinical_genes_service()
        cg_syms = cg_svc.get_all_symbols() if cg_svc.is_available() else set()

        # LoF + missense + synonymous consequence sets
        coding_csqs = set(
            CONSEQUENCE_CATEGORIES.get("pLoF", [])
            + CONSEQUENCE_CATEGORIES.get("missense", [])
            + CONSEQUENCE_CATEGORIES.get("synonymous", [])
        )
        splice_csqs = set(CONSEQUENCE_CATEGORIES.get("splice_region", []))

        def _sf(v) -> Optional[float]:
            if v is None:
                return None
            try:
                f = float(v)
                return None if (math.isnan(f) or math.isinf(f)) else f
            except (TypeError, ValueError):
                return None

        pgs_rows: list = []

        for pgs_id, gdf in df.groupby("PRS_ID"):
            n = len(gdf)
            if n < min_variants:
                continue

            # Genetic load  (sum |BETA|)
            load = 0.0
            if "BETA" in gdf.columns:
                raw = pd.to_numeric(gdf["BETA"], errors="coerce").abs().sum()
                load = float(raw) if not math.isnan(raw) else 0.0

            # Consequence counts
            coding_count = 0
            splice_count = 0
            if "consequence" in gdf.columns:
                csq_series   = gdf["consequence"].fillna("")
                coding_count = int(csq_series.isin(coding_csqs).sum())
                splice_count = int(csq_series.isin(splice_csqs).sum())

            # CADD
            mean_cadd = None
            if "cadd_phred" in gdf.columns:
                mean_cadd = _sf(
                    pd.to_numeric(gdf["cadd_phred"], errors="coerce").mean()
                )

            # Ranking scores
            scores     = _sdf(gdf, cg_syms)
            max_score  = float(scores.max())
            mean_score = float(scores.mean())

            pgs_rows.append({
                "pgs_id":             str(pgs_id),
                "n_variants":         n,
                "coding_count":       coding_count,
                "splice_count":       splice_count,
                "sum_abs_beta":       round(load, 6),
                "mean_cadd":          round(mean_cadd, 2) if mean_cadd is not None else None,
                "max_ranking_score":  round(max_score,  4),
                "mean_ranking_score": round(mean_score, 4),
            })

        if not pgs_rows:
            return _empty

        # Sort
        sort_key_map = {
            "sum_abs_beta":      lambda r: r["sum_abs_beta"],
            "max_ranking_score": lambda r: r["max_ranking_score"],
            "variant_count":     lambda r: r["n_variants"],
        }
        pgs_rows.sort(key=sort_key_map.get(sort_by, sort_key_map["sum_abs_beta"]), reverse=True)

        # Heatmap matrix: x = PGS IDs, y = metrics
        pgs_ids    = [r["pgs_id"] for r in pgs_rows]
        metric_key = ["n_variants", "coding_count", "splice_count",
                      "sum_abs_beta", "mean_cadd", "max_ranking_score", "mean_ranking_score"]

        def _norm_row(vals):
            clean = [0.0 if v is None else float(v) for v in vals]
            mx    = max(clean) if clean else 1.0
            return [round(v / mx, 4) if mx else 0.0 for v in clean]

        z_matrix = [_norm_row([r[k] for r in pgs_rows]) for k in metric_key]

        return {
            "gene_symbol": gene_symbol,
            "total_pgs":   len(pgs_rows),
            "pgs_rows":    pgs_rows,
            "heatmap": {
                "pgs_ids": pgs_ids,
                "metrics": metric_key,
                "z":       z_matrix,
            },
        }

    def _load_gene_variants_from_path(
        self,
        parquet_path: Path,
        tsv_path: Path,
        gene_upper: str,
        usecols: Optional[List[str]] = None,
    ):
        """
        Fast per-file gene filter: DuckDB WHERE clause when parquet exists,
        otherwise fall back to full read + pandas filter.
        Filter: gene_name == gene_upper only (all_overlapping_genes excluded).
        """
        import pandas as pd

        if parquet_path.exists():
            try:
                import duckdb
                con = duckdb.connect()
                try:
                    # Detect available columns via LIMIT 0 (DuckDB-version agnostic)
                    avail = set(
                        con.execute(
                            "SELECT * FROM read_parquet(?) LIMIT 0",
                            [str(parquet_path)]
                        ).df().columns
                    )
                    cols = ([c for c in (usecols or []) if c in avail] or list(avail))
                    if "gene_name" in avail and "gene_name" not in cols:
                        cols.append("gene_name")

                    col_sql   = ", ".join(f'"{c}"' for c in cols)
                    safe_gene = gene_upper.replace("'", "''")
                    df = con.execute(
                        f"SELECT {col_sql} FROM read_parquet(?) WHERE upper(coalesce(gene_name,'')) = '{safe_gene}'",
                        [str(parquet_path)]
                    ).df()
                finally:
                    con.close()
                return df if not df.empty else None
            except Exception as exc:
                logger.debug("duckdb gene filter failed %s: %s", parquet_path.name, exc)

        # Fallback: full read + pandas filter
        df_full = self._read_file(parquet_path, tsv_path, usecols=usecols)
        if df_full is None or df_full.empty:
            return None
        filtered = self._pandas_filter(df_full, gene_upper)
        return filtered if not filtered.empty else None

    # ── PRS × Clinical genes panel ────────────────────────────────────────────

    def get_prs_clinical_panel(self) -> Dict[str, Any]:
        """
        Intersect PRS variants with clinical genes and return genetic load metrics.

        Returns:
        {
          "n_clinical_genes_with_variants": int,
          "genes": [
            {
              "gene_name":          str,
              "n_variants":         int,
              "genetic_load":       float,   # sum |BETA|
              "cumulative_score":   float,   # sum ranking_scores
              "mean_ranking_score": float,
              "confidence":         str | None,
              "moi":                str | None,
              "disease":            str | None,
            }, ...
          ],
          "total_genetic_load": float,
        }
        """
        import pandas as pd
        import math

        from services.variant_ranking import score_dataframe as _sdf

        cg_svc  = get_clinical_genes_service()
        cg_syms = cg_svc.get_all_symbols() if cg_svc.is_available() else set()

        if not cg_syms:
            return {
                "n_clinical_genes_with_variants": 0,
                "genes": [],
                "total_genetic_load": 0.0,
                "error": "Clinical genes dataset not available.",
            }

        frames: list = []
        if not self.annotations_dir.exists():
            return {"n_clinical_genes_with_variants": 0, "genes": [], "total_genetic_load": 0.0}

        for subdir in sorted(self.annotations_dir.iterdir()):
            if not subdir.is_dir():
                continue
            parquet_path = subdir / f"{subdir.name}_variants_annotated.parquet"
            tsv_path     = subdir / f"{subdir.name}_variants_annotated.tsv.gz"
            df = self._read_file(
                parquet_path, tsv_path,
                usecols=["gene_name", "CHROM", "POS", "EFFECT_ALLELE",
                         "BETA", "af_global", "consequence",
                         "cadd_phred", "revel_score", "rarity_class", "is_lof"],
            )
            if df is not None and not df.empty:
                frames.append(df)

        if not frames:
            return {"n_clinical_genes_with_variants": 0, "genes": [], "total_genetic_load": 0.0}

        df = pd.concat(frames, ignore_index=True)

        # Keep only clinical gene rows
        if "gene_name" not in df.columns:
            return {"n_clinical_genes_with_variants": 0, "genes": [], "total_genetic_load": 0.0}

        df = df[df["gene_name"].str.upper().isin(cg_syms)].copy()
        if df.empty:
            return {"n_clinical_genes_with_variants": 0, "genes": [], "total_genetic_load": 0.0}

        df["_score"] = _sdf(df, cg_syms)

        gene_rows: list = []
        total_load = 0.0

        for gene, gdf in df.groupby("gene_name", sort=False):
            gene_upper = str(gene).upper()
            cg_entry   = cg_svc.get_gene_info(gene_upper) or {}

            load = float(gdf["BETA"].abs().sum()) if "BETA" in gdf.columns else 0.0
            if math.isnan(load):
                load = 0.0

            cum_score  = float(gdf["_score"].sum())
            mean_score = float(gdf["_score"].mean())
            total_load += load

            _dedup_cols = [c for c in ("CHROM", "POS", "EFFECT_ALLELE") if c in gdf.columns]
            n_unique = (
                len(gdf.drop_duplicates(subset=_dedup_cols))
                if _dedup_cols else len(gdf)
            )

            gene_rows.append({
                "gene_name":          str(gene),
                "n_variants":         n_unique,           # unique positions (CHROM,POS,EFFECT_ALLELE)
                "n_variant_records":  len(gdf),           # raw rows across all PGS files
                "genetic_load":       round(load, 6),
                "cumulative_score":   round(cum_score, 4),
                "mean_ranking_score": round(mean_score, 4),
                "confidence":         cg_entry.get("confidence"),
                "moi":                cg_entry.get("moi"),
                "disease":            cg_entry.get("disease"),
            })

        gene_rows.sort(key=lambda x: x["genetic_load"], reverse=True)

        return {
            "n_clinical_genes_with_variants": len(gene_rows),
            "genes":               gene_rows,
            "total_genetic_load":  round(total_load, 6),
        }
