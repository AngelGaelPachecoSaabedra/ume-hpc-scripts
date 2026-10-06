"""
Variant Classification Logic
================================
Assigns biological region class, consequence, and boolean flags to a variant
based on the list of GFF3Feature objects overlapping its position.

Priority order (most specific wins):
  1  CDS              → coding           / coding_sequence_variant
  2  five_prime_UTR   → UTR_5prime       / 5_prime_UTR_variant
  3  three_prime_UTR  → UTR_3prime       / 3_prime_UTR_variant
  4  exon (non-CDS)   → non_coding_exon  / non_coding_exon_variant
  5  transcript        → intronic         / intron_variant
  6  gene (no tx hit)  → intronic         / intron_variant  (non-basic transcripts)
  7  within UPSTREAM_WINDOW bp of a TSS → near_gene / upstream_gene_variant
  8  none              → intergenic       / intergenic_variant

Splice site classification:
  GFF3-only (coordinate-based, v1.1):
    distance ≤ SPLICE_CANONICAL_DIST  → splice_site_variant   (±1-2 bp)
    distance ≤ SPLICE_REGION_DIST     → splice_region_variant  (±3-8 bp)
  FASTA-backed (GT/AG inspection, v1.2):
    splice_type='donor'   → splice_donor_variant
    splice_type='acceptor' → splice_acceptor_variant

Coding consequence model (v1.2, requires FASTAEngine):
  CDS variants are further classified via codon context:
    missense_variant | synonymous_variant | stop_gained | stop_lost |
    start_lost | frameshift_variant | inframe_insertion | inframe_deletion

Multi-gene overlap model:
  All overlapping genes are enumerated; each gets a best-feature score.
  The highest-priority gene (protein_coding preferred at ties) becomes the
  "best hit" and fills the flat output columns for backward compatibility.
  All gene names, transcript IDs, and region classes are stored as
  comma-separated strings in all_overlapping_genes / all_overlapping_transcripts
  / all_region_classes.

Regulatory annotation:
  is_regulatory = True when:
    • The variant overlaps a feature in the RegulatoryBEDIndex (if loaded).
      regulatory_source → "bed:<source_name>"
    • OR the variant is within UPSTREAM_WINDOW bp of a gene boundary (heuristic).
      regulatory_source → "heuristic"
"""
from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, List, Optional, TYPE_CHECKING

if TYPE_CHECKING:
    from gff3_parser import GFF3Feature

# ── Constants ─────────────────────────────────────────────────────────────────

# bp upstream of a gene's start (strand-aware TSS) counted as "near_gene"
UPSTREAM_WINDOW: int = 2_000

# Splice-site distance thresholds (bp from nearest exon boundary)
SPLICE_CANONICAL_DIST: int = 2   # ±1-2 bp → canonical splice site
SPLICE_REGION_DIST:    int = 8   # ±3-8 bp → splice region

# Numeric priority: lower = more specific / higher biological impact
_FEATURE_PRIORITY: Dict[str, int] = {
    "CDS":            0,
    "five_prime_UTR": 1,
    "three_prime_UTR": 2,
    "exon":           3,
    "transcript":     4,
    "gene":           5,
}

# Consequence priority (lower = higher impact) — used for consequence_priority column
CONSEQUENCE_PRIORITY: Dict[str, int] = {
    # High-impact coding consequences (v1.2 FASTA-backed)
    "stop_gained":                 50,
    "frameshift_variant":          60,
    "stop_lost":                   70,
    "start_lost":                  80,
    # Splice site (canonical GT/AG, FASTA-backed in v1.2)
    "splice_donor_variant":       100,
    "splice_acceptor_variant":    110,
    # Splice site (coordinate-based fallback, v1.1)
    "splice_site_variant":        120,
    # Moderate-impact coding
    "inframe_insertion":          150,
    "inframe_deletion":           160,
    "splice_region_variant":      200,
    "missense_variant":           250,
    "synonymous_variant":         280,
    "coding_sequence_variant":    300,
    # Non-coding
    "5_prime_UTR_variant":        400,
    "3_prime_UTR_variant":        500,
    "non_coding_exon_variant":    600,
    "intron_variant":             700,
    "upstream_gene_variant":      800,
    "regulatory_region_variant":  850,
    "intergenic_variant":         900,
}

_REGION_CLASS: Dict[Optional[str], str] = {
    "CDS":            "coding",
    "five_prime_UTR": "UTR_5prime",
    "three_prime_UTR": "UTR_3prime",
    "exon":           "non_coding_exon",
    "transcript":     "intronic",
    "gene":           "intronic",
    None:             "intergenic",
}

_CONSEQUENCE: Dict[str, str] = {
    "coding":          "coding_sequence_variant",
    "UTR_5prime":      "5_prime_UTR_variant",
    "UTR_3prime":      "3_prime_UTR_variant",
    "non_coding_exon": "non_coding_exon_variant",
    "intronic":        "intron_variant",
    "splice_region":   "splice_region_variant",
    "near_gene":       "upstream_gene_variant",
    "regulatory":      "regulatory_region_variant",
    "intergenic":      "intergenic_variant",
}

# Feature types that represent exonic positions (exon or sub-exon features)
_EXON_FEATURES = frozenset({"CDS", "exon", "five_prime_UTR", "three_prime_UTR"})

# Feature types that represent intronic positions (inside gene/transcript but no exon)
_INTRONIC_FEATURES = frozenset({"transcript", "gene"})


# ── Data model ────────────────────────────────────────────────────────────────

@dataclass
class AnnotationResult:
    """
    Full per-variant annotation result.
    All fields must be JSON-serialisable.

    Best-hit fields (backward-compatible with v1.0 schema):
        gene_name, gene_id, gene_type, transcript_id, transcript_type,
        feature_type, region_class, consequence, is_coding, is_regulatory,
        is_intergenic, n_overlapping_genes, strand, distance_nearest_gene

    New fields (v1.1+):
        distance_to_splice_site   — bp to nearest exon boundary (all variants)
        consequence_priority      — numeric impact rank (lower = higher impact)
        all_overlapping_genes     — comma-separated gene names of all hits
        all_overlapping_transcripts — comma-separated transcript IDs
        all_region_classes        — comma-separated region_class per gene
        regulatory_source         — "heuristic", "bed:<name>", or ""
    """
    # ── Best-hit (v1.0-compatible) ────────────────────────────────────────────
    gene_name:           str   = ""
    gene_id:             str   = ""
    gene_type:           str   = ""
    transcript_id:       str   = ""
    transcript_type:     str   = ""
    feature_type:        str   = ""          # raw GFF3 feature driving the classification
    region_class:        str   = "intergenic"
    consequence:         str   = "intergenic_variant"
    is_coding:           bool  = False
    is_regulatory:       bool  = False
    is_intergenic:       bool  = True
    n_overlapping_genes: int   = 0
    strand:              str   = "."
    distance_nearest_gene: Optional[int] = None  # None for non-intergenic

    # ── New fields (v1.1+) ───────────────────────────────────────────────────
    distance_to_splice_site:    Optional[int] = None  # bp to nearest exon boundary
    consequence_priority:       int           = 900   # lower = higher impact
    all_overlapping_genes:      str           = ""    # "GENE1,GENE2,..."
    all_overlapping_transcripts: str          = ""    # "ENST1,ENST2,..."
    all_region_classes:         str           = ""    # "coding,intronic,..."
    regulatory_source:          str           = ""    # "heuristic" | "bed:name" | ""

    # ── New fields (v1.5, regulatory element detail) ─────────────────────────
    regulatory_element_type: str           = ""    # "promoter" | "enhancer" | "insulator" | "TFBS" | "open_chromatin" | …
    regulatory_element_id:   str           = ""    # TF name, OREG ID, or element name
    distance_to_regulatory:  Optional[int] = None  # 0 if overlapping, >0 bp to nearest element

    # ── New fields (v1.2, FASTA-backed) ──────────────────────────────────────
    # Codon/AA consequence
    codon_ref:          str           = ""    # reference codon, e.g. "ATG"
    codon_alt:          str           = ""    # alternate codon, e.g. "TTG"
    aa_ref:             str           = ""    # 1-letter reference amino acid
    aa_alt:             str           = ""    # 1-letter alternate amino acid
    aa_ref_3:           str           = ""    # 3-letter reference amino acid
    aa_alt_3:           str           = ""    # 3-letter alternate amino acid
    # Splice site direction (FASTA GT/AG inspection)
    splice_type:        str           = ""    # "donor" | "acceptor" | ""
    # dbNSFP5 functional scores
    cadd_phred:         Optional[float] = None
    revel_score:        Optional[float] = None
    sift_pred:          str           = ""    # "D" | "T" | ""
    polyphen2_pred:     str           = ""    # "D" | "P" | "B" | ""
    clinvar_clnsig:     str           = ""
    # Boolean consequence flags
    is_missense:        bool          = False
    is_synonymous:      bool          = False
    is_lof:             bool          = False  # stop_gained/lost/start_lost/frameshift/splice_donor/acceptor

    def to_dict(self) -> dict:
        return {
            # v1.0 fields
            "gene_name":            self.gene_name,
            "gene_id":              self.gene_id,
            "gene_type":            self.gene_type,
            "transcript_id":        self.transcript_id,
            "transcript_type":      self.transcript_type,
            "feature_type":         self.feature_type,
            "region_class":         self.region_class,
            "consequence":          self.consequence,
            "is_coding":            self.is_coding,
            "is_regulatory":        self.is_regulatory,
            "is_intergenic":        self.is_intergenic,
            "n_overlapping_genes":  self.n_overlapping_genes,
            "strand":               self.strand,
            "distance_nearest_gene": self.distance_nearest_gene,
            # v1.1 fields
            "distance_to_splice_site":     self.distance_to_splice_site,
            "consequence_priority":        self.consequence_priority,
            "all_overlapping_genes":       self.all_overlapping_genes,
            "all_overlapping_transcripts": self.all_overlapping_transcripts,
            "all_region_classes":          self.all_region_classes,
            "regulatory_source":           self.regulatory_source,
            # v1.5 fields
            "regulatory_element_type":     self.regulatory_element_type,
            "regulatory_element_id":       self.regulatory_element_id,
            "distance_to_regulatory":      self.distance_to_regulatory,
            # v1.2 fields
            "codon_ref":        self.codon_ref,
            "codon_alt":        self.codon_alt,
            "aa_ref":           self.aa_ref,
            "aa_alt":           self.aa_alt,
            "aa_ref_3":         self.aa_ref_3,
            "aa_alt_3":         self.aa_alt_3,
            "splice_type":      self.splice_type,
            "cadd_phred":       self.cadd_phred,
            "revel_score":      self.revel_score,
            "sift_pred":        self.sift_pred,
            "polyphen2_pred":   self.polyphen2_pred,
            "clinvar_clnsig":   self.clinvar_clnsig,
            "is_missense":      self.is_missense,
            "is_synonymous":    self.is_synonymous,
            "is_lof":           self.is_lof,
        }


# ── Classification ────────────────────────────────────────────────────────────

def _best_feature_for_gene(
    features: List["GFF3Feature"],
) -> "GFF3Feature":
    """
    From a list of features for a single gene, return the one with
    highest specificity.  Tie-break: protein_coding before other biotypes.
    """
    best = features[0]
    best_prio = _FEATURE_PRIORITY.get(best.feature, 999)
    for f in features[1:]:
        prio = _FEATURE_PRIORITY.get(f.feature, 999)
        if prio < best_prio:
            best, best_prio = f, prio
        elif prio == best_prio:
            if best.gene_type != "protein_coding" and f.gene_type == "protein_coding":
                best = f
    return best


def _apply_splice_classification(
    result: AnnotationResult,
    nearest_exon_dist: Optional[int],
    splice_type: Optional[str] = None,
) -> None:
    """
    Override consequence/region_class for intronic variants near exon boundaries.
    Mutates *result* in place.

    FASTA-backed (v1.2, splice_type provided):
      splice_type='donor'   → splice_donor_variant
      splice_type='acceptor' → splice_acceptor_variant

    Coordinate-based fallback (v1.1):
      Canonical splice site: distance ≤ SPLICE_CANONICAL_DIST (1-2 bp) → splice_site_variant
      Splice region:         distance ≤ SPLICE_REGION_DIST    (3-8 bp) → splice_region_variant

    Only reclassifies intronic variants (feature_type in _INTRONIC_FEATURES).
    For exonic/CDS variants the distance field is set but the consequence
    is not overridden (the exon consequence takes priority).
    """
    if nearest_exon_dist is None:
        return

    result.distance_to_splice_site = nearest_exon_dist

    # Only reclassify intronic variants
    if result.feature_type not in _INTRONIC_FEATURES:
        return

    if nearest_exon_dist <= SPLICE_CANONICAL_DIST:
        # FASTA-backed: refine to donor/acceptor if available
        if splice_type == "donor":
            result.splice_type          = "donor"
            result.consequence          = "splice_donor_variant"
            result.consequence_priority = CONSEQUENCE_PRIORITY["splice_donor_variant"]
            result.is_lof               = True
        elif splice_type == "acceptor":
            result.splice_type          = "acceptor"
            result.consequence          = "splice_acceptor_variant"
            result.consequence_priority = CONSEQUENCE_PRIORITY["splice_acceptor_variant"]
            result.is_lof               = True
        else:
            # Coordinate-based fallback
            result.consequence          = "splice_site_variant"
            result.consequence_priority = CONSEQUENCE_PRIORITY["splice_site_variant"]

    elif nearest_exon_dist <= SPLICE_REGION_DIST:
        result.region_class       = "splice_region"
        result.consequence        = "splice_region_variant"
        result.consequence_priority = CONSEQUENCE_PRIORITY["splice_region_variant"]


def classify(
    features: List["GFF3Feature"],
    nearest_gene_dist: Optional[int] = None,
    nearest_exon_dist: Optional[int] = None,
    regulatory_hits: Optional[List[dict]] = None,
    codon_result: Optional[dict] = None,
    splice_type: Optional[str] = None,
    dbnsfp_scores: Optional[dict] = None,
) -> AnnotationResult:
    """
    Classify a variant given its overlapping GFF3 features.

    Args:
        features:           All GFF3Feature objects overlapping the variant position.
        nearest_gene_dist:  Distance to nearest gene boundary (for intergenic variants).
        nearest_exon_dist:  Distance to nearest exon boundary (for splice classification).
                            Computed by GFF3Index.nearest_exon_boundary_distance().
        regulatory_hits:    List of regulatory element dicts from RegulatoryBEDIndex.query().
                            Each dict has at least: {"source": str, "name": str}.
        codon_result:       Dict from FASTAEngine.get_codon_consequence() — applied for CDS
                            variants to refine consequence (missense, synonymous, stop_gained…).
        splice_type:        'donor' | 'acceptor' from FASTAEngine.get_splice_type() — applied
                            for near-splice intronic variants to distinguish donor/acceptor.
        dbnsfp_scores:      Dict from DbNSFP5Fetcher.fetch() with functional scores.

    Returns:
        AnnotationResult with all annotation fields populated.
    """
    result = AnnotationResult()

    # ── Intergenic branch ─────────────────────────────────────────────────────
    if not features:
        result.is_intergenic = True
        result.region_class  = "intergenic"
        result.distance_nearest_gene = nearest_gene_dist

        # Regulatory BED takes precedence over heuristic upstream window
        if regulatory_hits:
            hit    = regulatory_hits[0]
            source = hit.get("source", "")
            result.region_class             = "regulatory"
            result.consequence              = "regulatory_region_variant"
            result.is_regulatory            = True
            result.is_intergenic            = False
            result.regulatory_source        = f"bed:{source}" if source else "bed"
            result.regulatory_element_type  = hit.get("element_type", "")
            result.regulatory_element_id    = hit.get("element_id", "") or hit.get("name", "")
            result.distance_to_regulatory   = hit.get("distance", 0)
            result.consequence_priority     = CONSEQUENCE_PRIORITY["regulatory_region_variant"]
        elif nearest_gene_dist is not None and nearest_gene_dist <= UPSTREAM_WINDOW:
            result.region_class      = "near_gene"
            result.is_regulatory     = True
            result.regulatory_source = "heuristic"
            result.consequence       = _CONSEQUENCE["near_gene"]
            result.consequence_priority = CONSEQUENCE_PRIORITY["upstream_gene_variant"]
        else:
            result.consequence        = "intergenic_variant"
            result.consequence_priority = CONSEQUENCE_PRIORITY["intergenic_variant"]

        result.distance_to_splice_site = nearest_exon_dist
        return result

    # ── Count distinct overlapping genes ──────────────────────────────────────
    # Group features by gene_id (fall back to gene_name if no gene_id)
    gene_buckets: Dict[str, List["GFF3Feature"]] = defaultdict(list)
    for f in features:
        key = f.gene_id or f.gene_name or "__unknown__"
        gene_buckets[key].append(f)

    # Per-gene best feature
    gene_hits: List[tuple] = []  # (feature_priority, not_protein_coding, best_feat, region_class)
    for gfeats in gene_buckets.values():
        best = _best_feature_for_gene(gfeats)
        prio = _FEATURE_PRIORITY.get(best.feature, 999)
        rc   = _REGION_CLASS.get(best.feature, "intergenic")
        gene_hits.append((prio, int(best.gene_type != "protein_coding"), best, rc))

    # Sort: most specific feature first, protein_coding preferred at ties
    gene_hits.sort(key=lambda x: (x[0], x[1]))

    # ── Populate best-hit fields ──────────────────────────────────────────────
    _prio, _not_pc, best_feat, best_rc = gene_hits[0]

    result.gene_name       = best_feat.gene_name
    result.gene_id         = best_feat.gene_id
    result.gene_type       = best_feat.gene_type
    result.transcript_id   = best_feat.transcript_id
    result.transcript_type = best_feat.transcript_type
    result.feature_type    = best_feat.feature
    result.strand          = best_feat.strand
    result.region_class    = best_rc
    result.consequence     = _CONSEQUENCE.get(best_rc, "intergenic_variant")
    result.is_coding       = (best_rc == "coding")
    result.is_intergenic   = False
    result.n_overlapping_genes = len(gene_buckets)
    result.consequence_priority = CONSEQUENCE_PRIORITY.get(result.consequence, 900)

    # Regulatory BED check for genic variants (e.g. promoter within a gene body)
    if regulatory_hits:
        hit    = regulatory_hits[0]
        source = hit.get("source", "")
        result.is_regulatory           = True
        result.regulatory_source       = f"bed:{source}" if source else "bed"
        result.regulatory_element_type = hit.get("element_type", "")
        result.regulatory_element_id   = hit.get("element_id", "") or hit.get("name", "")
        result.distance_to_regulatory  = hit.get("distance", 0)

    # ── All-hits lists ────────────────────────────────────────────────────────
    all_gene_names = [h[2].gene_name       for h in gene_hits]
    all_tx_ids     = [h[2].transcript_id   for h in gene_hits]
    all_rcs        = [h[3]                  for h in gene_hits]

    result.all_overlapping_genes      = ",".join(g for g in all_gene_names if g)
    result.all_overlapping_transcripts = ",".join(t for t in all_tx_ids   if t)
    result.all_region_classes         = ",".join(all_rcs)

    # ── Splice site classification ────────────────────────────────────────────
    _apply_splice_classification(result, nearest_exon_dist, splice_type)

    # ── Codon consequence (v1.2, CDS variants only) ───────────────────────────
    if codon_result and result.is_coding and "error" not in codon_result:
        csq = codon_result.get("consequence", "")
        if csq and csq != "coding_sequence_variant":
            result.consequence          = csq
            result.consequence_priority = CONSEQUENCE_PRIORITY.get(csq, 300)
        result.codon_ref = codon_result.get("codon_ref", "")
        result.codon_alt = codon_result.get("codon_alt", "")
        result.aa_ref    = codon_result.get("aa_ref", "")
        result.aa_alt    = codon_result.get("aa_alt", "")
        result.aa_ref_3  = codon_result.get("aa_ref_3", "")
        result.aa_alt_3  = codon_result.get("aa_alt_3", "")

        # Boolean flags
        result.is_missense   = (csq == "missense_variant")
        result.is_synonymous = (csq == "synonymous_variant")
        result.is_lof        = csq in {
            "stop_gained", "stop_lost", "start_lost", "frameshift_variant"
        }

    # ── dbNSFP5 scores (v1.2) ─────────────────────────────────────────────────
    if dbnsfp_scores:
        result.cadd_phred    = dbnsfp_scores.get("cadd_phred")
        result.revel_score   = dbnsfp_scores.get("revel_score")
        result.sift_pred     = dbnsfp_scores.get("sift_pred") or ""
        result.polyphen2_pred = dbnsfp_scores.get("polyphen2_pred") or ""
        result.clinvar_clnsig = dbnsfp_scores.get("clinvar_clnsig") or ""

    return result


# ── Summary ───────────────────────────────────────────────────────────────────

def summarise(annotations: List[AnnotationResult]) -> dict:
    """
    Build a summary dict counting variants per region class and consequence.
    Suitable for direct JSON serialisation.
    """
    from collections import Counter

    total = len(annotations)
    region_counts:     Counter = Counter(a.region_class for a in annotations)
    consequence_counts: Counter = Counter(a.consequence for a in annotations)
    gene_type_counts:  Counter = Counter(
        a.gene_type for a in annotations if a.gene_type
    )

    n_coding     = sum(1 for a in annotations if a.is_coding)
    n_regulatory = sum(1 for a in annotations if a.is_regulatory)

    # v1.5: regulatory element type counts
    from collections import Counter as _Counter
    reg_type_counts: _Counter = _Counter(
        a.regulatory_element_type
        for a in annotations
        if a.is_regulatory and a.regulatory_element_type
    )
    n_intergenic = sum(1 for a in annotations if a.is_intergenic)
    n_splice_site   = sum(1 for a in annotations if a.consequence == "splice_site_variant")
    n_splice_region = sum(1 for a in annotations if a.consequence == "splice_region_variant")
    # v1.2 counts
    n_splice_donor    = sum(1 for a in annotations if a.consequence == "splice_donor_variant")
    n_splice_acceptor = sum(1 for a in annotations if a.consequence == "splice_acceptor_variant")
    n_missense    = sum(1 for a in annotations if a.is_missense)
    n_synonymous  = sum(1 for a in annotations if a.is_synonymous)
    n_lof         = sum(1 for a in annotations if a.is_lof)
    n_stop_gained = sum(1 for a in annotations if a.consequence == "stop_gained")
    n_frameshift  = sum(1 for a in annotations if a.consequence == "frameshift_variant")

    # Splice distance distribution (for variants with a computed distance)
    splice_dists = [
        a.distance_to_splice_site
        for a in annotations
        if a.distance_to_splice_site is not None
    ]

    return {
        "total_variants":  total,
        "n_coding":        n_coding,
        "n_regulatory":    n_regulatory,
        "n_intergenic":    n_intergenic,
        "n_genic":         total - n_intergenic,
        "n_splice_site":   n_splice_site,
        "n_splice_region": n_splice_region,
        # v1.2
        "n_splice_donor":    n_splice_donor,
        "n_splice_acceptor": n_splice_acceptor,
        "n_missense":        n_missense,
        "n_synonymous":      n_synonymous,
        "n_lof":             n_lof,
        "n_stop_gained":     n_stop_gained,
        "n_frameshift":      n_frameshift,
        "pct_coding":      round(100 * n_coding / total, 2) if total else 0,
        "pct_intergenic":  round(100 * n_intergenic / total, 2) if total else 0,
        "pct_missense":    round(100 * n_missense / total, 2) if total else 0,
        "pct_lof":         round(100 * n_lof / total, 2) if total else 0,
        "region_class_counts":              dict(region_counts.most_common()),
        "consequence_counts":               dict(consequence_counts.most_common()),
        "gene_type_counts":                 dict(gene_type_counts.most_common(20)),
        "regulatory_element_type_counts":   dict(reg_type_counts.most_common()),
        "splice_distance_stats": {
            "n_with_distance": len(splice_dists),
            "n_canonical":     sum(1 for d in splice_dists if d <= SPLICE_CANONICAL_DIST),
            "n_region":        sum(1 for d in splice_dists if SPLICE_CANONICAL_DIST < d <= SPLICE_REGION_DIST),
        },
    }
