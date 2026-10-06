"""
Regulatory BED Index  (v2)
==========================
Loads one or more BED files containing regulatory elements (ENCODE TFBS/DHS,
ORegAnno, Ensembl Regulatory Build, custom) and builds per-chromosome
interval trees for fast overlap and nearest-element queries.

v2 additions over v1:
  • element_type inference from TF/element name (CTCF → insulator, H3K27ac → enhancer, …)
  • distance_to_regulatory_element — nearest element distance even when no overlap
  • query() returns dicts with: name, source, element_type, element_id
  • nearest() returns (hit_or_None, distance) for non-overlapping queries

Supported input formats
-----------------------
Standard BED (default):
    col 1: chrom        e.g. chr1 or 1
    col 2: chromStart   0-based inclusive
    col 3: chromEnd     0-based exclusive (standard BED half-open)
    col 4: name         optional
    Lines beginning with '#', 'track', or 'browser' are skipped.

UCSC table format (load_ucsc_table()):
    Same as BED but with a leading integer "bin" column that is ignored.
    Used for files from /mnt/cephfs/hot_nvme/ucsc/:
      encRegTfbsClustered.txt.gz     (ENCODE TFBS clusters)
      wgEncodeRegDnaseClustered.txt.gz (ENCODE DNase clusters)
      oreganno.txt.gz                (ORegAnno curated elements)

Both plain and .gz-compressed files are supported.

Element type inference
----------------------
Inferred from the name column using keyword heuristics:
  CTCF                    → insulator
  H3K4me3                 → promoter
  H3K27ac, H3K4me1        → enhancer
  EP300 / p300 / H3K4me1  → enhancer
  DNase, ATAC, open_chrom → open_chromatin
  contains "promoter"     → promoter
  contains "enhancer"     → enhancer
  OREG…                   → regulatory_region (ORegAnno)
  other named TF          → TFBS
  unnamed                 → regulatory_element

Usage example
-------------
    from pathlib import Path
    from regulatory import RegulatoryBEDIndex

    reg = RegulatoryBEDIndex()
    # Standard BED
    reg.load_bed(Path("/ref/encode_cCRE.bed.gz"), "ENCODE_cCRE")
    # UCSC table format (strip bin column)
    reg.load_ucsc_table(Path("/ucsc/encRegTfbsClustered.txt.gz"), "ENCODE_TFBS")
    reg.load_ucsc_table(Path("/ucsc/wgEncodeRegDnaseClustered.txt.gz"), "ENCODE_DHS",
                        name_col=None, default_name="open_chromatin")

    hits = reg.query("1", 12345)
    # → [{"name": "CTCF", "source": "ENCODE_TFBS", "element_type": "insulator",
    #      "element_id": "CTCF", "start": 12300, "end": 12400}, ...]

    # Nearest element (even if no overlap)
    hit, dist = reg.nearest("1", 12345)
    # → ({"element_type": "enhancer", ...}, 120)

When no BED files are loaded (reg.is_empty), all queries return [] / (None, None).
"""
import gzip
import logging
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from intervaltree import IntervalTree, Interval

logger = logging.getLogger(__name__)

# ── Chromosome normalisation ──────────────────────────────────────────────────

_CHROM_REMAP = {"M": "MT", "chrM": "MT"}


def _strip_chr(name: str) -> str:
    bare = name[3:] if name.startswith("chr") else name
    return _CHROM_REMAP.get(bare, bare)


def _open_file(path: Path):
    p = str(path)
    return gzip.open(p, "rt", encoding="utf-8", errors="replace") if p.endswith(".gz") else open(p, "r", encoding="utf-8")


# ── Source priority (lower = higher priority) ─────────────────────────────────

_SOURCE_PRIORITY: Dict[str, int] = {
    "ENCODE_cCRE_V3":         1,
    "cCRE_V3":                1,
    "encodeCcreCombined":     1,
    "ENCODE_TFBS":            2,
    "encRegTfbsClustered":    2,
    "ENCODE_DHS":             3,
    "wgEncodeRegDnaseClustered": 3,
    "ORegAnno":               4,
    "oreganno":               4,
}


# ── cCRE V3 element type mapping ──────────────────────────────────────────────

# ENCODE cCRE encodeLabel values → element type
_CCRE_TYPE_MAP: Dict[str, str] = {
    "PLS":           "promoter",
    "pELS":          "proximal_enhancer",
    "dELS":          "distal_enhancer",
    "CTCF-only":     "insulator",
    "DNase-H3K4me3": "open_chromatin",
}


# ── Element type inference ────────────────────────────────────────────────────

# TF → canonical regulatory element type
_TF_TO_TYPE: Dict[str, str] = {
    # Insulators
    "CTCF":     "insulator",
    "CTCFL":    "insulator",
    "RAD21":    "insulator",
    "SMC3":     "insulator",
    # Active promoter marks
    "H3K4ME3":  "promoter",
    "H3K9AC":   "promoter",
    "H3K14AC":  "promoter",
    # Active enhancer marks
    "H3K27AC":  "enhancer",
    "H3K4ME1":  "enhancer",
    "H3K36ME3": "enhancer",
    "EP300":    "enhancer",
    "P300":     "enhancer",
    # Open chromatin
    "DNASE":       "open_chromatin",
    "OPEN_CHROMATIN": "open_chromatin",
    "ATAC":        "open_chromatin",
}

_NAME_KEYWORDS_TO_TYPE = [
    ("promoter",      "promoter"),
    ("enhancer",      "enhancer"),
    ("insulator",     "insulator"),
    ("ctcf",          "insulator"),
    ("dnase",         "open_chromatin"),
    ("open_chrom",    "open_chromatin"),
    ("atac",          "open_chromatin"),
    ("h3k27ac",       "enhancer"),
    ("h3k4me1",       "enhancer"),
    ("h3k4me3",       "promoter"),
    ("h3k9ac",        "promoter"),
    ("ep300",         "enhancer"),
    ("oreg",          "regulatory_region"),    # ORegAnno IDs start with OREG
    ("silencer",      "silencer"),
    ("repressor",     "repressor"),
]


def infer_element_type(name: str, source: str = "") -> str:
    """
    Infer regulatory element type from element name and source label.

    Returns one of: promoter | proximal_enhancer | distal_enhancer | enhancer |
                    insulator | TFBS | open_chromatin | regulatory_region |
                    silencer | repressor | regulatory_element
    """
    if not name:
        name_l = ""
    else:
        name_l = name.upper()

    # cCRE V3 encodeLabel exact match (highest priority — unambiguous)
    if name in _CCRE_TYPE_MAP:
        return _CCRE_TYPE_MAP[name]

    # Direct TF name lookup
    if name_l in _TF_TO_TYPE:
        return _TF_TO_TYPE[name_l]

    # Keyword scan on lowercase
    name_lc = name.lower()
    for kw, etype in _NAME_KEYWORDS_TO_TYPE:
        if kw in name_lc:
            return etype

    # If name is non-empty and doesn't match known marks, treat as generic TFBS
    if name:
        return "TFBS"

    return "regulatory_element"


# ── Main class ────────────────────────────────────────────────────────────────

class RegulatoryBEDIndex:
    """
    Interval index built from one or more BED or UCSC-table regulatory files.

    Per-chromosome IntervalTree stores feature dicts with:
        name, source, element_type, element_id, start (0-based), end (exclusive)

    Key methods:
        load_bed()          — load standard BED file
        load_ucsc_table()   — load UCSC table (strip leading bin column)
        query()             — return all overlapping elements at a 1-based position
        nearest()           — return (nearest_feature, distance); 0 if overlapping
        is_empty            — True if no features loaded
        stats()             — summary by chromosome and source
    """

    def __init__(self) -> None:
        self.trees: Dict[str, IntervalTree] = defaultdict(IntervalTree)
        self.n_features: int = 0
        self.sources: List[str] = []

    # ── Internal loader ───────────────────────────────────────────────────────

    def _add_record(
        self,
        chrom_raw: str,
        start: int,
        end: int,
        name: str,
        source: str,
        element_id: Optional[str] = None,
    ) -> bool:
        """Validate and add one record to the index. Returns True if added."""
        if start < 0 or start >= end:
            return False
        chrom = _strip_chr(chrom_raw)
        etype = infer_element_type(name, source)
        feature = {
            "name":         name,
            "source":       source,
            "element_type": etype,
            "element_id":   element_id if element_id is not None else (name or source),
            "priority":     _SOURCE_PRIORITY.get(source, 5),
            "start":        start,
            "end":          end,
        }
        self.trees[chrom][start:end] = feature
        return True

    # ── Public loaders ────────────────────────────────────────────────────────

    def load_bed(
        self,
        bed_path: Path,
        source_name: str = "",
        canonical_only: bool = True,
    ) -> int:
        """
        Load a standard BED file (plain or .gz) into the index.

        Columns: chrom start end [name ...]
        Lines beginning with '#', 'track', or 'browser' are skipped.

        Args:
            bed_path:       Path to the BED file.
            source_name:    Human-readable source label (defaults to file basename).
            canonical_only: Skip alt scaffolds; only load chr1-22, X, Y, MT.

        Returns:
            Number of features loaded.
        """
        if not source_name:
            source_name = bed_path.name
            for suffix in (".bed.gz", ".bed", ".tsv.gz", ".tsv"):
                if source_name.endswith(suffix):
                    source_name = source_name[: -len(suffix)]
                    break

        logger.info("Loading regulatory BED: %s  source=%s", bed_path, source_name)
        n_loaded = 0
        n_skipped = 0
        _canonical = _canonical_chroms()

        with _open_file(bed_path) as fh:
            for raw_line in fh:
                line = raw_line.strip()
                if not line or line.startswith(("#", "track", "browser")):
                    continue
                parts = line.split("\t")
                if len(parts) < 3:
                    n_skipped += 1
                    continue
                chrom = parts[0]
                if canonical_only and _strip_chr(chrom) not in _canonical:
                    continue
                try:
                    start = int(parts[1])
                    end   = int(parts[2])
                except ValueError:
                    n_skipped += 1
                    continue
                name = parts[3].strip() if len(parts) > 3 else ""
                if self._add_record(chrom, start, end, name, source_name):
                    n_loaded += 1
                else:
                    n_skipped += 1

        self._register_source(source_name, n_loaded, n_skipped, bed_path)
        return n_loaded

    def load_ucsc_table(
        self,
        table_path: Path,
        source_name: str = "",
        name_col: Optional[int] = 4,
        default_name: str = "",
        canonical_only: bool = True,
    ) -> int:
        """
        Load a UCSC table file (leading bin column → strip it).

        UCSC table layout (1-indexed columns):
          col1=bin  col2=chrom  col3=chromStart  col4=chromEnd  col5=name  ...

        Args:
            table_path:   Path to the UCSC table file (.txt.gz or .txt).
            source_name:  Human-readable source label.
            name_col:     0-based Python index of the *name* column (after the
                          bin column is stripped). Default 4 → col5 in UCSC =
                          col4 after stripping bin.  Set None to use default_name.
            default_name: Name to assign when name_col is None or the column is
                          empty (e.g. "open_chromatin" for DNase DHS files).
            canonical_only: Skip alt scaffolds.

        Returns:
            Number of features loaded.
        """
        if not source_name:
            source_name = table_path.name
            for suffix in (".txt.gz", ".txt", ".tsv.gz", ".tsv"):
                if source_name.endswith(suffix):
                    source_name = source_name[: -len(suffix)]
                    break

        logger.info("Loading UCSC table: %s  source=%s", table_path, source_name)
        n_loaded = 0
        n_skipped = 0
        _canonical = _canonical_chroms()

        with _open_file(table_path) as fh:
            for raw_line in fh:
                line = raw_line.strip()
                if not line or line.startswith(("#", "track", "browser")):
                    continue
                parts = line.split("\t")
                # UCSC: bin chrom chromStart chromEnd name ...
                # minimum: bin + chrom + start + end = 4 columns
                if len(parts) < 4:
                    n_skipped += 1
                    continue
                chrom = parts[1]    # col[1] after bin col
                if canonical_only and _strip_chr(chrom) not in _canonical:
                    continue
                try:
                    start = int(parts[2])
                    end   = int(parts[3])
                except ValueError:
                    n_skipped += 1
                    continue

                name = default_name
                if name_col is not None and name_col < len(parts):
                    name = parts[name_col].strip() or default_name

                if self._add_record(chrom, start, end, name, source_name):
                    n_loaded += 1
                else:
                    n_skipped += 1

        self._register_source(source_name, n_loaded, n_skipped, table_path)
        return n_loaded

    def load_ccre_bed(
        self,
        bed_path: Path,
        source_name: str = "ENCODE_cCRE_V3",
        canonical_only: bool = True,
    ) -> int:
        """
        Load ENCODE cCRE V3 BED file (5 columns: chrom start end ccreType accession).

        Generated by build_regulatory_beds.py --ccre from encodeCcreCombined.bb.

        Column layout (0-indexed):
          col[0] = chrom
          col[1] = chromStart  (0-based)
          col[2] = chromEnd
          col[3] = ccreType    (PLS / pELS / dELS / CTCF-only / DNase-H3K4me3)
          col[4] = accession   (EH38E…)

        The ccreType is the element name passed to infer_element_type().
        The accession becomes element_id for downstream display.

        Returns: number of features loaded.
        """
        logger.info("Loading ENCODE cCRE V3 BED: %s  source=%s", bed_path, source_name)
        n_loaded = 0
        n_skipped = 0
        _canonical = _canonical_chroms()

        with _open_file(bed_path) as fh:
            for raw_line in fh:
                line = raw_line.strip()
                if not line or line.startswith(("#", "track", "browser")):
                    continue
                parts = line.split("\t")
                if len(parts) < 3:
                    n_skipped += 1
                    continue
                chrom = parts[0]
                if canonical_only and _strip_chr(chrom) not in _canonical:
                    continue
                try:
                    start = int(parts[1])
                    end   = int(parts[2])
                except ValueError:
                    n_skipped += 1
                    continue
                ccre_type = parts[3].strip() if len(parts) > 3 else "cCRE"
                accession = parts[4].strip() if len(parts) > 4 else ""
                if self._add_record(chrom, start, end, ccre_type, source_name,
                                    element_id=accession or ccre_type):
                    n_loaded += 1
                else:
                    n_skipped += 1

        self._register_source(source_name, n_loaded, n_skipped, bed_path)
        return n_loaded

    def _register_source(
        self, source_name: str, n_loaded: int, n_skipped: int, path: Path
    ) -> None:
        self.n_features += n_loaded
        if source_name not in self.sources:
            self.sources.append(source_name)
        logger.info(
            "  → %d features loaded from %s (%d skipped)", n_loaded, source_name, n_skipped
        )

    # ── Queries ───────────────────────────────────────────────────────────────

    def query(self, chrom: str, pos: int) -> List[dict]:
        """
        Return all regulatory feature dicts overlapping 1-based position pos,
        sorted by source priority (cCRE V3 first, then TFBS, DHS, ORegAnno).

        Each dict: {name, source, element_type, element_id, priority, start, end}
        Returns [] when no overlap or index is empty.
        """
        if self.is_empty:
            return []
        bare = _strip_chr(str(chrom))
        if bare not in self.trees:
            return []
        qp = pos - 1  # 1-based → 0-based for IntervalTree
        hits = [iv.data for iv in self.trees[bare][qp]]
        hits.sort(key=lambda h: h.get("priority", 5))
        return hits

    def nearest(
        self,
        chrom: str,
        pos: int,
        max_distance: int = 100_000,
    ) -> Tuple[Optional[dict], Optional[int]]:
        """
        Return (feature_dict, distance) for the nearest regulatory element.

        If the position overlaps an element, distance = 0.
        If the nearest element is farther than max_distance, returns (None, None).

        Distance is measured in bp from the query position to the nearest
        edge of the interval (0-based coords internally, 1-based pos input).

        Returns:
            (feature_dict, 0) for overlapping elements.
            (feature_dict, distance_bp) for non-overlapping closest element.
            (None, None) if chrom not in index or nothing within max_distance.
        """
        if self.is_empty:
            return (None, None)
        bare = _strip_chr(str(chrom))
        if bare not in self.trees:
            return (None, None)

        tree = self.trees[bare]
        qp = pos - 1  # 1-based → 0-based

        # Overlapping: distance = 0; pick highest-priority source
        overlapping = sorted(list(tree[qp]), key=lambda iv: iv.data.get("priority", 5))
        if overlapping:
            return (overlapping[0].data, 0)

        # Search expanding windows for nearest element
        # IntervalTree doesn't natively support nearest-neighbour, so we use
        # incremental window expansion — efficient for common cases.
        best_iv:   Optional[Interval] = None
        best_dist: int = max_distance + 1

        # Check in a series of expanding windows
        for radius in [500, 2_000, 10_000, 50_000, max_distance]:
            lo = max(0, qp - radius)
            hi = qp + radius
            candidates = list(tree[lo:hi])
            if candidates:
                for iv in candidates:
                    # Distance to nearest edge
                    if iv.end <= qp:
                        d = qp - iv.end + 1  # position is to the right of interval
                    elif iv.begin > qp:
                        d = iv.begin - qp    # position is to the left
                    else:
                        d = 0
                    if d < best_dist:
                        best_dist, best_iv = d, iv
                break  # found at least one candidate in this radius

        if best_iv is None or best_dist > max_distance:
            return (None, None)

        return (best_iv.data, best_dist)

    # ── Properties / helpers ──────────────────────────────────────────────────

    @property
    def is_empty(self) -> bool:
        return self.n_features == 0

    def stats(self) -> dict:
        return {
            "total_features": self.n_features,
            "sources":        self.sources,
            "chromosomes": {
                chrom: len(tree) for chrom, tree in self.trees.items()
            },
        }


# ── Helpers ───────────────────────────────────────────────────────────────────

def _canonical_chroms() -> frozenset:
    return frozenset(
        [str(i) for i in range(1, 23)] + ["X", "Y", "MT"]
    )


# ── Factory: build from env var REGULATORY_BED ───────────────────────────────

def build_from_env(env_var: str = "REGULATORY_BED") -> "RegulatoryBEDIndex":
    """
    Build a RegulatoryBEDIndex from the REGULATORY_BED environment variable.

    REGULATORY_BED may be a comma- or space-separated list of BED file paths.
    Files that end in .ucsc.gz or are in the ucsc/ directory are loaded as
    UCSC tables; all others as standard BED.

    Example:
        export REGULATORY_BED=/work/regulatory_derived/encode_tfbs_GRCh38.bed.gz,\
                               /work/regulatory_derived/encode_dhs_GRCh38.bed.gz

    Returns an empty index if the env var is not set.
    """
    import os
    raw = os.environ.get(env_var, "").strip()
    idx = RegulatoryBEDIndex()
    if not raw:
        return idx

    # Support comma or space separated
    paths = [p.strip() for p in raw.replace(",", " ").split() if p.strip()]
    for p_str in paths:
        path = Path(p_str)
        if not path.exists():
            logger.warning("REGULATORY_BED: file not found — skipping: %s", path)
            continue
        try:
            idx.load_bed(path)
        except Exception as exc:
            logger.warning("REGULATORY_BED: failed to load %s: %s", path, exc)

    return idx
