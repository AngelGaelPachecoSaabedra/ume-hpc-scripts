"""
GFF3 Parser and Interval Index Builder
========================================
Stream-parses a GENCODE GFF3 file and builds per-chromosome, per-feature-type
interval trees for fast variant-to-genome overlap queries.

Features retained: gene, transcript, exon, CDS, five_prime_UTR, three_prime_UTR

Chromosome naming:
  GFF3 uses "chr1", "chrX", "chrM" → normalised to bare: "1", "X", "MT"
  Betamap CHROM column uses bare integers/letters: 1, 2, ..., 22, X, Y, MT

Interval convention:
  GFF3 is 1-based inclusive [start, end].
  IntervalTree uses 0-based half-open [begin, end).
  Stored as: begin = gff_start - 1, end = gff_end
  Query for 1-based POS: query_point = POS - 1

GFF3Feature.begin / .end:
  Each feature stores its own interval coordinates (0-based half-open) so that
  callers (e.g. splice-distance computation) can compute boundary distances without
  a second tree query.
"""
import gzip
import logging
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from intervaltree import IntervalTree

logger = logging.getLogger(__name__)

# Only these GFF3 feature types are loaded into the index.
KEPT_FEATURES = frozenset({
    "gene",
    "transcript",
    "exon",
    "CDS",
    "five_prime_UTR",
    "three_prime_UTR",
})

# Chromosome name normalisation
_CHROM_REMAP = {"M": "MT", "chrM": "MT"}


def _strip_chr(name: str) -> str:
    """Remove 'chr' prefix and remap M→MT."""
    bare = name[3:] if name.startswith("chr") else name
    return _CHROM_REMAP.get(bare, bare)


def _parse_attrs(raw: str) -> dict:
    """Parse semicolon-separated GFF3 attribute string into a dict."""
    result: dict = {}
    for part in raw.split(";"):
        part = part.strip()
        if "=" in part:
            k, v = part.split("=", 1)
            result[k.strip()] = v.strip()
    return result


def _strip_version(ensembl_id: str) -> str:
    """Remove Ensembl version suffix: ENSG00000223972.5 → ENSG00000223972."""
    return ensembl_id.split(".")[0] if ensembl_id else ""


def _open_gff3(path: Path):
    """Open GFF3 (plain or .gz) as a text stream."""
    p = str(path)
    return gzip.open(p, "rt", encoding="utf-8") if p.endswith(".gz") else open(p, "r", encoding="utf-8")


class GFF3Feature:
    """
    Lightweight, slot-allocated record stored as interval data.
    All string fields default to "" to keep JSON serialisation simple.

    begin / end: 0-based half-open interval [begin, end) matching IntervalTree
    convention.  They mirror the iv.begin / iv.end of the containing interval
    and are used by distance calculations without a second tree query.
    """
    __slots__ = (
        "feature", "gene_name", "gene_id", "gene_type",
        "transcript_id", "transcript_type", "strand",
        "begin", "end",  # 0-based half-open coordinates
    )

    def __init__(
        self,
        feature: str,
        gene_name: str = "",
        gene_id: str = "",
        gene_type: str = "",
        transcript_id: str = "",
        transcript_type: str = "",
        strand: str = ".",
        begin: int = 0,
        end: int = 0,
    ) -> None:
        self.feature = feature
        self.gene_name = gene_name
        self.gene_id = gene_id
        self.gene_type = gene_type
        self.transcript_id = transcript_id
        self.transcript_type = transcript_type
        self.strand = strand
        self.begin = begin
        self.end = end

    def to_dict(self) -> dict:
        """Return a JSON-serialisable dict (coordinates excluded)."""
        return {
            "feature":          self.feature,
            "gene_name":        self.gene_name,
            "gene_id":          self.gene_id,
            "gene_type":        self.gene_type,
            "transcript_id":    self.transcript_id,
            "transcript_type":  self.transcript_type,
            "strand":           self.strand,
        }


class GFF3Index:
    """
    Hierarchical interval index built from a GENCODE GFF3 file.

    Structure:
        trees[chrom_bare][feature_type] → IntervalTree[GFF3Feature]

    Usage:
        idx = GFF3Index()
        idx.build_from_gff3(Path("/path/to/gencode.v49.basic.annotation.gff3"))
        hits = idx.query("1", 12345)   # returns List[GFF3Feature]
    """

    def __init__(self) -> None:
        # trees[chrom][feature] = IntervalTree
        self.trees: Dict[str, Dict[str, IntervalTree]] = defaultdict(
            lambda: defaultdict(IntervalTree)
        )
        self.n_features: int = 0
        self.chromosomes_seen: set = set()

        # transcript_cds_map[transcript_id] = [(chrom, begin, end, strand), ...]
        # begin/end are 0-based half-open coordinates (same as IntervalTree).
        # Populated during build_from_gff3 for CDS features only.
        self.transcript_cds_map: Dict[str, List[Tuple[str, int, int, str]]] = defaultdict(list)

    # ── Build ─────────────────────────────────────────────────────────────────

    def build_from_gff3(self, gff3_path: Path) -> None:
        """
        Stream-parse GFF3 and populate interval trees.
        Reads .gz or plain files automatically.
        Logs progress every 500k lines.
        """
        logger.info("Building GFF3 index from: %s", gff3_path)
        n_lines = 0
        n_kept = 0
        n_skipped = 0

        with _open_gff3(gff3_path) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue

                n_lines += 1
                if n_lines % 500_000 == 0:
                    logger.info(
                        "  GFF3 progress: %d lines, %d features kept, %d skipped",
                        n_lines, n_kept, n_skipped,
                    )

                parts = line.rstrip("\n").split("\t")
                if len(parts) < 9:
                    continue

                feature = parts[2]
                if feature not in KEPT_FEATURES:
                    n_skipped += 1
                    continue

                seqname = parts[0]
                try:
                    gff_start = int(parts[3])  # 1-based inclusive
                    gff_end   = int(parts[4])  # 1-based inclusive
                except ValueError:
                    n_skipped += 1
                    continue

                if gff_start > gff_end:
                    n_skipped += 1
                    continue

                strand = parts[6]
                attrs  = _parse_attrs(parts[8])

                iv_begin = gff_start - 1  # 0-based
                iv_end   = gff_end         # 0-based exclusive

                feat = GFF3Feature(
                    feature=feature,
                    gene_name=attrs.get("gene_name", ""),
                    gene_id=_strip_version(attrs.get("gene_id", "")),
                    gene_type=attrs.get("gene_type", ""),
                    transcript_id=_strip_version(attrs.get("transcript_id", "")),
                    transcript_type=attrs.get("transcript_type", ""),
                    strand=strand,
                    begin=iv_begin,
                    end=iv_end,
                )

                chrom = _strip_chr(seqname)

                # Store as 0-based half-open: [gff_start-1, gff_end)
                self.trees[chrom][feature][iv_begin : iv_end] = feat
                self.chromosomes_seen.add(chrom)
                n_kept += 1

                # Populate CDS map for codon-offset lookups (FASTAEngine)
                if feature == "CDS":
                    tid = _strip_version(attrs.get("transcript_id", ""))
                    if tid:
                        self.transcript_cds_map[tid].append(
                            (chrom, iv_begin, iv_end, strand)
                        )

        self.n_features = n_kept
        logger.info(
            "GFF3 index ready: %d features across %d chromosomes (skipped %d)",
            self.n_features, len(self.chromosomes_seen), n_skipped,
        )

    # ── Query ─────────────────────────────────────────────────────────────────

    def query(self, chrom: str, pos: int) -> List[GFF3Feature]:
        """
        Return all GFF3Feature objects overlapping 1-based position pos.

        Args:
            chrom: chromosome without 'chr' prefix (e.g. "1", "X", "MT")
            pos:   variant position, 1-based

        Returns:
            List of GFF3Feature; empty list if no overlap or unknown chrom.
        """
        chrom = _strip_chr(str(chrom))
        if chrom not in self.trees:
            return []

        query_pt = pos - 1  # convert 1-based → 0-based for IntervalTree query
        results: List[GFF3Feature] = []
        for tree in self.trees[chrom].values():
            for iv in tree[query_pt]:
                results.append(iv.data)
        return results

    def nearest_exon_boundary_distance(
        self, chrom: str, pos: int
    ) -> Optional[int]:
        """
        Return distance (bp) from pos to the nearest exon boundary (start or end).

        Uses an expanding search window so the query remains O(log n + k) for
        nearby exons while still finding boundaries in long introns.

        Window sequence: 200 bp → 2000 bp → 20 000 bp.
        Stops expanding once a candidate is found within the current window.

        Returns:
            0  if the variant lies exactly on an exon boundary.
            >0 distance to the nearest exon start/end.
            None if no exon data is available for this chromosome.

        Args:
            chrom: bare chromosome name (e.g. "1", "X", "MT")
            pos:   variant position, 1-based
        """
        chrom = _strip_chr(str(chrom))
        if chrom not in self.trees or "exon" not in self.trees[chrom]:
            return None

        exon_tree = self.trees[chrom]["exon"]
        if not exon_tree:
            return None

        qp = pos - 1  # convert to 0-based

        min_dist: Optional[int] = None
        for window in (200, 2_000, 20_000):
            lo = max(0, qp - window)
            hi = qp + window + 1
            hits = exon_tree[lo:hi]
            if hits:
                for iv in hits:
                    # Exon boundaries in 0-based coords:
                    #   start boundary: iv.begin  (first nt of exon)
                    #   end boundary:   iv.end - 1 (last nt of exon)
                    d = min(abs(qp - iv.begin), abs(qp - (iv.end - 1)))
                    if min_dist is None or d < min_dist:
                        min_dist = d
                # Any boundary found in this window is ≤ window bp away —
                # no exon outside the window can be closer.
                if min_dist is not None and min_dist <= window:
                    break

        return int(min_dist) if min_dist is not None else None

    def nearest_gene_distance(self, chrom: str, pos: int) -> Optional[int]:
        """
        For intergenic variants: return the distance (bp) to the nearest gene boundary.
        Returns None if no gene data available for this chromosome.
        """
        chrom = _strip_chr(str(chrom))
        if chrom not in self.trees or "gene" not in self.trees[chrom]:
            return None

        gene_tree = self.trees[chrom]["gene"]
        if not gene_tree:
            return None

        # Scan gene intervals for nearest boundary (O(n) — only called for intergenic)
        min_dist = None
        qp = pos - 1  # 0-based query point
        for iv in gene_tree.all_intervals:
            dist = max(0, iv.begin - qp) if qp < iv.begin else max(0, qp - iv.end + 1)
            if min_dist is None or dist < min_dist:
                min_dist = dist
        return int(min_dist) if min_dist is not None else None

    def get_cds_offset(
        self, transcript_id: str, chrom: str, pos: int
    ) -> Optional[dict]:
        """
        For a CDS variant at 1-based pos on chrom, return its CDS offset.

        Returns:
            dict with keys:
                cds_pos      — 0-based nucleotide position within the CDS
                codon_index  — 0-based codon index (cds_pos // 3)
                pos_in_codon — position within codon (0, 1, or 2)
                strand       — '+' or '-'
            or None if transcript_id not in map or pos not inside any CDS interval.

        Handles multi-exon transcripts correctly:
            + strand: CDS intervals sorted ascending by begin.
            - strand: CDS intervals sorted descending by end (reverse transcription order).
        """
        chrom_bare = _strip_chr(str(chrom))
        intervals = self.transcript_cds_map.get(transcript_id, [])
        if not intervals:
            return None

        strand = intervals[0][3]
        if strand == "+":
            sorted_ivs = sorted(intervals, key=lambda x: x[1])   # ascending begin
        else:
            sorted_ivs = sorted(intervals, key=lambda x: -x[2])  # descending end

        qp = pos - 1  # convert to 0-based
        cumulative = 0

        for (iv_chrom, begin, end, _) in sorted_ivs:
            if iv_chrom != chrom_bare:
                cumulative += (end - begin)
                continue
            exon_len = end - begin
            if begin <= qp < end:
                if strand == "+":
                    offset_in_exon = qp - begin
                else:
                    offset_in_exon = end - 1 - qp
                cds_pos = cumulative + offset_in_exon
                return {
                    "cds_pos":      cds_pos,
                    "codon_index":  cds_pos // 3,
                    "pos_in_codon": cds_pos % 3,
                    "strand":       strand,
                }
            cumulative += exon_len

        return None

    def stats(self) -> dict:
        """Return summary counts by chromosome and feature type."""
        out: dict = {"total_features": self.n_features, "chromosomes": {}}
        for chrom, feat_dict in self.trees.items():
            out["chromosomes"][chrom] = {
                ft: len(tree) for ft, tree in feat_dict.items()
            }
        return out
