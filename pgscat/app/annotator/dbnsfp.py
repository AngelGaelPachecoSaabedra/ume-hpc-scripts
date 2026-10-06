"""
dbNSFP5 Score Fetcher
=====================
Queries the tabix-indexed dbNSFP5 file for variant functional scores.

Fetches per-variant:
  - CADD_phred        CADD pathogenicity score (phred-scaled)
  - REVEL_score       REVEL ensemble missense pathogenicity score
  - SIFT_pred         SIFT prediction (D=deleterious, T=tolerated)
  - Polyphen2_HDIV_pred  PolyPhen-2 HumDiv prediction (D/P/B)
  - clinvar_clnsig    ClinVar clinical significance

dbNSFP5 is 1-based, tab-separated, bgzip-compressed and tabix-indexed.
Column indices are resolved once from the header line at open time.
Multi-transcript fields (semicolon-separated) are collapsed to the first
non-missing value.

Requires:
  pysam >= 0.22
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional, TYPE_CHECKING

logger = logging.getLogger(__name__)

# Fields we want to extract (dbNSFP5 column names, exact match)
_WANTED_COLS = {
    "CADD_phred",
    "REVEL_score",
    "SIFT_pred",
    "Polyphen2_HDIV_pred",
    "clinvar_clnsig",
}

# Chromosome name remapping for dbNSFP5 tabix query (uses bare names, same as betamap)
# dbNSFP5_grch38 has contigs named "1", "2", ..., "X", "Y", "MT" (no "chr")
_CHROM_REMAP_IN: dict[str, str] = {"MT": "MT", "M": "MT"}


def _strip_chr(name: str) -> str:
    bare = name[3:] if name.startswith("chr") else name
    return _CHROM_REMAP_IN.get(bare, bare)


def _first_nonmissing(value: str) -> str:
    """
    dbNSFP5 stores multi-transcript values as 'val1;val2;...'.
    Returns the first entry that is not '.' or empty, else ''.
    """
    if not value or value == ".":
        return ""
    for part in value.split(";"):
        part = part.strip()
        if part and part != ".":
            return part
    return ""


class DbNSFP5Fetcher:
    """
    Random-access fetcher for dbNSFP5 tabix-indexed scores.

    Usage:
        fetcher = DbNSFP5Fetcher("/path/to/dbNSFP5.0a_grch38.gz")
        scores  = fetcher.fetch("1", 149934520, "T", "C")
        # → {"cadd_phred": 12.3, "revel_score": 0.45, ...}
    """

    def __init__(self, dbnsfp_path: Path) -> None:
        import pysam
        self.path = Path(dbnsfp_path)
        logger.info("Opening dbNSFP5: %s", self.path)
        self._tbx = pysam.TabixFile(str(self.path))

        # Resolve column indices from header
        self._col_idx: dict[str, int] = {}
        try:
            header_line = next(
                line for line in self._tbx.header
                if line.startswith("#chr\t") or line.startswith("chr\t")
                or "\tref\t" in line
            )
            # Strip leading '#' if present
            header_line = header_line.lstrip("#")
            cols = header_line.rstrip("\n").split("\t")
            for wanted in _WANTED_COLS:
                try:
                    self._col_idx[wanted] = cols.index(wanted)
                except ValueError:
                    logger.warning("dbNSFP5 column not found: %s", wanted)
        except StopIteration:
            logger.warning("dbNSFP5 header not found — falling back to hardcoded indices")
            # Fallback: dbNSFP5.0a_grch38 known column indices (0-based)
            # Verified against dbNSFP5.0a header
            self._col_idx = {
                "CADD_phred":         133,   # column 134 (1-based)
                "REVEL_score":         74,   # column 75
                "SIFT_pred":           37,   # column 38
                "Polyphen2_HDIV_pred": 43,   # column 44
                "clinvar_clnsig":     370,   # column 371
            }

        # Fixed columns: chromosome=0, pos=1, ref=2, alt=3
        self._chr_col = 0
        self._pos_col = 1
        self._ref_col = 2
        self._alt_col = 3

        logger.info(
            "dbNSFP5 ready. Columns mapped: %s",
            {k: v for k, v in self._col_idx.items()},
        )

    def fetch(
        self,
        chrom: str,
        pos: int,
        ref: str,
        alt: str,
    ) -> dict:
        """
        Fetch scores for a specific variant.

        Args:
            chrom:  Bare chromosome name (e.g. "1", "X", "MT")
            pos:    1-based genomic position
            ref:    Reference allele (uppercase)
            alt:    Alternate allele (uppercase)

        Returns:
            dict with keys: cadd_phred, revel_score, sift_pred,
                            polyphen2_pred, clinvar_clnsig
            Values are float | str | None.  None = not found.
        """
        empty = {
            "cadd_phred":    None,
            "revel_score":   None,
            "sift_pred":     None,
            "polyphen2_pred": None,
            "clinvar_clnsig": None,
        }

        chrom_q = _strip_chr(str(chrom))
        ref_up  = ref.upper()
        alt_up  = alt.upper()

        try:
            rows = list(self._tbx.fetch(chrom_q, pos - 1, pos))
        except (ValueError, KeyError):
            return empty

        if not rows:
            return empty

        for row in rows:
            fields = row.rstrip("\n").split("\t")
            # Match exact position, ref, alt
            try:
                if (int(fields[self._pos_col]) != pos
                        or fields[self._ref_col].upper() != ref_up
                        or fields[self._alt_col].upper() != alt_up):
                    continue
            except (ValueError, IndexError):
                continue

            result = {}

            # CADD_phred → float
            raw = self._get_field(fields, "CADD_phred")
            v = _first_nonmissing(raw)
            result["cadd_phred"] = float(v) if v else None

            # REVEL_score → float
            raw = self._get_field(fields, "REVEL_score")
            v = _first_nonmissing(raw)
            result["revel_score"] = float(v) if v else None

            # SIFT_pred → string (D/T)
            raw = self._get_field(fields, "SIFT_pred")
            result["sift_pred"] = _first_nonmissing(raw) or None

            # Polyphen2_HDIV_pred → string (D/P/B)
            raw = self._get_field(fields, "Polyphen2_HDIV_pred")
            result["polyphen2_pred"] = _first_nonmissing(raw) or None

            # clinvar_clnsig → string
            raw = self._get_field(fields, "clinvar_clnsig")
            result["clinvar_clnsig"] = _first_nonmissing(raw) or None

            return result

        return empty

    def _get_field(self, fields: list, col_name: str) -> str:
        idx = self._col_idx.get(col_name)
        if idx is None or idx >= len(fields):
            return "."
        return fields[idx]

    def close(self) -> None:
        try:
            self._tbx.close()
        except Exception:
            pass

    def __del__(self) -> None:
        self.close()
