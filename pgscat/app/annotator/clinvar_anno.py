"""
clinvar_anno.py — ClinVar GRCh38 Variant Annotator (Annotator)
===============================================================
Streams through the local ClinVar VCF (GRCh38) to annotate each variant with:
  - clinvar_clnsig       : germline classification (Pathogenic, Benign, VUS, etc.)
  - clinvar_clndn        : disease name
  - clinvar_clnrevstat   : review status (criteria_provided_single_submitter, etc.)
  - clinvar_alleleid     : ClinVar Allele ID (integer)
  - clinvar_clnvc        : variant type (SNV, Deletion, etc.)

This replaces the dbNSFP5-derived clinvar_clnsig with the authoritative
ClinVar VCF source, which is more complete and up-to-date.

Source:
    /mnt/cephfs/hot_nvme/clinvar/GRCh38/clinvar.vcf.gz  (BGZF, no .tbi)
    Updated: April 2026

ClinVar chromosome naming: bare numbers (1, 2, ..., X, Y) — no chr prefix.

Strategy (no tabix index):
    Pre-index approach — build a dict keyed by (chrom, pos, ref, alt) at
    annotator startup. Memory: ~1.2M variants × ~200 bytes ≈ ~240MB RAM.
    Acceptable for HPC annotation jobs.

Usage in annotate_variants.py:
    from clinvar_anno import ClinVarAnnotator
    cv = ClinVarAnnotator("/mnt/cephfs/hot_nvme/clinvar/GRCh38/clinvar.vcf.gz")
    result = cv.lookup("1", 66926, "AG", "A")
"""
from __future__ import annotations

import logging
import time
from pathlib import Path
from typing import Dict, Optional, Tuple

logger = logging.getLogger(__name__)

ClinVarKey = Tuple[str, int, str, str]  # (chrom, pos, ref, alt)

ClinVarRecord = Dict[str, Optional[object]]


class ClinVarAnnotator:
    """
    In-memory ClinVar annotation index (build once, query many times).

    Loads the full ClinVar VCF into a dict at construction time (~30-90s).
    After loading, lookups are O(1).
    """

    def __init__(self, vcf_path: str) -> None:
        self._path = Path(vcf_path)
        self._index: Dict[ClinVarKey, ClinVarRecord] = {}
        self._loaded = False
        self._load()

    def _load(self) -> None:
        if not self._path.exists():
            logger.warning("ClinVar VCF not found: %s", self._path)
            return

        try:
            import pysam
        except ImportError:
            logger.warning("pysam not available — ClinVar annotation disabled")
            return

        logger.info("Loading ClinVar index from %s …", self._path)
        t0 = time.time()
        n = 0

        try:
            f = pysam.BGZFile(str(self._path), "r")
            for raw in f:
                line = raw.decode("utf-8", errors="replace") if isinstance(raw, bytes) else raw
                if line.startswith("#"):
                    continue

                parts = line.rstrip("\n").split("\t")
                if len(parts) < 8:
                    continue

                chrom = parts[0]
                try:
                    pos = int(parts[1])
                except ValueError:
                    continue

                ref = parts[3].upper()
                alt = parts[4].upper()
                info_str = parts[7]

                info: Dict[str, str] = {}
                for kv in info_str.split(";"):
                    if "=" in kv:
                        k, v = kv.split("=", 1)
                        info[k] = v

                record: ClinVarRecord = {
                    "clinvar_clnsig":     self._clean(info.get("CLNSIG")),
                    "clinvar_clndn":      self._clean(info.get("CLNDN")),
                    "clinvar_clnrevstat": self._clean(info.get("CLNREVSTAT")),
                    "clinvar_alleleid":   self._safe_int(info.get("ALLELEID")),
                    "clinvar_clnvc":      self._clean(info.get("CLNVC")),
                }

                # Multi-allelic: split by comma in ALT
                for a in alt.split(","):
                    a = a.strip()
                    if a:
                        self._index[(chrom, pos, ref, a)] = record
                        # Also store without chr prefix and with chr prefix
                        bare = chrom[3:] if chrom.startswith("chr") else chrom
                        self._index[(bare, pos, ref, a)] = record

                n += 1

            f.close()
            elapsed = time.time() - t0
            logger.info("ClinVar index loaded: %d variants in %.1fs", n, elapsed)
            self._loaded = True

        except Exception as exc:
            logger.error("ClinVar index load failed: %s", exc)

    def lookup(
        self,
        chrom: str,
        pos: int,
        ref: str,
        alt: str,
    ) -> Optional[ClinVarRecord]:
        """
        Lookup a variant in the ClinVar index.

        Returns ClinVarRecord dict or None if not found.
        """
        if not self._loaded:
            return None

        ref_u = ref.upper()
        alt_u = alt.upper()
        bare = chrom[3:] if chrom.startswith("chr") else chrom

        return (
            self._index.get((bare, pos, ref_u, alt_u))
            or self._index.get((f"chr{bare}", pos, ref_u, alt_u))
            or self._index.get((chrom, pos, ref_u, alt_u))
        )

    @staticmethod
    def _clean(v: Optional[str]) -> Optional[str]:
        if v is None or v in (".", ""):
            return None
        return v.replace("_", " ").replace("|", "; ")

    @staticmethod
    def _safe_int(v: Optional[str]) -> Optional[int]:
        try:
            return int(v)
        except (TypeError, ValueError):
            return None

    def is_loaded(self) -> bool:
        return self._loaded

    def __len__(self) -> int:
        return len(self._index)
