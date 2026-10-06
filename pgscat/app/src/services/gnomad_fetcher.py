"""
gnomad_fetcher.py — gnomAD 4.1.1 Population Frequency Fetcher
==============================================================
Queries tabix-indexed gnomAD 4.1.1 VCF (exome + genome) for:
  - af_gnomad          : global allele frequency (exome preferred, genome fallback)
  - af_gnomad_afr      : African/African-American AF
  - af_gnomad_amr      : Admixed American AF (closest proxy for MEX)
  - af_gnomad_eas      : East Asian AF
  - af_gnomad_nfe      : Non-Finnish European AF
  - af_gnomad_sas      : South Asian AF
  - af_gnomad_mid      : Middle Eastern AF
  - an_gnomad          : Total allele number
  - ac_gnomad          : Alt allele count
  - nhomalt_gnomad     : Homozygous count
  - rarity_class_gnomad: rarity classification based on gnomAD AF

Data:
    /mnt/cephfs/hot_nvme/gnomad_4.1.1/variants/exome/   (preferred)
    /mnt/cephfs/hot_nvme/gnomad_4.1.1/variants/genome/  (fallback)

Files: gnomad.exomes.v4.1.1.sites.chrN.vcf.bgz + .tbi
Chromosomes use chr-prefixed names (chr1, chrX, etc.)

Usage (web service):
    fetcher = GnomadFetcher()
    result = fetcher.fetch("1", 925952, "G", "A")
    # {"af_gnomad": 5.48e-6, "af_gnomad_amr": 0.0, ...}
"""
from __future__ import annotations

import logging
import os
from functools import lru_cache
from pathlib import Path
from typing import Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)

# ── Paths — configurable via env vars ─────────────────────────────────────────
# GNOMAD_DIR overrides the base directory (default: /mnt/cephfs/hot_nvme/gnomad_4.1.1/variants)
# GNOMAD_EXOME_DIR / GNOMAD_GENOME_DIR override individual sub-dirs.
_GNOMAD_BASE = Path(os.environ.get("GNOMAD_DIR", "/mnt/cephfs/hot_nvme/gnomad_4.1.1/variants"))
_EXOME_DIR   = Path(os.environ.get("GNOMAD_EXOME_DIR", str(_GNOMAD_BASE / "exome")))
_GENOME_DIR  = Path(os.environ.get("GNOMAD_GENOME_DIR", str(_GNOMAD_BASE / "genome")))

# ── Rarity thresholds (gnomAD recommended) ────────────────────────────────────
_RARITY_THRESHOLDS = [
    ("common",        0.01),
    ("low_frequency", 0.001),
    ("rare",          0.0001),
]


def classify_rarity(af: Optional[float]) -> str:
    if af is None or af <= 0.0:
        return "ultra_rare"
    for label, threshold in _RARITY_THRESHOLDS:
        if af >= threshold:
            return label
    return "ultra_rare"


def _chrom_to_vcf(chrom: str) -> str:
    """Normalise chromosome to chr-prefixed (gnomAD convention)."""
    bare = chrom[3:] if chrom.startswith("chr") else chrom
    return f"chr{bare}"


def _exome_path(chrom_vcf: str) -> Path:
    return _EXOME_DIR / f"gnomad.exomes.v4.1.1.sites.{chrom_vcf}.vcf.bgz"


def _genome_path(chrom_vcf: str) -> Path:
    return _GENOME_DIR / f"gnomad.genomes.v4.1.1.sites.{chrom_vcf}.vcf.bgz"


def _parse_info(info_str: str) -> Dict[str, str]:
    """Parse VCF INFO field to dict."""
    result = {}
    for part in info_str.split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            result[k] = v
        else:
            result[part] = "True"
    return result


def _extract_af(info: Dict[str, str], key: str, alt_idx: int = 0) -> Optional[float]:
    """Extract float AF from INFO, handling multi-allelic comma-split."""
    v = info.get(key)
    if v is None or v in (".", ""):
        return None
    parts = v.split(",")
    try:
        val = float(parts[alt_idx] if alt_idx < len(parts) else parts[0])
        return None if val < 0 else val
    except (ValueError, IndexError):
        return None


def _extract_int(info: Dict[str, str], key: str, alt_idx: int = 0) -> Optional[int]:
    """Extract int from INFO, handling multi-allelic comma-split."""
    v = info.get(key)
    if v is None or v in (".", ""):
        return None
    parts = v.split(",")
    try:
        return int(parts[alt_idx] if alt_idx < len(parts) else parts[0])
    except (ValueError, IndexError):
        return None


class GnomadFetcher:
    """
    Per-variant gnomAD 4.1.1 frequency fetcher using pysam tabix.

    Maintains open TabixFile handles per chromosome (lazy-loaded) to avoid
    repeated open/close overhead when annotating many variants on the same
    chromosome.

    Thread safety: Not thread-safe by default. Use one instance per worker
    process (gunicorn pre-fork model is safe).
    """

    def __init__(self) -> None:
        self._exome_handles: Dict[str, object] = {}   # chrom_vcf → TabixFile
        self._genome_handles: Dict[str, object] = {}
        self._available: Optional[bool] = None
        self._check_availability()

    def _check_availability(self) -> None:
        """
        Check if gnomAD data is available.

        Strategy (in priority order):
          1. Any .vcf.bgz in GNOMAD_EXOME_DIR
          2. Any .vcf.bgz in GNOMAD_GENOME_DIR
        Does NOT require chr1 specifically or a .tbi — those are validated
        lazily per-chromosome in _get_handle().
        """
        for source, directory in [("exome", _EXOME_DIR), ("genome", _GENOME_DIR)]:
            if directory.exists():
                try:
                    files = list(directory.iterdir())
                except OSError:
                    continue
                if any(f.suffix == ".bgz" and f.stem.endswith(".vcf") for f in files):
                    self._available = True
                    logger.info(
                        "GnomadFetcher ready: found %s files in %s (%s)",
                        sum(1 for f in files if f.suffix == ".bgz"),
                        directory,
                        source,
                    )
                    return
        self._available = False
        logger.warning(
            "GnomadFetcher: no .vcf.bgz files found in %s or %s",
            _EXOME_DIR, _GENOME_DIR,
        )

    def is_available(self) -> bool:
        return bool(self._available)

    def _get_handle(self, chrom_vcf: str, source: str):
        """Lazily open and cache a TabixFile handle."""
        import pysam

        if source == "exome":
            handles = self._exome_handles
            path = _exome_path(chrom_vcf)
        else:
            handles = self._genome_handles
            path = _genome_path(chrom_vcf)

        if chrom_vcf not in handles:
            if not path.exists():
                return None
            tbi = Path(str(path) + ".tbi")
            if not tbi.exists():
                logger.debug("No .tbi for %s", path)
                return None
            try:
                handles[chrom_vcf] = pysam.TabixFile(str(path))
                logger.debug("Opened gnomAD %s %s", source, chrom_vcf)
            except Exception as exc:
                logger.warning("Cannot open gnomAD %s %s: %s", source, chrom_vcf, exc)
                return None

        return handles.get(chrom_vcf)

    def fetch(
        self,
        chrom: str,
        pos: int,
        ref: str,
        alt: str,
    ) -> Optional[Dict]:
        """
        Query gnomAD for allele frequencies at chrom:pos REF>ALT.

        Tries exome first (higher sample size for coding variants),
        falls back to genome if not found.

        Parameters
        ----------
        chrom : str  bare chromosome ("1", "chr1", "X")
        pos   : int  1-based position
        ref   : str  reference allele
        alt   : str  alternate allele

        Returns
        -------
        dict with gnomAD fields, or None if not found.
        """
        if not self._available:
            return None

        chrom_vcf = _chrom_to_vcf(chrom)
        ref_u = ref.upper()
        alt_u = alt.upper()

        for source in ("exome", "genome"):
            handle = self._get_handle(chrom_vcf, source)
            if handle is None:
                continue

            try:
                rows = list(handle.fetch(chrom_vcf, pos - 1, pos))
            except (ValueError, KeyError):
                continue

            for row in rows:
                if isinstance(row, bytes):
                    row = row.decode("utf-8", errors="replace")
                parts = row.split("\t")
                if len(parts) < 8:
                    continue
                try:
                    row_pos = int(parts[1])
                except ValueError:
                    continue
                if row_pos != pos:
                    continue

                row_ref = parts[3].upper()
                vcf_alts = [a.upper() for a in parts[4].split(",")]

                # Exact match
                if row_ref == ref_u and alt_u in vcf_alts:
                    alt_idx = vcf_alts.index(alt_u)
                    info = _parse_info(parts[7])
                    return self._build_result(info, alt_idx, source)

                # Swapped alleles
                if row_ref == alt_u and ref_u in vcf_alts:
                    alt_idx = vcf_alts.index(ref_u)
                    info = _parse_info(parts[7])
                    return self._build_result(info, alt_idx, source)

        return None

    def _build_result(
        self,
        info: Dict[str, str],
        alt_idx: int,
        source: str,
    ) -> Dict:
        """Build result dict from parsed INFO field."""
        af = _extract_af(info, "AF", alt_idx)
        ac = _extract_int(info, "AC", alt_idx)
        an = _extract_int(info, "AN", 0)
        nhomalt = _extract_int(info, "nhomalt", alt_idx)

        return {
            "af_gnomad":       af,
            "af_gnomad_afr":   _extract_af(info, "AF_afr", alt_idx),
            "af_gnomad_amr":   _extract_af(info, "AF_amr", alt_idx),
            "af_gnomad_eas":   _extract_af(info, "AF_eas", alt_idx),
            "af_gnomad_nfe":   _extract_af(info, "AF_nfe", alt_idx),
            "af_gnomad_sas":   _extract_af(info, "AF_sas", alt_idx),
            "af_gnomad_mid":   _extract_af(info, "AF_mid", alt_idx),
            "an_gnomad":       an,
            "ac_gnomad":       ac,
            "nhomalt_gnomad":  nhomalt,
            "rarity_class_gnomad": classify_rarity(af),
            "gnomad_source":   source,  # "exome" or "genome"
        }

    def fetch_region(
        self,
        chrom: str,
        pos_min: int,
        pos_max: int,
    ) -> Dict[Tuple[int, str, str], Dict]:
        """
        Fetch all gnomAD variants in a genomic region with a single tabix query.

        Dramatically faster than N individual fetch() calls when annotating
        a gene or genomic window (one tabix seek vs. N seeks).

        Returns
        -------
        dict mapping (pos, ref, alt) → gnomAD result dict.
        Only includes variants found in gnomAD; missing ones map to None.

        Usage for gene browser enrichment:
            region = fetcher.fetch_region("13", 32315226, 32400081)
            af = region.get((pos, ref, alt), {}).get("af_gnomad")
        """
        if not self._available:
            return {}

        chrom_vcf = _chrom_to_vcf(chrom)
        hits: Dict[Tuple[int, str, str], Dict] = {}

        for source in ("exome", "genome"):
            handle = self._get_handle(chrom_vcf, source)
            if handle is None:
                continue
            try:
                rows = list(handle.fetch(chrom_vcf, pos_min - 1, pos_max))
            except (ValueError, KeyError):
                continue

            for row in rows:
                if isinstance(row, bytes):
                    row = row.decode("utf-8", errors="replace")
                parts = row.split("\t")
                if len(parts) < 8:
                    continue
                try:
                    row_pos = int(parts[1])
                except ValueError:
                    continue

                row_ref = parts[3].upper()
                vcf_alts = [a.upper() for a in parts[4].split(",")]
                info = _parse_info(parts[7])

                for alt_idx, alt_u in enumerate(vcf_alts):
                    key = (row_pos, row_ref, alt_u)
                    if key not in hits:
                        hits[key] = self._build_result(info, alt_idx, source)
                    # Also add swapped-allele key for IS_FLIP variants
                    key_flip = (row_pos, alt_u, row_ref)
                    if key_flip not in hits:
                        hits[key_flip] = self._build_result(info, alt_idx, source)

        return hits

    def fetch_batch(
        self,
        variants: List[Tuple[str, int, str, str]],
    ) -> Dict[Tuple, Optional[Dict]]:
        """
        Batch fetch for a list of (chrom, pos, ref, alt) tuples.

        When all variants are on the same chromosome and span a compact region,
        uses fetch_region internally for a single tabix query.
        Falls back to individual fetch() calls for multi-chromosome batches.

        Returns
        -------
        dict mapping (chrom, pos, ref, alt) → result dict or None.
        """
        if not variants:
            return {}

        # Group by chromosome
        by_chrom: Dict[str, List[Tuple[str, int, str, str]]] = {}
        for v in variants:
            by_chrom.setdefault(v[0], []).append(v)

        results: Dict[Tuple, Optional[Dict]] = {}

        for chrom, chrom_variants in by_chrom.items():
            positions = [v[1] for v in chrom_variants]
            pos_min, pos_max = min(positions), max(positions)
            # Use region fetch when range is compact (< 10 Mb) — avoids N tabix seeks
            if pos_max - pos_min < 10_000_000:
                region = self.fetch_region(chrom, pos_min, pos_max)
                for v in chrom_variants:
                    _, pos, ref, alt = v
                    key = (pos, ref.upper(), alt.upper())
                    results[v] = region.get(key)
            else:
                for v in chrom_variants:
                    results[v] = self.fetch(*v)

        return results

    def close(self) -> None:
        """Close all open TabixFile handles."""
        for h in self._exome_handles.values():
            try:
                h.close()
            except Exception:
                pass
        for h in self._genome_handles.values():
            try:
                h.close()
            except Exception:
                pass
        self._exome_handles.clear()
        self._genome_handles.clear()

    def __del__(self):
        try:
            self.close()
        except Exception:
            pass


# ── Module-level singleton (lazy init) ────────────────────────────────────────

_singleton: Optional[GnomadFetcher] = None


def get_gnomad_fetcher() -> GnomadFetcher:
    """Return module-level GnomadFetcher singleton (thread-unsafe, gunicorn pre-fork safe)."""
    global _singleton
    if _singleton is None:
        _singleton = GnomadFetcher()
    return _singleton
