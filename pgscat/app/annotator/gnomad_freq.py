"""
gnomad_freq.py — gnomAD 4.1.1 Population Frequency Fetcher (Annotator)
========================================================================
Drop-in annotator-side companion to src/services/gnomad_fetcher.py.

Used inside Apptainer during the annotation pipeline to enrich each variant
with gnomAD 4.1.1 allele frequencies.

Provides richer population AFs than dbSNP GRAF-pop:
  - af_gnomad          global AF (exome + genome combined when both present)
  - af_gnomad_afr      African/African-American
  - af_gnomad_amr      Admixed American  ← best proxy for MEX/LATAM
  - af_gnomad_eas      East Asian
  - af_gnomad_nfe      Non-Finnish European
  - af_gnomad_sas      South Asian
  - af_gnomad_mid      Middle Eastern
  - an_gnomad          total allele number (quality indicator)
  - nhomalt_gnomad     homozygous count
  - rarity_class_gnomad  "common" | "low_frequency" | "rare" | "ultra_rare" | "novel"

Integration in annotate_variants.py:
    from gnomad_freq import GnomadAnnotator
    gnomad = GnomadAnnotator("--gnomad-exome-dir", "--gnomad-genome-dir")
    # per-variant:
    result = gnomad.fetch(chrom, pos, ref, alt)  → dict or None

Usage CLI flag:
    --gnomad-dir  /mnt/cephfs/hot_nvme/gnomad_4.1.1/variants
    (auto-discovers exome/ and genome/ subdirectories)
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)

_RARITY_THRESHOLDS = [
    ("common",        0.01),
    ("low_frequency", 0.001),
    ("rare",          0.0001),
]


def _classify_rarity(af: Optional[float]) -> str:
    if af is None or af <= 0.0:
        return "novel"
    for label, thr in _RARITY_THRESHOLDS:
        if af >= thr:
            return label
    return "ultra_rare"


def _to_chr(chrom: str) -> str:
    return chrom if chrom.startswith("chr") else f"chr{chrom}"


def _parse_info(info_str: str) -> Dict[str, str]:
    d: Dict[str, str] = {}
    for part in info_str.split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            d[k] = v
    return d


def _getf(info: Dict[str, str], key: str, idx: int = 0) -> Optional[float]:
    v = info.get(key)
    if v is None or v in (".", ""):
        return None
    parts = v.split(",")
    try:
        return float(parts[idx] if idx < len(parts) else parts[0])
    except (ValueError, IndexError):
        return None


def _geti(info: Dict[str, str], key: str, idx: int = 0) -> Optional[int]:
    v = info.get(key)
    if v is None or v in (".", ""):
        return None
    parts = v.split(",")
    try:
        return int(parts[idx] if idx < len(parts) else parts[0])
    except (ValueError, IndexError):
        return None


class GnomadAnnotator:
    """
    Annotator-side gnomAD 4.1.1 frequency fetcher.

    Fast path (preferred): per-chromosome parquet index queried via DuckDB.
      - Build index once with gnomad_index.py; point GNOMAD_PARQUET_DIR at it.
      - One DuckDB query per chromosome batch → eliminates CephFS seek latency.

    Slow fallback: per-variant pysam tabix + range queries.
      - Active when parquet_dir is not set or parquet file is absent.
      - Acceptable for small PGS (<100 variants); very slow for large sets.
    """

    def __init__(
        self,
        exome_dir: Optional[str] = None,
        genome_dir: Optional[str] = None,
        gnomad_dir: Optional[str] = None,
        parquet_dir: Optional[str] = None,
    ) -> None:
        import pysam  # noqa: F401 — ensure pysam is available

        if gnomad_dir:
            base = Path(gnomad_dir)
            self._exome_dir = base / "exome"
            self._genome_dir = base / "genome"
        else:
            self._exome_dir = Path(exome_dir) if exome_dir else None
            self._genome_dir = Path(genome_dir) if genome_dir else None

        self._parquet_dir: Optional[Path] = Path(parquet_dir) if parquet_dir else None
        self._duckdb_conn = None   # lazy DuckDB connection for parquet queries

        self._handles_e: Dict[str, object] = {}
        self._handles_g: Dict[str, object] = {}
        self.n_matched: int = 0
        self.n_novel:   int = 0

        if self._parquet_dir and self._parquet_dir.exists():
            logger.info(
                "GnomadAnnotator: parquet index at %s [fast path]", self._parquet_dir
            )
        else:
            logger.info(
                "GnomadAnnotator: tabix fallback — exome=%s  genome=%s",
                self._exome_dir, self._genome_dir,
            )

    # ── Parquet / DuckDB fast path ────────────────────────────────────────────

    def _get_duckdb(self):
        """Lazy persistent in-memory DuckDB connection for parquet reads."""
        if self._duckdb_conn is not None:
            try:
                self._duckdb_conn.execute("SELECT 1")
                return self._duckdb_conn
            except Exception:
                self._duckdb_conn = None
        try:
            import duckdb
            self._duckdb_conn = duckdb.connect()
        except Exception as exc:
            logger.warning("gnomAD DuckDB unavailable: %s", exc)
        return self._duckdb_conn

    def has_parquet(self, chrom: str) -> bool:
        """True if at least one parquet file exists for this chromosome."""
        if not self._parquet_dir:
            return False
        bare = chrom.lstrip("chr")
        return (
            (self._parquet_dir / f"chr{bare}_exome.parquet").exists()
            or (self._parquet_dir / f"chr{bare}_genome.parquet").exists()
        )

    def fetch_batch_parquet(
        self,
        chrom: str,
        variants: List[Tuple[int, str, str]],   # (pos, ref, alt)
    ) -> Dict[Tuple[int, str, str], Dict]:
        """
        Batch lookup for all (pos, ref, alt) on one chromosome via DuckDB.

        Queries exome parquet first; genome fills missing positions.
        Returns dict keyed by (pos, ref, alt) → result dict.
        """
        if not variants:
            return {}
        conn = self._get_duckdb()
        if conn is None:
            return {}

        bare = chrom.lstrip("chr")
        exome_par  = self._parquet_dir / f"chr{bare}_exome.parquet"
        genome_par = self._parquet_dir / f"chr{bare}_genome.parquet"

        results: Dict[Tuple[int, str, str], Dict] = {}

        for which, par in [("exome", exome_par), ("genome", genome_par)]:
            if not par.exists():
                continue
            # Only look up variants not yet matched
            remaining = [v for v in variants if v not in results]
            # Also try swapped alleles (handle ref/alt orientation mismatches)
            remaining_swapped = [(p, a, r) for p, r, a in remaining if (p, a, r) not in results]
            lookups = list({*remaining, *remaining_swapped})
            if not lookups:
                break

            values_sql = ", ".join(f"({p}, '{r}', '{a}')" for p, r, a in lookups)
            try:
                rows = conn.execute(f"""
                    SELECT g.pos, g.ref, g.alt,
                           g.af_raw, g.af_afr, g.af_amr, g.af_eas,
                           g.af_nfe, g.af_sas, g.af_mid, g.an, g.nhomalt
                    FROM read_parquet('{par!s}') g
                    JOIN (VALUES {values_sql}) AS v(pos, ref, alt)
                      ON g.pos = v.pos AND g.ref = v.ref AND g.alt = v.alt
                """).fetchall()
            except Exception as exc:
                logger.warning("gnomAD DuckDB query chr%s %s: %s", bare, which, exc)
                continue

            for row in rows:
                pos, ref, alt = row[0], row[1], row[2]
                af = row[3]
                entry = {
                    "af_gnomad":           af,
                    "af_gnomad_afr":       row[4],
                    "af_gnomad_amr":       row[5],
                    "af_gnomad_eas":       row[6],
                    "af_gnomad_nfe":       row[7],
                    "af_gnomad_sas":       row[8],
                    "af_gnomad_mid":       row[9],
                    "an_gnomad":           row[10],
                    "nhomalt_gnomad":      row[11],
                    "rarity_class_gnomad": _classify_rarity(af),
                    "gnomad_source":       which,
                }
                # Store under canonical key AND swapped key so caller lookup works
                results[(pos, ref, alt)] = entry
                results[(pos, alt, ref)] = entry   # also register swapped

        return results

    # ── Tabix fallback ────────────────────────────────────────────────────────

    def _open(self, chrom_vcf: str, which: str) -> Optional[object]:
        import pysam

        handles = self._handles_e if which == "exome" else self._handles_g
        src_dir = self._exome_dir if which == "exome" else self._genome_dir

        if chrom_vcf in handles:
            return handles[chrom_vcf]
        if src_dir is None:
            return None

        pfx = "exomes" if which == "exome" else "genomes"
        vcf = src_dir / f"gnomad.{pfx}.v4.1.1.sites.{chrom_vcf}.vcf.bgz"
        tbi = Path(str(vcf) + ".tbi")
        if not vcf.exists() or not tbi.exists():
            logger.debug("gnomAD: missing %s", vcf)
            handles[chrom_vcf] = None
            return None

        try:
            h = pysam.TabixFile(str(vcf))
            handles[chrom_vcf] = h
            logger.debug("gnomAD opened: %s %s", which, chrom_vcf)
            return h
        except Exception as exc:
            logger.warning("gnomAD open failed %s %s: %s", which, chrom_vcf, exc)
            handles[chrom_vcf] = None
            return None

    def fetch_range(
        self,
        chrom: str,
        start: int,
        end: int,
    ) -> Dict[Tuple[int, str, str], Dict]:
        """
        Fetch all gnomAD records overlapping [start, end] (1-based, inclusive).

        Returns a dict keyed by (pos, ref, alt) → result dict.
        Exome records take priority; genome fills positions absent from exome.

        Note: htslib "index file is older than data file" warnings are non-fatal
        and indicate only that the .tbi predates the .vcf.bgz by mtime.
        """
        chrom_vcf = _to_chr(chrom.replace("chr", "").strip())
        results: Dict[Tuple[int, str, str], Dict] = {}

        for which in ("exome", "genome"):
            handle = self._open(chrom_vcf, which)
            if handle is None:
                continue
            try:
                rows = list(handle.fetch(chrom_vcf, start - 1, end))
            except (ValueError, KeyError):
                continue
            for row in rows:
                if isinstance(row, bytes):
                    row = row.decode("utf-8", errors="replace")
                parts = row.split("\t")
                if len(parts) < 8:
                    continue
                try:
                    pos = int(parts[1])
                except ValueError:
                    continue
                row_ref  = parts[3].upper()
                vcf_alts = [a.upper() for a in parts[4].split(",")]
                info     = _parse_info(parts[7])
                for match_idx, vcf_alt in enumerate(vcf_alts):
                    key = (pos, row_ref, vcf_alt)
                    if key not in results:   # exome takes priority over genome
                        af = _getf(info, "AF", match_idx)
                        results[key] = {
                            "af_gnomad":           af,
                            "af_gnomad_afr":       _getf(info, "AF_afr", match_idx),
                            "af_gnomad_amr":       _getf(info, "AF_amr", match_idx),
                            "af_gnomad_eas":       _getf(info, "AF_eas", match_idx),
                            "af_gnomad_nfe":       _getf(info, "AF_nfe", match_idx),
                            "af_gnomad_sas":       _getf(info, "AF_sas", match_idx),
                            "af_gnomad_mid":       _getf(info, "AF_mid", match_idx),
                            "an_gnomad":           _geti(info, "AN", 0),
                            "nhomalt_gnomad":      _geti(info, "nhomalt", match_idx),
                            "rarity_class_gnomad": _classify_rarity(af),
                            "gnomad_source":       which,
                        }
        return results

    def fetch(
        self,
        chrom: str,
        pos: int,
        ref: str,
        alt: str,
    ) -> Optional[Dict]:
        chrom_vcf = _to_chr(chrom.replace("chr", "").strip())
        ref_u = ref.upper()
        alt_u = alt.upper()

        for which in ("exome", "genome"):
            handle = self._open(chrom_vcf, which)
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
                    if int(parts[1]) != pos:
                        continue
                except ValueError:
                    continue

                row_ref = parts[3].upper()
                vcf_alts = [a.upper() for a in parts[4].split(",")]

                match_idx: Optional[int] = None
                if row_ref == ref_u and alt_u in vcf_alts:
                    match_idx = vcf_alts.index(alt_u)
                elif row_ref == alt_u and ref_u in vcf_alts:
                    match_idx = vcf_alts.index(ref_u)

                if match_idx is not None:
                    info = _parse_info(parts[7])
                    af = _getf(info, "AF", match_idx)
                    result = {
                        "af_gnomad":          af,
                        "af_gnomad_afr":      _getf(info, "AF_afr", match_idx),
                        "af_gnomad_amr":      _getf(info, "AF_amr", match_idx),
                        "af_gnomad_eas":      _getf(info, "AF_eas", match_idx),
                        "af_gnomad_nfe":      _getf(info, "AF_nfe", match_idx),
                        "af_gnomad_sas":      _getf(info, "AF_sas", match_idx),
                        "af_gnomad_mid":      _getf(info, "AF_mid", match_idx),
                        "an_gnomad":          _geti(info, "AN", 0),
                        "nhomalt_gnomad":     _geti(info, "nhomalt", match_idx),
                        "rarity_class_gnomad": _classify_rarity(af),
                        "gnomad_source":      which,
                    }
                    self.n_matched += 1
                    return result

        self.n_novel += 1
        return None

    def log_summary(self) -> None:
        logger.info(
            "gnomAD summary — matched: %d | novel: %d",
            self.n_matched, self.n_novel,
        )

    def close(self) -> None:
        for d in (self._handles_e, self._handles_g):
            for h in d.values():
                if h is not None:
                    try:
                        h.close()
                    except Exception:
                        pass
            d.clear()

    def __del__(self):
        try:
            self.close()
        except Exception:
            pass
