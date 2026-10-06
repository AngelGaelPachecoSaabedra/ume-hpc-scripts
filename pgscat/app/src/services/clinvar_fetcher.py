"""
clinvar_fetcher.py — ClinVar GRCh38 Variant Fetcher
=====================================================
Provides per-variant ClinVar lookups (CLNSIG, CLNDN, CLNREVSTAT, ALLELEID)
from the local ClinVar VCF (GRCh38, April 2026).

Strategy: Build a derived Parquet index once (at first use) in the work
directory, then use DuckDB for O(log n) lookups.

Source:
    /mnt/cephfs/hot_nvme/clinvar/GRCh38/clinvar.vcf.gz  (BGZF, no .tbi)

Derived index:
    /mnt/cephfs/hot/pgscat/work/clinvar_idx/clinvar_grch38.parquet
    (~50-150MB, built once, rebuilt if source is newer)

ClinVar chromosome format: bare numbers (1, 2, ..., X, Y)

Usage:
    fetcher = ClinVarFetcher(work_dir=cfg.WORK_DIR)
    result = fetcher.fetch("1", 66926, "AG", "A")
    # {"clnsig": "Pathogenic", "clndn": "Retinitis_pigmentosa", ...}
"""
from __future__ import annotations

import logging
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)

_SOURCE_VCF = Path("/mnt/cephfs/hot_nvme/clinvar/GRCh38/clinvar.vcf.gz")

# INFO fields to extract
_INFO_KEYS = ["CLNSIG", "CLNDN", "CLNREVSTAT", "ALLELEID", "CLNVC", "RS"]


def _bare_chrom(chrom: str) -> str:
    return chrom[3:] if chrom.startswith("chr") else chrom


class ClinVarFetcher:
    """
    ClinVar variant lookup backed by a derived DuckDB/Parquet index.

    The index is built on first access from the source BGZ file via pysam.
    Subsequent lookups use DuckDB for sub-millisecond queries.
    """

    def __init__(self, work_dir: Path) -> None:
        self._source = _SOURCE_VCF
        self._idx_dir = work_dir / "clinvar_idx"
        self._parquet = self._idx_dir / "clinvar_grch38.parquet"
        self._available: Optional[bool] = None
        self._conn = None   # duckdb connection (lazy)

    def is_available(self) -> bool:
        if self._available is None:
            self._available = self._source.exists()
        return bool(self._available)

    def _ensure_index(self) -> bool:
        """Build Parquet index if missing or stale. Returns True if ready."""
        if not self.is_available():
            return False

        self._idx_dir.mkdir(parents=True, exist_ok=True)

        # Check if index is up-to-date
        if self._parquet.exists():
            src_mtime = self._source.stat().st_mtime
            idx_mtime = self._parquet.stat().st_mtime
            if idx_mtime >= src_mtime:
                logger.debug("ClinVar index up-to-date: %s", self._parquet)
                return True
            logger.info("ClinVar index stale, rebuilding...")

        return self._build_index()

    def _build_index(self) -> bool:
        """Stream through ClinVar BGZ and build Parquet index."""
        try:
            import pysam
            import pyarrow as pa
            import pyarrow.parquet as pq
        except ImportError as e:
            logger.warning("ClinVar index build requires pysam+pyarrow: %s", e)
            return False

        logger.info("Building ClinVar GRCh38 index from %s …", self._source)
        t0 = time.time()

        schema = pa.schema([
            ("chrom", pa.string()),
            ("pos",   pa.int32()),
            ("ref",   pa.string()),
            ("alt",   pa.string()),
            ("clnsig",      pa.string()),
            ("clndn",       pa.string()),
            ("clnrevstat",  pa.string()),
            ("alleleid",    pa.int32()),
            ("clnvc",       pa.string()),
            ("rsid",        pa.string()),
        ])

        tmp_path = self._parquet.with_suffix(".parquet.tmp")
        writer = None
        batch_rows: list = []
        BATCH_SIZE = 50_000
        n_written = 0

        try:
            f = pysam.BGZFile(str(self._source), "r")
            writer = pq.ParquetWriter(str(tmp_path), schema, compression="snappy")

            for raw_line in f:
                if isinstance(raw_line, bytes):
                    line = raw_line.decode("utf-8", errors="replace")
                else:
                    line = raw_line

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
                ref = parts[3]
                alt = parts[4]
                info_str = parts[7]

                # Parse INFO
                info: Dict[str, str] = {}
                for kv in info_str.split(";"):
                    if "=" in kv:
                        k, v = kv.split("=", 1)
                        info[k] = v

                clnsig     = info.get("CLNSIG", "")
                clndn      = info.get("CLNDN", "")
                clnrevstat = info.get("CLNREVSTAT", "")
                alleleid_s = info.get("ALLELEID", "")
                clnvc      = info.get("CLNVC", "")
                rs         = info.get("RS", "")

                try:
                    alleleid = int(alleleid_s)
                except (ValueError, TypeError):
                    alleleid = -1

                batch_rows.append({
                    "chrom": chrom,
                    "pos":   pos,
                    "ref":   ref,
                    "alt":   alt,
                    "clnsig":     clnsig or None,
                    "clndn":      clndn or None,
                    "clnrevstat": clnrevstat or None,
                    "alleleid":   alleleid,
                    "clnvc":      clnvc or None,
                    "rsid":       rs or None,
                })

                if len(batch_rows) >= BATCH_SIZE:
                    batch = pa.RecordBatch.from_pylist(batch_rows, schema=schema)
                    writer.write_batch(batch)
                    n_written += len(batch_rows)
                    batch_rows.clear()

            f.close()

            if batch_rows:
                batch = pa.RecordBatch.from_pylist(batch_rows, schema=schema)
                writer.write_batch(batch)
                n_written += len(batch_rows)

            writer.close()
            writer = None

            tmp_path.rename(self._parquet)
            elapsed = time.time() - t0
            logger.info(
                "ClinVar index built: %d variants, %.1fs → %s",
                n_written, elapsed, self._parquet,
            )
            return True

        except Exception as exc:
            logger.error("ClinVar index build failed: %s", exc)
            if writer:
                try:
                    writer.close()
                except Exception:
                    pass
            tmp_path.unlink(missing_ok=True)
            return False

    def _get_conn(self):
        """Return (or create) a DuckDB connection to the index."""
        if self._conn is not None:
            return self._conn
        try:
            import duckdb
            self._conn = duckdb.connect(database=":memory:")
            # Register the parquet as a view
            self._conn.execute(
                f"CREATE VIEW clinvar AS SELECT * FROM read_parquet('{self._parquet}')"
            )
            logger.info("ClinVar DuckDB view ready: %s", self._parquet)
        except Exception as exc:
            logger.warning("ClinVar DuckDB init failed: %s", exc)
            self._conn = None
        return self._conn

    def fetch(
        self,
        chrom: str,
        pos: int,
        ref: str,
        alt: str,
    ) -> Optional[Dict[str, Any]]:
        """
        Lookup ClinVar annotation for a variant.

        Parameters
        ----------
        chrom : str  chromosome (bare or chr-prefixed)
        pos   : int  1-based genomic position
        ref   : str  reference allele
        alt   : str  alternate allele

        Returns
        -------
        dict with clnsig, clndn, clnrevstat, alleleid, clnvc, rsid
        or None if not found.
        """
        if not self.is_available():
            return None

        if not self._parquet.exists():
            ok = self._ensure_index()
            if not ok:
                return None

        conn = self._get_conn()
        if conn is None:
            # Fall back to streaming if DuckDB unavailable
            return self._fetch_streaming(chrom, pos, ref, alt)

        bare = _bare_chrom(chrom)
        ref_u = ref.upper()
        alt_u = alt.upper()

        try:
            row = conn.execute(
                """
                SELECT clnsig, clndn, clnrevstat, alleleid, clnvc, rsid
                FROM clinvar
                WHERE chrom = ? AND pos = ?
                  AND upper(ref) = ? AND upper(alt) = ?
                LIMIT 1
                """,
                [bare, pos, ref_u, alt_u],
            ).fetchone()

            if row is None:
                # Try with chrom prefixed (sometimes ClinVar has "chr1")
                row = conn.execute(
                    """
                    SELECT clnsig, clndn, clnrevstat, alleleid, clnvc, rsid
                    FROM clinvar
                    WHERE chrom = ? AND pos = ?
                      AND upper(ref) = ? AND upper(alt) = ?
                    LIMIT 1
                    """,
                    [f"chr{bare}", pos, ref_u, alt_u],
                ).fetchone()

            if row:
                return {
                    "clnsig":     _clean(row[0]),
                    "clndn":      _clean(row[1]),
                    "clnrevstat": _clean(row[2]),
                    "alleleid":   row[3] if row[3] and row[3] >= 0 else None,
                    "clnvc":      _clean(row[4]),
                    "clinvar_rsid": _clean(row[5]),
                }
        except Exception as exc:
            logger.warning("ClinVar DuckDB query failed: %s", exc)

        return None

    def fetch_batch(
        self,
        variants: List[Tuple[str, int, str, str]],
    ) -> Dict[Tuple, Optional[Dict]]:
        """
        Batch ClinVar lookup using a single DuckDB VALUES join.

        ~50× faster than N individual fetch() calls for gene-sized batches.

        Returns
        -------
        dict mapping (chrom, pos, ref, alt) → result dict or None.
        """
        if not variants:
            return {}

        results: Dict[Tuple, Optional[Dict]] = {v: None for v in variants}

        if not self.is_available() or not self._parquet.exists():
            return results

        conn = self._get_conn()
        if conn is None:
            # Fallback: serial fetch
            for v in variants:
                results[v] = self.fetch(*v)
            return results

        try:
            placeholders = ", ".join(
                f"('{_bare_chrom(v[0])}', {int(v[1])}, '{v[2].upper()}', '{v[3].upper()}')"
                for v in variants
            )
            rows = conn.execute(f"""
                SELECT q.chrom, q.pos, q.ref, q.alt,
                       c.clnsig, c.clndn, c.clnrevstat, c.alleleid, c.clnvc, c.rsid
                FROM (VALUES {placeholders}) AS q(chrom, pos, ref, alt)
                LEFT JOIN clinvar AS c
                  ON (c.chrom = q.chrom OR c.chrom = 'chr' || q.chrom)
                 AND c.pos = q.pos
                 AND upper(c.ref) = q.ref
                 AND upper(c.alt) = q.alt
            """).fetchall()

            for row in rows:
                q_chrom, q_pos, q_ref, q_alt = row[0], row[1], row[2], row[3]
                if row[4] is None:
                    continue  # no ClinVar hit
                hit = {
                    "clnsig":       _clean(row[4]),
                    "clndn":        _clean(row[5]),
                    "clnrevstat":   _clean(row[6]),
                    "alleleid":     row[7] if row[7] and row[7] >= 0 else None,
                    "clnvc":        _clean(row[8]),
                    "clinvar_rsid": _clean(row[9]),
                }
                # Match back to original variant tuple
                for v in variants:
                    if (
                        _bare_chrom(v[0]) == q_chrom
                        and int(v[1]) == q_pos
                        and v[2].upper() == q_ref
                        and v[3].upper() == q_alt
                    ):
                        results[v] = hit

        except Exception as exc:
            logger.warning("ClinVar batch fetch failed: %s — falling back to serial", exc)
            for v in variants:
                results[v] = self.fetch(*v)

        return results

    def _fetch_streaming(
        self,
        chrom: str,
        pos: int,
        ref: str,
        alt: str,
    ) -> Optional[Dict[str, Any]]:
        """Last-resort: stream through BGZ looking for the variant (slow!)."""
        try:
            import pysam
            bare = _bare_chrom(chrom)
            ref_u = ref.upper()
            alt_u = alt.upper()
            f = pysam.BGZFile(str(self._source), "r")
            for raw_line in f:
                if isinstance(raw_line, bytes):
                    line = raw_line.decode("utf-8", errors="replace")
                else:
                    line = raw_line
                if line.startswith("#"):
                    continue
                parts = line.split("\t", 8)
                if len(parts) < 8:
                    continue
                try:
                    lpos = int(parts[1])
                except ValueError:
                    continue
                if lpos != pos:
                    continue
                if parts[0] not in (bare, f"chr{bare}"):
                    continue
                if parts[3].upper() == ref_u and parts[4].upper() == alt_u:
                    info: Dict[str, str] = {}
                    for kv in parts[7].split(";"):
                        if "=" in kv:
                            k, v = kv.split("=", 1)
                            info[k] = v
                    f.close()
                    return {
                        "clnsig":       _clean(info.get("CLNSIG")),
                        "clndn":        _clean(info.get("CLNDN")),
                        "clnrevstat":   _clean(info.get("CLNREVSTAT")),
                        "alleleid":     _safe_int(info.get("ALLELEID")),
                        "clnvc":        _clean(info.get("CLNVC")),
                        "clinvar_rsid": _clean(info.get("RS")),
                    }
            f.close()
        except Exception as exc:
            logger.warning("ClinVar streaming fetch failed: %s", exc)
        return None

    def build_index_if_needed(self) -> bool:
        """Public method to trigger index build (e.g. at startup)."""
        return self._ensure_index()

    def close(self) -> None:
        if self._conn:
            try:
                self._conn.close()
            except Exception:
                pass
            self._conn = None


def _clean(v: Optional[str]) -> Optional[str]:
    if v is None or v in (".", ""):
        return None
    return v.replace("_", " ").replace("|", "; ")


def _safe_int(v: Optional[str]) -> Optional[int]:
    try:
        return int(v)
    except (TypeError, ValueError):
        return None


# ── Module-level singleton ────────────────────────────────────────────────────

_singleton: Optional[ClinVarFetcher] = None


def get_clinvar_fetcher(work_dir: Optional[Path] = None) -> ClinVarFetcher:
    global _singleton
    if _singleton is None:
        from pathlib import Path as _Path
        _work = work_dir or _Path("/mnt/cephfs/hot/pgscat/work")
        _singleton = ClinVarFetcher(work_dir=_work)
    return _singleton
