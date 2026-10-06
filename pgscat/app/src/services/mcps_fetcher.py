"""
mcps_fetcher.py — MCPS Population Allele Frequency Fetcher (v2)
================================================================
Fast per-variant MCPS allele frequency lookups using per-chromosome
Parquet files built via DuckDB's native CSV reader.

Build strategy (v2):
    DuckDB read_csv() on chrN.freq.tsv.gz  →  chrN.parquet (Snappy)
    ~1-3 s per chromosome, ~2 min total for all 23 chromosomes.
    Old Python gzip+executemany approach: ~5h per chromosome.

Source:
    /mnt/cephfs/hot_nvme/mcps/mcps-variant-browser-afs/chrN.freq.tsv.gz

Output (per-chromosome Parquet, incremental):
    {work_dir}/mcps_idx/parquet/chrN.parquet
    {work_dir}/mcps_idx/done/chrN.done          ← sentinel

Format (TSV, header):
    ID  CPRA  SOURCE  STON  GHOM  AHOM
    AN_RAW  AC_RAW  AF_RAW
    AN_AFR  AC_AFR  AF_AFR
    AN_EUR  AC_EUR  AF_EUR
    AN_MEX  AC_MEX  AF_MEX

Key: CPRA = "CHROM:POS:REF:ALT"  (bare chrom, no chr prefix)

IMPORTANT — Privacy:
    Aggregated allele frequencies only — no individual-level data.
    Classification: INTERNAL_AGGREGATED_ONLY
"""
from __future__ import annotations

import logging
import os
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)

_SOURCE_DIR = Path(os.environ.get("MCPS_DIR", "/mnt/cephfs/hot_nvme/mcps/mcps-variant-browser-afs"))
_CHROMS = [str(i) for i in range(1, 23)] + ["X"]


def _bare(chrom: str) -> str:
    """Strip 'chr' prefix if present."""
    return chrom[3:] if chrom.startswith("chr") else chrom


class McpsFetcher:
    """
    MCPS population frequency fetcher backed by per-chromosome Parquet files.

    Build once with build_index(); subsequent fetches are fast (<5ms each).
    """

    def __init__(self, work_dir: Path) -> None:
        self._source_dir = _SOURCE_DIR
        self._idx_dir    = work_dir / "mcps_idx"
        self._par_dir    = self._idx_dir / "parquet"
        self._done_dir   = self._idx_dir / "done"
        self._available: Optional[bool] = None
        self._conn = None       # persistent in-memory DuckDB for reads

    # ── Availability ──────────────────────────────────────────────────────────

    def is_available(self) -> bool:
        if self._available is None:
            self._available = self._source_dir.exists() and any(
                (self._source_dir / f"chr{c}.freq.tsv.gz").exists()
                for c in ["1", "X"]
            )
        return bool(self._available)

    def index_ready(self) -> bool:
        """True if at least one chromosome Parquet file has been built."""
        if not self._par_dir.exists():
            return False
        return any(self._par_dir.glob("chr*.parquet"))

    def built_chroms(self) -> List[str]:
        """Return list of chromosome names whose Parquet files exist."""
        if not self._par_dir.exists():
            return []
        return sorted(
            p.stem[3:]                          # "chr1" → "1"
            for p in self._par_dir.glob("chr*.parquet")
        )

    def parquet_path(self, chrom: str) -> Path:
        return self._par_dir / f"chr{_bare(chrom)}.parquet"

    def _done_path(self, chrom: str) -> Path:
        return self._done_dir / f"chr{_bare(chrom)}.done"

    # ── Read connection ───────────────────────────────────────────────────────

    def _get_conn(self):
        """Persistent in-memory DuckDB connection for read queries."""
        if self._conn is not None:
            try:
                self._conn.execute("SELECT 1")
                return self._conn
            except Exception:
                self._conn = None
        try:
            import duckdb
            self._conn = duckdb.connect()           # in-memory, read-only semantics
            logger.debug("MCPS: opened in-memory DuckDB connection")
        except Exception as exc:
            logger.warning("MCPS: could not open DuckDB: %s", exc)
            self._conn = None
        return self._conn

    # ── Build (fast path via DuckDB native CSV reader) ────────────────────────

    def build_index(
        self,
        chroms: Optional[List[str]] = None,
        force: bool = False,
    ) -> bool:
        """
        Build per-chromosome Parquet files from source TSV.GZ.

        Uses DuckDB's native gzip CSV reader — no Python-level decompression
        or row-by-row parsing.  ~1-3 s per chromosome.

        Parameters
        ----------
        chroms : chromosomes to build (default: all 1-22, X)
        force  : rebuild even if done sentinel exists
        """
        if not self.is_available():
            logger.error("MCPS source not found: %s", self._source_dir)
            return False

        self._par_dir.mkdir(parents=True, exist_ok=True)
        self._done_dir.mkdir(parents=True, exist_ok=True)

        target_chroms = [_bare(c) for c in (chroms or _CHROMS)]
        t_total = time.time()
        ok_all = True

        import duckdb

        for chrom in target_chroms:
            tsv_gz = self._source_dir / f"chr{chrom}.freq.tsv.gz"
            if not tsv_gz.exists():
                logger.warning("MCPS: source missing — skipping chr%s (%s)", chrom, tsv_gz)
                continue

            done = self._done_path(chrom)
            out_par = self.parquet_path(chrom)

            if done.exists() and out_par.exists() and not force:
                logger.info("MCPS: chr%s already built — skipping (use --force to rebuild)", chrom)
                continue

            src_mb = tsv_gz.stat().st_size / 1e6
            logger.info(
                "MCPS: building chr%s  source=%.0fMB  →  %s",
                chrom, src_mb, out_par,
            )
            t0 = time.time()

            # Remove stale output if forcing
            if force:
                out_par.unlink(missing_ok=True)
                done.unlink(missing_ok=True)

            try:
                # Each chromosome uses a fresh in-memory connection so we
                # never hold a write lock on any persistent file.
                conn = duckdb.connect()

                conn.execute(f"""
                    COPY (
                        SELECT
                            string_split(CPRA, ':')[1]                       AS chrom,
                            TRY_CAST(string_split(CPRA, ':')[2] AS INTEGER)  AS pos,
                            upper(string_split(CPRA, ':')[3])                AS ref,
                            upper(string_split(CPRA, ':')[4])                AS alt,
                            SOURCE                                           AS source,
                            TRY_CAST(AF_RAW AS DOUBLE)                       AS af_raw,
                            TRY_CAST(AF_AFR AS DOUBLE)                       AS af_afr,
                            TRY_CAST(AF_EUR AS DOUBLE)                       AS af_eur,
                            TRY_CAST(AF_MEX AS DOUBLE)                       AS af_mex,
                            TRY_CAST(AN_RAW AS DOUBLE)                       AS an_raw,
                            TRY_CAST(AC_RAW AS DOUBLE)                       AS ac_raw
                        FROM read_csv(
                            '{tsv_gz!s}',
                            delim        = '\t',
                            header       = true,
                            compression  = 'gzip',
                            all_varchar  = true,
                            ignore_errors = true
                        )
                        WHERE length(CPRA) > 0
                          AND len(string_split(CPRA, ':')) = 4
                    ) TO '{out_par!s}' (
                        FORMAT      PARQUET,
                        COMPRESSION SNAPPY,
                        ROW_GROUP_SIZE 100000
                    )
                """)

                conn.close()

                elapsed = time.time() - t0
                out_mb = out_par.stat().st_size / 1e6
                # Row count from parquet metadata
                n_rows = duckdb.connect().execute(
                    f"SELECT count(*) FROM read_parquet('{out_par!s}')"
                ).fetchone()[0]
                mb_s = src_mb / elapsed if elapsed > 0 else 0

                logger.info(
                    "MCPS: chr%s done — %s rows, %.0fMB parquet, %.1fs (%.0fMB/s source)",
                    chrom, f"{n_rows:,}", out_mb, elapsed, mb_s,
                )

                # Write done sentinel
                done.write_text(
                    f"chrom={chrom}\nrows={n_rows}\nelapsed={elapsed:.1f}\n"
                    f"parquet_mb={out_mb:.1f}\nbuild_time={time.strftime('%Y-%m-%dT%H:%M:%S')}\n"
                )

            except Exception as exc:
                logger.error("MCPS: chr%s FAILED — %s", chrom, exc)
                out_par.unlink(missing_ok=True)
                done.unlink(missing_ok=True)
                ok_all = False

        total_elapsed = time.time() - t_total
        built = self.built_chroms()
        logger.info(
            "MCPS: build complete — %d chromosomes ready, %.1fs total",
            len(built), total_elapsed,
        )
        return ok_all

    # ── Lookups ───────────────────────────────────────────────────────────────

    def fetch(
        self,
        chrom: str,
        pos: int,
        ref: str,
        alt: str,
    ) -> Optional[Dict[str, Any]]:
        """
        Lookup MCPS allele frequency for a single variant.

        Returns dict with af_mcps, af_mcps_afr, af_mcps_eur, af_mcps_mex,
        an_mcps, ac_mcps, mcps_source — or None if not found.
        """
        if not self.is_available():
            return None

        bare = _bare(chrom)
        par = self.parquet_path(bare)
        if not par.exists():
            return None

        conn = self._get_conn()
        if conn is None:
            return None

        try:
            row = conn.execute(
                f"""
                SELECT af_raw, af_afr, af_eur, af_mex, an_raw, ac_raw, source
                FROM read_parquet('{par!s}')
                WHERE chrom = ? AND pos = ? AND ref = ? AND alt = ?
                LIMIT 1
                """,
                [bare, int(pos), ref.upper(), alt.upper()],
            ).fetchone()

            if row:
                return {
                    "af_mcps":     row[0],
                    "af_mcps_afr": row[1],
                    "af_mcps_eur": row[2],
                    "af_mcps_mex": row[3],
                    "an_mcps":     row[4],
                    "ac_mcps":     row[5],
                    "mcps_source": row[6],
                }
        except Exception as exc:
            logger.warning("MCPS fetch failed: %s", exc)

        return None

    def fetch_batch(
        self,
        variants: List[Tuple[str, int, str, str]],
    ) -> Dict[Tuple, Optional[Dict]]:
        """
        Batch fetch for a list of (chrom, pos, ref, alt) tuples.

        Groups variants by chromosome to minimise parquet file opens.
        """
        if not variants:
            return {}

        results: Dict[Tuple, Optional[Dict]] = {v: None for v in variants}

        if not self.is_available():
            return results

        conn = self._get_conn()
        if conn is None:
            return results

        # Group by bare chromosome
        by_chrom: Dict[str, List[Tuple]] = {}
        for v in variants:
            b = _bare(v[0])
            by_chrom.setdefault(b, []).append(v)

        for bare, chrom_vars in by_chrom.items():
            par = self.parquet_path(bare)
            if not par.exists():
                continue

            try:
                placeholders = ",".join(
                    f"('{bare}', {int(v[1])}, '{v[2].upper()}', '{v[3].upper()}')"
                    for v in chrom_vars
                )
                rows = conn.execute(f"""
                    SELECT m.chrom, m.pos, m.ref, m.alt,
                           m.af_raw, m.af_afr, m.af_eur, m.af_mex,
                           m.an_raw, m.ac_raw, m.source
                    FROM read_parquet('{par!s}') AS m
                    JOIN (VALUES {placeholders}) AS q(chrom, pos, ref, alt)
                      ON m.chrom = q.chrom AND m.pos = q.pos
                     AND m.ref   = q.ref   AND m.alt = q.alt
                """).fetchall()

                for row in rows:
                    key_chrom, key_pos, key_ref, key_alt = row[0], row[1], row[2], row[3]
                    for orig in chrom_vars:
                        if (
                            _bare(orig[0]) == key_chrom
                            and orig[1] == key_pos
                            and orig[2].upper() == key_ref
                            and orig[3].upper() == key_alt
                        ):
                            results[orig] = {
                                "af_mcps":     row[4],
                                "af_mcps_afr": row[5],
                                "af_mcps_eur": row[6],
                                "af_mcps_mex": row[7],
                                "an_mcps":     row[8],
                                "ac_mcps":     row[9],
                                "mcps_source": row[10],
                            }
                            break

            except Exception as exc:
                logger.warning("MCPS batch fetch chr%s failed: %s", bare, exc)

        return results

    def status(self) -> Dict[str, Any]:
        """Return build status summary."""
        built = self.built_chroms()
        total = len(_CHROMS)
        done_info = {}
        for c in built:
            dp = self._done_path(c)
            if dp.exists():
                done_info[c] = dp.read_text().strip()
        return {
            "built_chroms":    built,
            "n_built":         len(built),
            "n_total":         total,
            "complete":        len(built) >= total,
            "parquet_dir":     str(self._par_dir),
        }

    def close(self) -> None:
        if self._conn:
            try:
                self._conn.close()
            except Exception:
                pass
            self._conn = None


# ── Module-level singleton ────────────────────────────────────────────────────

_singleton: Optional[McpsFetcher] = None


def get_mcps_fetcher(work_dir: Optional[Path] = None) -> McpsFetcher:
    global _singleton
    if _singleton is None:
        _work = work_dir or Path("/mnt/cephfs/hot/pgscat/work")
        _singleton = McpsFetcher(work_dir=_work)
    return _singleton
