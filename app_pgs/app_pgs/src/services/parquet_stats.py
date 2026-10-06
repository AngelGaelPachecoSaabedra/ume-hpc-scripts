"""
ParquetStats service
─────────────────────
Reads PRS results files (parquet preferred, TSV fallback) using DuckDB
and returns descriptive statistics + histogram data for the web dashboard.

Caching strategy (two layers):
  1. FileCache (/work/cache/) — file-based, TTL + mtime check, survives
     process restart within a gunicorn worker lifetime.
  2. DBPool (PostgreSQL stats_cache) — persistent across restarts, used
     as a secondary fallback when file cache is cold.

DuckDB usage:
  - New in-memory connection per request.
  - Only the score column is read (projection pushdown on parquet).
  - Paths are SQL-escaped before embedding (single-quote doubling).
  - Never writes to /data.

Column detection:
  NOT assumed to always be called "PRS".
  detect_score_column() priority:
    1. Exact match against canonical candidates.
    2. Starts with "PRS" (case-insensitive).
    3. Contains "score" (case-insensitive).
    4. Second non-ID column as last resort.
"""
import logging
from pathlib import Path
from typing import Optional

import duckdb

from config import Config
from services.local_catalog import validate_pgs_id

logger = logging.getLogger(__name__)

_SCORE_CANDIDATES = [
    "PRS_total", "PRS", "prs_total", "prs",
    "SCORE", "Score", "score", "PGS_score", "weighted_sum",
]


def detect_score_column(columns: list[str]) -> Optional[str]:
    col_set = set(columns)
    for c in _SCORE_CANDIDATES:
        if c in col_set:
            return c
    for c in columns:
        if c.upper().startswith("PRS"):
            return c
    for c in columns:
        if "score" in c.lower():
            return c
    id_like = {"sample_id", "id", "iid", "patid", "fid"}
    non_id = [c for c in columns if c.lower() not in id_like]
    if non_id:
        logger.warning("Score column guessed as %r from columns %s", non_id[0], columns)
        return non_id[0]
    return None


def _safe_path(path: Path) -> str:
    return str(path).replace("'", "''")


class ParquetStats:
    def __init__(self, cfg: Config, cache=None, db=None) -> None:
        self.cfg = cfg
        self._cache = cache   # FileCache | None
        self._db = db         # DBPool | None

    # ── Public API ───────────────────────────────────────────────────────────

    def get_stats(self, pgs_id: str) -> dict:
        """
        Return stats dict with n, mean, stddev, min, max, median,
        p05/p25/p75/p95, histogram, score_column, data_format.
        Returns {"error": "..."} on failure.
        """
        if not validate_pgs_id(pgs_id):
            return {"error": "Invalid PGS ID format", "pgs_id": pgs_id}

        path, fmt = self._find_data_file(pgs_id)
        if not path:
            return {
                "error": f"No results file found for {pgs_id} "
                         "(expected .parquet or _PRS_total.tsv)",
                "pgs_id": pgs_id,
            }

        # ── Cache lookup ──────────────────────────────────────────────────
        try:
            source_mtime = int(path.stat().st_mtime)
        except OSError:
            source_mtime = 0

        cached = self._cache_get(pgs_id, source_mtime)
        if cached:
            cached["_from_cache"] = True
            return cached

        # ── Compute via DuckDB ────────────────────────────────────────────
        columns = self._get_columns(path, fmt)
        if columns is None:
            return {"error": "Could not read column list from data file", "pgs_id": pgs_id}

        score_col = detect_score_column(columns)
        if not score_col:
            return {"error": "Could not detect score column", "pgs_id": pgs_id, "columns": columns}

        con = duckdb.connect(":memory:")
        try:
            con.execute(self._create_view_sql(path, fmt))

            row = con.execute(f"""
                SELECT
                    count(*)                                                     AS n,
                    avg("{score_col}")                                           AS mean,
                    stddev("{score_col}")                                        AS stddev,
                    min("{score_col}")                                           AS min,
                    max("{score_col}")                                           AS max,
                    median("{score_col}")                                        AS median,
                    percentile_cont(0.25) WITHIN GROUP (ORDER BY "{score_col}") AS p25,
                    percentile_cont(0.75) WITHIN GROUP (ORDER BY "{score_col}") AS p75,
                    percentile_cont(0.05) WITHIN GROUP (ORDER BY "{score_col}") AS p05,
                    percentile_cont(0.95) WITHIN GROUP (ORDER BY "{score_col}") AS p95
                FROM _data
                WHERE "{score_col}" IS NOT NULL
            """).fetchone()

            n, mean, std, mn, mx, med, p25, p75, p05, p95 = row
            hist = self._histogram(con, score_col, mn, mx)

            result = {
                "pgs_id": pgs_id,
                "score_column": score_col,
                "all_columns": columns,
                "data_format": fmt,
                "data_file": path.name,
                "n": int(n) if n is not None else None,
                "mean":   _r(mean),
                "stddev": _r(std),
                "min":    _r(mn),
                "max":    _r(mx),
                "median": _r(med),
                "p25":    _r(p25),
                "p75":    _r(p75),
                "p05":    _r(p05),
                "p95":    _r(p95),
                "histogram": hist,
                "_source_mtime": source_mtime,
                "_from_cache": False,
            }

            self._cache_set(pgs_id, result, source_mtime)
            return result

        except Exception as exc:
            logger.exception("DuckDB error for %s", pgs_id)
            return {"error": str(exc), "pgs_id": pgs_id}
        finally:
            con.close()

    # ── Cache helpers ─────────────────────────────────────────────────────────

    def _cache_get(self, pgs_id: str, source_mtime: int) -> Optional[dict]:
        # L1: file cache
        if self._cache:
            data = self._cache.get_with_mtime_check(
                self._cache.stats_key(pgs_id), source_mtime
            )
            if data:
                logger.debug("Stats file-cache hit: %s", pgs_id)
                return data
        # L2: DB cache
        if self._db and self._db.available:
            data = self._db.get_stats_cache(pgs_id, source_mtime)
            if data:
                logger.debug("Stats DB-cache hit: %s", pgs_id)
                # Warm file cache from DB
                if self._cache:
                    self._cache.set(self._cache.stats_key(pgs_id), data)
                return data
        return None

    def _cache_set(self, pgs_id: str, result: dict, source_mtime: int) -> None:
        if self._cache:
            self._cache.set(self._cache.stats_key(pgs_id), result)
        if self._db and self._db.available:
            self._db.upsert_stats_cache(pgs_id, result, source_mtime)

    # ── Internal helpers ──────────────────────────────────────────────────────

    def _find_data_file(self, pgs_id: str) -> tuple[Optional[Path], Optional[str]]:
        pgs_dir = self.cfg.scores_dir / pgs_id
        parquet = pgs_dir / f"{pgs_id}_PRS_total.parquet"
        tsv = pgs_dir / f"{pgs_id}_PRS_total.tsv"
        if parquet.exists():
            return parquet, "parquet"
        if tsv.exists():
            return tsv, "tsv"
        return None, None

    def _get_columns(self, path: Path, fmt: str) -> Optional[list[str]]:
        con = duckdb.connect(":memory:")
        try:
            con.execute(self._create_view_sql(path, fmt))
            return [r[0] for r in con.execute("DESCRIBE _data").fetchall()]
        except Exception as exc:
            logger.error("Cannot read columns from %s: %s", path, exc)
            return None
        finally:
            con.close()

    def _create_view_sql(self, path: Path, fmt: str) -> str:
        sp = _safe_path(path)
        if fmt == "parquet":
            return f"CREATE OR REPLACE VIEW _data AS SELECT * FROM read_parquet('{sp}')"
        return (
            f"CREATE OR REPLACE VIEW _data AS "
            f"SELECT * FROM read_csv_auto('{sp}', delim='\\t', header=true)"
        )

    def _histogram(
        self,
        con: duckdb.DuckDBPyConnection,
        score_col: str,
        min_v: Optional[float],
        max_v: Optional[float],
        bins: int = 50,
    ) -> dict:
        if min_v is None or max_v is None or min_v >= max_v:
            return {"edges": [], "counts": [], "bins": 0}
        try:
            width = (max_v - min_v) / bins
            rows = con.execute(f"""
                SELECT
                    CAST(floor(("{score_col}" - {min_v}) / {width}) AS INTEGER) AS bucket,
                    count(*) AS cnt
                FROM _data
                WHERE "{score_col}" IS NOT NULL
                GROUP BY bucket ORDER BY bucket
            """).fetchall()
            edges = [round(min_v + i * width, 10) for i in range(bins + 1)]
            counts = [0] * bins
            for bucket, cnt in rows:
                idx = int(min(max(bucket, 0), bins - 1))
                counts[idx] = int(cnt)
            return {"edges": edges, "counts": counts, "bins": bins}
        except Exception as exc:
            logger.error("Histogram error: %s", exc)
            return {"edges": [], "counts": [], "bins": 0, "error": str(exc)}


def _r(v, digits: int = 8) -> Optional[float]:
    return round(float(v), digits) if v is not None else None
