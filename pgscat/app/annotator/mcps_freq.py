"""
mcps_freq.py — MCPS Population AF Annotator (Annotation Pipeline)
==================================================================
Reads from MCPS per-chromosome TSV.GZ source files, builds a lightweight
per-run DuckDB parquet index in the annotation output directory, and provides
per-variant lookups.

Used by annotate_variants.py when --mcps-dir is supplied.

Source layout:
    {source_dir}/chrN.freq.tsv.gz

Index output (written to outdir/mcps_idx/):
    chrN.parquet

TSV format (tab-separated, has header):
    ID  CPRA  SOURCE  STON  GHOM  AHOM
    AN_RAW  AC_RAW  AF_RAW
    AN_AFR  AC_AFR  AF_AFR
    AN_EUR  AC_EUR  AF_EUR
    AN_MEX  AC_MEX  AF_MEX

Key: CPRA = "CHROM:POS:REF:ALT"  (bare chrom, no chr prefix)
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Dict, List, Optional

logger = logging.getLogger(__name__)

_RARITY_THRESHOLDS = [
    ("common",        0.01),
    ("low_frequency", 0.001),
    ("rare",          0.0001),
]


def classify_rarity_mcps(af: Optional[float]) -> str:
    if af is None or af <= 0.0:
        return "novel"
    for label, thr in _RARITY_THRESHOLDS:
        if af >= thr:
            return label
    return "ultra_rare"


class McpsAnnotator:
    """
    MCPS population frequency annotator for the annotation pipeline.

    Workflow:
      1. Call build_for_chroms(chroms) with the chromosomes present in the betamap.
         This builds a parquet index for each chromosome (~1-3 s each) using DuckDB.
      2. Call fetch(chrom, pos, ref, alt) per variant.
    """

    def __init__(self, source_dir: str, index_dir: str) -> None:
        self._source_dir = Path(source_dir)
        self._index_dir  = Path(index_dir)
        self._index_dir.mkdir(parents=True, exist_ok=True)
        self._par_cache: Dict[str, Optional[Path]] = {}   # bare_chrom → parquet path or None
        self._conn = None   # persistent DuckDB connection for reads
        self.n_matched = 0
        self.n_novel   = 0

    def _open_conn(self):
        """Lazy open a persistent in-memory DuckDB connection."""
        if self._conn is not None:
            try:
                self._conn.execute("SELECT 1")
                return self._conn
            except Exception:
                self._conn = None
        try:
            import duckdb
            self._conn = duckdb.connect()
        except Exception as exc:
            logger.warning("MCPS: cannot open DuckDB: %s", exc)
        return self._conn

    def build_for_chroms(self, chroms: List[str]) -> None:
        """
        Build parquet index for the requested chromosomes.

        Only builds chromosomes not yet indexed.  Safe to call multiple times.
        Chroms may use bare ("1") or prefixed ("chr1") form.
        """
        try:
            import duckdb
        except ImportError:
            logger.error("MCPS: duckdb not available — cannot build index")
            return

        for chrom in chroms:
            bare = chrom.lstrip("chr")
            if bare in self._par_cache:
                continue   # already known

            par = self._index_dir / f"chr{bare}.parquet"

            if par.exists():
                logger.debug("MCPS: chr%s index already exists: %s", bare, par)
                self._par_cache[bare] = par
                continue

            tsv = self._source_dir / f"chr{bare}.freq.tsv.gz"
            if not tsv.exists():
                logger.warning("MCPS: chr%s source not found: %s — skipping", bare, tsv)
                self._par_cache[bare] = None
                continue

            logger.info("MCPS: building parquet index for chr%s …", bare)
            try:
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
                            TRY_CAST(AN_RAW AS BIGINT)                       AS an_raw,
                            TRY_CAST(AC_RAW AS BIGINT)                       AS ac_raw
                        FROM read_csv(
                            '{tsv!s}',
                            delim         = '\t',
                            header        = true,
                            compression   = 'gzip',
                            all_varchar   = true,
                            ignore_errors = true
                        )
                        WHERE length(CPRA) > 0
                          AND len(string_split(CPRA, ':')) = 4
                    ) TO '{par!s}' (
                        FORMAT       PARQUET,
                        COMPRESSION  SNAPPY,
                        ROW_GROUP_SIZE 100000
                    )
                """)
                conn.close()
                self._par_cache[bare] = par
                logger.info("MCPS: chr%s index ready", bare)
            except Exception as exc:
                logger.error("MCPS: chr%s build failed: %s", bare, exc)
                self._par_cache[bare] = None

    def fetch(
        self,
        chrom: str,
        pos: int,
        ref: str,
        alt: str,
    ) -> Optional[Dict]:
        """Look up MCPS allele frequencies for a single variant."""
        bare = chrom.lstrip("chr")
        par = self._par_cache.get(bare)
        if par is None:
            self.n_novel += 1
            return None

        conn = self._open_conn()
        if conn is None:
            self.n_novel += 1
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
                self.n_matched += 1
                return {
                    "af_mcps":          row[0],
                    "af_mcps_afr":      row[1],
                    "af_mcps_eur":      row[2],
                    "af_mcps_mex":      row[3],
                    "an_mcps":          int(row[4]) if row[4] is not None else None,
                    "ac_mcps":          int(row[5]) if row[5] is not None else None,
                    "mcps_source":      row[6],
                    "rarity_class_mcps": classify_rarity_mcps(row[0]),
                }
        except Exception as exc:
            logger.warning("MCPS fetch failed %s:%d: %s", chrom, pos, exc)

        self.n_novel += 1
        return None

    def log_summary(self) -> None:
        logger.info(
            "MCPS summary — matched: %d | novel: %d",
            self.n_matched, self.n_novel,
        )

    def close(self) -> None:
        if self._conn is not None:
            try:
                self._conn.close()
            except Exception:
                pass
            self._conn = None

    def __del__(self):
        try:
            self.close()
        except Exception:
            pass
