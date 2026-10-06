"""
Database connection pool – PostgreSQL via psycopg2
════════════════════════════════════════════════════

Design goals:
  • Graceful degradation: if PostgreSQL is unavailable, all operations
    return safe empty/None values; the app continues in filesystem-only mode.
  • Thin wrapper: no ORM, no migrations in Python — schema lives in sql/schema.sql.
  • Thread-safe pool: one ThreadedConnectionPool per worker process (gunicorn
    pre-fork model means each worker has its own pool, which is correct).
  • Password loaded from file if APP_DB_PASSWORD_FILE exists; env var as fallback.

Usage:
    db = DBPool(cfg)
    if db.available:
        row = db.fetchone("SELECT * FROM local_pgs_index WHERE pgs_id = %s", (pid,))
    # Or with cursor context manager:
    try:
        with db.cursor() as cur:
            cur.execute("INSERT ...", params)
    except DBNotAvailable:
        pass  # degraded mode
"""
import logging
import threading
from contextlib import contextmanager
from typing import Any, Generator, Optional

logger = logging.getLogger(__name__)

# ── Import guard ──────────────────────────────────────────────────────────────
try:
    import psycopg2
    import psycopg2.extras
    import psycopg2.pool
    _PSYCOPG2_OK = True
except ImportError:
    psycopg2 = None  # type: ignore
    _PSYCOPG2_OK = False
    logger.warning(
        "psycopg2 not installed (pip install psycopg2-binary). "
        "Running in filesystem-only mode."
    )


class DBNotAvailable(Exception):
    """Raised when DB pool is not available (degraded mode)."""


class DBPool:
    """
    Thread-safe PostgreSQL connection pool with graceful degradation.

    Attributes:
        available  – True when the pool is connected and operational.
    """

    def __init__(self, cfg) -> None:
        self._cfg = cfg
        self._pool = None
        self._available = False
        self._lock = threading.Lock()
        self._connect_attempts = 0

        if not _PSYCOPG2_OK:
            return
        self._try_connect()

    # ── Connection management ─────────────────────────────────────────────────

    def _try_connect(self) -> bool:
        """Attempt to create the pool. Returns True on success."""
        self._connect_attempts += 1
        try:
            dsn = self._build_dsn()
            pool = psycopg2.pool.ThreadedConnectionPool(
                minconn=1,
                maxconn=5,
                dsn=dsn,
            )
            # Validate with a quick test
            conn = pool.getconn()
            conn.cursor().execute("SELECT 1")
            conn.rollback()
            pool.putconn(conn)

            with self._lock:
                if self._pool:
                    self._pool.closeall()
                self._pool = pool
                self._available = True

            logger.info(
                "PostgreSQL connected: %s:%s/%s (attempt %d)",
                self._cfg.DB_HOST, self._cfg.DB_PORT,
                self._cfg.DB_NAME, self._connect_attempts,
            )
            return True

        except Exception as exc:
            logger.warning(
                "PostgreSQL not available (attempt %d): %s",
                self._connect_attempts, exc,
            )
            self._available = False
            return False

    def retry_connect(self) -> bool:
        """
        Try to reconnect if currently unavailable.
        Safe to call from any thread or health check.
        """
        if self._available:
            return True
        if not _PSYCOPG2_OK:
            return False
        return self._try_connect()

    def _build_dsn(self) -> str:
        cfg = self._cfg
        password = cfg.db_password  # property on Config handles file vs env
        return (
            f"host={cfg.DB_HOST} "
            f"port={cfg.DB_PORT} "
            f"dbname={cfg.DB_NAME} "
            f"user={cfg.DB_USER} "
            f"password={password} "
            f"connect_timeout=5 "
            f"application_name=pgs-dashboard"
        )

    @property
    def available(self) -> bool:
        return self._available and self._pool is not None

    def close(self) -> None:
        with self._lock:
            if self._pool:
                try:
                    self._pool.closeall()
                except Exception:
                    pass
                self._pool = None
                self._available = False

    # ── Context manager ───────────────────────────────────────────────────────

    @contextmanager
    def cursor(
        self, dict_cursor: bool = False
    ) -> Generator[Any, None, None]:
        """
        Yield a cursor within a transaction that auto-commits on exit.
        Raises DBNotAvailable if pool is down.
        Rolls back automatically on exception.
        """
        if not self.available:
            raise DBNotAvailable("PostgreSQL pool is not available")

        conn = None
        try:
            conn = self._pool.getconn()  # type: ignore[union-attr]
            factory = (
                psycopg2.extras.RealDictCursor  # type: ignore[attr-defined]
                if dict_cursor else None
            )
            kwargs = {"cursor_factory": factory} if factory else {}
            with conn.cursor(**kwargs) as cur:
                yield cur
            conn.commit()
        except DBNotAvailable:
            raise
        except Exception:
            if conn:
                try:
                    conn.rollback()
                except Exception:
                    pass
            raise
        finally:
            if conn and self._pool:
                try:
                    self._pool.putconn(conn)
                except Exception:
                    pass

    # ── Convenience helpers ───────────────────────────────────────────────────

    def execute(self, sql: str, params=None) -> bool:
        """
        Execute a DML statement.
        Returns True on success, False if DB unavailable or error.
        Never raises.
        """
        try:
            with self.cursor() as cur:
                cur.execute(sql, params)
            return True
        except DBNotAvailable:
            return False
        except Exception as exc:
            logger.error("DB execute error: %s | sql=%s", exc, sql[:120])
            return False

    def fetchone(self, sql: str, params=None) -> Optional[dict]:
        """
        Fetch a single row as dict.
        Returns None if DB unavailable, not found, or on error.
        """
        try:
            with self.cursor(dict_cursor=True) as cur:
                cur.execute(sql, params)
                row = cur.fetchone()
                return dict(row) if row else None
        except DBNotAvailable:
            return None
        except Exception as exc:
            logger.error("DB fetchone error: %s | sql=%s", exc, sql[:120])
            return None

    def fetchall(self, sql: str, params=None) -> list[dict]:
        """
        Fetch all rows as list of dicts.
        Returns [] if DB unavailable or on error.
        """
        try:
            with self.cursor(dict_cursor=True) as cur:
                cur.execute(sql, params)
                return [dict(r) for r in cur.fetchall()]
        except DBNotAvailable:
            return []
        except Exception as exc:
            logger.error("DB fetchall error: %s | sql=%s", exc, sql[:120])
            return []

    def upsert_local_pgs(self, info: dict) -> bool:
        """
        Upsert a single PGS entry into local_pgs_index.
        Safe to call repeatedly (idempotent).
        """
        import json as _json
        return self.execute(
            """
            INSERT INTO local_pgs_index
                (pgs_id, trait_name, trait_efo, n_variants,
                 has_parquet, has_tsv, has_metadata,
                 chromosomes, n_chromosomes, last_synced_at)
            VALUES
                (%s, %s, %s::jsonb, %s,
                 %s, %s, %s,
                 %s::jsonb, %s, NOW())
            ON CONFLICT (pgs_id) DO UPDATE SET
                trait_name      = EXCLUDED.trait_name,
                trait_efo       = EXCLUDED.trait_efo,
                n_variants      = EXCLUDED.n_variants,
                has_parquet     = EXCLUDED.has_parquet,
                has_tsv         = EXCLUDED.has_tsv,
                has_metadata    = EXCLUDED.has_metadata,
                chromosomes     = EXCLUDED.chromosomes,
                n_chromosomes   = EXCLUDED.n_chromosomes,
                last_synced_at  = NOW()
            """,
            (
                info["pgs_id"],
                info.get("trait_name"),
                _json.dumps(info.get("trait_efo", [])),
                info.get("n_variants"),
                info.get("has_parquet", False),
                info.get("has_tsv", False),
                info.get("meta_available", False),
                _json.dumps(info.get("chromosomes", [])),
                info.get("n_chromosomes", 0),
            ),
        )

    def upsert_remote_cache(self, pgs_id: str, data: dict, ttl_seconds: int = 86400) -> bool:
        """Cache a normalised pgscat response. TTL default = 24 hours."""
        import json as _json
        from datetime import datetime, timezone, timedelta
        expires = datetime.now(timezone.utc) + timedelta(seconds=ttl_seconds)
        pub = data.get("publication") or {}
        return self.execute(
            """
            INSERT INTO remote_pgs_cache
                (pgs_id, name, trait_reported, trait_efo, variants_number,
                 is_harmonized, ftp_scoring_file, publication, raw_response,
                 fetched_at, ttl_expires_at)
            VALUES
                (%s, %s, %s, %s::jsonb, %s,
                 %s, %s, %s::jsonb, %s::jsonb,
                 NOW(), %s)
            ON CONFLICT (pgs_id) DO UPDATE SET
                name              = EXCLUDED.name,
                trait_reported    = EXCLUDED.trait_reported,
                trait_efo         = EXCLUDED.trait_efo,
                variants_number   = EXCLUDED.variants_number,
                is_harmonized     = EXCLUDED.is_harmonized,
                ftp_scoring_file  = EXCLUDED.ftp_scoring_file,
                publication       = EXCLUDED.publication,
                raw_response      = EXCLUDED.raw_response,
                fetched_at        = NOW(),
                ttl_expires_at    = EXCLUDED.ttl_expires_at
            """,
            (
                pgs_id,
                data.get("name"),
                data.get("trait_reported"),
                _json.dumps(data.get("trait_efo", [])),
                data.get("variants_number"),
                data.get("is_harmonized", False),
                data.get("ftp_scoring_file"),
                _json.dumps(pub),
                _json.dumps(data.get("_raw") or data),
                expires,
            ),
        )

    def get_remote_cache(self, pgs_id: str) -> Optional[dict]:
        """Return cached remote metadata if not expired, else None."""
        return self.fetchone(
            """
            SELECT pgs_id, name, trait_reported, trait_efo,
                   variants_number, is_harmonized, ftp_scoring_file,
                   publication, fetched_at
            FROM remote_pgs_cache
            WHERE pgs_id = %s AND ttl_expires_at > NOW()
            """,
            (pgs_id,),
        )

    def upsert_stats_cache(self, pgs_id: str, stats: dict, source_mtime: int) -> bool:
        """Cache DuckDB-computed stats for a PGS dataset."""
        import json as _json
        return self.execute(
            """
            INSERT INTO stats_cache (pgs_id, score_column, data_format, stats_json, source_mtime, computed_at)
            VALUES (%s, %s, %s, %s::jsonb, %s, NOW())
            ON CONFLICT (pgs_id) DO UPDATE SET
                score_column = EXCLUDED.score_column,
                data_format  = EXCLUDED.data_format,
                stats_json   = EXCLUDED.stats_json,
                source_mtime = EXCLUDED.source_mtime,
                computed_at  = NOW()
            """,
            (
                pgs_id,
                stats.get("score_column"),
                stats.get("data_format"),
                _json.dumps(stats),
                source_mtime,
            ),
        )

    def get_stats_cache(self, pgs_id: str, source_mtime: int) -> Optional[dict]:
        """Return cached stats if source file mtime matches."""
        import json as _json
        row = self.fetchone(
            "SELECT stats_json FROM stats_cache WHERE pgs_id = %s AND source_mtime = %s",
            (pgs_id, source_mtime),
        )
        if row and row.get("stats_json"):
            v = row["stats_json"]
            return v if isinstance(v, dict) else _json.loads(v)
        return None

    def log_search(
        self,
        query: str,
        search_type: str,
        n_results: int,
        from_cache: bool,
        error: Optional[str],
        duration_ms: int,
        client_ip: Optional[str] = None,
    ) -> None:
        """Append a row to search_audit. Silently ignores failures."""
        self.execute(
            """
            INSERT INTO search_audit
                (query, search_type, n_results, from_cache, error, duration_ms, client_ip)
            VALUES (%s, %s, %s, %s, %s, %s, %s::inet)
            """,
            (query, search_type, n_results, from_cache, error, duration_ms, client_ip),
        )
