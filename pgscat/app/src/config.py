"""
Central configuration – loaded once from environment variables.
All paths, credentials, and feature flags live here.
No other module reads os.environ directly.

Environment variables (full reference):
  APP_DATA_DIR              /data           scores directory (read-only mount)
  APP_WORK_DIR              /work           writable temp/cache directory
  APP_SCRIPTS_DIR           /pipeline_scripts  pipeline scripts (read-only mount)
  APP_ANNOTATIONS_DIR       /annotations    variant annotation outputs (read-only from web)

  APP_PGSCAT_MODE           python          "python" | "cli"
  APP_PGSCAT_BIN            pgscat          pgscat executable name/path
  APP_PGSCAT_TIMEOUT        30              seconds before remote query is killed
  APP_REMOTE_ENABLED        true            enable/disable remote PGS Catalog queries
  APP_REMOTE_CACHE_TTL      86400           seconds to cache remote metadata (default 24h)

  APP_DB_HOST               postgres        hostname (use "postgres" inside container network,
                                           "127.0.0.1" when running on host)
  APP_DB_PORT               5432
  APP_DB_NAME               pgs_dashboard
  APP_DB_USER               pgs_user
  APP_DB_PASSWORD_FILE                      path to file containing DB password (preferred)
  APP_DB_PASSWORD                           fallback env var for DB password

  APP_STATS_CACHE_TTL       3600            seconds for DuckDB stats file cache
  APP_TRUSTED_PROXY         true            trust X-Forwarded-* headers from reverse proxy
  APP_HOST                  0.0.0.0         bind host (dev only; gunicorn uses --bind)
  APP_PORT                  8080            bind port (dev only)
"""
import logging
import os
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)


class Config:
    # ── Paths ────────────────────────────────────────────────────────────────
    DATA_DIR: Path
    WORK_DIR: Path
    SCRIPTS_DIR: Path
    ANNOTATIONS_DIR: Path   # variant annotation outputs (written by Apptainer pipeline)

    # ── pgscat ───────────────────────────────────────────────────────────────
    PGSCAT_MODE: str
    PGSCAT_BIN: str
    PGSCAT_TIMEOUT: int
    REMOTE_ENABLED: bool
    REMOTE_CACHE_TTL: int     # seconds to cache remote pgscat responses in DB

    # ── PostgreSQL ───────────────────────────────────────────────────────────
    DB_HOST: str
    DB_PORT: int
    DB_NAME: str
    DB_USER: str
    _db_password: Optional[str]      # resolved once at startup

    # ── Cache ────────────────────────────────────────────────────────────────
    STATS_CACHE_TTL: int      # seconds for file-based DuckDB stats cache

    # ── Server ───────────────────────────────────────────────────────────────
    HOST: str
    PORT: int
    TRUSTED_PROXY: bool       # trust X-Forwarded-* headers

    def __init__(self) -> None:
        self.DATA_DIR = Path(os.environ.get("APP_DATA_DIR", "/data"))
        self.WORK_DIR = Path(os.environ.get("APP_WORK_DIR", "/work"))
        self.SCRIPTS_DIR = Path(os.environ.get(
            "APP_SCRIPTS_DIR",
            "/pipeline_scripts",
        ))
        self.ANNOTATIONS_DIR = Path(os.environ.get(
            "APP_ANNOTATIONS_DIR",
            "/annotations",
        ))

        self.PGSCAT_MODE = os.environ.get("APP_PGSCAT_MODE", "python").lower()
        self.PGSCAT_BIN = os.environ.get("APP_PGSCAT_BIN", "pgscat")
        self.PGSCAT_TIMEOUT = int(os.environ.get("APP_PGSCAT_TIMEOUT", "30"))
        self.REMOTE_ENABLED = os.environ.get("APP_REMOTE_ENABLED", "true").lower() == "true"
        self.REMOTE_CACHE_TTL = int(os.environ.get("APP_REMOTE_CACHE_TTL", "86400"))

        self.DB_HOST = os.environ.get("APP_DB_HOST", "postgres")
        self.DB_PORT = int(os.environ.get("APP_DB_PORT", "5432"))
        self.DB_NAME = os.environ.get("APP_DB_NAME", "pgs_dashboard")
        self.DB_USER = os.environ.get("APP_DB_USER", "pgs_user")
        self._db_password = self._resolve_db_password()

        self.STATS_CACHE_TTL = int(os.environ.get("APP_STATS_CACHE_TTL", "3600"))
        self.TRUSTED_PROXY = os.environ.get("APP_TRUSTED_PROXY", "true").lower() == "true"
        self.HOST = os.environ.get("APP_HOST", "0.0.0.0")
        self.PORT = int(os.environ.get("APP_PORT", "8080"))

    # ── Properties ────────────────────────────────────────────────────────────

    @property
    def scores_dir(self) -> Path:
        return self.DATA_DIR

    @property
    def db_password(self) -> str:
        return self._db_password or ""

    # ── Helpers ───────────────────────────────────────────────────────────────

    def _resolve_db_password(self) -> Optional[str]:
        """
        Load DB password from file if APP_DB_PASSWORD_FILE is set and readable;
        otherwise fall back to APP_DB_PASSWORD env var.
        Logs a warning (not an error) if neither is set — DB will fail to connect.
        """
        pw_file = os.environ.get("APP_DB_PASSWORD_FILE", "").strip()
        if pw_file:
            try:
                pw = Path(pw_file).read_text(encoding="utf-8").strip()
                logger.info("DB password loaded from file: %s", pw_file)
                return pw
            except Exception as exc:
                logger.warning("Cannot read DB password file %s: %s", pw_file, exc)

        pw_env = os.environ.get("APP_DB_PASSWORD", "").strip()
        if pw_env:
            logger.info("DB password loaded from APP_DB_PASSWORD env var")
            return pw_env

        logger.warning(
            "No DB password configured (APP_DB_PASSWORD_FILE or APP_DB_PASSWORD). "
            "PostgreSQL connection will likely fail."
        )
        return None

    def __repr__(self) -> str:
        return (
            f"Config("
            f"DATA_DIR={self.DATA_DIR}, "
            f"WORK_DIR={self.WORK_DIR}, "
            f"ANNOTATIONS_DIR={self.ANNOTATIONS_DIR}, "
            f"DB={self.DB_USER}@{self.DB_HOST}:{self.DB_PORT}/{self.DB_NAME}, "
            f"PGSCAT_MODE={self.PGSCAT_MODE}, "
            f"REMOTE_ENABLED={self.REMOTE_ENABLED}, "
            f"TRUSTED_PROXY={self.TRUSTED_PROXY}"
            f")"
        )
