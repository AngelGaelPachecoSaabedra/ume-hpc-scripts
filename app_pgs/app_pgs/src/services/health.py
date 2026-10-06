"""
Health and readiness checks
════════════════════════════

/healthz  – liveness:  is the process alive? (almost always 200)
/readyz   – readiness: can the app serve requests?

Both return JSON:
  {
    "status": "ok" | "degraded" | "error",
    "components": { "data_dir": "ok", "database": "degraded", "pgscat": "ok" },
    "degraded": ["database"],
    "version": "1.0"
  }

Degraded ≠ error: a degraded app continues serving dashboards from
parquet files.  Only "error" status (HTTP 503) means nothing works.

HTTP codes:
  /healthz → 200 always (liveness; if the process is dead it can't respond)
  /readyz  → 200 (ok/degraded) | 503 (error — data dir missing)
"""
import logging
import time
from typing import Optional

from config import Config
from services.db import DBPool

logger = logging.getLogger(__name__)

APP_VERSION = "1.1.0"


class HealthService:
    def __init__(self, cfg: Config, db: DBPool) -> None:
        self.cfg = cfg
        self.db = db

    def liveness(self) -> dict:
        """
        Minimal liveness check.
        Returns 200 unless the process can't even answer (it can't).
        """
        return {
            "status": "ok",
            "version": APP_VERSION,
            "timestamp": _now(),
        }

    def readiness(self) -> tuple[dict, int]:
        """
        Full readiness check.
        Returns (response_dict, http_status_code).
        """
        components: dict[str, str] = {}
        degraded: list[str] = []

        # 1. Data directory (critical – without this, nothing works)
        data_ok = self._check_data_dir()
        components["data_dir"] = "ok" if data_ok else "error"
        if not data_ok:
            return {
                "status": "error",
                "components": components,
                "degraded": [],
                "error": f"Data directory not accessible: {self.cfg.DATA_DIR}",
                "version": APP_VERSION,
                "timestamp": _now(),
            }, 503

        # 2. PostgreSQL (optional)
        db_ok = self._check_db()
        components["database"] = "ok" if db_ok else "degraded"
        if not db_ok:
            degraded.append("database")

        # 3. Work directory (optional)
        work_ok = self._check_work_dir()
        components["work_dir"] = "ok" if work_ok else "degraded"
        if not work_ok:
            degraded.append("work_dir")

        status = "ok" if not degraded else "degraded"
        return {
            "status": status,
            "components": components,
            "degraded": degraded,
            "version": APP_VERSION,
            "timestamp": _now(),
        }, 200

    # ── Component checks ──────────────────────────────────────────────────────

    def _check_data_dir(self) -> bool:
        try:
            return self.cfg.DATA_DIR.is_dir()
        except Exception:
            return False

    def _check_db(self) -> bool:
        if not self.db.available:
            # Try a lightweight reconnect probe
            return self.db.retry_connect()
        try:
            row = self.db.fetchone("SELECT 1 AS alive")
            return row is not None
        except Exception:
            return False

    def _check_work_dir(self) -> bool:
        try:
            if not self.cfg.WORK_DIR.exists():
                return False
            # Try writing a probe file
            probe = self.cfg.WORK_DIR / ".health_probe"
            probe.write_text("ok")
            probe.unlink(missing_ok=True)
            return True
        except Exception:
            return False


def _now() -> str:
    from datetime import datetime, timezone
    return datetime.now(timezone.utc).isoformat()
