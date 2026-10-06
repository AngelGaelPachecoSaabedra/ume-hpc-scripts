"""
File-based JSON cache with TTL
══════════════════════════════

Stores computed stats as JSON files in $APP_WORK_DIR/cache/.
This is the primary cache layer — fast, no DB required.
PostgreSQL stats_cache is the secondary, persistent layer.

Cache key examples:
  stats_PGS000004   → /work/cache/stats_PGS000004.json

Invalidation strategy:
  Each cache entry stores the source file's mtime.
  On read, if the source file's mtime differs → cache is stale.

Thread safety:
  Write uses atomic rename (write to .tmp then rename), so concurrent
  gunicorn workers cannot corrupt each other's cache files.
"""
import json
import logging
import os
import tempfile
import time
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)


class FileCache:
    def __init__(self, work_dir: Path, ttl: int = 3600) -> None:
        self.cache_dir = work_dir / "cache"
        self.ttl = ttl
        self._ensure_dir()

    def _ensure_dir(self) -> None:
        try:
            self.cache_dir.mkdir(parents=True, exist_ok=True)
        except Exception as exc:
            logger.warning("Cannot create cache dir %s: %s", self.cache_dir, exc)

    def _path(self, key: str) -> Path:
        # Sanitise: only allow alphanumeric + underscore
        safe = "".join(c if c.isalnum() or c == "_" else "_" for c in key)
        return self.cache_dir / f"{safe}.json"

    def get(self, key: str) -> Optional[dict]:
        """
        Return cached data if it exists and is within TTL.
        Returns None on miss, expiry, or read error.
        """
        path = self._path(key)
        try:
            if not path.exists():
                return None
            age = time.time() - path.stat().st_mtime
            if age > self.ttl:
                logger.debug("Cache expired for %s (age=%.0fs)", key, age)
                return None
            data = json.loads(path.read_text(encoding="utf-8"))
            logger.debug("Cache hit: %s", key)
            return data
        except Exception as exc:
            logger.debug("Cache read error for %s: %s", key, exc)
            return None

    def get_with_mtime_check(self, key: str, source_mtime: int) -> Optional[dict]:
        """
        Return cached data only if _source_mtime matches the given value.
        Useful for file-backed data where content changes mean re-compute.
        """
        data = self.get(key)
        if data and data.get("_source_mtime") == source_mtime:
            return data
        if data and data.get("_source_mtime") != source_mtime:
            logger.debug("Cache mtime mismatch for %s (cached=%s, current=%s)",
                         key, data.get("_source_mtime"), source_mtime)
        return None

    def set(self, key: str, data: dict) -> bool:
        """
        Write data to cache atomically using a temp file + rename.
        Returns True on success.
        """
        path = self._path(key)
        try:
            self._ensure_dir()
            with tempfile.NamedTemporaryFile(
                mode="w",
                encoding="utf-8",
                dir=self.cache_dir,
                suffix=".tmp",
                delete=False,
            ) as tmp:
                json.dump(data, tmp, default=str)
                tmp_path = tmp.name
            os.replace(tmp_path, path)  # atomic on same filesystem
            logger.debug("Cache set: %s", key)
            return True
        except Exception as exc:
            logger.warning("Cache write error for %s: %s", key, exc)
            return False

    def invalidate(self, key: str) -> None:
        try:
            self._path(key).unlink(missing_ok=True)
        except Exception:
            pass

    def stats_key(self, pgs_id: str) -> str:
        return f"stats_{pgs_id}"
