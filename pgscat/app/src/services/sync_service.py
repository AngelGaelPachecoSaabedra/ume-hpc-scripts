"""
SyncService – local filesystem → PostgreSQL catalog sync
═════════════════════════════════════════════════════════

Reads the scores directory and upserts every valid PGS directory
into local_pgs_index + file_inventory.  Completely non-destructive:
uses INSERT … ON CONFLICT DO UPDATE everywhere.

Runs:
  • Automatically in a background daemon thread at startup (5s delay).
  • On demand via POST /api/admin/sync (no authentication in Phase 1;
    add auth before exposing externally).

If PostgreSQL is unavailable the sync logs a warning and exits immediately.
The app continues in filesystem-only mode.

Thread safety:
  The sync is idempotent. Multiple gunicorn workers may run it in parallel
  — the upserts are safe for concurrent execution.
"""
import json
import logging
import threading
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

from config import Config
from services.db import DBPool
from services.local_catalog import LocalCatalog, CHROMOSOMES

logger = logging.getLogger(__name__)

# File type classification (order matters — checked in order)
_FILE_TYPE_RULES: list[tuple[str, str]] = [
    ("_PRS_total.parquet", "parquet"),
    ("_PRS_total.tsv", "tsv_total"),
    ("_PRS_total_metadata.json", "metadata"),
    ("_hmPOS_GRCh38.betamap.tsv.gz", "betamap"),
    (".json", "json"),
    ("_metadata.json", "metadata"),
]


def _classify_file(filename: str) -> tuple[str, Optional[str]]:
    """Return (file_type, chrom_or_None) for a filename."""
    # per-chrom score file: PGS000004_chr1_scores.tsv
    for c in CHROMOSOMES:
        if f"_chr{c}_scores.tsv" in filename:
            return "tsv_chrom", c
        if f"_chr{c}_metadata.json" in filename:
            return "metadata", c
    for suffix, ftype in _FILE_TYPE_RULES:
        if filename.endswith(suffix):
            return ftype, None
    return "other", None


class SyncService:
    def __init__(self, cfg: Config, db: DBPool, catalog: LocalCatalog) -> None:
        self.cfg = cfg
        self.db = db
        self.catalog = catalog
        self._thread: Optional[threading.Thread] = None

    # ── Public API ─────────────────────────────────────────────────────────────

    def start_background(self, delay: int = 8) -> None:
        """
        Launch sync in a daemon background thread after `delay` seconds.
        The delay lets gunicorn finish forking before the sync starts.
        """
        if self._thread and self._thread.is_alive():
            logger.debug("Sync thread already running, skipping launch")
            return

        def _run() -> None:
            time.sleep(delay)
            try:
                result = self.run_sync()
                logger.info("Background sync complete: %s", result)
            except Exception as exc:
                logger.error("Background sync error: %s", exc)

        self._thread = threading.Thread(
            target=_run, name="catalog-sync", daemon=True
        )
        self._thread.start()
        logger.info("Catalog sync scheduled (delay=%ds)", delay)

    def run_sync(self) -> dict:
        """
        Scan scores directory and upsert all PGS entries into PostgreSQL.
        Returns a summary dict.  Never raises; errors are logged and returned.
        """
        if not self.db.available:
            logger.debug("Sync skipped: DB not available")
            return {"status": "skipped", "reason": "db_unavailable"}

        started = time.monotonic()
        run_id = self._start_run()
        added = updated = scanned = 0

        try:
            pgs_list = self.catalog.list_pgs()
            scanned = len(pgs_list)

            for quick_info in pgs_list:
                pgs_id = quick_info["pgs_id"]
                # Check if already indexed
                existing = self.db.fetchone(
                    "SELECT last_synced_at FROM local_pgs_index WHERE pgs_id = %s",
                    (pgs_id,),
                )
                was_new = existing is None

                ok = self.db.upsert_local_pgs(quick_info)
                if ok:
                    if was_new:
                        added += 1
                    else:
                        updated += 1
                    self._sync_file_inventory(pgs_id)

        except Exception as exc:
            logger.error("Sync error after %d/%d: %s", added + updated, scanned, exc)
            self._finish_run(run_id, scanned, added, updated, status="failed", error=str(exc))
            return {"status": "failed", "error": str(exc), "scanned": scanned}

        elapsed = round(time.monotonic() - started, 2)
        self._finish_run(run_id, scanned, added, updated, status="done")
        return {
            "status": "done",
            "scanned": scanned,
            "added": added,
            "updated": updated,
            "elapsed_s": elapsed,
        }

    # ── Internal ──────────────────────────────────────────────────────────────

    def _sync_file_inventory(self, pgs_id: str) -> None:
        pgs_dir = self.cfg.scores_dir / pgs_id
        if not pgs_dir.is_dir():
            return
        for f in pgs_dir.iterdir():
            if not f.is_file():
                continue
            ftype, chrom = _classify_file(f.name)
            try:
                size = f.stat().st_size
            except OSError:
                size = None
            self.db.execute(
                """
                INSERT INTO file_inventory (pgs_id, filename, file_type, chrom, size_bytes, last_seen_at)
                VALUES (%s, %s, %s, %s, %s, NOW())
                ON CONFLICT (pgs_id, filename) DO UPDATE SET
                    file_type    = EXCLUDED.file_type,
                    chrom        = EXCLUDED.chrom,
                    size_bytes   = EXCLUDED.size_bytes,
                    last_seen_at = NOW()
                """,
                (pgs_id, f.name, ftype, chrom, size),
            )

    def _start_run(self) -> Optional[int]:
        """Insert a sync_runs row and return its id."""
        try:
            with self.db.cursor() as cur:
                cur.execute(
                    "INSERT INTO sync_runs (status) VALUES ('running') RETURNING id"
                )
                row = cur.fetchone()
                return row[0] if row else None
        except Exception:
            return None

    def _finish_run(
        self,
        run_id: Optional[int],
        scanned: int,
        added: int,
        updated: int,
        status: str,
        error: Optional[str] = None,
    ) -> None:
        if not run_id:
            return
        self.db.execute(
            """
            UPDATE sync_runs
            SET finished_at = NOW(), pgs_scanned = %s,
                pgs_added = %s, pgs_updated = %s,
                status = %s, error = %s
            WHERE id = %s
            """,
            (scanned, added, updated, status, error, run_id),
        )
