"""
PGS/PRS Dashboard  –  main Bottle application
══════════════════════════════════════════════
Production URL: https://pgscat.mcps-epcm.org/

Routes:
  GET  /                              – local catalog
  GET  /dashboard/<pgs_id>            – per-PGS dashboard
  GET  /api/data/<pgs_id>             – JSON stats (DuckDB)
  GET  /search                        – remote search UI
  POST /search                        – execute remote search
  GET  /pgs/<pgs_id>/source           – local vs remote status
  GET  /api/pipeline/<pgs_id>/plan    – pipeline plan JSON (no execution)
  POST /api/admin/sync                – trigger catalog→DB sync manually
  GET  /healthz                       – liveness probe
  GET  /readyz                        – readiness probe (JSON)
  GET  /static/<path>                 – static assets
"""
import json
import logging
import time
from functools import wraps
from pathlib import Path

import bottle
from bottle import Bottle, abort, request, response, static_file, template

from config import Config
from services.cache import FileCache
from services.db import DBPool
from services.health import HealthService
from services.local_catalog import LocalCatalog, validate_pgs_id
from services.parquet_stats import ParquetStats
from services.pgscat_client import PGSCatClient
from services.pipeline_inspector import PipelineInspector
from services.sync_service import SyncService

# ── Logging ──────────────────────────────────────────────────────────────────
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(name)s: %(message)s",
)
logger = logging.getLogger(__name__)

# ── Template path ─────────────────────────────────────────────────────────────
HERE = Path(__file__).parent
bottle.TEMPLATE_PATH.insert(0, str(HERE / "views"))

# ── Services ──────────────────────────────────────────────────────────────────
cfg = Config()
db = DBPool(cfg)
cache = FileCache(cfg.WORK_DIR, ttl=cfg.STATS_CACHE_TTL)
catalog = LocalCatalog(cfg)
stats_svc = ParquetStats(cfg, cache=cache, db=db)
pgscat_client = PGSCatClient(cfg, db=db)
inspector = PipelineInspector(cfg)
sync_svc = SyncService(cfg, db, catalog)
health_svc = HealthService(cfg, db)

logger.info("PGS Dashboard starting — %s", cfg)

# ── Bottle app ────────────────────────────────────────────────────────────────
_bottle_app = Bottle()


# ── Proxy middleware (WSGI-level) ─────────────────────────────────────────────

class _ProxyMiddleware:
    """
    Honour X-Forwarded-* headers set by the upstream reverse proxy (Nginx/Caddy).
    Only applied when APP_TRUSTED_PROXY=true.

    Handles:
      X-Forwarded-For   → REMOTE_ADDR  (first IP in the chain)
      X-Forwarded-Proto → wsgi.url_scheme (http/https)
      X-Forwarded-Host  → HTTP_HOST
      X-Forwarded-Prefix → SCRIPT_NAME (for sub-path deployments)
    """
    def __init__(self, wsgi_app):
        self._app = wsgi_app

    def __call__(self, environ, start_response):
        xff = environ.get("HTTP_X_FORWARDED_FOR", "")
        if xff:
            environ["REMOTE_ADDR"] = xff.split(",")[0].strip()

        proto = environ.get("HTTP_X_FORWARDED_PROTO", "")
        if proto in ("http", "https"):
            environ["wsgi.url_scheme"] = proto

        host = environ.get("HTTP_X_FORWARDED_HOST", "")
        if host:
            environ["HTTP_HOST"] = host.split(",")[0].strip()

        prefix = environ.get("HTTP_X_FORWARDED_PREFIX", "")
        if prefix:
            environ["SCRIPT_NAME"] = prefix.rstrip("/")

        return self._app(environ, start_response)


# The exported WSGI app (used by gunicorn)
app = _ProxyMiddleware(_bottle_app) if cfg.TRUSTED_PROXY else _bottle_app


# ── Background sync ───────────────────────────────────────────────────────────
# Start catalog→DB sync after 8s so gunicorn finishes forking first.
# Each worker will run its own sync (idempotent upserts make this safe).
sync_svc.start_background(delay=8)


# ── Helpers ───────────────────────────────────────────────────────────────────

def _json(data: dict, status: int = 200) -> str:
    response.content_type = "application/json"
    response.status = status
    return json.dumps(data, default=str)


def _require_pgs_id(fn):
    @wraps(fn)
    def _wrapper(pgs_id: str, *args, **kwargs):
        if not validate_pgs_id(pgs_id):
            abort(
                400,
                f"Invalid PGS ID {pgs_id!r}. "
                "Expected 'PGS' followed by exactly 6 digits (e.g. PGS000004).",
            )
        return fn(pgs_id, *args, **kwargs)
    return _wrapper


def _client_ip() -> str:
    """Best-effort client IP, respecting X-Forwarded-For if trusted."""
    return request.environ.get("REMOTE_ADDR", "")


# ── Static files ──────────────────────────────────────────────────────────────

@_bottle_app.route("/static/<filename:path>")
def static(filename: str):
    return static_file(filename, root=str(HERE / "static"))


# ── Health & readiness ────────────────────────────────────────────────────────

@_bottle_app.route("/healthz")
def healthz():
    data = health_svc.liveness()
    return _json(data)


@_bottle_app.route("/readyz")
def readyz():
    data, status = health_svc.readiness()
    return _json(data, status=status)


# ── Home / Catalog ────────────────────────────────────────────────────────────

@_bottle_app.route("/")
def index():
    # Try DB index first (faster), fall back to filesystem scan
    if db.available:
        rows = db.fetchall(
            """
            SELECT pgs_id, trait_name, n_variants,
                   has_parquet, has_tsv, (has_parquet OR has_tsv) AS has_results,
                   n_chromosomes, last_synced_at
            FROM local_pgs_index
            ORDER BY pgs_id
            """
        )
        if rows:
            return template("index", pgs_list=rows, cfg=cfg, from_db=True)

    pgs_list = catalog.list_pgs()
    return template("index", pgs_list=pgs_list, cfg=cfg, from_db=False)


# ── Dashboard ─────────────────────────────────────────────────────────────────

@_bottle_app.route("/dashboard/<pgs_id>")
@_require_pgs_id
def dashboard(pgs_id: str):
    info = catalog.get_pgs_info(pgs_id)
    if info is None:
        abort(404, f"PGS {pgs_id} not found in local catalog.")
    return template("dashboard", pgs_id=pgs_id, info=info, cfg=cfg)


# ── API: stats ────────────────────────────────────────────────────────────────

@_bottle_app.route("/api/data/<pgs_id>")
@_require_pgs_id
def api_data(pgs_id: str):
    stats = stats_svc.get_stats(pgs_id)
    if "error" in stats:
        code = 404 if "not found" in stats["error"].lower() else 500
        return _json(stats, status=code)
    return _json(stats)


# ── Search ────────────────────────────────────────────────────────────────────

@_bottle_app.route("/search", method=["GET", "POST"])
def search():
    query = ""
    search_type = "text"
    results = None
    error = None

    if request.method == "POST":
        query = (request.forms.get("query") or "").strip()
        search_type = (request.forms.get("search_type") or "text").strip()
        if search_type not in ("id", "trait", "pmid", "text"):
            search_type = "text"

        if not query:
            error = "Please enter a search term."
        else:
            t0 = time.monotonic()
            raw_results, from_cache = pgscat_client.search_with_cache(query, search_type)
            duration_ms = int((time.monotonic() - t0) * 1000)

            search_error = None
            results = []
            for r in raw_results:
                pid = r.get("pgs_id", "")
                r["exists_locally"] = catalog.exists_locally(pid) if pid else False
                r["has_results"] = catalog.has_results(pid) if pid else False
                results.append(r)

            # Collapse single-error result into page-level error
            if len(results) == 1 and "error" in results[0] and not results[0].get("pgs_id"):
                search_error = results[0]["error"]
                results = []
                error = search_error

            # Audit log (best-effort, non-blocking)
            db.log_search(
                query=query,
                search_type=search_type,
                n_results=len(results),
                from_cache=from_cache,
                error=search_error,
                duration_ms=duration_ms,
                client_ip=_client_ip(),
            )

    return template(
        "search",
        query=query,
        search_type=search_type,
        results=results,
        error=error,
        remote_enabled=cfg.REMOTE_ENABLED,
        cfg=cfg,
    )


# ── PGS Source (local vs remote) ─────────────────────────────────────────────

@_bottle_app.route("/pgs/<pgs_id>/source")
@_require_pgs_id
def pgs_source(pgs_id: str):
    local_info = catalog.get_pgs_info(pgs_id) if catalog.exists_locally(pgs_id) else None
    remote_info = None
    remote_error = None

    if cfg.REMOTE_ENABLED:
        result, _ = pgscat_client.get_score_with_cache(pgs_id)
        if "error" in result:
            remote_error = result["error"]
        else:
            remote_info = result

    if local_info and remote_info:
        source = "both"
    elif local_info:
        source = "local_only"
    elif remote_info:
        source = "remote_only"
    else:
        source = "not_found"

    plan = inspector.get_pipeline_plan(pgs_id) if source == "remote_only" else None

    return template(
        "pgs_remote",
        pgs_id=pgs_id,
        source=source,
        local_info=local_info,
        remote_info=remote_info,
        remote_error=remote_error,
        plan=plan,
        cfg=cfg,
    )


# ── API: Pipeline plan ────────────────────────────────────────────────────────

@_bottle_app.route("/api/pipeline/<pgs_id>/plan")
@_require_pgs_id
def api_pipeline_plan(pgs_id: str):
    plan = inspector.get_pipeline_plan(pgs_id)
    if "error" in plan:
        return _json(plan, status=400)
    return _json(plan)


# ── Admin: manual sync ────────────────────────────────────────────────────────

@_bottle_app.route("/api/admin/sync", method="POST")
def admin_sync():
    """
    Trigger a synchronous catalog→DB sync.
    Not authenticated in Phase 1 — protect via reverse proxy auth or network policy.
    """
    if not db.available:
        return _json({"status": "skipped", "reason": "db_unavailable"}, status=503)
    result = sync_svc.run_sync()
    return _json(result)


# ── Error handlers ────────────────────────────────────────────────────────────

@_bottle_app.error(400)
def error_400(err):
    return template("error", code=400, message=err.body, cfg=cfg)


@_bottle_app.error(404)
def error_404(err):
    return template("error", code=404, message=err.body, cfg=cfg)


@_bottle_app.error(500)
def error_500(err):
    logger.error("Unhandled 500: %s", err)
    return template("error", code=500, message="Internal server error. Check logs.", cfg=cfg)


# ── Dev runner ────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    bottle.run(
        app=_bottle_app,   # dev: use Bottle directly (no proxy middleware needed)
        host=cfg.HOST,
        port=cfg.PORT,
        debug=True,
        reloader=True,
    )
