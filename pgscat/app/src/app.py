"""
PGS/PRS Platform  –  main Bottle application
═════════════════════════════════════════════
Platform modes:
  A  Local catalog    – list, visualise, and inspect locally computed PGS/PRS
  B  Remote catalog   – search PGS Catalog, view prep plan for missing scores
  C  Custom PRS       – define a custom score and get a pipeline plan
  D  Variant Annotator– view GFF3-based variant annotation results (HPC pipeline)

Routes:
  GET  /                              – local catalog (Mode A)
  GET  /dashboard/<pgs_id>            – per-PGS stats dashboard (Mode A)
  GET  /api/data/<pgs_id>             – JSON stats (DuckDB)
  GET  /search                        – remote search UI (Mode B)
  POST /search                        – execute remote search
  GET  /pgs/<pgs_id>/source           – local vs remote status + prep plan
  GET  /pipeline/<pgs_id>/plan        – HTML pipeline plan view
  GET  /api/pipeline/<pgs_id>/plan    – pipeline plan JSON
  GET  /api/pipeline/<pgs_id>/manifest – preparation manifest JSON
  GET  /custom-prs/new                – custom PRS form (Mode C)
  POST /api/custom-prs/plan           – custom PRS plan JSON
  POST /api/pgs/<pgs_id>/prepare      – execute score preparation (download + betamap)
  POST /api/admin/sync                – trigger catalog→DB sync manually
  GET  /variants/<pgs_id>             – variant annotation view (Mode D)
  GET  /api/variants/<pgs_id>         – annotated variants JSON (paginated)
  GET  /api/variants/<pgs_id>/summary – annotation summary JSON
  GET  /api/variants/<pgs_id>/ideogram – ideogram data JSON (future visualisation)
  GET  /api/variants/<pgs_id>/charts  – pre-aggregated analytical chart data (G1–G8)
  GET  /healthz                       – liveness probe
  GET  /readyz                        – readiness probe (JSON)
  GET  /static/<path>                 – static assets

Mode E — Gene Browser (gene-centric interpretation):
  GET  /genes                             – gene list / search UI
  GET  /gene/<gene_symbol>               – gene browser HTML (Plotly tracks)
  GET  /api/gene/<gene_symbol>           – gene info JSON
  GET  /api/gene/<gene_symbol>/tracks    – Plotly track data JSON
  GET  /api/gene/<gene_symbol>/variants  – variants in gene (paginated JSON)
  GET  /api/genes                        – annotated gene list JSON (autocomplete)
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
from services.normalizer import normalize_local_info, normalize_remote_info, normalize_search_result
from services.parquet_stats import ParquetStats
from services.pgscat_client import PGSCatClient
from services.pipeline_inspector import PipelineInspector
from services.score_preparer import ScorePreparer
from services.sync_service import SyncService
from services.variant_annotation import VariantAnnotationService
from services.gene_browser import GeneBrowserService
from services.variant_ranking import WEIGHTS as _RANKING_WEIGHTS, CONSEQUENCE_SCORES as _CSQ_SCORES

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
preparer = ScorePreparer(cfg)
sync_svc = SyncService(cfg, db, catalog)
health_svc = HealthService(cfg, db)
variant_svc = VariantAnnotationService(cfg)
gene_svc    = GeneBrowserService(cfg)

logger.info("PGS Platform starting — %s", cfg)

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


# ── Home / Catalog (Mode A) ───────────────────────────────────────────────────

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


# ── Dashboard (Mode A) ────────────────────────────────────────────────────────

@_bottle_app.route("/dashboard/<pgs_id>")
@_require_pgs_id
def dashboard(pgs_id: str):
    info = catalog.get_pgs_info(pgs_id)
    if info is None:
        abort(404, f"PGS {pgs_id} not found in local catalog.")
    return template("dashboard", pgs_id=pgs_id, info=normalize_local_info(info), cfg=cfg)


# ── API: stats ────────────────────────────────────────────────────────────────

@_bottle_app.route("/api/data/<pgs_id>")
@_require_pgs_id
def api_data(pgs_id: str):
    stats = stats_svc.get_stats(pgs_id)
    if "error" in stats:
        code = 404 if "not found" in stats["error"].lower() else 500
        return _json(stats, status=code)
    return _json(stats)


# ── Search (Mode B) ───────────────────────────────────────────────────────────

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
                results.append(normalize_search_result(r))

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


# ── PGS Source (Mode A + B) ───────────────────────────────────────────────────

@_bottle_app.route("/pgs/<pgs_id>/source")
@_require_pgs_id
def pgs_source(pgs_id: str):
    raw_local = catalog.get_pgs_info(pgs_id) if catalog.exists_locally(pgs_id) else None
    local_info = normalize_local_info(raw_local) if raw_local else None
    remote_info = None
    remote_error = None

    if cfg.REMOTE_ENABLED:
        result, _ = pgscat_client.get_score_with_cache(pgs_id)
        if "error" in result:
            remote_error = result["error"]
        else:
            remote_info = normalize_remote_info(result)

    if local_info and remote_info:
        source = "both"
    elif local_info:
        source = "local_only"
    elif remote_info:
        source = "remote_only"
    else:
        source = "not_found"

    # For remote-only: also produce manifest for preparation UI
    plan = None
    manifest = None
    if source in ("remote_only", "not_found"):
        plan = inspector.get_pipeline_plan(
            pgs_id,
            flow_mode="catalog_remote_prepare",
            local_exists=False,
            remote_info=remote_info,
        )
        manifest = preparer.build_manifest(pgs_id, remote_info=remote_info)

    return template(
        "pgs_remote",
        pgs_id=pgs_id,
        source=source,
        local_info=local_info,
        remote_info=remote_info,
        remote_error=remote_error,
        plan=plan,
        manifest=manifest,
        cfg=cfg,
    )


# ── Pipeline plan HTML view ───────────────────────────────────────────────────

@_bottle_app.route("/pipeline/<pgs_id>/plan")
@_require_pgs_id
def pipeline_plan_html(pgs_id: str):
    local_exists = catalog.exists_locally(pgs_id)
    remote_info = None

    if cfg.REMOTE_ENABLED:
        result, _ = pgscat_client.get_score_with_cache(pgs_id)
        if "error" not in result:
            remote_info = normalize_remote_info(result)

    flow_mode = "local_existing" if local_exists else "catalog_remote_prepare"
    plan = inspector.get_pipeline_plan(
        pgs_id,
        flow_mode=flow_mode,
        local_exists=local_exists,
        remote_info=remote_info,
    )
    manifest = preparer.build_manifest(pgs_id, remote_info=remote_info)
    local_info = normalize_local_info(catalog.get_pgs_info(pgs_id)) if local_exists else None

    return template(
        "pipeline_plan",
        pgs_id=pgs_id,
        plan=plan,
        manifest=manifest,
        local_info=local_info,
        remote_info=remote_info,
        local_exists=local_exists,
        cfg=cfg,
    )


# ── API: Pipeline plan JSON ───────────────────────────────────────────────────

@_bottle_app.route("/api/pipeline/<pgs_id>/plan")
@_require_pgs_id
def api_pipeline_plan(pgs_id: str):
    local_exists = catalog.exists_locally(pgs_id)
    remote_info = None
    if cfg.REMOTE_ENABLED:
        result, _ = pgscat_client.get_score_with_cache(pgs_id)
        if "error" not in result:
            remote_info = result

    flow_mode = "local_existing" if local_exists else "catalog_remote_prepare"
    plan = inspector.get_pipeline_plan(
        pgs_id,
        flow_mode=flow_mode,
        local_exists=local_exists,
        remote_info=remote_info,
    )
    if "error" in plan:
        return _json(plan, status=400)
    return _json(plan)


# ── API: Preparation manifest JSON ───────────────────────────────────────────

@_bottle_app.route("/api/pipeline/<pgs_id>/manifest")
@_require_pgs_id
def api_pipeline_manifest(pgs_id: str):
    remote_info = None
    if cfg.REMOTE_ENABLED:
        result, _ = pgscat_client.get_score_with_cache(pgs_id)
        if "error" not in result:
            remote_info = result

    manifest = preparer.build_manifest(pgs_id, remote_info=remote_info)
    if "error" in manifest:
        return _json(manifest, status=400)
    return _json(manifest)


# ── Custom PRS (Mode C) ───────────────────────────────────────────────────────

@_bottle_app.route("/custom-prs/new")
def custom_prs_new():
    return template("custom_prs_new", cfg=cfg)


@_bottle_app.route("/api/custom-prs/plan", method="POST")
def api_custom_prs_plan():
    """
    Generate a pipeline plan for a user-defined PRS.
    Accepts JSON or form-encoded body.
    No jobs are submitted or files created.
    """
    ct = request.content_type or ""
    if "application/json" in ct:
        try:
            body = request.json or {}
        except Exception:
            body = {}
    else:
        body = dict(request.forms)

    # Validate required fields
    analysis_name = (body.get("analysis_name") or "").strip()
    pgs_id_raw = (body.get("pgs_id") or "").strip().upper()
    input_type = (body.get("input_type") or "scoring_file").strip()
    genome_build = (body.get("genome_build") or "GRCh38").strip()
    source_file = (body.get("source_file") or "").strip()
    needs_betamap = str(body.get("needs_betamap", "true")).lower() in ("true", "1", "yes")
    chromosomes = body.get("chromosomes") or list(range(1, 23))
    output_dir = (body.get("output_dir") or "").strip()
    notes = (body.get("notes") or "").strip()

    if not analysis_name:
        return _json({"error": "analysis_name is required"}, status=400)

    # Use provided pgs_id or synthesise a placeholder
    if pgs_id_raw and validate_pgs_id(pgs_id_raw):
        pgs_id = pgs_id_raw
    else:
        # For custom analyses, require a valid PGS ID format or default to PGS000000
        pgs_id = pgs_id_raw if pgs_id_raw else "PGS000000"
        if not validate_pgs_id(pgs_id):
            pgs_id = "PGS000000"

    scores_dir = str(cfg.scores_dir)
    scripts_dir = str(cfg.SCRIPTS_DIR)
    out_dir = output_dir or f"{scores_dir}/{pgs_id}"

    steps = inspector.get_pipeline_plan(
        pgs_id, flow_mode="custom_prs"
    ).get("steps", [])

    # Compute suggested resources
    suggested_resources = {
        "partition": "gpu",
        "account": "researchers",
        "qos": "vip",
        "nodes": 1,
        "cpus_per_task": 16,
        "mem": "120G",
        "gres": "gpu:1",
        "time": "06:00:00",
        "note": "Use CPU fallback (compute_prs.py) if no GPU is available.",
    }

    return _json({
        "analysis_name": analysis_name,
        "pgs_id": pgs_id,
        "input_type": input_type,
        "genome_build": genome_build,
        "source_file": source_file or f"{out_dir}/{pgs_id}_hmPOS_GRCh38.txt.gz",
        "needs_betamap": needs_betamap,
        "chromosomes": chromosomes,
        "output_dir": out_dir,
        "betamap_columns": ["PRS_ID", "CHROM", "POS", "ID",
                             "EFFECT_ALLELE", "OTHER_ALLELE", "BETA", "IS_FLIP"],
        "expected_files": [
            f"{pgs_id}_hmPOS_GRCh38.txt.gz",
            f"{pgs_id}_hmPOS_GRCh38.columns.json",
            f"{pgs_id}_hmPOS_GRCh38.betamap.tsv.gz",
            f"{pgs_id}_hmPOS_GRCh38.betamap.tsv.gz.meta.json",
            *[f"{pgs_id}_chr{c}_scores.tsv" for c in chromosomes],
            *[f"{pgs_id}_chr{c}_metadata.json" for c in chromosomes],
            f"{pgs_id}_PRS_total.tsv",
            f"{pgs_id}_PRS_total_metadata.json",
            f"{pgs_id}_PRS_total.parquet",
        ],
        "zarr_base": "{GENOTYPES_DIR}/chr{N}.zarr",
        "zarr_index": "{GENOTYPES_DIR}/zarr_index.duckdb",
        "scripts_dir": scripts_dir,
        "suggested_resources": suggested_resources,
        "steps": steps,
        "notes": notes,
        "note": (
            "This plan is informational only. "
            "No files have been created and no jobs have been submitted."
        ),
    })


# ── API: Execute score preparation ───────────────────────────────────────────

@_bottle_app.route("/api/pgs/<pgs_id>/prepare", method="POST")
@_require_pgs_id
def api_prepare_score(pgs_id: str):
    """
    Execute input-preparation steps for pgs_id:
      1. Create output directory
      2. Download catalog metadata  → {PGS_ID}.json
      3. Download harmonised file   → {PGS_ID}_hmPOS_GRCh38.txt.gz
      4. Build betamap              → .betamap.tsv.gz + .columns.json + .meta.json

    Returns JSON:
      { status: "ok"|"partial"|"error",
        steps:  [{name, status, message, elapsed_s}, ...],
        current_stage, ready_for_compute }

    Already-completed steps are skipped (idempotent).
    Gunicorn timeout is 120 s; betamap build is pure-Python so it completes
    well within that limit for typical scoring files.
    """
    logger.info("Preparation requested for %s by %s", pgs_id, _client_ip())
    try:
        result = preparer.prepare_score(pgs_id)
    except Exception as exc:
        logger.exception("Unexpected error preparing %s", pgs_id)
        return _json({"status": "error", "error": str(exc), "steps": []}, status=500)

    http_status = 200 if result.get("status") == "ok" else 207
    return _json(result, status=http_status)


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


# ── Variant Annotator (Mode D) ────────────────────────────────────────────────

@_bottle_app.route("/variants/<pgs_id>")
@_require_pgs_id
def variants_view(pgs_id: str):
    """
    HTML view for variant annotation results.
    Shows annotation status (annotated / running / not_annotated) and,
    when available, summary statistics and an interactive variant table.
    Does NOT trigger annotation — the pipeline runs externally via Apptainer.
    """
    status_info = variant_svc.get_status(pgs_id)
    return template(
        "variants",
        pgs_id=pgs_id,
        status_info=status_info,
        cfg=cfg,
    )


@_bottle_app.route("/api/variants/<pgs_id>")
@_require_pgs_id
def api_variants(pgs_id: str):
    """
    Paginated JSON table of annotated variants.

    Query parameters:
      page        int  (default 1)
      page_size   int  (default 500, max 5000)
      chrom       str  filter by chromosome (bare, e.g. "1", "X", "MT")
      region_class str  filter by region class
      gene_name   str  filter by gene symbol (case-insensitive exact match)
      only_coding bool filter to coding variants only
    """
    page               = int(request.query.get("page",        1))
    page_size          = int(request.query.get("page_size", 500))
    chrom              = (request.query.get("chrom")              or "").strip() or None
    region_class       = (request.query.get("region_class")       or "").strip() or None
    gene_name          = (request.query.get("gene_name")          or "").strip() or None
    only_coding        = request.query.get("only_coding", "").lower() in ("1", "true", "yes")
    clinical_confidence= (request.query.get("clinical_confidence") or "").strip().lower() or None
    add_ranking        = request.query.get("add_ranking", "").lower() in ("1", "true", "yes")

    if clinical_confidence and clinical_confidence not in ("high", "medium", "low"):
        clinical_confidence = None

    result = variant_svc.get_variants(
        pgs_id,
        page=page,
        page_size=page_size,
        chrom=chrom,
        region_class=region_class,
        gene_name=gene_name,
        only_coding=only_coding,
        clinical_confidence=clinical_confidence,
        add_ranking=add_ranking,
    )
    if "error" in result:
        code = 404 if "not found" in result["error"].lower() else 500
        return _json(result, status=code)
    return _json(result)


@_bottle_app.route("/api/variants/<pgs_id>/summary")
@_require_pgs_id
def api_variants_summary(pgs_id: str):
    """
    Full annotation summary JSON including per-chromosome breakdown,
    region-class counts, and run metadata.
    """
    result = variant_svc.get_summary(pgs_id)
    if "error" in result:
        return _json(result, status=404)
    return _json(result)


@_bottle_app.route("/api/variants/<pgs_id>/ideogram")
@_require_pgs_id
def api_variants_ideogram(pgs_id: str):
    """
    Chromosome-level variant positions and counts for ideogram visualisation.

    Contract (stable; implements the ideogram data schema):
    {
      "pgs_id": str,
      "chromosomes": [
        {
          "chrom": "1",
          "n_variants": int,
          "n_coding": int,
          "n_intergenic": int,
          "n_regulatory": int,
          "region_class_counts": {"coding": n, ...},
          "positions_coding":    [int, ...],   # capped at 2000
          "positions_intronic":  [int, ...],
          "positions_intergenic":[int, ...],
        }
      ],
      "chrom_order": ["1",...,"22","X","Y","MT"],
      "note": str
    }
    """
    result = variant_svc.get_ideogram_data(pgs_id)
    if "error" in result:
        return _json(result, status=404)
    return _json(result)


@_bottle_app.route("/api/variants/<pgs_id>/charts")
@_require_pgs_id
def api_variants_charts(pgs_id: str):
    """
    Pre-aggregated analytical chart data for the /variants dashboard.

    Reads the annotation parquet once and returns all 6 chart datasets:
      cadd_by_consequence  G1 – CADD phred arrays by consequence category
      af_spectrum          G2 – log₁₀(AF) histogram bins
      beta_vs_ranking      G3 – |BETA| × ranking_score scatter
      gnomad_vs_mcps       G4 – gnomAD AF × MCPS AF scatter
      cumulative_beta      G7 – cumulative |BETA| fraction curve
      manhattan            G8 – genomic position × ranking_score

    No query parameters required.
    """
    result = variant_svc.get_charts_data(pgs_id)
    if "error" in result:
        code = 404 if "not found" in result["error"].lower() else 500
        return _json(result, status=code)
    return _json(result)


@_bottle_app.route("/api/variants/<pgs_id>/run-info")
@_require_pgs_id
def api_variants_run_info(pgs_id: str):
    """
    Return the Apptainer run command and paths needed to annotate pgs_id.
    Useful for the UI 'how to run' panel.
    """
    return _json(variant_svc.get_run_command(pgs_id))


@_bottle_app.route("/api/variants/ranked")
def api_variants_ranked():
    """
    Return variants ranked by clinical priority score (descending).

    Ranking formula:
        score = 0.25*rarity + 0.30*consequence + 0.20*cadd_norm
              + 0.10*revel_norm + 0.10*clinical_bonus + 0.05*lof_bonus

    All components normalised to [0, 1].

    Query parameters:
        pgs_id        str   required — PGS ID to rank
        page          int   (default 1)
        page_size     int   (default 500, max 5000)
        min_score     float filter: only variants with score >= min_score
        clinical_only bool  filter: only variants in clinical genes
    """
    pgs_id        = (request.query.get("pgs_id") or "").strip().upper() or None
    page          = int(request.query.get("page",      1))
    page_size     = int(request.query.get("page_size", 500))
    min_score_raw = request.query.get("min_score", "0")
    clinical_only = request.query.get("clinical_only", "").lower() in ("1", "true", "yes")

    try:
        min_score = float(min_score_raw)
    except (TypeError, ValueError):
        min_score = 0.0

    if not pgs_id:
        return _json({"error": "pgs_id query parameter is required."}, status=400)

    if not validate_pgs_id(pgs_id):
        return _json({"error": f"Invalid PGS ID {pgs_id!r}."}, status=400)

    result = variant_svc.get_ranked_variants(
        pgs_id,
        page=page,
        page_size=page_size,
        min_score=min_score,
        clinical_only=clinical_only,
    )
    if "error" in result:
        code = 404 if "not found" in result.get("error", "").lower() else 500
        return _json(result, status=code)
    return _json(result)


@_bottle_app.route("/api/genes/summary")
def api_genes_summary():
    """
    Gene heatmap data (Heatmap 1: genes × metrics).

    Returns the top-N filtered/sorted genes plus a Plotly-ready heatmap matrix.
    Used by the heatmap panel in /genes.

    Query parameters
    ----------------
    top_n               int   max genes in heatmap (default 60, max 200)
    pgs_id              str   restrict to one PGS
    sort_by             str   (same options as /api/genes)
    sort_dir            str
    clinical_only       bool
    clinical_confidence str
    min_variants        int
    """
    top_n               = int(request.query.get("top_n",        60))
    top_n               = max(1, min(top_n, 200))
    pgs_id_filter       = (request.query.get("pgs_id")              or "").strip().upper() or None
    sort_by             = (request.query.get("sort_by")             or "sum_abs_beta").strip()
    sort_dir            = (request.query.get("sort_dir")            or "desc").strip()
    clinical_confidence = (request.query.get("clinical_confidence") or "").strip().lower() or None
    clinical_only       = request.query.get("clinical_only", "").lower() in ("1", "true", "yes")
    min_variants        = int(request.query.get("min_variants", 1))

    if pgs_id_filter and not validate_pgs_id(pgs_id_filter):
        return _json({"error": f"Invalid PGS ID {pgs_id_filter!r}."}, status=400)
    if clinical_confidence and clinical_confidence not in ("high", "medium", "low"):
        clinical_confidence = None
    if sort_by not in ("ranking_score", "variant_count", "sum_abs_beta", "mean_cadd", "alphabetical"):
        sort_by = "sum_abs_beta"

    result = gene_svc.get_gene_summary(
        limit=top_n,
        pgs_id=pgs_id_filter,
        sort_by=sort_by,
        sort_dir=sort_dir,
        clinical_only=clinical_only,
        clinical_confidence=clinical_confidence,
        min_variants=min_variants,
    )
    return _json(result)


@_bottle_app.route("/api/gene/<gene_symbol>/pgs_matrix")
def api_gene_pgs_matrix(gene_symbol: str):
    """
    Gene × PGS matrix (Heatmap 2): per-PGS metrics for a given gene.

    For each PGS that contains variants in this gene, returns:
      n_variants, coding_count, splice_count, sum_abs_beta,
      mean_cadd, max_ranking_score, mean_ranking_score

    Used by the Gene × PGS comparison panel in /gene/<symbol>.

    Query parameters
    ----------------
    min_variants  int  minimum variant count to include a PGS (default 1)
    sort_by       str  sum_abs_beta | max_ranking_score | variant_count
    """
    gene_symbol  = gene_symbol.strip().upper()
    min_variants = int(request.query.get("min_variants", 1))
    sort_by      = (request.query.get("sort_by") or "sum_abs_beta").strip()

    if sort_by not in ("sum_abs_beta", "max_ranking_score", "variant_count"):
        sort_by = "sum_abs_beta"

    try:
        result = gene_svc.get_gene_pgs_matrix(
            gene_symbol,
            min_variants=min_variants,
            sort_by=sort_by,
        )
    except Exception as exc:
        return _json({"error": str(exc), "gene_symbol": gene_symbol, "total_pgs": 0}, status=500)
    if result["total_pgs"] == 0:
        return _json(result, status=404)
    return _json(result)


@_bottle_app.route("/api/prs/clinical")
def api_prs_clinical():
    """
    PRS × Clinical genes: genetic load and ranking scores for clinical genes
    that have at least one PRS variant annotated.

    No query parameters required.

    Response:
        n_clinical_genes_with_variants  int
        total_genetic_load              float  (sum |BETA| across all clinical variants)
        genes[]                         sorted by genetic_load descending:
          gene_name, n_variants, genetic_load, cumulative_score,
          mean_ranking_score, confidence, moi, disease
    """
    result = gene_svc.get_prs_clinical_panel()
    return _json(result)


# ── Gene Browser (Mode E) ─────────────────────────────────────────────────────
# Gene-centric interpretation browser — Plotly stacked genomic tracks.
# Uses annotations produced by the Variant Annotator pipeline (parquet).
# MTR (Missense Tolerance Ratio w31) track is intentionally NOT implemented.

@_bottle_app.route("/genes")
def genes_list():
    """
    Gene list / search page.
    GET /genes                        → list all annotated genes (global)
    GET /genes?pgs_id=PGS000004       → list genes only from that PGS
    GET /genes?q=GENE                 → search within current scope
    GET /genes?clinical_only=1        → filter to clinical genes only
    """
    query         = (request.query.get("q")          or "").strip().upper()
    clinical_only = request.query.get("clinical_only", "") in ("1", "true", "yes")
    pgs_id_filter = (request.query.get("pgs_id")     or "").strip().upper() or None

    # Validate PGS ID if provided
    if pgs_id_filter and not validate_pgs_id(pgs_id_filter):
        pgs_id_filter = None

    all_genes = gene_svc.list_annotated_genes(pgs_id=pgs_id_filter)
    available_pgs_ids = gene_svc.list_available_pgs_ids()

    genes = []
    if query:
        genes = [g for g in all_genes if query in g["gene_symbol"].upper()]
        # If exact match, redirect directly to gene browser
        exact = [g for g in genes if g["gene_symbol"].upper() == query]
        if exact:
            from bottle import redirect
            redirect(f"/gene/{exact[0]['gene_symbol']}")

    if clinical_only:
        if genes:
            genes = [g for g in genes if g.get("is_clinical_gene")]
        else:
            genes = [g for g in all_genes if g.get("is_clinical_gene")]

    return template(
        "genes",
        query=query,
        genes=genes,
        all_genes=all_genes,
        clinical_only=clinical_only,
        pgs_id_filter=pgs_id_filter,
        available_pgs_ids=available_pgs_ids,
        cfg=cfg,
    )


@_bottle_app.route("/gene/<gene_symbol>")
def gene_browser(gene_symbol: str):
    """
    Gene browser HTML view.
    Shows Plotly-based stacked genomic tracks for the requested gene.
    Tracks: Variants | Gene Model | Functional Summary
    Excluded: MTR (Missense Tolerance Ratio w31)
    """
    # Normalise gene symbol: uppercase
    gene_symbol = gene_symbol.strip().upper()
    gene_info   = gene_svc.get_gene_info(gene_symbol)
    return template("gene", gene_symbol=gene_symbol, gene_info=gene_info, cfg=cfg)


@_bottle_app.route("/api/gene/<gene_symbol>")
def api_gene_info(gene_symbol: str):
    """
    Gene summary JSON: genomic range, variant counts, PGS IDs, rarity distribution.
    """
    gene_symbol = gene_symbol.strip().upper()
    result = gene_svc.get_gene_info(gene_symbol)
    if not result.get("found"):
        return _json(result, status=404)
    return _json(result)


@_bottle_app.route("/api/gene/<gene_symbol>/tracks")
def api_gene_tracks(gene_symbol: str):
    """
    Track data JSON for Plotly subplots.

    Returns:
      - Track 1: Variants scatter (consequence priority vs. position)
      - Track 2: Gene model (exon blocks from coding variant clusters)
      - Track 3: Functional summary (consequence distribution by genomic bin)
      - Excluded: MTR (Missense Tolerance Ratio w31) — not implemented by design
    """
    gene_symbol = gene_symbol.strip().upper()
    result = gene_svc.get_tracks(gene_symbol)
    if not result.get("found"):
        return _json(result, status=404)
    return _json(result)


@_bottle_app.route("/api/gene/<gene_symbol>/variants")
def api_gene_variants(gene_symbol: str):
    """
    Paginated JSON table of annotated variants in a gene.
    Includes dbSNP population frequency fields (rsid, af_global, rarity_class)
    when available in the annotation parquets.

    Query parameters:
      page       int  (default 1)
      page_size  int  (default 200, max 2000)
      format     str  "json" (default) | "tsv" (download)
    """
    gene_symbol = gene_symbol.strip().upper()
    page      = int(request.query.get("page", 1))
    page_size = int(request.query.get("page_size", 200))
    fmt       = (request.query.get("format") or "json").strip().lower()
    page_size = max(1, min(page_size, 2000))

    result = variant_svc.get_variants(
        pgs_id=None,      # signal: gene-mode query
        page=page,
        page_size=page_size,
        gene_name=gene_symbol,
        _scan_all=True,   # scan across all annotation dirs
    )

    if "error" in result:
        code = 404 if "not found" in result["error"].lower() else 500
        return _json(result, status=code)

    if fmt == "tsv":
        import io
        response.content_type = "text/tab-separated-values"
        response.headers["Content-Disposition"] = (
            f'attachment; filename="{gene_symbol}_variants.tsv"'
        )
        cols = result.get("columns", [])
        rows = result.get("rows", [])
        lines = ["\t".join(str(c) for c in cols)]
        for row in rows:
            lines.append("\t".join(str(row.get(c, "")) for c in cols))
        return "\n".join(lines)

    return _json(result)


@_bottle_app.route("/api/genes")
def api_genes_list():
    """
    Paginated, filterable list of annotated genes with per-gene metrics.

    Query parameters
    ----------------
    q                   str   symbol substring filter (case-insensitive)
    pgs_id              str   restrict to one PGS (optional; global if omitted)
    page                int   (default 1)
    page_size           int   (default 50, max 500)
    sort_by             str   ranking_score | variant_count | sum_abs_beta |
                              mean_cadd | alphabetical  (default: sum_abs_beta)
    sort_dir            str   asc | desc  (default: desc)
    clinical_only       bool  only clinical genes
    clinical_confidence str   high | medium | low
    min_variants        int   (default 1)

    Response
    --------
    {
      pgs_id, total_genes, total_filtered, page, page_size, total_pages,
      genes[]:  [{ gene_name, n_variants, sum_abs_beta, mean_cadd,
                   mean_af, is_clinical_gene, consequence_counts,
                   ranking_score_mean }]
    }
    """
    q                   = (request.query.get("q")                   or "").strip().upper()
    pgs_id_filter       = (request.query.get("pgs_id")              or "").strip().upper() or None
    sort_by             = (request.query.get("sort_by")             or "sum_abs_beta").strip()
    sort_dir            = (request.query.get("sort_dir")            or "desc").strip()
    clinical_confidence = (request.query.get("clinical_confidence") or "").strip().lower() or None
    clinical_only       = request.query.get("clinical_only", "").lower() in ("1", "true", "yes")
    page                = int(request.query.get("page",        1))
    page_size           = int(request.query.get("page_size",  50))
    min_variants        = int(request.query.get("min_variants", 1))
    page_size           = max(1, min(page_size, 500))

    if pgs_id_filter and not validate_pgs_id(pgs_id_filter):
        pgs_id_filter = None
    if clinical_confidence and clinical_confidence not in ("high", "medium", "low"):
        clinical_confidence = None
    if sort_by not in ("ranking_score", "variant_count", "sum_abs_beta", "mean_cadd", "alphabetical"):
        sort_by = "sum_abs_beta"

    result = gene_svc.get_gene_summary(
        pgs_id=pgs_id_filter,
        sort_by=sort_by,
        sort_dir=sort_dir,
        clinical_only=clinical_only,
        clinical_confidence=clinical_confidence,
        min_variants=min_variants,
        page=page,
        page_size=page_size,
        q=q,
    )

    # Drop heatmap from gene-list endpoint (keep response lean)
    result.pop("heatmap", None)
    return _json(result)


# ── External source lookup endpoints ─────────────────────────────────────────

@_bottle_app.route("/api/variant/lookup")
def api_variant_lookup():
    """
    Live per-variant lookup from gnomAD 4.1.1, ClinVar, and MCPS.

    Query parameters
    ----------------
    chrom   str  chromosome (1-22, X, Y — with or without chr prefix)
    pos     int  1-based genomic position
    ref     str  reference allele
    alt     str  alternate allele

    Response
    --------
    {
      "gnomad":   { af_gnomad, af_gnomad_amr, ... } | null,
      "clinvar":  { clnsig, clndn, clnrevstat, alleleid } | null,
      "mcps":     { af_mcps, af_mcps_mex, ... } | null,
    }
    """
    chrom = (request.query.get("chrom") or "").strip()
    pos_s = (request.query.get("pos")   or "").strip()
    ref   = (request.query.get("ref")   or "").strip().upper()
    alt   = (request.query.get("alt")   or "").strip().upper()

    if not chrom or not pos_s or not ref or not alt:
        response.status = 400
        return _json({"error": "chrom, pos, ref, alt are required."})

    try:
        pos = int(pos_s)
    except ValueError:
        response.status = 400
        return _json({"error": "pos must be an integer."})

    result: dict = {"chrom": chrom, "pos": pos, "ref": ref, "alt": alt}

    # gnomAD
    try:
        from services.gnomad_fetcher import get_gnomad_fetcher
        gf = get_gnomad_fetcher()
        result["gnomad"] = gf.fetch(chrom, pos, ref, alt) if gf.is_available() else None
    except Exception as exc:
        logger.warning("variant/lookup gnomAD error: %s", exc)
        result["gnomad"] = None

    # ClinVar
    try:
        from services.clinvar_fetcher import get_clinvar_fetcher
        cf = get_clinvar_fetcher(work_dir=cfg.WORK_DIR)
        result["clinvar"] = cf.fetch(chrom, pos, ref, alt) if (cf.is_available() and cf._parquet.exists()) else None
    except Exception as exc:
        logger.warning("variant/lookup ClinVar error: %s", exc)
        result["clinvar"] = None

    # MCPS
    try:
        from services.mcps_fetcher import get_mcps_fetcher
        mf = get_mcps_fetcher(work_dir=cfg.WORK_DIR)
        result["mcps"] = mf.fetch(chrom, pos, ref, alt) if (mf.is_available() and mf.index_ready()) else None
    except Exception as exc:
        logger.warning("variant/lookup MCPS error: %s", exc)
        result["mcps"] = None

    return _json(result)


@_bottle_app.route("/api/admin/build-indexes", method="POST")
def api_build_indexes():
    """
    Trigger one-time build of derived indexes (ClinVar Parquet, MCPS DuckDB).

    POST /api/admin/build-indexes
    Body (JSON, optional): {"sources": ["clinvar", "mcps"]}

    Runs synchronously — may take 5-30 minutes. Call from HPC or admin only.
    """
    import json as _json_mod
    try:
        body = _json_mod.loads(request.body.read() or b"{}") or {}
    except Exception:
        body = {}
    sources = body.get("sources", ["clinvar", "mcps"])

    results: dict = {}

    if "clinvar" in sources:
        try:
            from services.clinvar_fetcher import get_clinvar_fetcher
            cf = get_clinvar_fetcher(work_dir=cfg.WORK_DIR)
            ok = cf.build_index_if_needed()
            results["clinvar"] = "built" if ok else "failed"
        except Exception as exc:
            results["clinvar"] = f"error: {exc}"

    if "mcps" in sources:
        try:
            from services.mcps_fetcher import get_mcps_fetcher
            mf = get_mcps_fetcher(work_dir=cfg.WORK_DIR)
            ok = mf.build_index()
            results["mcps"] = "built" if ok else "failed"
        except Exception as exc:
            results["mcps"] = f"error: {exc}"

    return _json({"status": "done", "results": results})


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
