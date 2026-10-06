# Architecture

## Data flow overview

```
Browser → Nginx/Caddy (TLS, auth) → Gunicorn:8080 (4 sync workers)
                                         │
                          ┌──────────────┼──────────────┐
                          │              │              │
                    Bottle routes   FileCache      DBPool
                          │         /work/cache/   psycopg2
                          │              │         PostgreSQL 16
                     DuckDB (in-mem)     │         10.89.0.11
                     /data/*.parquet     │         (mcps-net)
                     /data/*.tsv         │
                          │              │
                     pgscat lib ──→ PGS Catalog API
                     (remote, optional)
```

## Storage boundaries

### What goes in DuckDB / Parquet (analytical, read-only)

| Data | Location |
|------|----------|
| Per-sample PRS scores (140k rows × 24 cols) | `/data/{PGS_ID}/{PGS_ID}_PRS_total.parquet` |
| Per-chromosome scores | `/data/{PGS_ID}/{PGS_ID}_chr{N}_scores.tsv` |
| Descriptive stats, histograms | Computed on-demand by DuckDB, cached in `/work/cache/` |
| Input weights (betamap) | `/data/{PGS_ID}/*.betamap.tsv.gz` |

DuckDB opens a new `:memory:` connection per request and closes it in a `finally` block. No persistent DuckDB database file. No writes to `/data`.

### What goes in PostgreSQL (metadata & operational)

| Table | Contents |
|-------|----------|
| `local_pgs_index` | Filesystem scan results: has_parquet, has_tsv, chromosomes, trait |
| `remote_pgs_cache` | Normalised pgscat API responses, TTL=24h |
| `stats_cache` | DuckDB-computed aggregates (NOT raw scores), mtime-keyed |
| `file_inventory` | Per-file inventory: name, type, size |
| `search_audit` | Query, search_type, n_results, duration_ms, IP |
| `pipeline_plans` | Generated plans (status: generated/manifest_written/submitted) |
| `sync_runs` | History of filesystem→DB sync operations |

### What NEVER goes in PostgreSQL

- `PRS_total` values per sample
- `call_DS` dosages
- Any genomic variant-level data
- Parquet file contents
- Betamap coefficients

---

## Component details

### `config.py`
Single source of truth. Every env var is read exactly once at startup and stored as a typed attribute. Services receive `cfg` by reference.

### `services/db.py` — DBPool
- `psycopg2.pool.ThreadedConnectionPool` (1–5 connections per gunicorn worker)
- Graceful degradation: if PostgreSQL is down at startup, `db.available = False` and all `execute/fetchone/fetchall` return empty/False without raising.
- `db.retry_connect()` is safe to call from health checks.
- All DML helpers use `%s` parameterisation — never string interpolation.

### `services/cache.py` — FileCache
- Writes JSON to `/work/cache/{key}.json` atomically (write to `.tmp`, then `os.replace`).
- `get_with_mtime_check(key, mtime)`: returns None if stored `_source_mtime` differs from current file mtime. Prevents serving stale stats after pipeline re-runs.
- TTL (wall-clock age of the cache file) is a secondary guard.
- Safe for concurrent gunicorn workers (`os.replace` is atomic on Linux).

### `services/sync_service.py` — SyncService
- Runs `catalog.list_pgs()` and upserts each entry into `local_pgs_index` + `file_inventory`.
- Launched in a **daemon background thread** 8s after worker startup.
- Idempotent: safe to run from multiple workers simultaneously (upserts on primary keys).
- If DB unavailable, returns immediately without error.

### `services/health.py` — HealthService
- `/healthz`: always 200 (liveness — process is running).
- `/readyz`: 200 if data dir accessible (even if DB is degraded), 503 if data dir missing.
- Tries `db.retry_connect()` on each readiness check to recover from transient DB outages.

### `services/parquet_stats.py` — ParquetStats
Cache hierarchy:
```
1. FileCache.get_with_mtime_check()      ← L1, fast, per-worker
2. DBPool.get_stats_cache(mtime)         ← L2, persistent, shared across restarts
3. DuckDB scan of parquet/tsv            ← compute (slow first time, ~20s on cold CephFS)
```
On compute, writes to both L1 and L2.

### `services/pgscat_client.py` — PGSCatClient
- `get_score_with_cache(pgs_id)`: checks `remote_pgs_cache` before hitting API.
- `search_with_cache(query, type)`: caches individual results by `pgs_id`.
- Audit logging via `db.log_search()` (fire-and-forget, never raises).
- pgscat free-text limitation documented inline: only EFO IDs, PGS IDs, and PMIDs work reliably.

### WSGI proxy middleware (`app.py`)
```python
class _ProxyMiddleware:
    """Reads X-Forwarded-For/Proto/Host/Prefix from Nginx/Caddy."""
```
Applied when `APP_TRUSTED_PROXY=true`. Sets `REMOTE_ADDR`, `wsgi.url_scheme`, `HTTP_HOST`, and `SCRIPT_NAME` from upstream headers. This makes `request.url` and `request.environ["REMOTE_ADDR"]` correct behind a reverse proxy.

---

## Startup sequence

```
1. Config.__init__()           – read all env vars, resolve DB password
2. DBPool.__init__()           – attempt PostgreSQL connection (5s timeout)
3. FileCache.__init__()        – create /work/cache/ if missing
4. LocalCatalog.__init__()     – no I/O (lazy)
5. ParquetStats.__init__()     – inject cache + db
6. PGSCatClient.__init__()     – import pgscat library (may warn if not installed)
7. SyncService.__init__()      – no I/O
8. gunicorn forks N workers    – each worker re-runs steps 1-7
9. sync_svc.start_background() – per-worker daemon thread, sleeps 8s then syncs
```

If PostgreSQL is unavailable at step 2, the app still starts. Dashboards serve from filesystem. Sync is skipped. Search audit is skipped. Cache layer 1 (file) still works.

---

## Deployment topology

```
CephFS (hot_nvme)
  └── /pgscatalog/scores/     ← mounted ro at /data inside container
  └── /postgres/data/         ← PostgreSQL volume (separate container)

Host (Linux, Podman)
  ├── pgs-dashboard container ─┐
  │   port 8080               │  mcps-net (Podman bridge)
  ├── postgres container ──────┘  10.89.0.11
  └── Nginx/Caddy             ← terminates TLS for pgscat.mcps-epcm.org
```

---

## Security model (Phase 1)

- All data access is read-only (`/data` mounted `:ro`).
- `pgs_id` validated against `^PGS\d{6}$` before any filesystem/DB access.
- `pipeline_inspector.read_script()` uses a strict allowlist — no path traversal.
- DuckDB paths are SQL-escaped (single-quote doubling) before embedding in strings.
- DB credentials loaded from file (preferred) or env var — never in code.
- `/api/admin/sync` is not authenticated in Phase 1; restrict via proxy ACL or firewall.
- Authentication for the entire app should be handled at the reverse proxy level (Nginx `auth_basic` or OAuth2 proxy).
