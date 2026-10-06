# Architecture

## Data flow overview

```
Browser â†’ Nginx/Caddy (TLS, auth) â†’ Gunicorn:8080 (4 sync workers)
                                         â”‚
                          â”Œâ”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”¼â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”
                          â”‚              â”‚              â”‚
                    Bottle routes   FileCache      DBPool
                          â”‚         /work/cache/   psycopg2
                          â”‚              â”‚         PostgreSQL 16
                     DuckDB (in-mem)     â”‚         10.x.x.x (red interna)
                     /data/*.parquet     â”‚         (mcps-net)
                     /data/*.tsv         â”‚
                          â”‚              â”‚
                     pgscat lib â”€â”€â†’ PGS Catalog API
                     (remote, optional)
```

## Storage boundaries

### What goes in DuckDB / Parquet (analytical, read-only)

| Data | Location |
|------|----------|
| Per-sample PRS scores (140k rows Ã— 24 cols) | `/data/{PGS_ID}/{PGS_ID}_PRS_total.parquet` |
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
| `sync_runs` | History of filesystemâ†’DB sync operations |

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

### `services/db.py` â€” DBPool
- `psycopg2.pool.ThreadedConnectionPool` (1â€“5 connections per gunicorn worker)
- Graceful degradation: if PostgreSQL is down at startup, `db.available = False` and all `execute/fetchone/fetchall` return empty/False without raising.
- `db.retry_connect()` is safe to call from health checks.
- All DML helpers use `%s` parameterisation â€” never string interpolation.

### `services/cache.py` â€” FileCache
- Writes JSON to `/work/cache/{key}.json` atomically (write to `.tmp`, then `os.replace`).
- `get_with_mtime_check(key, mtime)`: returns None if stored `_source_mtime` differs from current file mtime. Prevents serving stale stats after pipeline re-runs.
- TTL (wall-clock age of the cache file) is a secondary guard.
- Safe for concurrent gunicorn workers (`os.replace` is atomic on Linux).

### `services/sync_service.py` â€” SyncService
- Runs `catalog.list_pgs()` and upserts each entry into `local_pgs_index` + `file_inventory`.
- Launched in a **daemon background thread** 8s after worker startup.
- Idempotent: safe to run from multiple workers simultaneously (upserts on primary keys).
- If DB unavailable, returns immediately without error.

### `services/health.py` â€” HealthService
- `/healthz`: always 200 (liveness â€” process is running).
- `/readyz`: 200 if data dir accessible (even if DB is degraded), 503 if data dir missing.
- Tries `db.retry_connect()` on each readiness check to recover from transient DB outages.

### `services/parquet_stats.py` â€” ParquetStats
Cache hierarchy:
```
1. FileCache.get_with_mtime_check()      â† L1, fast, per-worker
2. DBPool.get_stats_cache(mtime)         â† L2, persistent, shared across restarts
3. DuckDB scan of parquet/tsv            â† compute (slow first time, ~20s on cold CephFS)
```
On compute, writes to both L1 and L2.

### `services/pgscat_client.py` â€” PGSCatClient
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
1. Config.__init__()           â€“ read all env vars, resolve DB password
2. DBPool.__init__()           â€“ attempt PostgreSQL connection (5s timeout)
3. FileCache.__init__()        â€“ create /work/cache/ if missing
4. LocalCatalog.__init__()     â€“ no I/O (lazy)
5. ParquetStats.__init__()     â€“ inject cache + db
6. PGSCatClient.__init__()     â€“ import pgscat library (may warn if not installed)
7. SyncService.__init__()      â€“ no I/O
8. gunicorn forks N workers    â€“ each worker re-runs steps 1-7
9. sync_svc.start_background() â€“ per-worker daemon thread, sleeps 8s then syncs
```

If PostgreSQL is unavailable at step 2, the app still starts. Dashboards serve from filesystem. Sync is skipped. Search audit is skipped. Cache layer 1 (file) still works.

---

## Deployment topology

```
CephFS (hot_nvme)
  â””â”€â”€ /pgscatalog/scores/     â† mounted ro at /data inside container
  â””â”€â”€ /postgres/data/         â† PostgreSQL volume (separate container)

Host (Linux, Podman)
  â”œâ”€â”€ pgs-dashboard container â”€â”
  â”‚   port 8080               â”‚  mcps-net (Podman bridge)
  â”œâ”€â”€ postgres container â”€â”€â”€â”€â”€â”€â”˜  10.x.x.x (red interna)
  â””â”€â”€ Nginx/Caddy             â† terminates TLS for pgscat.mcps-epcm.org
```

---

## Security model (Phase 1)

- All data access is read-only (`/data` mounted `:ro`).
- `pgs_id` validated against `^PGS\d{6}$` before any filesystem/DB access.
- `pipeline_inspector.read_script()` uses a strict allowlist â€” no path traversal.
- DuckDB paths are SQL-escaped (single-quote doubling) before embedding in strings.
- DB credentials loaded from file (preferred) or env var â€” never in code.
- `/api/admin/sync` is not authenticated in Phase 1; restrict via proxy ACL or firewall.
- Authentication for the entire app should be handled at the reverse proxy level (Nginx `auth_basic` or OAuth2 proxy).
