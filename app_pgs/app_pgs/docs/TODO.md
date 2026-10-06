# TODO / Roadmap

## Phase 1 – Completed ✓

- [x] Local catalog scanner (filesystem)
- [x] DuckDB stats + histogram (in-memory, no pandas)
- [x] Robust score column detection
- [x] pgscat Python library + CLI fallback
- [x] Dashboard with Plotly histogram
- [x] Pipeline plan generator (read-only, no sbatch)
- [x] Containerfile for Podman
- [x] `pgs/<id>/source` local vs remote page
- [x] PostgreSQL schema (`sql/schema.sql`)
- [x] DBPool with graceful degradation
- [x] FileCache with mtime-based invalidation (2-layer: file + DB)
- [x] SyncService (filesystem → PostgreSQL, background thread)
- [x] HealthService (`/healthz`, `/readyz`)
- [x] Proxy middleware (`X-Forwarded-*`)
- [x] Search audit logging
- [x] Remote metadata caching in PostgreSQL
- [x] `/api/admin/sync` manual trigger

---

## Phase 2 – Near term

### Parquet export in pipeline
- [ ] Add DuckDB COPY step to `compute_prs_spark_gpu.py` after aggregation
  ```python
  # After PRS_total.tsv is written:
  import duckdb
  duckdb.execute("""
      COPY (SELECT * FROM read_csv_auto(?, delim='\t', header=true))
      TO ? (FORMAT PARQUET, COMPRESSION ZSTD)
  """, [str(tsv_path), str(parquet_path)])
  ```
- [ ] Run one-off TSV→Parquet conversion for all existing PGS (see README)

### Per-chromosome dashboard
- [ ] `/dashboard/<pgs_id>/chrom` page: chromosome-level stats table
- [ ] Completion badge (N/22 chromosomes computed)

### Pipeline manifest generation
- [ ] `POST /api/pipeline/<pgs_id>/manifest`:
  - Write `/work/manifests/{pgs_id}_manifest.json`
  - Save to `pipeline_plans` table with `status='manifest_written'`
  - Return manifest content + path
- [ ] UI button: "Prepare for processing" → POST manifest
- [ ] Manifest format:
  ```json
  {
    "pgs_id": "PGS000004",
    "betamap_path": "...",
    "zarr_base": "...",
    "output_dir": "...",
    "script": "run_prs_spark_gpu.sbatch",
    "slurm_account": "researchers",
    "slurm_partition": "gpu"
  }
  ```

### Sample lookup
- [ ] `/api/data/<pgs_id>/sample/<sample_id>` — DuckDB point query for one sample
- [ ] Requires parquet (TSV point query is too slow at 140k rows without index)

### Stats cache management
- [ ] `/api/admin/cache/invalidate/<pgs_id>` — force re-compute (POST)
- [ ] Display cache status in dashboard footer ("Stats computed X min ago")

---

## Phase 3 – Slurm integration

### Architecture (preferred)
```
Web → write manifest → /work/manifests/ → operator reviews → sbatch
```
The web app NEVER calls `sbatch`. An external agent or cron job reads manifests.

### Endpoints to add
- [ ] `GET /api/pipeline/<pgs_id>/status` — check Slurm job status via `sacct`
  (read-only, queries Slurm accounting DB or output files)
- [ ] `GET /api/pipeline/<pgs_id>/logs` — tail sbatch `--output` file (read-only)

### Security gate
- [ ] Add auth to `/api/admin/*` (API key in header, validated from env var)
- [ ] Rate limit remote search endpoint (pgscat API is shared/public)

---

## Phase 4 – Production hardening

- [ ] Structured JSON logging for log aggregation (Loki/ELK)
- [ ] `/metrics` endpoint (Prometheus format) with:
  - `pgs_count_total`, `pgs_parquet_count`, `pgs_tsv_only_count`
  - `search_duration_seconds` (histogram)
  - `db_pool_available` (gauge)
  - `cache_hit_total`, `cache_miss_total`
- [ ] Nginx auth_basic or oauth2-proxy in front for access control
- [ ] Gunicorn `--max-requests 1000` to prevent memory leaks from DuckDB connections
- [ ] systemd unit file for non-container deployment on host

---

## Known limitations

| # | Issue | Impact | Workaround |
|---|-------|--------|------------|
| 1 | pgscat no free-text search | "text" mode only works with EFO IDs, PMIDs, PGS IDs | Document in UI; point users to pgscatalog.org |
| 2 | No Parquet in pipeline yet | Dashboard falls back to TSV (slower) | Run one-off conversion; add to pipeline |
| 3 | DuckDB cold cache on CephFS | First stats request ~20s after restart | File cache warms on first hit; subsequent requests are fast |
| 4 | gunicorn workers all sync on startup | 4 workers × sync = 4 parallel upserts | Idempotent upserts are safe; slight DB overhead |
| 5 | `/api/admin/sync` unauthenticated | Anyone on the network can trigger sync | Restrict via Nginx `allow/deny` or network policy |
| 6 | pgscat heavy deps (zarr/pandas/numpy) | ~300 MB extra in container image | Use `APP_PGSCAT_MODE=cli` + host-installed pgscat |
| 7 | psycopg2-binary vs libpq | `psycopg2-binary` works on most Linux; wheels may differ on arm64 | Use `psycopg2` + `libpq-dev` on arm64 |
