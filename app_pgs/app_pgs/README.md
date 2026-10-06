# PGS/PRS Dashboard

**Production URL**: https://pgscat.mcps-epcm.org/

Internal web platform for exploring Polygenic Score (PGS/PRS) results from the MCPS cohort.

Stack: **Bottle · DuckDB · PostgreSQL 16 · pgscat · Gunicorn · Podman**

---

## Architecture summary

| Layer | Technology | Role |
|-------|-----------|------|
| Web framework | Bottle | Routes, templates, WSGI |
| Analytical queries | DuckDB (in-memory) | Parquet/TSV stats, histograms |
| Metadata persistence | PostgreSQL 16 | Catalog index, remote cache, audit |
| File cache | JSON in `/work/cache/` | Stats TTL cache (no DB hit on hot data) |
| Remote catalog | pgscat Python library | PGS Catalog API queries |
| WSGI server | Gunicorn (4 sync workers) | Production serving |
| Reverse proxy | Nginx/Caddy (external) | TLS, auth, rate limiting |

**Key rule**: sample-level scores (`PRS_total` values) never enter PostgreSQL. They stay in Parquet files on CephFS and are queried only by DuckDB.

---

## Quick start

### Direct (no container)

```bash
# 1. Install Python deps
pip install bottle duckdb gunicorn requests psycopg2-binary
pip install "pgscat @ git+https://github.com/fenandosr/pgscat.git"

# 2. Configure (minimum for filesystem-only mode)
export APP_DATA_DIR=/mnt/cephfs/hot_nvme/pgscatalog/scores
export APP_SCRIPTS_DIR=/mnt/cephfs/orgs/home/angel.pacheco/Ordanamiento_Bases/scripts

# 3. Configure PostgreSQL (optional — app degrades gracefully if absent)
export APP_DB_HOST=127.0.0.1
export APP_DB_NAME=pgs_dashboard
export APP_DB_USER=pgs_user
export APP_DB_PASSWORD=yourpassword   # or APP_DB_PASSWORD_FILE=/path/to/file

# 4. Create DB schema (one-time)
psql -h 127.0.0.1 -U pgs_user -d pgs_dashboard -f sql/schema.sql

# 5. Run (dev)
cd src && python3 app.py

# 5. Run (production)
cd src && gunicorn --bind 0.0.0.0:8080 --workers 4 --timeout 120 app:app
```

### Container (attached to mcps-net)

```bash
# Build
podman build -t pgs-dashboard:1.1 -f Containerfile .

# Run with DB on mcps-net
podman run --rm -p 8080:8080 \
  --network mcps-net \
  -v /mnt/cephfs/hot_nvme/pgscatalog/scores:/data:ro \
  -v /mnt/cephfs/orgs/home/angel.pacheco/Ordanamiento_Bases/scripts:/scripts:ro \
  -v /tmp/pgs-work:/work \
  -e APP_DB_HOST=postgres \
  -e APP_DB_NAME=pgs_dashboard \
  -e APP_DB_USER=pgs_user \
  -e APP_DB_PASSWORD_FILE=/run/secrets/db_password \
  -v /path/to/db_password_file:/run/secrets/db_password:ro \
  pgs-dashboard:1.1
```

---

## Environment variables

| Variable | Default | Description |
|----------|---------|-------------|
| `APP_DATA_DIR` | `/data` | Scores directory (read-only mount) |
| `APP_WORK_DIR` | `/work` | Writable temp/cache directory |
| `APP_SCRIPTS_DIR` | `/scripts` | Pipeline scripts (read-only, optional) |
| `APP_PGSCAT_MODE` | `python` | `python` = library; `cli` = subprocess |
| `APP_PGSCAT_BIN` | `pgscat` | pgscat binary path (CLI mode) |
| `APP_PGSCAT_TIMEOUT` | `30` | Remote query timeout (seconds) |
| `APP_REMOTE_ENABLED` | `true` | Enable/disable remote PGS Catalog |
| `APP_REMOTE_CACHE_TTL` | `86400` | Remote metadata cache TTL (seconds) |
| `APP_DB_HOST` | `postgres` | PostgreSQL host (`postgres` in container, `127.0.0.1` on host) |
| `APP_DB_PORT` | `5432` | PostgreSQL port |
| `APP_DB_NAME` | `pgs_dashboard` | Database name |
| `APP_DB_USER` | `pgs_user` | Database user |
| `APP_DB_PASSWORD_FILE` | _(none)_ | Path to file with DB password (preferred) |
| `APP_DB_PASSWORD` | _(none)_ | Fallback password env var |
| `APP_STATS_CACHE_TTL` | `3600` | DuckDB stats file-cache TTL (seconds) |
| `APP_TRUSTED_PROXY` | `true` | Trust `X-Forwarded-*` headers from proxy |
| `APP_HOST` | `0.0.0.0` | Bind host (dev only) |
| `APP_PORT` | `8080` | Bind port (dev only) |

---

## PostgreSQL setup

```sql
-- Run as postgres superuser once:
CREATE USER pgs_user WITH PASSWORD 'yourpassword';
CREATE DATABASE pgs_dashboard OWNER pgs_user;
\c pgs_dashboard
-- Apply schema:
\i /path/to/sql/schema.sql
```

Or from shell:
```bash
psql -h 127.0.0.1 -U postgres -c "CREATE USER pgs_user WITH PASSWORD 'yourpassword';"
psql -h 127.0.0.1 -U postgres -c "CREATE DATABASE pgs_dashboard OWNER pgs_user;"
psql -h 127.0.0.1 -U pgs_user -d pgs_dashboard -f sql/schema.sql
```

### PostgreSQL connection: host vs container

| Scenario | `APP_DB_HOST` value |
|----------|---------------------|
| App runs directly on host | `127.0.0.1` |
| App in container, DB in separate container on `mcps-net` | `postgres` (container alias) |
| App in container, DB on host | host IP on Podman bridge (not `localhost`) |

---

## Reverse proxy (Nginx example)

```nginx
server {
    listen 443 ssl;
    server_name pgscat.mcps-epcm.org;

    # TLS config (certificates managed externally)
    ssl_certificate     /etc/ssl/certs/pgscat.crt;
    ssl_certificate_key /etc/ssl/private/pgscat.key;

    location / {
        proxy_pass         http://127.0.0.1:8080;
        proxy_set_header   Host              $host;
        proxy_set_header   X-Real-IP         $remote_addr;
        proxy_set_header   X-Forwarded-For   $proxy_add_x_forwarded_for;
        proxy_set_header   X-Forwarded-Proto $scheme;
        proxy_set_header   X-Forwarded-Host  $host;
        proxy_read_timeout 120s;
    }

    # Optional: restrict admin endpoints at proxy level
    location /api/admin/ {
        allow 10.0.0.0/8;
        deny all;
        proxy_pass http://127.0.0.1:8080;
    }
}
```

`APP_TRUSTED_PROXY=true` makes Bottle honour `X-Forwarded-*` headers.

---

## Routes

| Method | Path | Description |
|--------|------|-------------|
| GET | `/` | Local catalog (DB index if available, filesystem fallback) |
| GET | `/dashboard/<pgs_id>` | Dashboard with stats + Plotly histogram |
| GET | `/api/data/<pgs_id>` | JSON: stats, histogram, column info |
| GET | `/search` | Remote PGS Catalog search UI |
| POST | `/search` | Execute remote search (cached + audit-logged) |
| GET | `/pgs/<pgs_id>/source` | Local vs remote status + pipeline plan |
| GET | `/api/pipeline/<pgs_id>/plan` | Pipeline plan JSON (no execution) |
| POST | `/api/admin/sync` | Trigger catalog→DB sync |
| GET | `/healthz` | Liveness probe (always 200) |
| GET | `/readyz` | Readiness probe (200/503 + JSON component status) |

---

## Parquet conversion (required for dashboard)

The pipeline currently produces TSV files. Convert to Parquet for faster web queries:

```bash
for pgs in PGS000004 PGS000363; do
  base="/mnt/cephfs/hot_nvme/pgscatalog/scores/${pgs}/${pgs}"
  duckdb -c "
    COPY (SELECT * FROM read_csv_auto('${base}_PRS_total.tsv', delim=chr(9), header=true))
    TO '${base}_PRS_total.parquet' (FORMAT PARQUET, COMPRESSION ZSTD)"
done
```

---

## Project structure

```
.
├── Containerfile
├── requirements.txt
├── README.md
├── sql/
│   └── schema.sql              ← PostgreSQL schema (apply once)
├── docs/
│   ├── ARCHITECTURE.md
│   ├── PIPELINE_INTEGRATION.md
│   └── TODO.md
└── src/
    ├── app.py                  ← Bottle routes + WSGI app (with proxy middleware)
    ├── config.py               ← All env-var config
    ├── services/
    │   ├── db.py               ← PostgreSQL pool + graceful degradation
    │   ├── cache.py            ← File-based JSON cache with TTL
    │   ├── sync_service.py     ← Catalog → PostgreSQL sync
    │   ├── health.py           ← /healthz + /readyz logic
    │   ├── local_catalog.py    ← Filesystem scanner
    │   ├── parquet_stats.py    ← DuckDB stats + histogram (with cache)
    │   ├── pgscat_client.py    ← Remote PGS Catalog (lib + CLI + DB cache)
    │   └── pipeline_inspector.py ← Pipeline plan generator (read-only)
    ├── views/                  ← Bottle SimpleTemplate (.tpl)
    └── static/                 ← CSS + JS
```
