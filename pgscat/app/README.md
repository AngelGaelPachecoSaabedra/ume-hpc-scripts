# PGS/PRS Dashboard

Web platform for exploring Polygenic Score (PGS/PRS) results computed from cohort genotype data.

**Stack:** Bottle · DuckDB · PostgreSQL 16 · pgscat · Gunicorn · Podman

---

## Features

- **Local catalog** — browse and visualize PGS/PRS scores stored as Parquet files
- **Remote catalog** — search [PGS Catalog](https://www.pgscatalog.org/) API and inspect score metadata
- **Variant Annotator viewer** — read-only view of HPC-computed variant annotations (GFF3-based)
- **Gene Browser** — gene-centric genomic track visualization (Plotly stacked subplots)
- **Custom PRS planner** — generate a pipeline plan for user-defined scoring files (no execution)
- **Health probes** — `/healthz` and `/readyz` for container orchestration

---

## Architecture summary

| Layer | Technology | Role |
|-------|-----------|------|
| Web framework | Bottle | Routes, templates, WSGI |
| Analytical queries | DuckDB (in-memory) | Parquet/TSV stats, histograms |
| Metadata persistence | PostgreSQL 16 | Catalog index, remote cache, audit |
| File cache | JSON (work dir) | Stats TTL cache |
| Remote catalog | pgscat Python library | PGS Catalog API queries |
| WSGI server | Gunicorn (4 sync workers) | Production serving |
| Reverse proxy | Nginx/Caddy (external) | TLS, auth, rate limiting |

> **Key rule:** per-sample score values never enter PostgreSQL. They stay in Parquet files and are queried only by DuckDB.

---

## Quick start

### Direct (no container)

```bash
# 1. Install Python deps
pip install -r requirements.txt
pip install "pgscat @ git+https://github.com/fenandosr/pgscat.git"

# 2. Set required paths
export APP_DATA_DIR=/path/to/scores       # directory containing PGS_ID subdirectories
export APP_WORK_DIR=/tmp/pgs-work         # writable cache directory

# 3. Configure PostgreSQL (optional — app degrades gracefully if absent)
export APP_DB_HOST=127.0.0.1
export APP_DB_NAME=pgs_dashboard
export APP_DB_USER=pgs_user
export APP_DB_PASSWORD=yourpassword

# 4. Create DB schema (one-time)
psql -h 127.0.0.1 -U pgs_user -d pgs_dashboard -f sql/schema.sql

# 5. Run (development)
cd src && python3 app.py

# 5. Run (production)
cd src && gunicorn --bind 0.0.0.0:8080 --workers 4 --timeout 120 app:app
```

### Container (Docker / Podman)

```bash
# Build
podman build -t pgs-dashboard:latest -f Containerfile .

# Run (minimal — filesystem only, no DB)
podman run --rm -p 8080:8080 \
  -v /path/to/scores:/data:ro \
  -v /tmp/pgs-work:/work \
  pgs-dashboard:latest

# Run (with PostgreSQL)
podman run --rm -p 8080:8080 \
  -v /path/to/scores:/data:ro \
  -v /tmp/pgs-work:/work \
  -e APP_DB_HOST=<db-host> \
  -e APP_DB_NAME=pgs_dashboard \
  -e APP_DB_USER=pgs_user \
  -e APP_DB_PASSWORD_FILE=/run/secrets/db_password \
  -v /path/to/db_password_file:/run/secrets/db_password:ro \
  pgs-dashboard:latest
```

See `deployment/docker-compose.yml` for a full stack example.

---

## Environment variables

| Variable | Default | Description |
|----------|---------|-------------|
| `APP_DATA_DIR` | `/data` | Scores directory (Parquet/TSV, read-only) |
| `APP_WORK_DIR` | `/work` | Writable cache directory |
| `APP_SCRIPTS_DIR` | `/pipeline_scripts` | Pipeline scripts mount (read-only, optional) |
| `APP_ANNOTATIONS_DIR` | `/annotations` | Variant annotation outputs (read-only) |
| `APP_PGSCAT_MODE` | `python` | `python` = library; `cli` = subprocess |
| `APP_PGSCAT_BIN` | `pgscat` | pgscat binary path (CLI mode) |
| `APP_PGSCAT_TIMEOUT` | `30` | Remote query timeout (seconds) |
| `APP_REMOTE_ENABLED` | `true` | Enable/disable remote PGS Catalog queries |
| `APP_REMOTE_CACHE_TTL` | `86400` | Remote metadata cache TTL (seconds) |
| `APP_DB_HOST` | `postgres` | PostgreSQL host |
| `APP_DB_PORT` | `5432` | PostgreSQL port |
| `APP_DB_NAME` | `pgs_dashboard` | Database name |
| `APP_DB_USER` | `pgs_user` | Database user |
| `APP_DB_PASSWORD_FILE` | _(none)_ | Path to file with DB password (preferred) |
| `APP_DB_PASSWORD` | _(none)_ | Fallback: DB password env var |
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

---

## Routes

| Method | Path | Description |
|--------|------|-------------|
| GET | `/` | Local catalog index |
| GET | `/dashboard/<pgs_id>` | Dashboard: stats + Plotly histogram |
| GET | `/api/data/<pgs_id>` | JSON: stats, histogram, column info |
| GET | `/search` | Remote PGS Catalog search UI |
| POST | `/search` | Execute remote search |
| GET | `/pgs/<pgs_id>/source` | Local vs remote status + pipeline plan |
| GET | `/api/pipeline/<pgs_id>/plan` | Pipeline plan JSON (read-only) |
| POST | `/api/admin/sync` | Trigger catalog→DB sync |
| GET | `/healthz` | Liveness probe |
| GET | `/readyz` | Readiness probe (200/503 + JSON) |

---

## Expected data layout

```
{APP_DATA_DIR}/{PGS_ID}/
  {PGS_ID}.json                          PGS Catalog metadata
  {PGS_ID}_hmPOS_GRCh38.txt.gz          Raw harmonized scoring file
  {PGS_ID}_hmPOS_GRCh38.betamap.tsv.gz  Pipeline input
  {PGS_ID}_PRS_total.tsv                Aggregated scores (sample_id, PRS_total, ...)
  {PGS_ID}_PRS_total.parquet            Parquet export (recommended for web)
```

See `examples/data/` for small dummy files in this format.

### Convert TSV → Parquet

```bash
duckdb -c "
  COPY (SELECT * FROM read_csv_auto('PGS000001_PRS_total.tsv', delim=chr(9), header=true))
  TO 'PGS000001_PRS_total.parquet' (FORMAT PARQUET, COMPRESSION ZSTD)"
```

---

## Project structure

```
.
├── Containerfile              Container image (Docker/Podman)
├── requirements.txt
├── README.md
├── LICENSE
├── sql/
│   └── schema.sql             PostgreSQL schema
├── docs/
│   ├── ARCHITECTURE.md        System architecture
│   ├── PIPELINE_INTEGRATION.md  HPC pipeline integration contract
│   └── VARIANT_ANNOTATION.md  Variant annotator specification
├── src/
│   ├── app.py                 Bottle routes + WSGI app
│   ├── config.py              All env-var configuration
│   ├── services/              Service layer (DB, cache, DuckDB, pgscat, ...)
│   ├── views/                 Bottle SimpleTemplate (.tpl)
│   └── static/                CSS + JS
├── resources/
│   ├── carrier_screen/        Carrier screening gene schema + loader
│   └── clinical_genes/        Clinical gene lists (public sources)
├── examples/
│   └── data/                  Dummy PGS datasets for testing
├── deployment/
│   └── docker-compose.yml     Full stack compose example
└── apptainer/
    └── variant_annotator.def  Apptainer definition for HPC annotation pipeline
```

---

## Reverse proxy (Nginx example)

```nginx
server {
    listen 443 ssl;
    server_name your-domain.example.com;

    ssl_certificate     /etc/ssl/certs/pgs-dashboard.crt;
    ssl_certificate_key /etc/ssl/private/pgs-dashboard.key;

    location / {
        proxy_pass         http://127.0.0.1:8080;
        proxy_set_header   Host              $host;
        proxy_set_header   X-Real-IP         $remote_addr;
        proxy_set_header   X-Forwarded-For   $proxy_add_x_forwarded_for;
        proxy_set_header   X-Forwarded-Proto $scheme;
        proxy_read_timeout 120s;
    }
}
```

Set `APP_TRUSTED_PROXY=true` to honour `X-Forwarded-*` headers.
