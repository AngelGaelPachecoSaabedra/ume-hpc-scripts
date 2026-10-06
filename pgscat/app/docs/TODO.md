# TODO / Roadmap

## Phase 1 – Completed ✓

- [x] Local catalog scanner (filesystem)
- [x] Per-PGS stats dashboard (DuckDB, Plotly histogram)
- [x] Remote search via pgscat (library + CLI fallback)
- [x] PostgreSQL integration (graceful degradation)
- [x] 2-layer cache (FileCache L1 + PostgreSQL L2)
- [x] Catalog→DB background sync
- [x] Health endpoints (/healthz, /readyz)
- [x] Reverse proxy support (X-Forwarded-* headers)
- [x] Search audit logging
- [x] Data normalizer (polymorphic pgscat types)
- [x] Pipeline plan generator (3 flow modes: local/catalog/custom)
- [x] Score preparer + local layout manifest
- [x] HTML pipeline plan view (/pipeline/<id>/plan)
- [x] Custom PRS form + plan API (/custom-prs/new, /api/custom-prs/plan)
- [x] Preparation manifest JSON API (/api/pipeline/<id>/manifest)
- [x] Updated nav (Mode A / B / C clearly surfaced)
- [x] Pipeline scripts source configurable via APP_SCRIPTS_DIR

## Phase 2 – Next priorities

### Preparation automation
- [ ] Button in UI to trigger pgscat download (needs job queue or direct subprocess)
- [ ] Button to trigger betamap preparation (sbrc_01_download_format.sbatch)
- [ ] Slurm job submission from web UI (auth required)
- [ ] Job status polling (Slurm sacct or squeue)
- [ ] Stage progress indicator (auto-refresh preparation status)

### Dashboard enhancements
- [ ] Per-chromosome score breakdown (individual chrom PRS distributions)
- [ ] Sample lookup endpoint (query by sample_id across scores)
- [ ] Comparison view (two PGS IDs side-by-side)
- [ ] Percentile calculator (input a PRS score → population percentile)

### Operations
- [ ] Parquet export automation (run DuckDB COPY after TSV appears)
- [ ] Add authentication to /api/admin/sync
- [ ] Prometheus metrics endpoint /metrics
- [ ] Alerting when new PGS directories appear without parquet

### Phase 3 – Future
- [ ] LDpred2 / SBayesRC integration plan
- [ ] Multi-ancestry PRS support
- [ ] Clinical validation report template
- [ ] Batch PGS preparation (queue multiple IDs)
