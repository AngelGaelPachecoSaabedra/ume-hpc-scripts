# PGS/PRS Platform — Architecture

## Platform modes

| Mode | Routes | Description |
|------|--------|-------------|
| A — Local catalog | `/`, `/dashboard/<id>`, `/api/data/<id>` | Visualise locally computed PGS/PRS results |
| B — Remote catalog | `/search`, `/pgs/<id>/source`, `/api/pipeline/<id>/plan` | Search PGS Catalog; inspect remote scores |
| C — Custom PRS | `/custom-prs/new`, `/api/custom-prs/plan` | User-defined score — generate pipeline plan (no execution) |
| D — Variant Annotator | `/variants/<id>`, `/api/variants/<id>/*` | View HPC-computed variant annotation results (read-only) |
| E — Gene Browser | `/genes`, `/gene/<symbol>`, `/api/gene/<symbol>/*` | Gene-centric interpretation browser (Plotly tracks) |

---

## Tech stack

| Component | Technology |
|-----------|-----------|
| Web framework | Bottle (sync, WSGI) |
| WSGI server | Gunicorn — 4 sync workers, 120 s timeout |
| Analytics queries | DuckDB (in-memory, per-request, read-only) |
| Operational DB | PostgreSQL 16 via psycopg2 ThreadedConnectionPool |
| Remote score API | pgscat Python library + CLI fallback |
| Templates | Bottle SimpleTemplate (`.tpl`) |
| CSS | TailwindCSS CDN |
| Charts | Plotly.js CDN |
| Container | Docker/Podman, non-root user (UID 1000) |

---

## Service layer

```
src/services/
  cache.py              FileCache — L1 atomic JSON cache (os.replace, TTL)
  db.py                 DBPool — PostgreSQL pool with graceful degradation
  health.py             HealthService — /healthz + /readyz
  local_catalog.py      LocalCatalog — read-only scan of scores directory
  normalizer.py         Type normalizer — polymorphic pgscat → clean Python types
  parquet_stats.py      ParquetStats — DuckDB stats + histogram (with file cache)
  pgscat_client.py      PGSCatClient — remote PGS Catalog queries (lib + CLI + DB cache)
  pipeline_inspector.py PipelineInspector — 3-flow plan generator (no execution)
  score_preparer.py     ScorePreparer — local layout manifest + preparation status
  sync_service.py       SyncService — catalog → DB background sync
  variant_annotation.py VariantAnnotationService — read-only annotation results
  gene_browser.py       GeneBrowserService — gene-centric track data (Mode E)
  clinical_genes.py     ClinicalGenesService — clinical gene list queries
```

---

## Gene Browser (Mode E)

### Track architecture

```
Track 1 — Variants (50% height)
  Scatter: x=POS, y=consequence_priority
  Color by consequence: pLoF=red, missense=orange, synonymous=green,
                        splice_region=dark_orange, non_coding=slate
  Point size: scaled by |BETA| (clamped 6–18 px)

Track 2 — Gene model (20% height)
  Exon blocks: rectangles from clustered coding variant positions
  Intron line:  horizontal line across genomic range
  Label:        gene name + strand + genomic coordinates
  Note:         approximate model derived from PRS variant positions

Track 3 — Functional summary (30% height)
  Stacked bar: consequence category counts per genomic bin
  Secondary axis: mean allele frequency per bin (when available)
```

**Visualization library:** Plotly `make_subplots` with shapes and scatter traces.
LocusZoom.js was evaluated and rejected — it requires GWAS-specific data adapters with no
benefit for this gene-centric, PRS-oriented use case.

---

## Data flow

```
External:  PGS Catalog API  ──► pgscat_client ──► DB cache (remote_pgs_cache)
                                                 └► search_audit_log

Local:     /data/{PGS_ID}/
           *.parquet / *.tsv  ──► DuckDB (in-memory) ──► parquet_stats ──► file cache
                                                       └► JSON response

HPC:       pipeline_scripts/  ──► writes *.parquet / annotation outputs
           annotator/         ──► writes annotation TSV/Parquet to /annotations/
                                  (read-only by web app via variant_annotation service)
```

---

## Data layout contract

```
{APP_DATA_DIR}/{PGS_ID}/
  {PGS_ID}.json                                  PGS Catalog metadata
  {PGS_ID}_hmPOS_GRCh38.txt.gz                  Raw harmonized scoring file
  {PGS_ID}_hmPOS_GRCh38.betamap.tsv.gz          Pipeline input
  {PGS_ID}_hmPOS_GRCh38.betamap.tsv.gz.meta.json  Betamap generation report
  {PGS_ID}_chr{N}_scores.tsv                    Per-chromosome PRS
  {PGS_ID}_chr{N}_metadata.json                 Per-chromosome run metadata
  {PGS_ID}_PRS_total.tsv                        Aggregated (sample_id, PRS_total, PRS_chr1..22)
  {PGS_ID}_PRS_total_metadata.json              Aggregation report
  {PGS_ID}_PRS_total.parquet                    Parquet export (ZSTD) — required for web dashboard
```

**Betamap columns:** `PRS_ID, CHROM, POS, ID, EFFECT_ALLELE, OTHER_ALLELE, BETA, IS_FLIP`

- `IS_FLIP=1` → effect allele on REF strand; dosage = `2 − call_DS`
- `IS_FLIP=0` → standard; dosage = `call_DS`

---

## Data separation rule

Per-sample PRS scores (~140k rows × 24 columns) are **never** stored in PostgreSQL.
They live only in Parquet/TSV files and are queried directly via DuckDB.

PostgreSQL stores only:

| Table | Contents |
|-------|----------|
| `local_pgs_index` | Catalog index (pgs_id, trait, counts, file flags) |
| `remote_pgs_cache` | PGS Catalog API metadata cache |
| `stats_cache` | Aggregated statistics (NOT raw sample scores) |
| `search_audit_log` | Search query log |

---

## Volume mounts (container)

| Mount | Mode | Contents |
|-------|------|----------|
| `/data` | read-only | PGS score files (Parquet/TSV) |
| `/work` | read-write | Cache, temp files |
| `/pipeline_scripts` | read-only | HPC pipeline scripts (optional) |
| `/annotations` | read-only | Variant annotation outputs |
| `/ref` | read-only | GFF3 reference (Apptainer only) |
| `/fasta` | read-only | Reference FASTA (Apptainer only) |
| `/dbnsfp` | read-only | dbNSFP5 scores (Apptainer only) |
| `/dbsnp_freq` | read-only | dbSNP population freq VCF (Apptainer only) |

---

## Gunicorn

```
gunicorn --workers 4 --timeout 120 --bind 0.0.0.0:8080 --chdir /app/src app:app
```

120 s timeout accommodates cold-read latency on large Parquet files over network storage.

---

## Variant Annotator (Apptainer pipeline, HPC)

The `apptainer/variant_annotator.def` defines a self-contained HPC container that:

1. Parses betamap TSV variants
2. Annotates against GENCODE GFF3 (gene/transcript overlap, coding consequence)
3. Adds dbNSFP5 scores (CADD, REVEL, SIFT, PolyPhen2, ClinVar)
4. Adds dbSNP population allele frequencies (rsid, af_global, rarity_class)
5. Writes `{PGS_ID}_variants_annotated.parquet` and `{PGS_ID}_annotation_summary.json`

The web app (Mode D) reads these outputs as read-only via `APP_ANNOTATIONS_DIR`.
The `annotator/` module (Python scripts) is **not included** in this repository.

See `docs/VARIANT_ANNOTATION.md` for full annotation column reference.
