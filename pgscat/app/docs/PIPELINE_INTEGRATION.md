# Pipeline Integration

This document describes how the platform connects to the HPC pipeline, what it can automate, and what requires manual Slurm submission.

## Pipeline scripts

Pipeline scripts are mounted read-only at `APP_SCRIPTS_DIR` (default: `/pipeline_scripts`).
The app reads script names and checks their existence but never executes them directly.

Key scripts:

| Script | Purpose |
|--------|---------|
| `sbrc_01_download_format.sbatch` | Download raw scoring file + build betamap |
| `compute_prs_spark_gpu.py` | Per-chrom scoring (Spark + CuPy GPU) |
| `prs_gpu_compute.sbatch` | Slurm array job wrapper for GPU scoring |
| `run_prs_spark_gpu.sbatch` | Single-chrom Slurm job wrapper |
| `compute_prs.py` | CPU/DuckDB fallback scoring |
| `prs_end2end.sbatch` | End-to-end pipeline (clumping + scoring) |

## Local file layout contract

Expected layout under `{APP_DATA_DIR}/{PGS_ID}/`.
This is the **definitive contract** — all services use this layout.

```
{APP_DATA_DIR}/{PGS_ID}/
  {PGS_ID}.json                                  PGS Catalog metadata JSON
  {PGS_ID}_hmPOS_GRCh38.txt.gz                  Raw GRCh38-harmonised scoring file
  {PGS_ID}_hmPOS_GRCh38.columns.json            Column mapping (auto-detected)
    Fields: file, delimiter, meta{...}, n_columns, header[], columns{CHROM, POS, ID,
            OTHER_ALLELE, EFFECT_ALLELE, EFFECT_WEIGHT}, missing[]
  {PGS_ID}_hmPOS_GRCh38.betamap.tsv.gz          Pipeline input betamap
    Columns: PRS_ID, CHROM, POS, ID, EFFECT_ALLELE, OTHER_ALLELE, BETA, IS_FLIP
    IS_FLIP=1 → effect allele on REF strand, dosage = 2 − call_DS
  {PGS_ID}_hmPOS_GRCh38.betamap.tsv.gz.meta.json  Betamap generation report
    Fields: prs_id, prs_path, colsinfo_path, betamap_path, weight_type,
            selected_columns[], n_total_rows, n_written_rows, n_parse_error,
            n_bad_allele, n_not_translated, n_skipped_ref_mismatch
  {PGS_ID}_chr{N}_scores.tsv                    Per-chrom PRS output
    Columns: sample_id, PRS
  {PGS_ID}_chr{N}_metadata.json                 Per-chrom run metadata
    Fields: pgs_id, chrom, run_timestamp, elapsed_seconds, missing_strategy,
            window_size, dosage_array_key, dosage_layout, index_path,
            n_samples, n_zarr_variants, n_variants_in_score_file,
            n_variants_matched, n_variants_excluded_allele, n_variants_excluded_missing,
            prs_mean, prs_std, prs_min, prs_max,
            excluded_allele_mismatch[], excluded_missing[]
  {PGS_ID}_PRS_total.tsv                        Aggregated PRS
    Columns: sample_id, PRS_total, PRS_chr1, PRS_chr2, ..., PRS_chr22
  {PGS_ID}_PRS_total_metadata.json              Aggregation report
    Fields: pgs_id, aggregation_timestamp, elapsed_seconds,
            chroms_requested[], chroms_found[], chroms_missing[],
            n_samples, n_chr_scores, prs_total_mean/std/min/max,
            per_chromosome{chrom → metadata}
  {PGS_ID}_PRS_total.parquet                    Parquet export (ZSTD)
    Same schema as PRS_total.tsv. Required for web dashboard.
```

## Genotype data

| Resource | Description |
|----------|-------------|
| Zarr genotype store | Per-chromosome zarr arrays: `{GENOTYPES_DIR}/chr{N}.zarr` |
| DuckDB zarr index | Sample + variant index: `{GENOTYPES_DIR}/zarr_index.duckdb` |
| Chromosomes | 1–22 |

Configure `GENOTYPES_DIR` in pipeline scripts; the web app does not access genotype data directly.

## Three pipeline flows

### Flow A: Local existing (Mode A — dashboard)
Score already computed. All files present. Dashboard reads Parquet via DuckDB.

### Flow B: Catalog remote prepare (Mode B — pipeline plan)
Score in PGS Catalog, not locally present.

1. **Download** — `pgscat download {PGS_ID} --build GRCh38` → `{PGS_ID}_hmPOS_GRCh38.txt.gz`
2. **Betamap** — `sbrc_01_download_format.sbatch` → `{PGS_ID}_hmPOS_GRCh38.betamap.tsv.gz`
3. **Score** — `sbatch --array=1-22 prs_gpu_compute.sbatch` → `{PGS_ID}_chr{N}_scores.tsv`
4. **Aggregate** → `{PGS_ID}_PRS_total.tsv`
5. **Parquet** — `duckdb -c "COPY ... TO ... (FORMAT PARQUET, COMPRESSION ZSTD)"` → `{PGS_ID}_PRS_total.parquet`

Steps 1–2 and 5 can be automated by the platform (no Slurm needed).
Steps 3–4 require Slurm submission.

### Flow C: Custom PRS (Mode C)
User provides their own scoring file. Same pipeline steps as Flow B after the user drops their file in place.

## Automation matrix

| Step | Status | Notes |
|------|--------|-------|
| pgscat metadata download | Automatable | Requires `REMOTE_ENABLED=true` |
| pgscat scoring file download | Automatable | Requires pgscat CLI |
| Betamap preparation | Automatable | sbrc_01_download_format.sbatch |
| Per-chrom GPU scoring | **Informational only** | Requires Slurm GPU allocation |
| PRS aggregation | **Informational only** | Depends on all chroms completing |
| Parquet export | Automatable | DuckDB, no Slurm needed |

## Phase 2 roadmap (not yet implemented)

- Submit Slurm jobs directly from the web UI (requires auth + job DB table)
- Poll Slurm job status and update preparation stage
- Email/webhook notification when a score becomes dashboard-ready
- Batch preparation of multiple PGS IDs
