# Pipeline Integration

This document describes how the web dashboard connects to the existing HPC pipeline,
and where to extend it for Phase 2 (Slurm job submission).

---

## Verified pipeline facts (from code inspection)

### Input: betamap format

File: `{PGS_ID}_hmPOS_GRCh38.betamap.tsv.gz`

| Column        | Type    | Description                                              |
|---------------|---------|----------------------------------------------------------|
| `PRS_ID`      | str     | PGS Catalog ID                                           |
| `CHROM`       | str     | Chromosome ("1"–"22", "X")                               |
| `POS`         | int     | GRCh38 1-based position                                  |
| `ID`          | str/NaN | rsID (may be NaN/null for many PGS)                      |
| `EFFECT_ALLELE`| str   | Allele with directional effect                           |
| `OTHER_ALLELE`| str     | Reference allele                                         |
| `BETA`        | float   | Effect size coefficient                                  |
| `IS_FLIP`     | int     | 0 = effect on ALT; 1 = effect on REF (dosage = 2−DS)    |

### Compute model

```
PRS[sample] = Σ_i  BETA_i × dosage_i

dosage_i =
  call_DS[i]       if IS_FLIP = 0
  2 − call_DS[i]   if IS_FLIP = 1
```

Missing dosages are imputed by the pipeline (`mean`, `zero`, or `skip` strategy).

### Genotype data layout

Path: `/mnt/cephfs/hot_nvme/mcps/imputed-topmed/zar_files/chr{N}.zarr`

Key arrays:
- `sample_id`: (140,831,) string
- `variant_position`: (n_variants,) int64 — GRCh38 1-based
- `variant_allele`: (n_variants, 2) bytes — REF/ALT
- `call_DS` or `call_genotype`: dosage or genotype arrays

### Per-chromosome output

File: `{PGS_ID}_chr{N}_scores.tsv` — tab-separated, 2 columns:

```
sample_id       PRS
MCPS_...        0.04512345
```

Metadata: `{PGS_ID}_chr{N}_metadata.json` — includes run timestamp, n_samples,
n_variants_matched, n_variants_excluded, prs_mean, prs_std, etc.

### Aggregated output

File: `{PGS_ID}_PRS_total.tsv` — tab-separated, 24 columns:

```
sample_id    PRS_total    PRS_chr1    PRS_chr2  …  PRS_chr22
MCPS_...     0.04491      0.00821     0.00195   …  0.00043
```

`PRS_total = Σ_chr PRS_chr` (sum across all 22 autosomes per sample).

---

## Recommended new step: TSV → Parquet export

The web dashboard prefers parquet for fast column-pruned DuckDB queries.
This step should be added **after aggregation** in the pipeline.

### Option A: standalone DuckDB command (add to sbatch)

```bash
# After PRS_total.tsv is produced:
duckdb -c "
COPY (
  SELECT * FROM read_csv_auto(
    '${SCORES_DIR}/${PGS_ID}/${PGS_ID}_PRS_total.tsv',
    delim=chr(9), header=true
  )
) TO '${SCORES_DIR}/${PGS_ID}/${PGS_ID}_PRS_total.parquet'
(FORMAT PARQUET, COMPRESSION ZSTD)"
```

Size comparison (PGS000363 example, 140k samples × 24 columns):
- TSV uncompressed: ~40 MB
- Parquet ZSTD:     ~12 MB  (3× reduction)
- Parquet SNAPPY:   ~18 MB

### Option B: Python one-liner in compute script

```python
import duckdb
duckdb.execute("""
    COPY (SELECT * FROM read_csv_auto(?, delim='\t', header=true))
    TO ? (FORMAT PARQUET, COMPRESSION ZSTD)
""", [tsv_path, parquet_path])
```

Add this at the end of `compute_prs_spark_gpu.py` (after the TSV write).

---

## Slurm scripts reference (Phase 1 review)

| Script                       | Purpose                          | Array?  | GPU?  | Spark? |
|------------------------------|----------------------------------|---------|-------|--------|
| `prs_gpu_compute.sbatch`     | Per-chrom GPU scoring            | 1–22    | Yes   | No     |
| `run_prs_spark_gpu.sbatch`   | Spark matching + CuPy scoring    | No      | Yes   | Yes    |
| `compute_prs.py`             | CPU/DuckDB scoring               | —       | No    | No     |
| `compute_prs_spark_gpu.py`   | Spark+GPU backend                | —       | Yes   | Yes    |
| `prs_end2end.sbatch`         | GWAS clumping + scoring v2       | No      | ?     | No     |
| `prs_end2end_v3.sbatch`      | GWAS clumping + scoring v3       | No      | ?     | No     |

Recommended for PGS Catalog scores: `run_prs_spark_gpu.sbatch` or `prs_gpu_compute.sbatch`.

Typical submission for a new PGS from catalog:
```bash
sbatch --export=PGS_ID=PGS000004 \
       --array=1-22 \
       /path/to/scripts/prs_gpu_compute.sbatch
```

---

## Phase 2: Web-triggered pipeline preparation

In Phase 2, the `/api/pipeline/<pgs_id>/plan` endpoint can be extended to:

1. **Generate a manifest JSON** (no side effects):
   ```json
   {
     "pgs_id": "PGS000004",
     "betamap_path": "/data/PGS000004/PGS000004_hmPOS_GRCh38.betamap.tsv.gz",
     "zarr_base": "/mnt/.../zar_files/chr{N}.zarr",
     "output_dir": "/data/PGS000004/",
     "script": "run_prs_spark_gpu.sbatch",
     "slurm_account": "researchers",
     "slurm_partition": "gpu"
   }
   ```

2. **Write manifest to `/work/{pgs_id}_manifest.json`** (writable volume, not `/data`).

3. **Phase 3 (future)**: An operator reads the manifest and submits:
   ```bash
   sbatch --export=$(jq -r 'to_entries|map("\(.key)=\(.value)")|join(",")' manifest.json) \
     run_prs_spark_gpu.sbatch
   ```

This keeps the web app from ever running `sbatch` itself while still automating the workflow.

---

## Dependency map

```
pgscat download PGS000004         →  betamap.tsv.gz
betamap.tsv.gz + chr{N}.zarr      →  chr{N}_scores.tsv  (22 jobs)
chr{1..22}_scores.tsv             →  PRS_total.tsv
PRS_total.tsv                     →  PRS_total.parquet   ← NEW (enables dashboard)
PRS_total.parquet                 →  /api/data/PGS000004  ← web reads this
```
