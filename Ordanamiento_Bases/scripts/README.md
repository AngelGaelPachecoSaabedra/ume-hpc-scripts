# Ordanamiento_Bases/scripts — Pipeline PRS end-to-end

Back-end de cálculo de **Polygenic Risk Scores** a partir de genotipos en Zarr,
con comparación contra el PGS Catalog y gráficas en R. Varios scripts aquí son
los que invoca el dashboard `pgscat`/`app_pgs` (vía `APP_SCRIPTS_DIR`).

## Cómputo del PRS

| Script | Qué hace |
|--------|----------|
| `compute_prs.py` | Cómputo del PRS por cromosoma (modelo aditivo `PRS = Σ βᵢ·dosageᵢ`), leyendo el Zarr por ventanas e indexando el matching posición-alelo con DuckDB. |
| `compute_prs_spark_gpu.py` | Misma interfaz de salida que `compute_prs.py`, pero con matching en PySpark y *scoring* vectorizado en GPU (CuPy). |
| `dask_runner.py` | Orquesta `compute_prs.py` sobre todos los cromosomas con Dask. |
| `run_prs_dask.sbatch` | Job Slurm que lanza el *runner* Dask (requiere `PGS_ID`). |
| `run_prs_spark_gpu.sbatch` | Job Slurm del PRS Spark+GPU (integración Spark↔Apptainer, modo `local[N]`). |
| `prs_gpu_compute.sbatch` | *Array* por cromosoma del scoring en GPU. |
| `prs_mcps_gwas.sbatch` | PRS a partir del GWAS propio de MCPS. |
| `check_chr22.py` | Diagnóstico: escanea `chr22.zarr` chunk a chunk para detectar bloques corruptos (fuerza la descompresión Blosc). |

## PRS de cáncer

| Script | Qué hace |
|--------|----------|
| `prs_cancer_scoring.sbatch` | Calcula el PRS oncológico por muestra (genera el `.tsv` que consume `prs_cancer/`). |
| `prs_cancer_deciles.sbatch` | Deciles de riesgo del PRS oncológico. |

## Comparación contra PGS Catalog (PGS000363)

| Script | Qué hace |
|--------|----------|
| `comparar_prs_pgs000363.py` | Compara el PRS nuevo (Spark+CuPy) contra el *baseline* previo; exporta `prs_comparativa_pgs000363.csv`. |
| `graficar_comparativa_prs_pgs000363.R` | 5 gráficas (PNG+PDF) de la comparación. |
| `diff_stats_pgs000363.R` | Estadísticos de la diferencia `PRS_NEW − PRS_OLD`. |
| `plot_pct_diff_pgs000363.R` | Distribución del % de diferencia entre ambos PRS. |

## Validación clínica y gráficas

| Script | Qué hace |
|--------|----------|
| `validate_prs_clinical_R.sbatch` | Validación clínica del PRS en R (asociación con fenotipo). |
| `val_prs_R_completo.sbatch` / `val_prs_R_completo_v4.sbatch` | Validación R completa (v1 / v4). |
| `prs_histogram.sbatch` | Histograma del PRS. |
| `plot_prs_mcps.sbatch` | Distribución del PRS en MCPS. |
| `plot_prs_amr_sex.sbatch` | PRS para AMR estratificado por sexo. |

## Otros pipelines (versiones concretas conservadas aquí)

| Script | Qué hace |
|--------|----------|
| `ldpred2_custom_mcps_v5.sbatch` | LDpred2 custom MCPS *pruned* v5 (chr 1-23), fix de *integer overflow* en `snp_cor`. |
| `ldpred2_v8_hapmap3.sbatch` | LDpred2-auto v8 HapMap3 (`thr_r2=0.005`). |
| `match_cadd_spark.sbatch` | Cross-match CADD con Spark distribuido. |
| `sbrc_01_download_format.sbatch`, `sbrc_02_ldm_sbayesrc_v6.sbatch`, `sbrc_fix2_*`, `sbrc_fix3_*` | Copias del pipeline SBayesRC (ver `prs_bayes_mcps/` para la versión canónica documentada). |
