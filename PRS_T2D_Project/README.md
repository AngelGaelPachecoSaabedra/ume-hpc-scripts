# PRS_T2D_Project — PRS de diabetes tipo 2 (end-to-end)

Proyecto PRS de diabetes tipo 2 con validación clínica en R y gráficas por
sexo. Incluye chequeo de concordancia de alelos (incluido chrX). Los scripts
viven en `1_scripts/` y usan `jupyter-biotools-1.2.sif`.

## Scripts

| Script | Qué hace |
|--------|----------|
| `check_alleles_gwas_zarr.sbatch` | Verifica la concordancia de alelos entre el GWAS y el Zarr de genotipos antes de calcular el PRS. |
| `check_alleles_chrX.sbatch` | Igual que el anterior pero para el cromosoma X (codificación especial de hombres). |
| `prs_end2end.sbatch` | Pipeline **end-to-end**: matching de variantes → scoring del PRS → salida por muestra. |
| `val_prs_R_completo.sbatch` | Validación clínica completa en R (asociación PRS ↔ fenotipo, métricas, curvas). |
| `plot_prs_mcps.sbatch` | Histograma/distribución del PRS en la cohorte MCPS. |
| `plot_prs_amr.sbatch` | Gráfica del PRS para la población AMR (global). |
| `plot_prs_amr_sex.sbatch` | Gráfica del PRS para AMR estratificada por sexo. |
