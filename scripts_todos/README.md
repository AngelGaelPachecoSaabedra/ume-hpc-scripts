# scripts_todos — Scripts de depuración / diagnóstico

Scripts sueltos de *debugging* y diagnóstico puntual (no forman un pipeline).
Útiles como referencia para reproducir problemas concretos de LDpred2/Zarr.

| Script | Qué hace |
|--------|----------|
| `debug_ldpred2_chr22_v8.sbatch` | Reproduce LDpred2 en chr22 aislado (contenedor `bigsnpr-1.12.21.sif`) para depurar convergencia. |
| `debug_minimal.sbatch` | Caso mínimo con datos simulados que *sí* convergen (control positivo de `bigsnpr`). |
| `debug_mle.sbatch` | Prueba el estimador MLE de `bigsnpr` con un GWAS simulado inline. |
| `diagnostico_chrX.sbatch` | Diagnóstico de la codificación de chrX y del software PLINK2 disponible. |
| `extract_variants.sbatch` | Extrae el mapa de variantes (`CHR,POS,REF,ALT`) de un Zarr a CSV. |
| `final_conversion.sbatch` | Conversión a Zarr con captura forense (variante de diagnóstico). |
| `clump_threshold_prs.sbatch` | PRS por *clumping + thresholding* (C+T) sobre el GWAS. |
