# prs_diabetes — PRS de diabetes (Spark + GPU)

Variante del PRS de diabetes tipo 2 (**PGS000014**) que combina **PySpark**
para el cruce de variantes y **CuPy (GPU)** para el *scoring* vectorizado,
leyendo los genotipos desde Zarr. Los contenedores son `jupyter-biotools-1.2.sif`.

## Scripts y orden de ejecución

| # | Script | Qué hace |
|---|--------|----------|
| 1 | `extract_variants_prs.sbatch` | Extrae el mapa de variantes (`CHR,POS,REF,ALT`) del Zarro de MCPS y limpia el archivo PGS Catalog (`PGS000014_hmPOS_GRCh38`) a hg38. |
| 2 | `match_prs_spark.sbatch` | Levanta un clúster **Spark distribuido** sobre Slurm y cruza el catálogo de variantes del Zarr contra los pesos PGS (coincidencia posición + alelo). |
| 3 | `prs_gpu_compute.sbatch` | *Array* por cromosoma: calcula el PRS parcial en **GPU** (CuPy) a partir de los dosages del Zarr y los pesos casados. |
| 4 | `prs_per_chromosome.sbatch` | Produce el PRS por cromosoma **sin sumar** (útil para diagnóstico/validación por cromosoma). |
| 5 | `prs_merge_final.sbatch` | Suma los PRS por cromosoma en el PRS total por muestra. |

> La versión CPU/secuencial de este mismo cómputo (`compute_prs.py`,
> `dask_runner.py`) y la lógica Spark de matching (`spark/`) viven en sus
> propias carpetas.
