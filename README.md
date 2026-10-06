# UME HPC — Scripts de bioinformática y cómputo de alto rendimiento

Colección curada de scripts de los pipelines de investigación genómica ejecutados en
el clúster HPC (Slurm) de la UME: cálculo de **Polygenic Risk Scores (PRS)**, **GWAS**,
construcción de **paneles de referencia / imputación**, anotación de variantes,
preparación de fenotipos y utilidades de cómputo distribuido (Spark + GPU, Dask, DDP).

> **Nota sobre datos.** Este repositorio contiene **únicamente código y documentación**.
> No incluye datos de sujetos, pesos por individuo, genotipos, imágenes, modelos ni
> credenciales. Las rutas a datos (`/mnt`, `/home`, NVMe, Ceph, etc.) son parámetros que
> deben ajustarse al entorno donde se ejecuten.

## Estructura

| Carpeta | Descripción |
|---|---|
| `LDPRED2_CUSTOM/` | Pipeline LDpred2 "custom" (v6): poda/pruning, cálculo por bloques LD, merge final. Jobs Slurm (launcher → prune → worker → merge). |
| `ldpredhapmap3/` | LDpred2 sobre el panel HapMap3: launcher paralelo, worker por bloque y merge final. |
| `prs_bayes_mcps/` | PRS bayesiano con **SBayesRC** para MCPS: descarga/formato de sumstats, generación del `.ma`, construcción de la matriz LD, ejecución y arreglos (orden de alelos, liftover hg19→hg38) y validación comparativa. |
| `PRS_T2D_Project/` | Proyecto PRS de diabetes tipo 2: chequeo de alelos (incl. chrX), scoring end-to-end, validación en R y gráficas (global y por sexo). |
| `prs_diabetes/` | Variante Spark+GPU del PRS de diabetes: extracción de variantes, match con Spark, cómputo en GPU, scoring por cromosoma y merge. |
| `prs_cancer/` | PRS de cáncer: cálculo de deciles de riesgo (modelo Fine–Gray). |
| `Ordanamiento_Bases/` | Preparación y ordenamiento de bases clínicas/genómicas: fases 1–5 (mortalidad, denominador genómico, consistencia de IDs, construcción de fenotipo, sitio de cáncer, exportación final) y una subcarpeta `scripts/` con el pipeline PRS end-to-end, comparación contra PGS Catalog y gráficas en R. |
| `gwas_miami/` | Generación de **Miami/Manhattan plots** de GWAS: fase 1 de cómputo en GPU, fase 2 de graficado, y versiones en R con formato tipo *Nature*. |
| `refpanel-pipeline/` | Pipeline de panel de referencia para imputación: concatenar/normalizar, QC, phasing (común y raro con SHAPEIT/Beagle), conversión de formatos y validación. |
| `pgscat/` | Aplicación de anotación de variantes y cliente de **PGS Catalog**: backend (Flask/gunicorn), motor de anotación (ClinVar, gnomAD, dbNSFP, dbSNP, regulatorios), *gene browser*, preparación de scores y despliegue (Apptainer / docker-compose). Incluye `docs/`. |
| `app_pgs/` | Dashboard PGS más ligero (subconjunto de `pgscat`): app Flask, servicios (DB, caché, estadísticas de Parquet, inspector de pipeline) y documentación de arquitectura. |
| `scripts_zarr/` | Conversión de VCF (WES/WGS) a **Zarr** con `bio2zarr`: jobs Slurm por cromosoma. |
| `spark/` | Utilidades Spark sobre Slurm: match de CADD/PRS y ejemplo de cluster Spark. |
| `training_ddp/` | Entrenamiento distribuido PyTorch **DDP** multi-GPU (job con `uv`). |
| `espirometrias/` | *Matching* de 4 fuentes de datos de espirometría: versiones simple, paralela, MPI y distribuida, con sus jobs Slurm. |
| `vep_annotation/` | Anotación funcional con **VEP**. |
| `Scripts_compartidos/` | Scripts compartidos entre proyectos (conversión a Zarr, PRS Spark+GPU, match CADD, cliente PGS Catalog). |
| `old_scripts/` | Versiones previas del pipeline LDpred2 (histórico). |
| `scripts_todos/` | Scripts sueltos de depuración/diagnóstico (LDpred2 chr22, chrX, conversiones, etc.). |
| `_root_working_scripts/` | Copias de trabajo que vivían sueltas en la raíz del respaldo (muchas son variantes/versiones de los scripts anteriores). Se conservan tal cual para no perder historial; pueden depurarse más adelante. |

## Tecnologías

- **Cómputo:** Slurm, Apptainer/Singularity (`.def`), MPI, Dask, Apache Spark, PyTorch DDP, GPU (CUDA).
- **Genómica:** PLINK2, bigsnpr/LDpred2, SBayesRC, PRS-CS/CSx, VEP, bcftools, SHAPEIT/Beagle, bio2zarr, IMPUTE5.
- **Lenguajes:** Python, R, Bash.
- **Apps:** Flask + gunicorn, PostgreSQL/DuckDB, lectura de Parquet.

## Convenciones

- Los `*.sbatch` / `*.slurm` son scripts de envío a Slurm; revisa cabeceras `#SBATCH`
  (particiones, GPUs, memoria, tiempo) antes de ejecutar.
- Hay archivos con sufijos de versión (`_v2`…`_v8`, `_fix`, `_recover`) que documentan la
  evolución de cada pipeline; se mantienen intencionalmente.
