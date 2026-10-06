# Scripts_compartidos — Scripts compartidos entre proyectos

Utilidades que usan varios pipelines y que no pertenecen a un único proyecto.

| Script | Qué hace |
|--------|----------|
| `pgs_catalog_client.py` | Cliente Python de la API REST del **PGS Catalog** (Scores, Publications, Traits/EFO, Performance Metrics, Sample Sets). Usable como menú interactivo o como módulo importable. |
| `convert_to_zarr/plink_update_to_zarr.sbatch` | Conversión PLINK→Zarr con captura forense de logs y aislamiento de `TMPDIR` en *scratch* (maneja los terabytes temporales en disco mecánico). |
| `convert_to_zarr/final_conversion.sbatch` | Conversión final a Zarr en "modo paranoico" (log forense en tiempo real + aislamiento de TMP). |
| `prs_zarr/prs_end2end_v3.sbatch` | PRS *end-to-end* v3 leyendo genotipos desde Zarr (clump+threshold → scoring). |

> Las versiones de matching Spark (CADD/PRS) que antes se duplicaban aquí viven
> ahora en `spark/` y `prs_diabetes/`.
