# pgscat — Anotador de variantes + dashboard PGS/PRS

Aplicación web (Bottle · DuckDB · PostgreSQL 16 · Gunicorn) para explorar
resultados de **Polygenic Scores** y un **motor de anotación de variantes**
contra bases externas (ClinVar, gnomAD, dbNSFP, dbSNP, elementos regulatorios).
Es la versión completa; `app_pgs/` es un dashboard más ligero derivado de ésta.

> Documentación detallada y *quick start*: **[`app/README.md`](app/README.md)**
> y **[`app/docs/`](app/docs/)** (`ARCHITECTURE.md`, `VARIANT_ANNOTATION.md`,
> `PIPELINE_INTEGRATION.md`, `HOT_NVME_DATA_INVENTORY.md`).
>
> **Regla clave:** los valores de *score* por muestra nunca entran a PostgreSQL;
> viven en Parquet y solo los consulta DuckDB.

## Estructura (`app/`)

| Carpeta | Contenido |
|---------|-----------|
| `src/` | App Bottle (`app.py`, `config.py`) y **servicios** (`services/`): catálogo local/remoto, caché, estadísticas Parquet, *gene browser*, anotación de variantes, *ranking*, normalizador, *fetchers* de ClinVar/gnomAD/MCPS, preparador de scores, inspector de pipeline, health probes. |
| `annotator/` | Motor de anotación: `annotate_variants.py`, clasificación, ClinVar/gnomAD/dbNSFP/dbSNP/MCPS freq, parser GFF3, motor FASTA, regulatorios, y `run_annotation.sh` (wrapper Apptainer, gnomAD 4.1.1 + frecuencias MCPS). |
| `scripts/` | Construcción de índices externos y BEDs regulatorios, lector BigBed, diagnóstico PGS, descarga de genes clínicos. |
| `resources/carrier_screen/` | *Carrier screening* (loader + datos). |
| `apptainer/` | Definición del contenedor del anotador (`variant_annotator.def`). |
| `deployment/` | `docker-compose.yml` para el despliegue. |
| `sql/` | Esquema PostgreSQL (`schema.sql`). |
| `docs/`, `examples/` | Documentación de arquitectura y ejemplos. |
