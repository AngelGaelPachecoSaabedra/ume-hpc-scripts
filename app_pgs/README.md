# app_pgs — Dashboard PGS/PRS (versión ligera)

Dashboard PGS/PRS más ligero (subconjunto de `pgscat/`), centrado en explorar
los scores de la cohorte MCPS. Mismo stack (Bottle · DuckDB · PostgreSQL 16 ·
Gunicorn · Podman) y comparte por diseño los módulos de servicio con `pgscat/`
(cada app mantiene su propia copia para ser autónoma; por eso **no** se
deduplicaron entre sí).

> Documentación y *quick start*: **[`app_pgs/README.md`](app_pgs/README.md)**
> y **[`app_pgs/docs/`](app_pgs/docs/)** (`ARCHITECTURE.md`,
> `PIPELINE_INTEGRATION.md`, `TODO.md`).

## Estructura (`app_pgs/`)

| Carpeta | Contenido |
|---------|-----------|
| `src/app.py`, `src/config.py` | App Bottle y configuración (variables `APP_*`). |
| `src/services/` | Catálogo local/remoto, caché, DB, estadísticas Parquet, cliente pgscat, inspector de pipeline, sincronización, health probes. |
| `sql/schema.sql` | Esquema PostgreSQL. |
| `docs/` | Arquitectura e integración con el pipeline. |
