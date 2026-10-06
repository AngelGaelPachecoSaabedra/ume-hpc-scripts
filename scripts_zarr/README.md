# scripts_zarr — Conversión de VCF a Zarr (bio2zarr)

Jobs Slurm para convertir VCF (WES/WGS, MCPS-TOPMed) a **Zarr** con `bio2zarr`,
usando un entorno `uv` en `uv-envs/vcf2zarr`. Los Zarr resultantes alimentan
los pipelines de PRS que leen dosages por cromosoma.

## Scripts (`jobs/`)

| Script | Qué hace |
|--------|----------|
| `run_bio2zarr.sbatch` | Conversión VCF→Zarr **multinodo** con bio2zarr (activa el venv `uv`). |
| `run_bio2zarr2.sbatch` | Igual que el anterior con `set -euo pipefail` estricto (segunda iteración). |
| `vcf2zarr_chr10.sbatch` | Conversión de un solo cromosoma (chr10) como ejemplo/validación puntual. |
