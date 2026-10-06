# vep_annotation — Anotación funcional con VEP

| Script | Qué hace |
|--------|----------|
| `run_vep.sbatch` | Ejecuta **Ensembl VEP 115** (contenedor `ensembl-vep-115.sif`) con los plugins REVEL, CADD, SpliceAI y PrimateAI. Entrada: VCF; salida: anotación en TSV. |

Requiere el `.sif` de VEP y los directorios `plugins/`, `input/` y `results/`
en el entorno de ejecución (no incluidos en el repo).
