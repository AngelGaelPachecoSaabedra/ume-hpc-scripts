# refpanel-pipeline — Panel de referencia para imputación (MCPS10k)

Pipeline por etapas para construir un **panel de referencia de imputación**
(MCPS10k, chr1-22 + chrX): concatenar/normalizar, QC, *phasing* (común y raro),
conversión de formatos y validación. Las etapas se encadenan con dependencias
Slurm (`afterok`).

## Scripts (`scripts/`)

| # | Script | Qué hace |
|---|--------|----------|
| 0 | `config.sh` | Configuración global: contenedores, rutas de entrada, directorios de salida y *scratch*. Lo *sourcean* los demás scripts. |
| 0 | `00_setup.sh` | Crea la estructura de directorios de salida (`msav`, `bref3`, `vcf`, `imp5_chunks`, `qc`, logs). |
| 1 | `01_concat_norm.sbatch` | Concatena y normaliza los VCF de entrada (bcftools). |
| 2 | `02_filter_qc.sbatch` | Filtrado y control de calidad de variantes. |
| 3 | `03_phase_common.sbatch` | *Phasing* de variantes comunes (SHAPEIT/Beagle). |
| 3 | `03_phase_common_version2.sbatch` | Variante alternativa del *phasing* común. |
| 4 | `04_phase_rare.sbatch` | *Phasing* de variantes raras. |
| 5 | `05_convert_formats.sbatch` | Conversión a los formatos del panel (`.msav`, `.bref3`, etc.). |
| 6 | `06_validate.sh` | Validación final: cuenta archivos y variantes, reporta PASS/FAIL. |
| — | `submit_all.sh` | **Orquestador:** envía todas las etapas con dependencias Slurm encadenadas. |
