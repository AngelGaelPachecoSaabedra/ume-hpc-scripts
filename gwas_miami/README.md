# gwas_miami — Miami / Manhattan plots de GWAS

Generación de **Miami plots** (Manhattan espejado: machos arriba, hembras
abajo) del GWAS sex-stratified de MCPS (sumstats REGENIE). Hay dos rutas: una
**acelerada por GPU** (CuPy + Datashader) en dos fases, y una **en R puro**
con formato tipo *Nature*.

## `miamiR/Gpu/` — ruta GPU (dos fases)

| # | Script | Qué hace |
|---|--------|----------|
| 1 | `miami_phase1_gpu.py` | **Fase 1 (cómputo):** calcula las coordenadas genómicas acumulativas en GPU (CuPy, `float64` para evitar overflow), ordena por CHR/POS y persiste los puntos. |
| 2 | `miami_phase2_plot.py` | **Fase 2 (plotting):** rasteriza con Datashader y dibuja el Miami plot (línea de umbral configurable, etiquetas de genes). |
| — | `miami_gpu.py` | Versión monolítica (cómputo GPU + Datashader + matplotlib en un solo script). |
| — | `miami_gpu_slurm.sh` | Job Slurm que lanza la ruta GPU dentro del contenedor `jupyter-biotools-1.2.sif`. |

## `miamiR/` — ruta R (Nature)

| Script | Qué hace |
|--------|----------|
| `miami_plot_nature.R` | Pipeline completo en R (sin GPU): carga GWAS F+M, coordenadas acumulativas, anotación GTF (gen más cercano por SNP significativo), selección de top SNPs por gen y Miami plot de panel único conforme a especificaciones Nature (183×150 mm, Wong 2011, sin gridlines). |
| `miami_nature_slurm.sh` | Job Slurm que ejecuta `miami_plot_nature.R`. |

## Raíz

| Script | Qué hace |
|--------|----------|
| `miami_plot.R` | Variante de `miami_plot_nature.R` que ensambla los paneles con `patchwork` (incluye *panel label* 8 pt). |
