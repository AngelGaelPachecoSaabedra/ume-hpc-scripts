# old_scripts — Versiones históricas del pipeline LDpred2

Versiones previas (monolíticas) del pipeline LDpred2-auto, conservadas como
referencia histórica. Fueron reemplazadas por las variantes paralelas por
cromosoma de `ldpredhapmap3/` y `LDPRED2_CUSTOM/`. **No usar en producción.**

| Script | Qué hace |
|--------|----------|
| `ldpred2_01_full_pipeline_v2.sbatch` | Pipeline B v2: LDpred2-auto procesando por cromosoma (sin merge PLINK de 390 GB), meta-análisis IVW F+M, chrX como pseudo-diploide. |
| `ldpred2_01_full_pipeline_v3.sbatch` | v3: pre-filtra variantes del GWAS, hace merge PLINK de ~15–20 M variantes y corre LDpred2-auto sobre el `.bed` combinado. |
| `ldpred2_v4_fix.sbatch` | v4 "maestro unificado": gestiona el `.bed` gigante, lo divide por cromosomas en *scratch* y lee secuencialmente para evitar el límite de 1.1 TB. |
| `ldpred2_v4_chrX.sbatch` | v4 con soporte de chrX (chr23): incluye el cromosoma X en el split y en el loop de R. |
