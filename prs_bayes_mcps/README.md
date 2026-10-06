# prs_bayes_mcps — PRS bayesiano con SBayesRC (Pipeline A)

PRS bayesiano con **SBayesRC** (GCTB, contenedor `sbayesrc-0.2.6.sif`) usando
el panel LD *in-house* de MCPS (`ldm_79k`, N=79,612 maximally unrelated).
Incluye la descarga/formateo de los sumstats, la construcción del `.ma`, los
arreglos de concordancia (orden de alelos, *liftover* hg19→hg38) y la
validación clínica comparativa.

## Scripts y orden de ejecución

| # | Script | Qué hace |
|---|--------|----------|
| 1 | `sbrc_01_download_format.sbatch` | Etapa 1: formatea los sumstats T2DGGI y los convierte a formato COJO `.ma` para GCTB/SBayesRC. |
| 2 | `sbrc_mcps_generate_ma.sbatch` | Genera el `.ma` de MCPS casando el meta-análisis sex-stratified (REGENIE F+M) contra el panel LD, con `data.table` merge (evita el segfault del hash de 43 M elementos). |
| 3 | `sbrc_fix2_snpid_allele_order.sbatch` | **Fix:** reconstruye el SNP ID probando ambas orientaciones de alelos (`REF:ALT` vs `ALT:REF`) contra el `.bim` para resolver los no-matches por orden de alelos. |
| 4 | `sbrc_fix3_liftover_hg19_to_hg38.sbatch` | **Fix:** *liftOver* de T2DGGI de hg19→hg38 (T2DGGI está en GRCh37, MCPS-TOPMed en GRCh38) y reconstruye el `.ma` con coordenadas hg38. |
| 5 | `sbrc_02_ldm_sbayesrc_v6.sbatch` | Etapa 2 (definitiva): construye los bloques LD (`LDstep3` bloque a bloque, arma `snp.info` con columna `Block`) y corre `sbayesrc()`. |
| 6 | `sbrc_mcps_run.sbatch` | Ejecuta SBayesRC con el GWAS MCPS sobre el panel LD MCPS; salida en `results_mcps/sbayesrc_mcps.{txt,par,rds,mcmcsamples}`. |
| 7 | `validate_sbrc_comparison.sbatch` | Validación clínica comparativa: Pipeline A (T2DGGI) vs A2 (GWAS MCPS) contra el fenotipo de diabetes (HbA1c ≥ 6.5% + diagnóstico + medicación). |

> Variantes *sex-specific* de este pipeline (Pipeline D) están en
> `_root_working_scripts/` (`sbrc_mcps_*_sexspecific.sbatch`), junto con las
> iteraciones v1–v5 de la etapa 2.
