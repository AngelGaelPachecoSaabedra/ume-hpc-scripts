# LDPRED2_CUSTOM — Pipeline C (LDpred2 "custom" v6/v7)

Pipeline propio de **LDpred2-auto** sobre el panel LD *in-house* de MCPS
(no HapMap3), calculado **por cromosoma en paralelo** con Slurm. Pensado para
la cohorte admixta MCPS, donde el LD es más extenso que en poblaciones EUR.

Todos los jobs usan el contenedor `jupyter-biotools-1.4.sif` (R + `bigsnpr`
1.10.8) y escriben los betas por cromosoma en
`prs_diabetes/ldpred2_custom/`.

## Scripts y orden de ejecución

| # | Script | Qué hace |
|---|--------|----------|
| 1 | `ldp2_custom_v6_prune.sbatch` | **Paso previo (1 sola vez, ~30 min).** Pruning por LD de los sumstats en MCPS; genera `sumstats_pruned.rds` y `pruned_ids.rds` que alimentan a los workers. |
| 2 | `ldp2_custom_v6_launcher.sbatch` | Lanza el *array* Slurm `1-23` (22 autosomas + chrX=23) con 3 tareas concurrentes. Valida que el pruning haya terminado antes de lanzar. |
| 3 | `ldp2_custom_v6_worker.sbatch` | **Worker por cromosoma.** Corre LDpred2-auto en un cromosoma: `n_eff` per-SNP calibrado con `sd_val` (Privé 2022), `h2_init=0.01`, 50 cadenas para cromosomas difíciles, `allow_jump_sign=FALSE`. Escribe `ldpred2_custom_betas_chr{N}.rds`. |
| 4 | `ldp2_custom_v6_merge.sbatch` | **Merge final** (correr solo cuando los 23 workers terminan): combina los `.rds` por cromosoma en el PRS total y produce un resumen (h2, nº de cadenas, variantes casadas). |

> Nota: el worker incluido corresponde a la revisión **v7** (fixes acumulados
> sobre v6). Las variantes sex-specific y genome-wide de este mismo pipeline
> viven en `_root_working_scripts/` (`custom_step*`, `ldpred2_gwide_*`).
