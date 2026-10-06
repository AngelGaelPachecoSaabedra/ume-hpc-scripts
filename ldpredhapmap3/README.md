# ldpredhapmap3 — LDpred2 sobre panel HapMap3

Variante del pipeline LDpred2-auto (`bigsnpr` 1.10.8) que usa el panel LD
**HapMap3** y procesa los 22 autosomas **en paralelo por cromosoma** con Slurm.
Reaprovecha matrices LD (`.sbk`/`.rds`) calculadas en una corrida previa.

## Scripts y orden de ejecución

| # | Script | Qué hace |
|---|--------|----------|
| 1 | `ldpred2_parallel_launcher.sbatch` | Verifica que existan los 22 `corr_chr*.sbk` + `chr*.rds` de la corrida previa y lanza el *array* `1-22` con 3 tareas concurrentes (24 cores c/u). |
| 2 | `ldpred2_worker.sbatch` | **Worker v4** por cromosoma: QC estricto `SD_ss` vs `SD_val`, SFBM fresco, `N_eff` corregido, LDpred2-auto con 50 cadenas y filtro `h2` relajado. Escribe `ldpred2_betas_chr{N}.rds`. |
| 3 | `ldpred2_merge_final.sbatch` | **Merge final**: combina los 22 `.rds` por cromosoma en el PRS total (correr solo cuando todos los workers terminen). |
