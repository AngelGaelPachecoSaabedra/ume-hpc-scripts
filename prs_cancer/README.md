# prs_cancer — PRS de cáncer (deciles de riesgo)

| Script | Qué hace |
|--------|----------|
| `fase_fg_deciles.py` | **Fases F+G** del PRS oncológico: une el PRS por muestra (`prs_cancer_mcps_all_samples.tsv`) con la cohorte caso-control oncológica construida en `Ordanamiento_Bases/` (`MCPS_Cohorte_Casos_Controles_Oncologicos_F145K.csv`), calcula deciles de riesgo, *odds ratios* por decil y genera las gráficas. Cubre 6 scores. |

> El *scoring* previo que genera el `.tsv` de entrada está en
> `Ordanamiento_Bases/scripts/` (`prs_cancer_scoring.sbatch`,
> `prs_cancer_deciles.sbatch`).
