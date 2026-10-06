# _root_working_scripts — Copias de trabajo (raíz del respaldo)

Scripts que vivían sueltos en la raíz del respaldo del clúster. Son
mayormente **iteraciones y variantes** de los pipelines documentados en las
carpetas de proyecto (muchas versiones `v2..v8`, *sex-specific*, *genome-wide*,
recuperación y diagnóstico). Se conservan por su valor histórico y porque
varias son versiones **únicas** que no están duplicadas en otro sitio.

> Para la versión "limpia" y documentada de cada pipeline, ver:
> `LDPRED2_CUSTOM/`, `ldpredhapmap3/`, `prs_bayes_mcps/`, `prs_diabetes/`.

## Contenedores Apptainer (`.def`)

| Archivo | Qué define |
|---------|-----------|
| `bigsnpr-ldpred2-v1108.def` | Imagen con `bigsnpr` 1.10.8 (LDpred2). |
| `bigsnpr_nativo.def` | Imagen con `bigsnpr` compilado nativo. |
| `bgenix.def` | Imagen con `bgenix` (indexado de ficheros BGEN). |

## SBayesRC (Pipeline A / D)

| Script | Qué hace |
|--------|----------|
| `sbrc_02_ldm_sbayesrc.sbatch` … `_v5.sbatch` | Iteraciones v1–v5 de la etapa 2 (bloques LD + MCMC) de SBayesRC. |
| `sbrc_fix2_snpid_allele_order.sbatch` | Fix de orden de alelos en el SNP ID. |
| `sbrc_fix3_liftover_hg19_to_hg38.sbatch` | *LiftOver* hg19→hg38 de T2DGGI. |
| `sbrc_mcps_generate_ma_sexspecific.sbatch` | Genera el `.ma` **por sexo** (Pipeline D). |
| `sbrc_mcps_run_sexspecific.sbatch` | Corre SBayesRC **por sexo**. |
| `compare_sbayesrc_runs.sbatch` | Compara los betas de dos modelos SBayesRC. |
| `validate_sbrc_fix.sbatch` | Verifica que los *scores* resultantes existan/sean válidos. |
| `verify_ld_panel.sh` | Comprueba qué panel LD se usó en un job SBayesRC. |

## LDpred2 custom (Pipeline C) — pasos y variantes

| Script | Qué hace |
|--------|----------|
| `custom_step0_pruning.sbatch` | Step 0: LD-pruning r²<0.9 en MCPS (una sola vez). |
| `custom_step1_blocks.sbatch` | Step 1: bloques LD por cromosoma (array paralelo). |
| `custom_step2_ldpred2.sbatch` | Step 2: LDpred2-auto genome-wide + PRS. |
| `custom_step0_sexspecific.sbatch` | Step 0 por sexo (QC + pruning r²<0.85). |
| `custom_step1_sexspecific.sbatch` | Step 1 por sexo. |
| `custom_step2_sexspecific.sbatch` | Step 2 por sexo (LDpred2-auto + PRS). |
| `custom_ldsc_bysex_rg.sbatch` | LDSC in-sample por sexo + correlación genética rg(M,F). |
| `ldpred2_custom_mcps_v6.sbatch` | LDpred2 custom MCPS v6 (chr 1-22 + X). |
| `ldpred2_custom_mcps_v7.sbatch` | LDpred2 custom MCPS v7 (genome-wide). |

## LDpred2 genome-wide / recuperación / HapMap3

| Script | Qué hace |
|--------|----------|
| `ldpred2_gwide.sbatch` | LDpred2 genome-wide con `snp_ldsplit` (bloques LD independientes). |
| `ldpred2_gwide_step1.sbatch` | Paso 1/2: bloques LD por cromosoma (paralelo). |
| `ldpred2_gwide_step2_genomewide.sbatch` | Paso 2/2: ensambla SFBM genome-wide + LDpred2-auto + PRS. |
| `ldpred2_gwide_step2_v3.sbatch` | Paso 2/2 v3: LDpred2-auto por cromosoma + suma de PRS. |
| `ldpred2_resplit_blocks.sbatch` | Re-particiona bloques LD grandes (>MAX_BLOCK SNPs). |
| `ldpred2_recover.sbatch` | Recuperación: LDpred2-auto sin recalcular las matrices LD. |
| `ldpred2_v8_hapmap3.sbatch` | LDpred2-auto v8 sobre HapMap3. |
| `ld_check_indels_v3.sbatch` | Chequeo de indels en el panel LD. |

## Diagnóstico y pruebas

| Script | Qué hace |
|--------|----------|
| `debug_ldpred2_chr22.sbatch` | chr22 con CERO paralelismo (evita el `assert_cores`). |
| `debug_ldpred2_chr22_v5.sbatch` | chr22 con el fix maestro `safe_as_SFBM`. |
| `diag_bigsnpr_slurm.sbatch` | Qué variables de SLURM ve `bigsnpr` dentro del contenedor. |
| `diag_neff_chr22.sbatch` | Qué `n_eff` hace converger LDpred2-auto. |
| `check_x.sbatch` | Comprobación del cromosoma X. |
| `explore_chrx.py` | Ubica hombres/mujeres a partir de covariables (chrX). |
| `test_exact_params.sh` | Test del fix de paralelismo BLAS. |
| `test_fieldswap_final.sh` | Confirma el *field-swap* con `.sbk` + recomputo de `p`/`first_i`. |
| `test_sfbm_reattach.sh` | Diagnóstico `bigsparser` 0.7.3 (evita crash en `externalptr`). |
| `test_dask.sbatch` | Prueba del entorno Dask. |

## Zarr / reindex / PRS / otros

| Script | Qué hace |
|--------|----------|
| `vcf2zarr_wes_wgs.sbatch` | Conversión VCF→Zarr (WES/WGS) con captura forense. |
| `reindex_all.sbatch` | Re-indexa todos los Zarr/ficheros. |
| `reindex_bulletproof.sbatch` | Re-indexado robusto (a prueba de fallos). |
| `reindex_fix.sbatch` | Re-indexado en modo estricto. |
| `index_script.sh` | Indexado de ficheros (bgenix/tabix). |
| `run_prs_dask.sbatch` | PRS con Dask (requiere `PGS_ID`). |
| `download_gtex_v8.sh` | Descarga de los QTLs de GTEx Analysis v8. |
| `vllm_cluster.slurm` | Levanta un clúster de inferencia LLM con **vLLM** sobre Slurm. |
