#!/bin/bash
#SBATCH --job-name=ldp2_custom
#SBATCH --account=researchers
#SBATCH --qos=vip
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=150G
#SBATCH --output=/mnt/cephfs/orgs/home/angel.pacheco/logs/ldp2_custom_%j.out
#SBATCH --error=/mnt/cephfs/orgs/home/angel.pacheco/logs/ldp2_custom_%j.err
#SBATCH --mail-user=angelpachecosaavedra@gmail.com
#SBATCH --mail-type=BEGIN,END,FAIL

set -euo pipefail

# ================================================================
# PIPELINE C — LDpred2 CUSTOM MCPS PRUNED v6 (Chr 1-22 + X)
#
# Cambios v6 (vs v5):
#   1. n_eff CORRECTO = 4/(1/ncas+1/ncon) = 78,304 (SIN ×K(1-K)).
#      La doc oficial de Privé define n_eff binario así; el ×K(1-K)
#      de v5 NO es una corrección reconocida (era parche numérico).
#   2. snp_ldsplit con retry escalado (thr_r2: 0.02→0.3) + re-partición
#      forzada a bloques ≤2,500 SNPs → evita divergencia del MCMC por
#      bloques gigantes (lección del genome-wide HapMap3).
#   3. Contenedor jupyter-biotools-1.5 (bigsnpr 1.10.8).
#   4. Imputación defensiva de NAs antes de snp_cor.
#   5. h2 se transformará a liability con Lee 2011 (post-proceso).
#   6. chrX (chr23) incluido: los bfiles ya tienen hombres pseudo-
#      diploidizados (0/2) — verificado vía HAP_REF_CT/HAP_ALT_CT en PLINK.
#      snp_readBed lo lee como 0/2; NO requiere --haploids ni recodificar.
# ================================================================

SIF="/mnt/cephfs/biocontainers/images/jupyter-biotools-1.5.sif"
PLINK_DIR="/mnt/cephfs/hot_nvme/mcps/imputed-topmed/plink_files/maximally_unrelated"
GWAS_DIR="/mnt/cephfs/orgs/home/angel.pacheco/gwas_miami"
OUT_DIR="/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes/ldpred2"
mkdir -p "$OUT_DIR"

TMPDIR="/mnt/cephfs/m.2/scratch/$USER/ldp2_custom_tmp/$SLURM_JOB_ID"
mkdir -p "$TMPDIR"; chmod 700 "$TMPDIR"
export TMPDIR TMP="$TMPDIR" TEMP="$TMPDIR"

echo "============================================================"
echo "  PIPELINE C — LDpred2 CUSTOM MCPS PRUNED v6"
echo "  Nodo: $SLURMD_NODENAME   RAM: $(free -g | awk 'NR==2{print $2}')G"
echo "  $(date)"
echo "============================================================"

R_SCRIPT="$TMPDIR/run_custom_v6.R"
cat << 'REOF' > "$R_SCRIPT"
Sys.setenv(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
           MKL_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1",
           BLIS_NUM_THREADS="1", NUMEXPR_NUM_THREADS="1")
suppressMessages({ library(bigsnpr); library(data.table); library(Matrix); library(methods) })
options(bigstatsr.check.parallel.blas=FALSE)
options(bigstatsr.ncores.max=24L)

bv <- packageVersion("bigsnpr")
cat(sprintf("  bigsnpr version: %s\n", as.character(bv)))
if (as.character(bv) == "1.12.21") stop("bigsnpr 1.12.21 roto. Usa 1.10.8.")

ncores    <- 12L
out_dir   <- "/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes/ldpred2"
gwas_dir  <- "/mnt/cephfs/orgs/home/angel.pacheco/gwas_miami"
plink_dir <- "/mnt/cephfs/hot_nvme/mcps/imputed-topmed/plink_files/maximally_unrelated"
tmp_dir   <- Sys.getenv("TMPDIR")

# Parámetros de re-partición (lección del genome-wide)
TARGET_BLOCK <- 2500L  # bloques ≤ esto → MCMC estable

# ── [1] GWAS meta-análisis IVW ─────────────────────────────────────
cat("\n[1] Cargando GWAS Meta-análisis...\n")
gf <- fread(file.path(gwas_dir, "combined-gwas-results-all-ctrl-females.txt.gz"),
            select=c("CHROM","GENPOS","ID","ALLELE0","ALLELE1","A1FREQ","INFO","N","BETA","SE","P"),
            showProgress=FALSE)
gm <- fread(file.path(gwas_dir, "combined-gwas-results-all-ctrl-males.txt.gz"),
            select=c("CHROM","GENPOS","ID","ALLELE0","ALLELE1","A1FREQ","INFO","N","BETA","SE","P"),
            showProgress=FALSE)
gwas <- merge(gf, gm, by="ID", suffixes=c(".f",".m"))
rm(gf, gm); gc()

w_f <- 1/gwas$SE.f^2; w_m <- 1/gwas$SE.m^2

# ══════════════════════════════════════════════════════════════════
# n_eff CORRECTO (doc oficial Privé): 4/(1/ncas + 1/ncon) por estrato.
# NO multiplicar por K(1-K). El n_eff binario ya maneja escala logística.
# La transformación a liability se hace DESPUÉS sobre h2 (Lee 2011).
# ══════════════════════════════════════════════════════════════════
neff_f <- 4 * 17300 * 54196 / (17300 + 54196)   # 52,456
neff_m <- 4 *  8522 * 26733 / ( 8522 + 26733)   # 25,848
neff_ldpred2 <- neff_f + neff_m                  # 78,304
cat(sprintf("  n_eff = %.0f (binario correcto, sin K(1-K))\n", neff_ldpred2))

sumstats <- data.frame(
    chr = gwas$CHROM.f, pos = gwas$GENPOS.f,
    a1 = gwas$ALLELE1.f, a0 = gwas$ALLELE0.f,
    freq = (gwas$A1FREQ.f*gwas$N.f + gwas$A1FREQ.m*gwas$N.m)/(gwas$N.f+gwas$N.m),
    info = pmin(gwas$INFO.f, gwas$INFO.m),
    beta = (gwas$BETA.f*w_f + gwas$BETA.m*w_m)/(w_f+w_m),
    beta_se = sqrt(1/(w_f+w_m)), n_eff = neff_ldpred2,
    stringsAsFactors = FALSE)
sumstats$p <- 2*pnorm(-abs(sumstats$beta/sumstats$beta_se))
rm(gwas); gc()

sumstats <- sumstats[sumstats$info>=0.8 & sumstats$freq>=0.01 & sumstats$freq<=0.99 &
                     !is.na(sumstats$beta) & !is.na(sumstats$beta_se) & sumstats$beta_se>0, ]
# Autosomas (1-22) + chrX (codificado internamente como 23).
# chrX en los bfiles ya está pseudo-diploidizado (hombres 0/2) — verificado
# previamente vía HAP_REF_CT/HAP_ALT_CT en PLINK. snp_readBed lo lee como
# 0/2 automáticamente, así que NO requiere --haploids ni recodificación.
chr_raw <- as.character(sumstats$chr)
chr_raw[chr_raw %in% c("X", "chrX", "23")] <- "23"
sumstats$chr <- suppressWarnings(as.integer(chr_raw))
sumstats <- sumstats[!is.na(sumstats$chr) & sumstats$chr >= 1 & sumstats$chr <= 23, ]
cat(sprintf("  GWAS tras QC: %s variantes (%s autosomas + %s chrX)\n",
    format(nrow(sumstats), big.mark=","),
    format(sum(sumstats$chr <= 22), big.mark=","),
    format(sum(sumstats$chr == 23), big.mark=",")))

# ── [2] LD-Pruning en MCPS (r²<0.9) ────────────────────────────────
cat("\n[2] LD-Pruning r²<0.9 en datos MCPS...\n")
gwas_ids_file <- file.path(tmp_dir, "gwas_ids_for_pruning.txt")
# Etiqueta de cromosoma para los IDs de PLINK: 1..22 literal, 23 -> "X"
chr_lab_vec <- ifelse(sumstats$chr == 23L, "X", as.character(sumstats$chr))
ids_fwd <- paste0("chr", chr_lab_vec, ":", sumstats$pos, ":", sumstats$a1, ":", sumstats$a0)
ids_rev <- paste0("chr", chr_lab_vec, ":", sumstats$pos, ":", sumstats$a0, ":", sumstats$a1)
writeLines(unique(c(ids_fwd, ids_rev)), gwas_ids_file)
rm(ids_fwd, ids_rev, chr_lab_vec); gc()

pruned_ids_all <- c()
for (chr in 1:23) {
    chr_label <- if (chr == 23L) "X" else as.character(chr)
    prefix <- file.path(plink_dir, sprintf("mcps-freeze150k_qcd_chr%s_ivs", chr_label))
    if (!file.exists(paste0(prefix, ".bed"))) next
    prune_out <- file.path(tmp_dir, sprintf("prune_chr%s", chr_label))
    system(sprintf(
        "plink2 --bfile %s --extract %s --indep-pairwise 500 50 0.9 --out %s --threads %d --memory 8000 2>/dev/null",
        prefix, gwas_ids_file, prune_out, ncores))
    prune_in <- paste0(prune_out, ".prune.in")
    if (file.exists(prune_in)) {
        ids <- readLines(prune_in)
        pruned_ids_all <- c(pruned_ids_all, ids)
        cat(sprintf("    chr%s: %s variantes\n", chr_label, format(length(ids), big.mark=",")))
    }
}
cat(sprintf("  Total filtro MCPS: %s variantes\n", format(length(pruned_ids_all), big.mark=",")))

# Parsear IDs: "chr1:..." o "chrX:..." → chr numérico (X -> 23)
pp <- strsplit(pruned_ids_all, ":")
pruned_chr_raw <- gsub("chr", "", sapply(pp, `[`, 1))
pruned_chr_raw[pruned_chr_raw == "X"] <- "23"
pruned_chr <- suppressWarnings(as.integer(pruned_chr_raw))
pruned_pos <- as.integer(sapply(pp, `[`, 2))
ok_parse   <- !is.na(pruned_chr) & !is.na(pruned_pos)
pruned_key <- paste(pruned_chr[ok_parse], pruned_pos[ok_parse], sep=":")
gwas_key   <- paste(sumstats$chr, sumstats$pos, sep=":")
sumstats   <- sumstats[gwas_key %in% pruned_key, ]
cat(sprintf("  GWAS filtrado: %s variantes (%s autosomas + %s chrX)\n",
    format(nrow(sumstats), big.mark=","),
    format(sum(sumstats$chr <= 22), big.mark=","),
    format(sum(sumstats$chr == 23), big.mark=",")))

# ── Helpers: re-partición y SFBM seguro ────────────────────────────
subdiv <- function(blk, target) {
  m <- nrow(blk)
  if (m <= target) return(list(blk))
  n_sub <- ceiling(m / target)
  sizes <- rep(floor(m/n_sub), n_sub)
  extra <- m - sum(sizes); if (extra > 0) sizes[seq_len(extra)] <- sizes[seq_len(extra)] + 1L
  starts <- c(1L, cumsum(sizes[-n_sub]) + 1L); ends <- cumsum(sizes)
  lapply(seq_len(n_sub), function(i) blk[starts[i]:ends[i], starts[i]:ends[i], drop=FALSE])
}

is_valid_split <- function(rs) {
  if (is.null(rs) || nrow(rs)==0) return(FALSE)
  if (!is.numeric(rs$cost) || !is.numeric(rs$n_block)) return(FALSE)
  any(!is.na(rs$cost) & !is.na(rs$n_block) & rs$n_block > 0)
}

# ── [3] LDpred2-auto por cromosoma (con ldsplit + re-partición) ────
cat("\n[3] MCMC y Scoring por Cromosoma...\n")
prs_total <- NULL; fam_global <- NULL; matched_total <- 0L; results_log <- list()

for (chr in 1:23) {
    chr_label <- if (chr == 23L) "X" else as.character(chr)
    prefix <- file.path(plink_dir, sprintf("mcps-freeze150k_qcd_chr%s_ivs", chr_label))
    bed_file <- paste0(prefix, ".bed")
    if (!file.exists(bed_file)) next
    sumstats_chr <- sumstats[sumstats$chr == chr, ]
    if (nrow(sumstats_chr) < 50) { cat(sprintf("  chr%s: <50 vars, skip\n", chr_label)); next }

    cat(sprintf("\n  >>> chr%s: %s variantes GWAS %s<<<\n", chr_label,
        format(nrow(sumstats_chr), big.mark=","),
        if (chr == 23L) "(chrX, hombres pseudo-diploides 0/2) " else ""))
    t0 <- Sys.time()

    bk_file <- file.path(tmp_dir, sprintf("chr%s", chr_label))
    snp_readBed(bed_file, backingfile=bk_file)
    obj <- snp_attach(paste0(bk_file, ".rds"))
    G <- obj$genotypes
    if (is.null(prs_total)) { prs_total <- rep(0, nrow(G)); fam_global <- obj$fam }

    # Imputación defensiva
    n_na <- big_counts(G, ind.col=seq_len(ncol(G)))[4, ] |> sum()
    if (n_na > 0) { cat(sprintf("    NAs=%s, imputando\n", format(n_na,big.mark=","))); G <- snp_fastImputeSimple(G, "mean0", ncores=ncores) }

    map_df <- data.frame(chr=obj$map$chromosome, pos=obj$map$physical.pos,
                         a1=obj$map$allele1, a0=obj$map$allele2, stringsAsFactors=FALSE)
    info_snp <- snp_match(sumstats_chr, map_df, strand_flip=FALSE, join_by_pos=TRUE, match.min.prop=0.01)
    cat(sprintf("    Matched: %s\n", format(nrow(info_snp), big.mark=",")))
    if (nrow(info_snp) < 50) { rm(obj, G); gc(); next }

    ind_col <- info_snp$`_NUM_ID_`
    sd_vals <- big_scale()(G, ind.col=ind_col)$scale
    keep_var <- which(sd_vals > 1e-6)
    if (length(keep_var) < nrow(info_snp)) {
        info_snp <- info_snp[keep_var, ]; ind_col <- info_snp$`_NUM_ID_`
    }
    if (nrow(info_snp) < 50) { rm(obj, G); gc(); next }
    matched_total <- matched_total + nrow(info_snp)

    cat(sprintf("    snp_cor (%s vars)...\n", format(length(ind_col), big.mark=",")))
    corr0 <- snp_cor(Gna=G, ind.col=ind_col, infos.pos=info_snp$pos/1000,
                     size=3000, ncores=ncores, thr_r2=0.005)
    m_chr <- ncol(corr0)

    # ── snp_ldsplit con retry escalado ──
    max_size_seq <- round(seq_log(m_chr/30, m_chr/5, length.out=20))
    max_size_seq <- unique(max_size_seq[max_size_seq >= 100 & max_size_seq <= m_chr])
    if (length(max_size_seq) == 0) max_size_seq <- min(m_chr, 5000)

    res_split <- NULL
    for (thr in c(0.02, 0.05, 0.1, 0.2, 0.3)) {
        rs <- tryCatch(snp_ldsplit(corr0, thr_r2=thr, min_size=50,
                       max_size=max_size_seq, max_K=500),
                       error=function(e) NULL)
        if (is_valid_split(rs)) { res_split <- rs; thr_used <- thr; break }
    }

    if (is.null(res_split)) {
        block_sizes <- m_chr; blocks0 <- list(corr0)
        cat("    ldsplit falló → 1 bloque (se re-particionará)\n")
    } else {
        vi <- which(!is.na(res_split$cost) & !is.na(res_split$n_block) & res_split$n_block>0)
        res_split <- res_split[vi, ]
        bi <- which.min(res_split$cost + 0.5*log(res_split$n_block))
        block_sizes <- res_split$all_size[[bi]]
        nb <- length(block_sizes)
        bs <- c(1L, cumsum(block_sizes[-nb])+1L); be <- cumsum(block_sizes)
        blocks0 <- lapply(seq_len(nb), function(b) corr0[bs[b]:be[b], bs[b]:be[b], drop=FALSE])
        cat(sprintf("    ldsplit: %d bloques (thr_r2=%.3f, max=%d)\n", nb, thr_used, max(block_sizes)))
    }

    # ── Re-partición forzada a ≤TARGET_BLOCK ──
    blocks <- list()
    for (b in blocks0) blocks <- c(blocks, subdiv(b, TARGET_BLOCK))
    new_sizes <- sapply(blocks, nrow)
    stopifnot(sum(new_sizes) == m_chr)
    cat(sprintf("    Re-particionado: %d bloques (max=%d)\n", length(blocks), max(new_sizes)))

    # SFBM
    corr_chr <- as(Matrix::bdiag(blocks), "symmetricMatrix")
    rm(corr0, blocks0, blocks); gc()
    sfbm <- as_SFBM(corr_chr, backingfile=file.path(tmp_dir, sprintf("corr_chr%s", chr_label)), compact=TRUE)
    rm(corr_chr); gc()

    # LDpred2-auto (n_eff correcto = 78,304)
    df_beta <- data.frame(beta=info_snp$beta, beta_se=info_snp$beta_se, n_eff=neff_ldpred2)
    ma <- snp_ldpred2_auto(corr=sfbm, df_beta=df_beta, h2_init=0.1,
        vec_p_init=seq_log(1e-4, 0.5, length.out=30),
        burn_in=500L, num_iter=500L, report_step=100L,
        ncores=ncores, sparse=TRUE,
        allow_jump_sign=FALSE, shrink_corr=0.95)

    h2_est <- sapply(ma, function(x) x$h2_est)
    p_est  <- sapply(ma, function(x) x$p_est)
    keep <- which(!is.na(h2_est) & h2_est>0.0001 & h2_est<0.8 & !is.na(p_est) & p_est>1e-6 & p_est<0.99)

    if (length(keep) < 3L) {
        cat(sprintf("    %d/30 cadenas — RETRY shrink_corr=0.5\n", length(keep)))
        ma <- snp_ldpred2_auto(corr=sfbm, df_beta=df_beta, h2_init=0.05,
            vec_p_init=seq_log(1e-5, 0.3, length.out=30),
            burn_in=1000L, num_iter=1000L, report_step=100L,
            ncores=ncores, sparse=TRUE, allow_jump_sign=FALSE, shrink_corr=0.5)
        h2_est <- sapply(ma, function(x) x$h2_est); p_est <- sapply(ma, function(x) x$p_est)
        keep <- which(!is.na(h2_est) & h2_est>0.0001 & h2_est<0.8 & !is.na(p_est) & p_est>1e-6 & p_est<0.99)
    }

    if (length(keep) < 1L) {
        cat(sprintf("    chr%s DIVERGIÓ (0/30) — beta=0\n", chr_label))
        beta_final <- rep(0, nrow(info_snp)); h2_med <- NA; p_med <- NA
    } else {
        bm <- sapply(ma[keep], function(x) x$beta_est)
        if (is.null(dim(bm))) bm <- matrix(bm, ncol=1)
        beta_final <- rowMeans(bm, na.rm=TRUE); beta_final[!is.finite(beta_final)] <- 0
        h2_med <- median(h2_est[keep]); p_med <- median(p_est[keep])
    }

    prs_chr <- big_prodVec(X=G, y.col=beta_final, ind.col=info_snp$`_NUM_ID_`, ncores=ncores)
    prs_total <- prs_total + prs_chr
    el <- as.numeric(difftime(Sys.time(), t0, units="mins"))
    cat(sprintf("    h2=%.5f, p=%.2e, cadenas=%d/30, nonzero=%s, %.1f min\n",
        ifelse(is.na(h2_med),0,h2_med), ifelse(is.na(p_med),0,p_med),
        length(keep), format(sum(beta_final!=0), big.mark=","), el))

    results_log[[chr_label]] <- list(n_matched=nrow(info_snp),
        h2=h2_med, p=p_med, n_valid=length(keep),
        n_nonzero=sum(beta_final!=0), elapsed_min=el)
    rm(obj, G, sfbm, ma, info_snp, beta_final, prs_chr); gc()
}

# ── [4] Exportar ───────────────────────────────────────────────────
cat(sprintf("\n[4] Exportando (%s variantes)...\n", format(matched_total, big.mark=",")))
cat("\n  Resumen por cromosoma:\n")
for (cn in names(results_log)) {
    r <- results_log[[cn]]
    cat(sprintf("    chr%s: %s vars, h2=%.5f, p=%.2e, %d/30 cadenas, %s nonzero, %.1f min\n",
        cn, format(r$n_matched, big.mark=","), ifelse(is.na(r$h2),0,r$h2),
        ifelse(is.na(r$p),0,r$p), r$n_valid, format(r$n_nonzero, big.mark=","), r$elapsed_min))
}

covar_file <- "/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes/MCPS BASELINE.csv"
link_file  <- "/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes/RGN_LINK_IID.csv"
df_link <- read.csv(link_file)
iid_to_patid <- setNames(as.character(df_link$PATID), as.character(df_link$IID))
df_baseline <- read.csv(covar_file); df_baseline <- df_baseline[!is.na(df_baseline$MALE), ]
patid_to_sex <- setNames(as.integer(df_baseline$MALE), as.character(df_baseline$PATID))
patids <- iid_to_patid[as.character(fam_global$sample.ID)]; patids[is.na(patids)] <- "UNKNOWN"
sexes <- patid_to_sex[patids]; sexes[is.na(sexes)] <- -1L

df_out <- data.frame(PATID=patids, IID=fam_global$sample.ID, SEX=sexes,
                     PRS_TOTAL=prs_total, stringsAsFactors=FALSE)
df_out <- df_out[df_out$SEX >= 0, ]
df_out$PRS_Z <- NA_real_; df_out$PERCENTILE <- NA_real_
for (sv in c(0,1)) {
    mask <- df_out$SEX==sv; vals <- df_out$PRS_TOTAL[mask]
    if (length(vals)>0 && sd(vals)>0) {
        df_out$PRS_Z[mask] <- (vals-mean(vals))/sd(vals)
        df_out$PERCENTILE[mask] <- rank(vals)/length(vals)*100
    }
}
df_out <- df_out[order(-df_out$PRS_TOTAL), ]

out_file <- file.path(out_dir, "prs_diabetes_ldpred2_custom_mcps_final.tsv")
write.table(df_out, out_file, sep="\t", row.names=FALSE, quote=FALSE)
saveRDS(list(results_log=results_log, n_eff=neff_ldpred2, matched_total=matched_total),
        file.path(out_dir, "ldpred2_custom_results_log.rds"))

cat(sprintf("\n ✓ Custom MCPS completado: %s indiv, %s variantes\n",
    format(nrow(df_out), big.mark=","), format(matched_total, big.mark=",")))
cat(sprintf("   Output: %s\n", out_file))
REOF

apptainer exec --bind /mnt/cephfs:/mnt/cephfs --bind "$TMPDIR:/tmp,$TMPDIR:/var/tmp" \
    --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 --env MKL_NUM_THREADS=1 \
    --env BLIS_NUM_THREADS=1 --env NUMEXPR_NUM_THREADS=1 \
    "$SIF" /opt/miniforge3/envs/r-bio/bin/Rscript "$R_SCRIPT"

rm -rf "$TMPDIR"
echo "============================================================"
echo "  PIPELINE C v6 COMPLETADO — $(date)"
echo "============================================================"
