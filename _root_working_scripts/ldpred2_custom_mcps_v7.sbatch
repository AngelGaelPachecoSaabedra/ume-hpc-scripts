#!/bin/bash
#SBATCH --job-name=ldp2_v7_gw
#SBATCH --account=researchers
#SBATCH --qos=vip
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=150G
#SBATCH --output=/mnt/cephfs/orgs/home/angel.pacheco/logs/ldp2_v7_%j.out
#SBATCH --error=/mnt/cephfs/orgs/home/angel.pacheco/logs/ldp2_v7_%j.err
#SBATCH --mail-user=angelpachecosaavedra@gmail.com
#SBATCH --mail-type=BEGIN,END,FAIL

set -euo pipefail

# ================================================================
# PIPELINE C — LDpred2 CUSTOM MCPS v7 (GENOME-WIDE)
#
# Cambios v7 (vs v6):
#   1. Arquitectura Genome-Wide: snp_ldpred2_auto se ejecuta SOBRE TODO 
#      el genoma a la vez, no por cromosoma. Evita sobreajuste local.
#   2. Pruning más estricto: PLINK a r2 < 0.8 para estabilizar MCMC.
#   3. MCMC iteraciones suben a 1000/2000 para convergencia global.
#   4. SFBM monolítico en el NVMe ensamblado mediante list() en RAM.
# ================================================================

SIF="/mnt/cephfs/biocontainers/images/jupyter-biotools-1.5.sif"
PLINK_DIR="/mnt/cephfs/hot_nvme/mcps/imputed-topmed/plink_files/maximally_unrelated"
GWAS_DIR="/mnt/cephfs/orgs/home/angel.pacheco/gwas_miami"
OUT_DIR="/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes/ldpred2"
mkdir -p "$OUT_DIR"

# Importante: Asegurar espacio en NVMe scratch
TMPDIR="/mnt/cephfs/m.2/scratch/$USER/ldp2_v7_tmp/$SLURM_JOB_ID"
mkdir -p "$TMPDIR"; chmod 700 "$TMPDIR"
export TMPDIR TMP="$TMPDIR" TEMP="$TMPDIR"

echo "============================================================"
echo "  PIPELINE C — LDpred2 CUSTOM MCPS v7 (GENOME-WIDE)"
echo "  Nodo: $SLURMD_NODENAME   RAM: $(free -g | awk 'NR==2{print $2}')G"
echo "  $(date)"
echo "============================================================"

R_SCRIPT="$TMPDIR/run_v7.R"
cat << 'REOF' > "$R_SCRIPT"
Sys.setenv(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
           MKL_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1",
           BLIS_NUM_THREADS="1", NUMEXPR_NUM_THREADS="1")
suppressMessages({ library(bigsnpr); library(data.table); library(Matrix); library(methods) })
options(bigstatsr.check.parallel.blas=FALSE)
options(bigstatsr.ncores.max=24L)

ncores    <- 32L
out_dir   <- "/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes/ldpred2"
gwas_dir  <- "/mnt/cephfs/orgs/home/angel.pacheco/gwas_miami"
plink_dir <- "/mnt/cephfs/hot_nvme/mcps/imputed-topmed/plink_files/maximally_unrelated"
tmp_dir   <- Sys.getenv("TMPDIR")

TARGET_BLOCK <- 2500L  

# ── [1] GWAS meta-análisis IVW ─────────────────────────────────────
cat("\n[1] Cargando GWAS Meta-análisis...\n")
gf <- fread(file.path(gwas_dir, "combined-gwas-results-all-ctrl-females.txt.gz"),
            select=c("CHROM","GENPOS","ID","ALLELE0","ALLELE1","A1FREQ","INFO","N","BETA","SE","P"), showProgress=FALSE)
gm <- fread(file.path(gwas_dir, "combined-gwas-results-all-ctrl-males.txt.gz"),
            select=c("CHROM","GENPOS","ID","ALLELE0","ALLELE1","A1FREQ","INFO","N","BETA","SE","P"), showProgress=FALSE)
gwas <- merge(gf, gm, by="ID", suffixes=c(".f",".m"))
rm(gf, gm); gc()

w_f <- 1/gwas$SE.f^2; w_m <- 1/gwas$SE.m^2
neff_f <- 4 * 17300 * 54196 / (17300 + 54196)
neff_m <- 4 * 8522 * 26733 / ( 8522 + 26733)
neff_ldpred2 <- neff_f + neff_m 
cat(sprintf("  n_eff = %.0f\n", neff_ldpred2))

sumstats <- data.frame(
    chr = gwas$CHROM.f, pos = gwas$GENPOS.f,
    a1 = gwas$ALLELE1.f, a0 = gwas$ALLELE0.f,
    freq = (gwas$A1FREQ.f*gwas$N.f + gwas$A1FREQ.m*gwas$N.m)/(gwas$N.f+gwas$N.m),
    info = pmin(gwas$INFO.f, gwas$INFO.m),
    beta = (gwas$BETA.f*w_f + gwas$BETA.m*w_m)/(w_f+w_m),
    beta_se = sqrt(1/(w_f+w_m)), n_eff = neff_ldpred2, stringsAsFactors = FALSE)
sumstats$p <- 2*pnorm(-abs(sumstats$beta/sumstats$beta_se))
rm(gwas); gc()

sumstats <- sumstats[sumstats$info>=0.8 & sumstats$freq>=0.01 & sumstats$freq<=0.99 &
                     !is.na(sumstats$beta) & !is.na(sumstats$beta_se) & sumstats$beta_se>0, ]

chr_raw <- as.character(sumstats$chr)
chr_raw[chr_raw %in% c("X", "chrX", "23")] <- "23"
sumstats$chr <- suppressWarnings(as.integer(chr_raw))
sumstats <- sumstats[!is.na(sumstats$chr) & sumstats$chr >= 1 & sumstats$chr <= 23, ]

# ── [2] LD-Pruning r²<0.8 ──────────────────────────────────────────
cat("\n[2] LD-Pruning r²<0.8 en datos MCPS...\n")
gwas_ids_file <- file.path(tmp_dir, "gwas_ids_for_pruning.txt")
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
    
    # UMBRAL AJUSTADO A 0.8
    system(sprintf(
        "plink2 --bfile %s --extract %s --indep-pairwise 500 50 0.8 --out %s --threads %d --memory 8000 2>/dev/null",
        prefix, gwas_ids_file, prune_out, ncores))
        
    prune_in <- paste0(prune_out, ".prune.in")
    if (file.exists(prune_in)) {
        ids <- readLines(prune_in)
        pruned_ids_all <- c(pruned_ids_all, ids)
    }
}
cat(sprintf("  Total filtro MCPS (0.8): %s variantes\n", format(length(pruned_ids_all), big.mark=",")))

pp <- strsplit(pruned_ids_all, ":")
pruned_chr_raw <- gsub("chr", "", sapply(pp, `[`, 1))
pruned_chr_raw[pruned_chr_raw == "X"] <- "23"
pruned_chr <- suppressWarnings(as.integer(pruned_chr_raw))
pruned_pos <- as.integer(sapply(pp, `[`, 2))
ok_parse   <- !is.na(pruned_chr) & !is.na(pruned_pos)
pruned_key <- paste(pruned_chr[ok_parse], pruned_pos[ok_parse], sep=":")
gwas_key   <- paste(sumstats$chr, sumstats$pos, sep=":")
sumstats   <- sumstats[gwas_key %in% pruned_key, ]

# Helpers
subdiv <- function(blk, target) {
  m <- nrow(blk); if (m <= target) return(list(blk))
  n_sub <- ceiling(m / target); sizes <- rep(floor(m/n_sub), n_sub)
  extra <- m - sum(sizes); if (extra > 0) sizes[seq_len(extra)] <- sizes[seq_len(extra)] + 1L
  starts <- c(1L, cumsum(sizes[-n_sub]) + 1L); ends <- cumsum(sizes)
  lapply(seq_len(n_sub), function(i) blk[starts[i]:ends[i], starts[i]:ends[i], drop=FALSE])
}
is_valid_split <- function(rs) {
  if (is.null(rs) || nrow(rs)==0) return(FALSE)
  any(!is.na(rs$cost) & !is.na(rs$n_block) & rs$n_block > 0)
}

# ── [3] Computar Matrices de Correlación (Iterativo) ───────────────
cat("\n[3] Calculando Correlaciones por Cromosoma...\n")
list_corr <- list()
df_beta_global <- data.frame()
info_snp_list <- list()

for (chr in 1:23) {
    chr_label <- if (chr == 23L) "X" else as.character(chr)
    prefix <- file.path(plink_dir, sprintf("mcps-freeze150k_qcd_chr%s_ivs", chr_label))
    bed_file <- paste0(prefix, ".bed")
    if (!file.exists(bed_file)) next
    sumstats_chr <- sumstats[sumstats$chr == chr, ]
    if (nrow(sumstats_chr) < 50) next

    bk_file <- file.path(tmp_dir, sprintf("chr%s", chr_label))
    snp_readBed(bed_file, backingfile=bk_file)
    obj <- snp_attach(paste0(bk_file, ".rds"))
    G <- obj$genotypes
    
    if (chr == 1) fam_global <- obj$fam 

    n_na <- big_counts(G, ind.col=seq_len(ncol(G)))[4, ] |> sum()
    if (n_na > 0) G <- snp_fastImputeSimple(G, "mean0", ncores=ncores)

    map_df <- data.frame(chr=obj$map$chromosome, pos=obj$map$physical.pos,
                         a1=obj$map$allele1, a0=obj$map$allele2, stringsAsFactors=FALSE)
    info_snp <- snp_match(sumstats_chr, map_df, strand_flip=FALSE, join_by_pos=TRUE, match.min.prop=0.01)
    
    ind_col <- info_snp$`_NUM_ID_`
    sd_vals <- big_scale()(G, ind.col=ind_col)$scale
    keep_var <- which(sd_vals > 1e-6)
    if (length(keep_var) < nrow(info_snp)) { info_snp <- info_snp[keep_var, ]; ind_col <- info_snp$`_NUM_ID_` }

    cat(sprintf("  chr%s: procesando %s vars...\n", chr_label, format(length(ind_col), big.mark=",")))
    
    corr0 <- snp_cor(Gna=G, ind.col=ind_col, infos.pos=info_snp$pos/1000, size=3000, ncores=ncores, thr_r2=0.005)
    m_chr <- ncol(corr0)

    max_size_seq <- unique(round(seq_log(m_chr/30, m_chr/5, length.out=20)))
    max_size_seq <- max_size_seq[max_size_seq >= 100 & max_size_seq <= m_chr]
    if (length(max_size_seq) == 0) max_size_seq <- min(m_chr, 5000)

    res_split <- NULL
    for (thr in c(0.02, 0.05, 0.1, 0.2, 0.3)) {
        rs <- tryCatch(snp_ldsplit(corr0, thr_r2=thr, min_size=50, max_size=max_size_seq, max_K=500), error=function(e) NULL)
        if (is_valid_split(rs)) { res_split <- rs; break }
    }

    if (!is.null(res_split)) {
        vi <- which(!is.na(res_split$cost) & !is.na(res_split$n_block) & res_split$n_block>0)
        bi <- which.min(res_split[vi, ]$cost + 0.5*log(res_split[vi, ]$n_block))
        block_sizes <- res_split[vi, ]$all_size[[bi]]
        nb <- length(block_sizes)
        bs <- c(1L, cumsum(block_sizes[-nb])+1L); be <- cumsum(block_sizes)
        blocks0 <- lapply(seq_len(nb), function(b) corr0[bs[b]:be[b], bs[b]:be[b], drop=FALSE])
    } else {
        blocks0 <- list(corr0)
    }

    blocks <- list()
    for (b in blocks0) blocks <- c(blocks, subdiv(b, TARGET_BLOCK))
    
    # Guardar en RAM como matriz dispersa global
    list_corr[[chr]] <- as(Matrix::bdiag(blocks), "symmetricMatrix")
    
    # Acumular sumstats globales
    df_beta_global <- rbind(df_beta_global, data.frame(beta=info_snp$beta, beta_se=info_snp$beta_se, n_eff=neff_ldpred2))
    info_snp_list[[chr_label]] <- info_snp

    # LIMPIEZA CRÍTICA DE RAM
    rm(obj, G, corr0, blocks0, blocks); gc()
}

# ── [4] Ensamblar SFBM Global en NVMe ──────────────────────────────
cat("\n[4] Ensamblando SFBM Genómico en disco...\n")
sfbm_global <- as_SFBM(list_corr, backingfile=file.path(tmp_dir, "corr_genomewide"), compact=TRUE)
rm(list_corr); gc() # Liberar las matrices en RAM

# ── [5] MCMC Genome-Wide ───────────────────────────────────────────
cat(sprintf("\n[5] MCMC Genome-Wide (%s variantes totales)...\n", format(nrow(df_beta_global), big.mark=",")))
t0 <- Sys.time()

ma <- snp_ldpred2_auto(corr=sfbm_global, df_beta=df_beta_global, h2_init=0.1,
                       vec_p_init=seq_log(1e-4, 0.5, length.out=30),
                       burn_in=1000L, num_iter=2000L, report_step=100L,
                       ncores=ncores, sparse=TRUE, allow_jump_sign=FALSE, shrink_corr=0.95)

h2_est <- sapply(ma, function(x) x$h2_est)
p_est  <- sapply(ma, function(x) x$p_est)
keep <- which(!is.na(h2_est) & h2_est>0.0001 & h2_est<0.8 & !is.na(p_est) & p_est>1e-6 & p_est<0.99)

if (length(keep) < 3L) {
    cat("  RETRY global con shrink_corr=0.5\n")
    ma <- snp_ldpred2_auto(corr=sfbm_global, df_beta=df_beta_global, h2_init=0.05,
                           vec_p_init=seq_log(1e-5, 0.3, length.out=30),
                           burn_in=1000L, num_iter=2000L, report_step=100L,
                           ncores=ncores, sparse=TRUE, allow_jump_sign=FALSE, shrink_corr=0.5)
    keep <- which(!is.na(sapply(ma, function(x) x$h2_est)))
}

bm <- sapply(ma[keep], function(x) x$beta_est)
if (is.null(dim(bm))) bm <- matrix(bm, ncol=1)
beta_final_gw <- rowMeans(bm, na.rm=TRUE); beta_final_gw[!is.finite(beta_final_gw)] <- 0

el <- as.numeric(difftime(Sys.time(), t0, units="mins"))
cat(sprintf("  MCMC Genómico finalizado en %.1f min. h2 global=%.5f\n", el, median(sapply(ma[keep], function(x) x$h2_est))))

# ── [6] Calcular PRS Total ─────────────────────────────────────────
cat("\n[6] Calculando PRS...\n")
prs_total <- rep(0, nrow(fam_global))

idx_start <- 1
for (chr in 1:23) {
    chr_label <- if (chr == 23L) "X" else as.character(chr)
    info_snp <- info_snp_list[[chr_label]]
    if (is.null(info_snp)) next
    
    n_vars <- nrow(info_snp)
    idx_end <- idx_start + n_vars - 1
    beta_chr <- beta_final_gw[idx_start:idx_end]
    idx_start <- idx_end + 1
    
    if (sum(beta_chr != 0) == 0) next
    
    bk_file <- file.path(tmp_dir, sprintf("chr%s.rds", chr_label))
    obj <- snp_attach(bk_file)
    prs_chr <- big_prodVec(X=obj$genotypes, y.col=beta_chr, ind.col=info_snp$`_NUM_ID_`, ncores=ncores)
    prs_total <- prs_total + prs_chr
    rm(obj); gc()
}

# ── [7] Exportar ───────────────────────────────────────────────────
cat("\n[7] Exportando resultados...\n")
covar_file <- "/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes/MCPS BASELINE.csv"
link_file  <- "/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes/RGN_LINK_IID.csv"
df_link <- read.csv(link_file)
iid_to_patid <- setNames(as.character(df_link$PATID), as.character(df_link$IID))
df_baseline <- read.csv(covar_file); df_baseline <- df_baseline[!is.na(df_baseline$MALE), ]
patid_to_sex <- setNames(as.integer(df_baseline$MALE), as.character(df_baseline$PATID))
patids <- iid_to_patid[as.character(fam_global$sample.ID)]; patids[is.na(patids)] <- "UNKNOWN"
sexes <- patid_to_sex[patids]; sexes[is.na(sexes)] <- -1L

df_out <- data.frame(PATID=patids, IID=fam_global$sample.ID, SEX=sexes, PRS_TOTAL=prs_total, stringsAsFactors=FALSE)
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

out_file <- file.path(out_dir, "prs_diabetes_ldpred2_custom_mcps_v7_gw.tsv")
write.table(df_out, out_file, sep="\t", row.names=FALSE, quote=FALSE)
saveRDS(list(h2_est=h2_est, p_est=p_est, keep=keep), file.path(out_dir, "ldpred2_v7_gw_mcmc.rds"))

cat(sprintf("\n ✓ Pipeline v7 completado. Output: %s\n", out_file))
REOF

apptainer exec --bind /mnt/cephfs:/mnt/cephfs --bind "$TMPDIR:/tmp,$TMPDIR:/var/tmp" \
    --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 --env MKL_NUM_THREADS=1 \
    --env BLIS_NUM_THREADS=1 --env NUMEXPR_NUM_THREADS=1 \
    "$SIF" /opt/miniforge3/envs/r-bio/bin/Rscript "$R_SCRIPT"

rm -rf "$TMPDIR"
echo "============================================================"
echo "  PIPELINE C v7 COMPLETADO — $(date)"
echo "============================================================"
