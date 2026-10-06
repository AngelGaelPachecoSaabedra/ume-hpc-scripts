#!/bin/bash
#SBATCH --job-name=ldp2_custom_chr
#SBATCH --account=researchers
#SBATCH --qos=vip
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=150G
#SBATCH --output=/mnt/cephfs/orgs/home/angel.pacheco/logs/ldp2_custom_%A_%a.out
#SBATCH --error=/mnt/cephfs/orgs/home/angel.pacheco/logs/ldp2_custom_%A_%a.err
#SBATCH --mail-user=angelpachecosaavedra@gmail.com
#SBATCH --mail-type=FAIL

set -euo pipefail

# ================================================================
# PIPELINE C v7 — WORKER: procesa UN cromosoma (1-23)
#
# Fixes (mismos que HapMap3 worker v3):
#   1. n_eff per-SNP calibrado con sd_val (Privé 2022)
#   2. h2_init=0.01 unificado, vec_p amplio
#   3. Filtro h2 relajado (>1e-5) para hard-chrs
#   4. 50 cadenas para hard-chrs
#   5. allow_jump_sign=FALSE (único param extra en v1.10.8)
#   6. Acepta >=1 cadena válida
#   7. SIN shrink_corr ni use_MLE (no existen en v1.10.8)
# ================================================================

CHR=$SLURM_ARRAY_TASK_ID

SIF="/mnt/cephfs/biocontainers/images/jupyter-biotools-1.4.sif"
PLINK_DIR="/mnt/cephfs/hot_nvme/mcps/imputed-topmed/plink_files/maximally_unrelated"
OUT_DIR="/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes/ldpred2_custom"

mkdir -p "$OUT_DIR"

TMPDIR="/mnt/cephfs/scratch/$USER/ldp2_custom_worker_tmp/${SLURM_ARRAY_JOB_ID}_${CHR}"
mkdir -p "$TMPDIR"; chmod 700 "$TMPDIR"
export TMPDIR TMP="$TMPDIR" TEMP="$TMPDIR"

if [[ "$CHR" == "23" ]]; then
  CHR_LABEL="X"
else
  CHR_LABEL="$CHR"
fi

echo "============================================================"
echo "  LDpred2 CUSTOM worker v7 — chr${CHR_LABEL}"
echo "  Job: ${SLURM_ARRAY_JOB_ID}_${CHR}   Nodo: ${SLURMD_NODENAME}"
echo "  $(date)"
echo "============================================================"

R_SCRIPT="$TMPDIR/worker_chr${CHR_LABEL}.R"

cat << REOF > "$R_SCRIPT"
Sys.setenv(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
           MKL_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1")
if (requireNamespace("RhpcBLASctl", quietly=TRUE)) {
  RhpcBLASctl::blas_set_num_threads(1L)
  RhpcBLASctl::omp_set_num_threads(1L)
}
suppressMessages({
  library(bigsnpr); library(bigsparser); library(data.table); library(Matrix)
})
options(bigstatsr.check.parallel.blas=FALSE)
options(bigstatsr.ncores.max=24L)

cat(sprintf("  bigsnpr version: %s\n", as.character(packageVersion("bigsnpr"))))

chr       <- ${CHR}L
chr_label <- "$CHR_LABEL"
ncores    <- 12L
out_dir   <- "$OUT_DIR"
plink_dir <- "$PLINK_DIR"
tmp_dir   <- Sys.getenv("TMPDIR")

cat(sprintf("\n==== chr%s (ncores=%d) ====\n", chr_label, ncores))

betas_file <- file.path(out_dir, sprintf("ldpred2_custom_betas_chr%s.rds", chr_label))
if (file.exists(betas_file)) {
  cat(sprintf("Ya existe %s — saltando.\n", basename(betas_file)))
  quit(save="no", status=0)
}

# ── safe_as_SFBM ──────────────────────────────────────────────────
safe_as_SFBM <- function(corr0, backingfile, compact=TRUE, eps=1e-10) {
  stopifnot(is(corr0, "dsCMatrix"))
  col_nnz <- diff(corr0@p); empty <- which(col_nnz == 0L)
  if (length(empty) == 0L) {
    cat("      [SFBM] Sin columnas vacías\n")
    return(as_SFBM(corr0, backingfile=backingfile, compact=compact))
  }
  n <- corr0@Dim[1]; n_add <- length(empty); old_p <- corr0@p
  old_nnz <- length(corr0@i); new_nnz <- old_nnz + n_add
  cat(sprintf("      [SFBM] %d columnas vacías\n", n_add))
  if (as.double(new_nnz) > .Machine\$integer.max)
    stop("Total nnz excedería límite 32-bit")
  cum_add <- cumsum(col_nnz == 0L)
  new_p   <- old_p + c(0L, cum_add)
  new_i   <- integer(new_nnz); new_x <- numeric(new_nnz)
  boundaries <- c(0L, empty); n_blocks <- length(boundaries)
  for (b in seq_len(n_blocks)) {
    block_first <- boundaries[b] + 1L
    block_last  <- if (b < n_blocks) boundaries[b+1L] - 1L else n
    shift <- b - 1L
    if (block_last >= block_first) {
      old_start <- old_p[block_first] + 1L; old_end <- old_p[block_last+1L]
      if (old_end >= old_start) {
        cnt <- old_end - old_start + 1L; ns <- old_start + shift
        new_i[ns:(ns+cnt-1L)] <- corr0@i[old_start:old_end]
        new_x[ns:(ns+cnt-1L)] <- corr0@x[old_start:old_end]
      }
    }
    if (b <= n_add) {
      ec <- empty[b]; dp <- new_p[ec] + 1L
      new_i[dp] <- ec - 1L; new_x[dp] <- eps
    }
  }
  corr_fixed <- new("dsCMatrix",
    i=new_i, p=new_p, x=new_x, Dim=corr0@Dim,
    Dimnames=corr0@Dimnames, uplo=corr0@uplo, factors=list())
  as_SFBM(corr_fixed, backingfile=backingfile, compact=compact)
}

# ── [1] Cargar sumstats pruneados ──────────────────────────────────
cat("\n[1] Cargando sumstats pruneados...\n")
sumstats_file <- file.path(out_dir, "sumstats_pruned.rds")
if (!file.exists(sumstats_file)) {
  stop("No existe sumstats_pruned.rds — primero correr prune")
}
sumstats_all <- readRDS(sumstats_file)
sumstats_chr <- sumstats_all[sumstats_all\$chr == chr, ]
rm(sumstats_all); gc()

cat(sprintf("  sumstats chr%s: %s\n", chr_label, format(nrow(sumstats_chr), big.mark=",")))

if (nrow(sumstats_chr) < 50L) {
  cat(sprintf("  Muy pocas variantes (%d<50), saltando\n", nrow(sumstats_chr)))
  saveRDS(list(chr=chr, n_matched=nrow(sumstats_chr), status="too_few_variants",
               beta_final=numeric(0), ind_col=integer(0),
               prs_chr=NULL, fam=NULL, info_snp=NULL,
               h2=NA_real_, p=NA_real_, n_valid_chains=0L,
               burn_in=NA_integer_, num_iter=NA_integer_, h2_init=NA_real_,
               t_cor_min=0, t_ldpred_min=0), betas_file)
  quit(save="no", status=0)
}

# ── [2] Leer .bed y matching ───────────────────────────────────────
prefix <- file.path(plink_dir, sprintf("mcps-freeze150k_qcd_chr%s_ivs", chr_label))
bed_file <- paste0(prefix, ".bed")
if (!file.exists(bed_file)) stop("No existe .bed: ", bed_file)

cat("\n[2] Leyendo .bed y matching...\n")
bk_file <- file.path(tmp_dir, sprintf("chr%s", chr_label))
snp_readBed(bed_file, backingfile=bk_file)
obj <- snp_attach(paste0(bk_file, ".rds"))
G   <- obj\$genotypes

map_df <- data.frame(
  chr=obj\$map\$chromosome, pos=obj\$map\$physical.pos,
  a1=obj\$map\$allele1, a0=obj\$map\$allele2,
  stringsAsFactors=FALSE)

if (chr == 23L) sumstats_chr\$chr <- "X"

info_snp <- snp_match(sumstats_chr, map_df, strand_flip=FALSE,
                      join_by_pos=TRUE, match.min.prop=0.01)
ind_col <- info_snp\$\`_NUM_ID_\`
cat(sprintf("  Matched: %s\n", format(nrow(info_snp), big.mark=",")))

sd_vals <- big_scale()(G, ind.col=ind_col)\$scale
keep_var <- which(sd_vals > 1e-6)
if (length(keep_var) < nrow(info_snp)) {
  cat(sprintf("  Removidas %d variantes con sd≈0\n", nrow(info_snp) - length(keep_var)))
  info_snp <- info_snp[keep_var, ]
  ind_col  <- info_snp\$\`_NUM_ID_\`
}
n_snp <- nrow(info_snp)
cat(sprintf("  Final: %s variantes\n", format(n_snp, big.mark=",")))

if (n_snp < 50L) {
  cat("  Muy pocas variantes post-filtro\n")
  saveRDS(list(chr=chr, n_matched=n_snp, status="too_few_variants",
               beta_final=numeric(0), ind_col=integer(0),
               prs_chr=NULL, fam=obj\$fam, info_snp=NULL,
               h2=NA_real_, p=NA_real_, n_valid_chains=0L,
               burn_in=NA_integer_, num_iter=NA_integer_, h2_init=NA_real_,
               t_cor_min=0, t_ldpred_min=0), betas_file)
  quit(save="no", status=0)
}

# ── [3] snp_cor ────────────────────────────────────────────────────
cat("\n[3] snp_cor (size=3000, thr_r2=0.005)...\n")
t0 <- Sys.time()
corr0 <- snp_cor(Gna=G, ind.col=ind_col, infos.pos=info_snp\$pos/1000,
                 size=3000, ncores=ncores, thr_r2=0.005)
t_cor <- as.numeric(difftime(Sys.time(), t0, units="mins"))
cat(sprintf("  snp_cor: %.1f min, dim=%dx%d, nnz=%s\n",
            t_cor, nrow(corr0), ncol(corr0),
            format(length(corr0@x), big.mark=",")))

# ── [4] SFBM fresco ───────────────────────────────────────────────
cat("\n[4] Creando SFBM...\n")
sfbm_file <- file.path(tmp_dir, sprintf("corr_chr%s", chr_label))
corr <- safe_as_SFBM(corr0, backingfile=sfbm_file, compact=TRUE)
cat(sprintf("  SFBM: %.1f MB\n",
            file.info(paste0(sfbm_file, ".sbk"))\$size / 1024^2))
rm(corr0); gc()

# ── [5] LDpred2-auto ───────────────────────────────────────────────
is_hard <- chr %in% c(18L, 19L, 21L, 22L)

h2_init_use <- 0.01
burn_in <- 1000L; num_iter <- 1000L
n_chains <- 30L

if (is_hard) {
  burn_in <- 2000L; num_iter <- 2000L
  n_chains <- 50L
  h2_lower <- 1e-5
  cat("\n[5] LDpred2-auto HARD-MODE\n")
} else {
  h2_lower <- 0.0005
  cat("\n[5] LDpred2-auto\n")
}
cat(sprintf("    burn_in=%d, num_iter=%d, h2_init=%.3f, chains=%d, h2_lower=%.1e\n",
            burn_in, num_iter, h2_init_use, n_chains, h2_lower))

# ── Calibrar n_eff per-SNP con sd_val (Privé 2022) ─────────────────
df_beta <- info_snp[, c("beta","beta_se","n_eff")]
sd_val <- sd_vals[keep_var]

sd_ss_orig <- with(df_beta, 2 / sqrt(n_eff * beta_se^2 + beta^2))
ratio_orig <- median(sd_ss_orig / sd_val)
cat(sprintf("  sd_ss/sd_val ratio ANTES: %.3f\n", ratio_orig))

df_beta\$n_eff <- (2 / (sd_val * df_beta\$beta_se))^2 - (df_beta\$beta / df_beta\$beta_se)^2
df_beta\$n_eff <- pmax(df_beta\$n_eff, 1000)

sd_ss_new <- with(df_beta, 2 / sqrt(n_eff * beta_se^2 + beta^2))
ratio_new <- median(sd_ss_new / sd_val)
cat(sprintf("  sd_ss/sd_val ratio DESPUÉS: %.3f\n", ratio_new))
cat(sprintf("  n_eff per-SNP: median=%.0f, range=[%.0f, %.0f]\n",
            median(df_beta\$n_eff), min(df_beta\$n_eff), max(df_beta\$n_eff)))

# ── Run ───────────────────────────────────────────────────────────
t0 <- Sys.time()
multi_auto <- snp_ldpred2_auto(
  corr=corr, df_beta=df_beta,
  h2_init=h2_init_use,
  vec_p_init=seq_log(1e-4, 0.5, length.out=n_chains),
  burn_in=burn_in, num_iter=num_iter, report_step=100L,
  ncores=ncores, sparse=TRUE,
  allow_jump_sign=FALSE)
t_ldp <- as.numeric(difftime(Sys.time(), t0, units="mins"))
cat(sprintf("  LDpred2-auto: %.1f min\n", t_ldp))

h2_est <- sapply(multi_auto, function(x) x\$h2_est)
p_est  <- sapply(multi_auto, function(x) x\$p_est)
keep <- which(!is.na(h2_est) & h2_est > h2_lower & h2_est < 0.5 &
              !is.na(p_est)  & p_est  > 1e-6   & p_est  < 0.99)
cat(sprintf("  Cadenas válidas: %d/%d\n", length(keep), length(multi_auto)))
cat("  h2_est:\n")
cat(sprintf("    %s\n", paste(round(h2_est, 6), collapse=", ")))

# RETRY
if (is_hard && length(keep) < 1L) {
  cat(sprintf("\n  [RETRY] chr%s\n", chr_label))
  multi_auto <- snp_ldpred2_auto(
    corr=corr, df_beta=df_beta,
    h2_init=0.001,
    vec_p_init=seq_log(1e-5, 0.3, length.out=50),
    burn_in=3000L, num_iter=3000L, report_step=100L,
    ncores=ncores, sparse=TRUE,
    allow_jump_sign=FALSE)
  h2_est <- sapply(multi_auto, function(x) x\$h2_est)
  p_est  <- sapply(multi_auto, function(x) x\$p_est)
  keep <- which(!is.na(h2_est) & h2_est > h2_lower & h2_est < 0.5 &
                !is.na(p_est)  & p_est  > 1e-6   & p_est  < 0.99)
  cat(sprintf("  RETRY chains: %d/%d\n", length(keep), length(multi_auto)))
  cat("  h2_est (retry):\n")
  cat(sprintf("    %s\n", paste(round(h2_est, 6), collapse=", ")))
}

if (length(keep) < 1L) {
  cat(sprintf("\n  chr%s divergió — beta=0\n", chr_label))
  beta_final <- rep(0.0, n_snp)
  prs_chr    <- rep(0.0, nrow(G))
  status     <- "diverged_zeroed"
  h2_med     <- NA_real_; p_med <- NA_real_
} else {
  beta_mat   <- sapply(multi_auto[keep], function(x) x\$beta_est)
  if (is.null(dim(beta_mat))) beta_mat <- matrix(beta_mat, ncol=1)
  beta_final <- rowMeans(beta_mat, na.rm=TRUE)
  beta_final[!is.finite(beta_final)] <- 0.0
  h2_med <- median(h2_est[keep])
  p_med  <- median(p_est[keep])
  cat(sprintf("  h2=%.6f, p=%.2e, nonzero=%s, valid=%d\n",
              h2_med, p_med, format(sum(beta_final != 0), big.mark=","), length(keep)))
  prs_chr <- big_prodVec(X=G, y.col=beta_final,
                         ind.col=info_snp\$\`_NUM_ID_\`, ncores=ncores)
  status <- "ok"
}

# ── [6] Guardar ────────────────────────────────────────────────────
saveRDS(list(
  chr=chr, chr_label=chr_label, n_matched=n_snp,
  h2=h2_med, p=p_med,
  n_valid_chains=length(keep), status=status,
  burn_in=burn_in, num_iter=num_iter, h2_init=h2_init_use,
  h2_lower=h2_lower,
  neff_persnp_median=median(df_beta\$n_eff),
  sd_ss_sd_val_ratio=ratio_new,
  t_cor_min=t_cor, t_ldpred_min=t_ldp,
  beta_final=beta_final,
  ind_col=info_snp\$\`_NUM_ID_\`,
  prs_chr=prs_chr, fam=obj\$fam,
  info_snp=info_snp[, c("chr","pos","a0","a1","_NUM_ID_")]
), betas_file)

cat(sprintf("\n ✓ chr%s (%s) — %.1f min total\n",
            chr_label, status, t_cor + t_ldp))
REOF

apptainer exec \
  --bind /mnt/cephfs:/mnt/cephfs \
  --bind "$TMPDIR:/tmp,$TMPDIR:/var/tmp" \
  --env OMP_NUM_THREADS=1 \
  --env OPENBLAS_NUM_THREADS=1 \
  --env MKL_NUM_THREADS=1 \
  "$SIF" /opt/miniforge3/envs/r-bio/bin/Rscript "$R_SCRIPT"

rm -rf "$TMPDIR"

echo ""
echo "============================================================"
echo "  worker chr${CHR_LABEL} FIN — $(date)"
echo "============================================================"
