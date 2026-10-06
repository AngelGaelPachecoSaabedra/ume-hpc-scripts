#!/bin/bash
#SBATCH --job-name=ldpred2_chr
#SBATCH --account=researchers
#SBATCH --qos=vip
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=150G
#SBATCH --output=/mnt/cephfs/orgs/home/angel.pacheco/logs/ldpred2_worker_%A_%a.out
#SBATCH --error=/mnt/cephfs/orgs/home/angel.pacheco/logs/ldpred2_worker_%A_%a.err
#SBATCH --mail-user=angelpachecosaavedra@gmail.com
#SBATCH --mail-type=FAIL

set -euo pipefail

# ================================================================
# Worker v4 — bigsnpr 1.10.8, SFBM fresco, filtro QC SD_ss vs SD_val
#
# Fixes acumulados:
#   1. N_eff corregido (sin doble K*(1-K))
#   2. n_eff per-SNP calibrado con sd_val (Privé 2022)
#   3. NUEVO: Filtro estricto SD_ss vs SD_val ANTES de snp_cor
#   4. h2_init=0.01 unificado, vec_p amplio
#   5. Filtro h2 relajado (>1e-5) para hard-chrs
#   6. 50 cadenas para hard-chrs, acepta >=1 válida
#   7. allow_jump_sign=FALSE
#   8. SFBM fresco (filtro QC cambia dimensiones → no field-swap)
# ================================================================

CHR=$SLURM_ARRAY_TASK_ID

SIF="/mnt/cephfs/biocontainers/images/jupyter-biotools-1.4.sif"
GWAS_DIR="/mnt/cephfs/orgs/home/angel.pacheco/gwas_miami"
OUT_DIR="/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes/ldpred2"
PLINK_DIR="/mnt/cephfs/hot_nvme/mcps/imputed-topmed/plink_files/maximally_unrelated"

mkdir -p "$OUT_DIR"

TMPDIR="/mnt/cephfs/scratch/$USER/ldpred2_worker_tmp/${SLURM_ARRAY_JOB_ID}_${CHR}"
mkdir -p "$TMPDIR"; chmod 700 "$TMPDIR"
export TMPDIR TMP="$TMPDIR" TEMP="$TMPDIR"

echo "============================================================"
echo "  LDpred2 worker v4 (QC+fresh SFBM, v1.10.8) — chr${CHR}"
echo "  Job: ${SLURM_ARRAY_JOB_ID}_${CHR}   Nodo: ${SLURMD_NODENAME}"
echo "  $(date)"
echo "============================================================"

R_SCRIPT="$TMPDIR/worker_chr${CHR}.R"

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
chr_label <- as.character(chr)
ncores    <- 12L
plink_dir <- "$PLINK_DIR"
out_dir   <- "$OUT_DIR"
gwas_dir  <- "$GWAS_DIR"
tmp_dir   <- Sys.getenv("TMPDIR")

cat(sprintf("\n==== chr%s (ncores=%d) ====\n", chr_label, ncores))

betas_file <- file.path(out_dir, sprintf("ldpred2_betas_chr%s.rds", chr_label))
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
  corr_fixed <- new("dsCMatrix", i=new_i, p=new_p, x=new_x,
    Dim=corr0@Dim, Dimnames=corr0@Dimnames, uplo=corr0@uplo, factors=list())
  as_SFBM(corr_fixed, backingfile=backingfile, compact=compact)
}

# ── [0] HapMap3 ────────────────────────────────────────────────────
map_hm3 <- readRDS(file.path(out_dir, "map_hm3_ldpred2.rds"))
if ("pos_hg38" %in% names(map_hm3)) map_hm3\$pos_use <- map_hm3\$pos_hg38

# ── [1] GWAS meta-análisis IVW ─────────────────────────────────────
cat(sprintf("\n[1] GWAS meta chr%s...\n", chr_label))
gf <- fread(file.path(gwas_dir, "combined-gwas-results-all-ctrl-females.txt.gz"),
  select=c("CHROM","GENPOS","ID","ALLELE0","ALLELE1","A1FREQ","INFO","N","BETA","SE","P"),
  showProgress=FALSE)
gm <- fread(file.path(gwas_dir, "combined-gwas-results-all-ctrl-males.txt.gz"),
  select=c("CHROM","GENPOS","ID","ALLELE0","ALLELE1","A1FREQ","INFO","N","BETA","SE","P"),
  showProgress=FALSE)
gf <- gf[gf\$CHROM == chr, ]; gm <- gm[gm\$CHROM == chr, ]
gwas <- merge(gf, gm, by="ID", suffixes=c(".f",".m"))
rm(gf, gm); gc()

w_f <- 1/gwas\$SE.f^2; w_m <- 1/gwas\$SE.m^2
neff_f <- 4*17300*54196/(17300+54196)
neff_m <- 4*8522*26733/(8522+26733)
neff_ldpred2 <- neff_f + neff_m

cat(sprintf("  N_eff (IVW sum): %.0f\n", neff_ldpred2))

sumstats <- data.frame(
  chr=gwas\$CHROM.f, pos=gwas\$GENPOS.f,
  a1=gwas\$ALLELE1.f, a0=gwas\$ALLELE0.f,
  freq=(gwas\$A1FREQ.f*gwas\$N.f + gwas\$A1FREQ.m*gwas\$N.m)/(gwas\$N.f+gwas\$N.m),
  info=pmin(gwas\$INFO.f, gwas\$INFO.m),
  beta=(gwas\$BETA.f*w_f + gwas\$BETA.m*w_m)/(w_f+w_m),
  beta_se=sqrt(1/(w_f+w_m)), n_eff=neff_ldpred2, stringsAsFactors=FALSE)
sumstats\$p <- 2*pnorm(-abs(sumstats\$beta/sumstats\$beta_se))
sumstats <- sumstats[sumstats\$info>=0.8 & sumstats\$freq>=0.01 & sumstats\$freq<=0.99 &
                     !is.na(sumstats\$beta) & !is.na(sumstats\$beta_se) & sumstats\$beta_se>0, ]
sumstats\$chr <- as.integer(sumstats\$chr)

key_hm3  <- paste(map_hm3\$chr, map_hm3\$pos_use, sep=":")
key_gwas <- paste(sumstats\$chr, sumstats\$pos, sep=":")
sumstats <- sumstats[key_gwas %in% key_hm3, ]
rm(gwas); gc()
cat(sprintf("  GWAS: %s\n", format(nrow(sumstats), big.mark=",")))


# ── [2] Leer .bed, Matching + Calibración n_eff + Filtro QC ────────
cat("\n[2] Leyendo .bed, Matching + QC...\n")

prefix <- file.path(plink_dir, sprintf("mcps-freeze150k_qcd_chr%s_ivs", chr_label))
bed_file <- paste0(prefix, ".bed")
if (!file.exists(bed_file)) stop("No existe .bed: ", bed_file)

# Leer directamente al TMPDIR del job actual
bk_file <- file.path(tmp_dir, sprintf("chr%s", chr_label))
snp_readBed(bed_file, backingfile=bk_file)
obj <- snp_attach(paste0(bk_file, ".rds"))
G   <- obj\$genotypes

map_df <- data.frame(chr=obj\$map\$chromosome, pos=obj\$map\$physical.pos,
                     a1=obj\$map\$allele1, a0=obj\$map\$allele2, stringsAsFactors=FALSE)

info_snp <- snp_match(sumstats, map_df, strand_flip=FALSE,
                      join_by_pos=TRUE, match.min.prop=0.01)
ind_col <- info_snp\$\`_NUM_ID_\`

# 2a. Remover variantes con sd ≈ 0
sd_vals <- big_scale()(G, ind.col=ind_col)\$scale
keep_var <- which(sd_vals > 1e-6)
if (length(keep_var) < nrow(info_snp)) {
  cat(sprintf("  Removidas %d variantes con sd≈0\n", nrow(info_snp) - length(keep_var)))
  info_snp <- info_snp[keep_var, ]
  ind_col  <- info_snp\$\`_NUM_ID_\`
  sd_vals  <- sd_vals[keep_var]
}

# ====================================================================
# CORRECCIÓN: Filtro QC PRIMERO, Calibración DESPUÉS
# ====================================================================
df_beta_tmp <- info_snp[, c("beta","beta_se","n_eff")]

# 2b. Calcular sd_ss ORIGINAL (usando N_eff base de 78,304)
sd_ss_orig <- with(df_beta_tmp, 2 / sqrt(n_eff * beta_se^2 + beta^2))
ratio_qc <- median(sd_ss_orig / sd_vals)
cat(sprintf("  sd_ss_orig/sd_val ratio: %.3f\n", ratio_qc))

# 2c. Filtro estricto de Privé: eliminar variantes irreconciliables
keep_qc <- which(sd_ss_orig > 0.5 * sd_vals &
                 sd_ss_orig < sd_vals + 0.1 &
                 sd_ss_orig > 0.1 &
                 sd_vals > 0.05)

n_removed_qc <- length(sd_ss_orig) - length(keep_qc)
if (n_removed_qc > 0) {
  cat(sprintf("  Filtro QC (SD_ss vs SD_val): removidas %d variantes (%.1f%%)\n",
              n_removed_qc, 100 * n_removed_qc / length(sd_ss_orig)))
  info_snp    <- info_snp[keep_qc, ]
  ind_col     <- info_snp\$\`_NUM_ID_\`
  sd_vals     <- sd_vals[keep_qc]
  df_beta_tmp <- df_beta_tmp[keep_qc, ]
} else {
  cat("  Filtro QC: 0 variantes removidas\n")
}

# 2d. Calibrar n_eff per-SNP SOLO para las variantes que sobrevivieron
df_beta_tmp\$n_eff <- (2 / (sd_vals * df_beta_tmp\$beta_se))^2 - (df_beta_tmp\$beta / df_beta_tmp\$beta_se)^2
info_snp\$n_eff <- pmax(df_beta_tmp\$n_eff, 1000)

rm(df_beta_tmp, sd_ss_orig)

n_snp <- nrow(info_snp)
cat(sprintf("  Matched & QC: %s variantes\n", format(n_snp, big.mark=",")))
cat(sprintf("  n_eff per-SNP: median=%.0f, range=[%.0f, %.0f]\n",
            median(info_snp\$n_eff), min(info_snp\$n_eff), max(info_snp\$n_eff)))


# ── [3] snp_cor ────────────────────────────────────────────────────
cat("\n[3] snp_cor (size=3000, thr_r2=0.005)...\n")
t0 <- Sys.time()
corr0 <- snp_cor(Gna=G, ind.col=ind_col, infos.pos=info_snp\$pos/1000,
                 size=3000, ncores=ncores, thr_r2=0.005)
t_cor <- as.numeric(difftime(Sys.time(), t0, units="mins"))
cat(sprintf("  snp_cor: %.1f min, dim=%dx%d, nnz=%s\n",
            t_cor, nrow(corr0), ncol(corr0),
            format(length(corr0@x), big.mark=",")))

# ── [4] SFBM fresco (dimensiones post-filtro QC) ───────────────────
cat("\n[4] Creando SFBM fresco...\n")
sfbm_file <- file.path(tmp_dir, sprintf("corr_chr%s", chr_label))
col_nnz <- diff(corr0@p)
n_empty <- sum(col_nnz == 0L)
if (n_empty > 0L) {
  cat(sprintf("  [SFBM] %d columnas vacías — eps en diagonal\n", n_empty))
  empty_cols <- which(col_nnz == 0L)
  for (ec in empty_cols) corr0[ec, ec] <- 1e-10
}
sfbm <- as_SFBM(corr0, backingfile=sfbm_file, compact=TRUE)
cat(sprintf("  SFBM: %.1f MB\n", file.info(paste0(sfbm_file, ".sbk"))\$size / 1024^2))
rm(corr0); gc()

# ── [5] LDpred2-auto ───────────────────────────────────────────────
is_hard <- chr %in% c(18L, 19L, 21L, 22L)

# Parámetros unificados
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

# df_beta ya viene limpio con n_eff calibrado desde el paso 2
df_beta <- info_snp[, c("beta","beta_se","n_eff")]

# ── Run ───────────────────────────────────────────────────────────
t0 <- Sys.time()
multi_auto <- snp_ldpred2_auto(
  corr=sfbm, df_beta=df_beta,
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
    corr=sfbm, df_beta=df_beta,
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
  prs_chr <- rep(0.0, nrow(G))
  status <- "diverged_zeroed"
  h2_med <- NA_real_; p_med <- NA_real_
} else {
  beta_mat <- sapply(multi_auto[keep], function(x) x\$beta_est)
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
  chr=chr, n_matched=n_snp, h2=h2_med, p=p_med,
  n_valid_chains=length(keep), status=status,
  burn_in=burn_in, num_iter=num_iter, h2_init=h2_init_use,
  h2_lower=h2_lower,
  neff_ivw=neff_ldpred2, neff_persnp_median=median(df_beta\$n_eff),
  n_removed_qc=n_removed_qc, sd_ss_sd_val_ratio=ratio_qc,
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
  --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 --env MKL_NUM_THREADS=1 \
  "$SIF" /opt/miniforge3/envs/r-bio/bin/Rscript "$R_SCRIPT"

rm -rf "$TMPDIR"

echo ""
echo "============================================================"
echo "  worker chr${CHR} FIN — $(date)"
echo "============================================================"
