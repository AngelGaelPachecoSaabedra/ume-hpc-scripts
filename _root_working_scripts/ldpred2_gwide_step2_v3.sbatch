#!/bin/bash
#SBATCH --job-name=ldp2_gw_run
#SBATCH --account=researchers
#SBATCH --qos=vip
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=150G
#SBATCH --output=/mnt/cephfs/orgs/home/angel.pacheco/logs/ldp2_gw_run_%j.out
#SBATCH --error=/mnt/cephfs/orgs/home/angel.pacheco/logs/ldp2_gw_run_%j.err
#SBATCH --mail-user=angelpachecosaavedra@gmail.com
#SBATCH --mail-type=END,FAIL

set -euo pipefail

# ================================================================
# PASO 2/2 v3: LDpred2-auto POR CROMOSOMA + suma de PRS
#
# Cambio vs v2: el LDpred2-auto genome-wide (1 SFBM con 150 bloques,
# 1M SNPs) DIVERGE — bloques admixed de hasta 16,716 SNPs colapsan el
# Gibbs sampler. El diagnóstico chr22 convergió 10/10 procesando un
# cromosoma a la vez. Esta versión replica el enfoque de Pipeline B v8
# (que funcionó): LDpred2-auto por cromosoma, betas por cromosoma,
# y suma de PRS. Reutiliza los bloques pre-calculados de Step 1.
# ================================================================

SIF="/mnt/cephfs/biocontainers/images/jupyter-biotools-1.5.sif"
PLINK_DIR="/mnt/cephfs/hot_nvme/mcps/imputed-topmed/plink_files/maximally_unrelated"
OUT_DIR="/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes/ldpred2_gwide"
BLOCKS_DIR="${OUT_DIR}/blocks_resplit"

TMPDIR="/mnt/cephfs/m.2/scratch/$USER/ldp2_gw_run_tmp/$SLURM_JOB_ID"
mkdir -p "$TMPDIR"; chmod 700 "$TMPDIR"
export TMPDIR TMP="$TMPDIR" TEMP="$TMPDIR"

echo "============================================================"
echo "  LDpred2 gwide STEP 2 v3 — POR CROMOSOMA + PRS"
echo "  Job: $SLURM_JOB_ID   Nodo: $SLURMD_NODENAME"
echo "  $(date)"
echo "============================================================"

PRESENT=0
for CHR in $(seq 1 22); do
  [[ -f "${BLOCKS_DIR}/chr${CHR}_blocks.rds" ]] && PRESENT=$((PRESENT+1)) || echo "  FALTA: chr${CHR}"
done
if [[ $PRESENT -lt 18 ]]; then
  echo "ERROR: solo $PRESENT/22 chr presentes"; exit 1
fi
echo "  $PRESENT/22 chr*_blocks.rds presentes"

R_SCRIPT="$TMPDIR/gwide_run.R"

cat << 'REOF' > "$R_SCRIPT"
Sys.setenv(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
           MKL_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1",
           BLIS_NUM_THREADS="1", NUMEXPR_NUM_THREADS="1")
if (requireNamespace("RhpcBLASctl", quietly=TRUE)) {
  RhpcBLASctl::blas_set_num_threads(1L); RhpcBLASctl::omp_set_num_threads(1L)
}
suppressMessages({
  library(bigsnpr); library(bigsparser); library(data.table); library(Matrix)
})
options(bigstatsr.check.parallel.blas=FALSE)
options(bigstatsr.ncores.max=48L)

bv <- packageVersion("bigsnpr")
cat(sprintf("  bigsnpr version: %s\n", as.character(bv)))
if (as.character(bv) == "1.12.21") stop("bigsnpr 1.12.21 roto. Usa 1.10.8.")

ncores     <- 24L
plink_dir  <- Sys.getenv("PLINK_DIR")
out_dir    <- Sys.getenv("OUT_DIR")
blocks_dir <- Sys.getenv("BLOCKS_DIR")
tmp_dir    <- Sys.getenv("TMPDIR")

betas_file <- file.path(out_dir, "ldpred2_gwide_results.rds")
if (file.exists(betas_file)) { cat("Ya existe resultado — saltando.\n"); quit(save="no", status=0) }

# ── n_eff: corrección K(1-K) (igual que Pipeline B v8) ─────────────
K_prev       <- 0.2419
neff_binary  <- 78304
neff_ldpred2 <- neff_binary * K_prev * (1 - K_prev)   # ≈ 14,359
cat(sprintf("  n_eff = %.0f (K(1-K), constante)\n", neff_ldpred2))

# ══════════════════════════════════════════════════════════════════
# LDpred2-auto POR CROMOSOMA
# ══════════════════════════════════════════════════════════════════
cat("\n[1] LDpred2-auto por cromosoma...\n")

n_chains <- 30L
burn_in  <- 500L
num_iter <- 500L

prs_total    <- NULL
fam_global   <- NULL
chr_ok       <- integer(0)
chr_failed   <- integer(0)
results_log  <- list()
total_snps   <- 0L
total_nonzero <- 0L
betas_all    <- list()

for (chr in 1:22) {
  f <- file.path(blocks_dir, sprintf("chr%d_blocks.rds", chr))
  if (!file.exists(f)) { chr_failed <- c(chr_failed, chr); next }
  d <- readRDS(f)
  if (is.null(d$status) || d$status != "ok") {
    cat(sprintf("  chr%d: status no-ok — saltando\n", chr))
    chr_failed <- c(chr_failed, chr); next
  }

  t0 <- Sys.time()
  info_snp <- as.data.frame(d$info_snp)
  m_chr <- nrow(info_snp)

  # SFBM del cromosoma (bdiag de sus bloques)
  corr_chr <- Matrix::bdiag(d$blocks)
  corr_chr <- as(corr_chr, "symmetricMatrix")
  sfbm <- as_SFBM(corr_chr, backingfile=file.path(tmp_dir, sprintf("corr_chr%d", chr)),
                  compact=TRUE)
  rm(corr_chr); gc()

  # df_beta con K(1-K)
  df_beta <- data.frame(beta=info_snp$beta, beta_se=info_snp$beta_se,
                        n_eff=neff_ldpred2)

  # LDpred2-auto
  ma <- tryCatch(
    snp_ldpred2_auto(corr=sfbm, df_beta=df_beta,
      h2_init=0.1, vec_p_init=seq_log(1e-4, 0.5, length.out=n_chains),
      burn_in=burn_in, num_iter=num_iter, report_step=100L,
      ncores=ncores, sparse=TRUE,
      allow_jump_sign=FALSE, shrink_corr=0.95),
    error=function(e) { cat(sprintf("  chr%d ERROR: %s\n", chr, conditionMessage(e))); NULL })

  if (is.null(ma)) { chr_failed <- c(chr_failed, chr); next }

  h2 <- sapply(ma, function(x) x$h2_est)
  p  <- sapply(ma, function(x) x$p_est)
  keep <- which(!is.na(h2) & h2 > 0.0001 & h2 < 0.8 &
                !is.na(p)  & p  > 1e-6   & p  < 0.99)

  # Retry por cromosoma si pocas cadenas
  if (length(keep) < 3L) {
    cat(sprintf("  chr%d: %d/%d cadenas — RETRY shrink_corr=0.5\n", chr, length(keep), n_chains))
    ma <- snp_ldpred2_auto(corr=sfbm, df_beta=df_beta,
      h2_init=0.05, vec_p_init=seq_log(1e-5, 0.3, length.out=n_chains),
      burn_in=1000L, num_iter=1000L, report_step=100L,
      ncores=ncores, sparse=TRUE,
      allow_jump_sign=FALSE, shrink_corr=0.5)
    h2 <- sapply(ma, function(x) x$h2_est)
    p  <- sapply(ma, function(x) x$p_est)
    keep <- which(!is.na(h2) & h2 > 0.0001 & h2 < 0.8 &
                  !is.na(p)  & p  > 1e-6   & p  < 0.99)
  }

  dt <- as.numeric(difftime(Sys.time(), t0, units="mins"))

  if (length(keep) < 1L) {
    cat(sprintf("  chr%d: DIVERGIÓ (0/%d) — beta=0, %.1f min\n", chr, n_chains, dt))
    beta_final <- rep(0, m_chr)
    h2_med <- NA_real_; p_med <- NA_real_
    chr_failed <- c(chr_failed, chr)
  } else {
    beta_mat <- sapply(ma[keep], function(x) x$beta_est)
    if (is.null(dim(beta_mat))) beta_mat <- matrix(beta_mat, ncol=1)
    beta_final <- rowMeans(beta_mat, na.rm=TRUE)
    beta_final[!is.finite(beta_final)] <- 0.0
    h2_med <- median(h2[keep]); p_med <- median(p[keep])
    chr_ok <- c(chr_ok, chr)
    cat(sprintf("  chr%d: %d/%d cadenas, h2=%.5f, p=%.2e, nonzero=%s, %.1f min\n",
        chr, length(keep), length(ma), h2_med, p_med,
        format(sum(beta_final != 0), big.mark=","), dt))
  }

  total_snps    <- total_snps + m_chr
  total_nonzero <- total_nonzero + sum(beta_final != 0)
  betas_all[[as.character(chr)]] <- data.frame(
    chr=info_snp$chr, pos=info_snp$pos, a0=info_snp$a0, a1=info_snp$a1,
    beta=beta_final, num_id=info_snp$`_NUM_ID_`)
  results_log[[as.character(chr)]] <- list(
    n_snp=m_chr, n_blocks=d$n_blocks, h2=h2_med, p=p_med,
    n_valid=length(keep), elapsed_min=dt)

  # ── PRS del cromosoma ──
  prefix <- file.path(plink_dir, sprintf("mcps-freeze150k_qcd_chr%d_ivs", chr))
  bk_file <- file.path(tmp_dir, sprintf("prs_chr%d", chr))
  snp_readBed(paste0(prefix, ".bed"), backingfile=bk_file)
  obj <- snp_attach(paste0(bk_file, ".rds"))
  G <- obj$genotypes
  n_na_g <- big_counts(G, ind.col=info_snp$`_NUM_ID_`)[4, ] |> sum()
  if (n_na_g > 0) G <- snp_fastImputeSimple(G, method="mean0", ncores=ncores)

  prs_chr <- big_prodVec(X=G, y.col=beta_final,
    ind.col=info_snp$`_NUM_ID_`, ncores=ncores)

  if (is.null(prs_total)) { prs_total <- prs_chr; fam_global <- obj$fam }
  else                    { prs_total <- prs_total + prs_chr }
  rm(obj, G, sfbm, ma); gc()
}

cat(sprintf("\n  Resumen: %d chr OK, %d chr fallidos\n", length(chr_ok), length(chr_failed)))
if (length(chr_failed) > 0)
  cat(sprintf("  Chr fallidos: %s\n", paste(chr_failed, collapse=",")))

# h2 total aproximado (suma de h2 por cromosoma)
h2_per_chr <- sapply(results_log, function(r) r$h2)
h2_total_approx <- sum(h2_per_chr, na.rm=TRUE)
cat(sprintf("  h2 total aproximado (suma por chr): %.4f\n", h2_total_approx))

# ══════════════════════════════════════════════════════════════════
# [2] Exportar PRS con z-score y percentil por sexo (como v8)
# ══════════════════════════════════════════════════════════════════
cat("\n[2] Exportando PRS...\n")

betas_df <- do.call(rbind, betas_all)

covar_file <- "/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes/MCPS BASELINE.csv"
link_file  <- "/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes/RGN_LINK_IID.csv"

df_link <- read.csv(link_file)
iid_to_patid <- setNames(as.character(df_link$PATID), as.character(df_link$IID))
df_baseline <- read.csv(covar_file)
df_baseline <- df_baseline[!is.na(df_baseline$MALE), ]
patid_to_sex <- setNames(as.integer(df_baseline$MALE), as.character(df_baseline$PATID))

patids <- iid_to_patid[as.character(fam_global$sample.ID)]
patids[is.na(patids)] <- "UNKNOWN"
sexes <- patid_to_sex[patids]; sexes[is.na(sexes)] <- -1L

df_out <- data.frame(PATID=patids, IID=fam_global$sample.ID,
                     SEX=sexes, PRS_TOTAL=prs_total, stringsAsFactors=FALSE)
df_out <- df_out[df_out$SEX >= 0, ]

df_out$PRS_Z <- NA_real_; df_out$PERCENTILE <- NA_real_
for (sv in c(0,1)) {
  mask <- df_out$SEX == sv
  vals <- df_out$PRS_TOTAL[mask]
  if (length(vals) > 0 && sd(vals) > 0) {
    df_out$PRS_Z[mask] <- (vals - mean(vals)) / sd(vals)
    df_out$PERCENTILE[mask] <- rank(vals) / length(vals) * 100
  }
}
df_out <- df_out[order(-df_out$PRS_TOTAL), ]

prs_file   <- file.path(out_dir, "prs_diabetes_ldpred2_gwide_final.tsv")
write.table(df_out, prs_file, sep="\t", row.names=FALSE, quote=FALSE)

saveRDS(list(
  status="ok",
  approach="per_chromosome",
  n_eff=neff_ldpred2,
  h2_per_chr=h2_per_chr,
  h2_total_approx=h2_total_approx,
  n_snp_total=total_snps, n_nonzero=total_nonzero,
  chr_ok=chr_ok, chr_failed=chr_failed,
  results_log=results_log,
  betas=betas_df,
  prs_total=prs_total, fam=fam_global,
  bigsnpr_version=as.character(bv)
), betas_file)

cat("\n  Resumen por cromosoma:\n")
cat("  chr | h2       | p        | cadenas | SNPs    | min\n")
for (cn in names(results_log)) {
  r <- results_log[[cn]]
  cat(sprintf("  %3s | %.6f | %.2e | %3d/%d | %7d | %.1f\n",
    cn, ifelse(is.na(r$h2),0,r$h2), ifelse(is.na(r$p),0,r$p),
    r$n_valid, n_chains, r$n_snp, r$elapsed_min))
}

cat(sprintf("\n ✓ PRS completado: %s individuos, %s SNPs nonzero\n",
    format(length(prs_total), big.mark=","), format(total_nonzero, big.mark=",")))
cat(sprintf("   Output: %s\n", prs_file))
if (length(chr_failed) > 0)
  cat(sprintf("   AVISO: PRS sin chr %s\n", paste(chr_failed, collapse=",")))
REOF

apptainer exec \
  --bind /mnt/cephfs:/mnt/cephfs \
  --bind "$TMPDIR:/tmp,$TMPDIR:/var/tmp" \
  --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 --env MKL_NUM_THREADS=1 \
  --env BLIS_NUM_THREADS=1 --env NUMEXPR_NUM_THREADS=1 \
  --env PLINK_DIR="$PLINK_DIR" --env OUT_DIR="$OUT_DIR" --env BLOCKS_DIR="$BLOCKS_DIR" \
  "$SIF" /opt/miniforge3/envs/r-bio/bin/Rscript "$R_SCRIPT"

rm -rf "$TMPDIR"

echo ""
echo "============================================================"
echo "  STEP 2 v3 FIN — $(date)"
echo "============================================================"
