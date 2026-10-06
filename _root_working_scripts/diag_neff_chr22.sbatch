#!/bin/bash
#SBATCH --job-name=diag_neff
#SBATCH --account=researchers
#SBATCH --qos=vip
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=90G
#SBATCH --time=01:00:00
#SBATCH --output=/mnt/cephfs/orgs/home/angel.pacheco/logs/diag_neff_%j.out
#SBATCH --error=/mnt/cephfs/orgs/home/angel.pacheco/logs/diag_neff_%j.err

set -uo pipefail

# ================================================================
# DIAGNÓSTICO: ¿qué n_eff hace converger LDpred2-auto?
# Usa solo chr22 (el más pequeño) → ~15 min en vez de 2h.
# Prueba 4 valores de n_eff y reporta cuántas cadenas convergen.
# ================================================================

SIF="/mnt/cephfs/biocontainers/images/jupyter-biotools-1.5.sif"
OUT_DIR="/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes/ldpred2_gwide"
BLOCKS_DIR="${OUT_DIR}/blocks"

TMPDIR="/mnt/cephfs/m.2/scratch/$USER/diag_neff_$SLURM_JOB_ID"
mkdir -p "$TMPDIR"; chmod 700 "$TMPDIR"
export TMPDIR TMP="$TMPDIR" TEMP="$TMPDIR"

echo "============================================================"
echo "  DIAG n_eff — chr22 — Job: $SLURM_JOB_ID  $(date)"
echo "============================================================"

R_SCRIPT="$TMPDIR/diag_neff.R"
cat << 'REOF' > "$R_SCRIPT"
Sys.setenv(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
           MKL_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1",
           BLIS_NUM_THREADS="1", NUMEXPR_NUM_THREADS="1")
suppressMessages({ library(bigsnpr); library(bigsparser); library(Matrix) })
options(bigstatsr.check.parallel.blas=FALSE)

ncores     <- 12L
blocks_dir <- Sys.getenv("BLOCKS_DIR")
tmp_dir    <- Sys.getenv("TMPDIR")

cat(sprintf("  bigsnpr: %s\n", as.character(packageVersion("bigsnpr"))))

# ── Cargar chr22 ───────────────────────────────────────────────────
d <- readRDS(file.path(blocks_dir, "chr22_blocks.rds"))
cat(sprintf("\n  chr22: %d SNPs, %d bloques, status=%s\n",
            d$n_snp, d$n_blocks, d$status))

info_snp <- as.data.frame(d$info_snp)
cat(sprintf("  Columnas info_snp: %s\n", paste(names(info_snp), collapse=", ")))

# Construir SFBM de chr22
corr22 <- Matrix::bdiag(d$blocks)
corr22 <- as(corr22, "symmetricMatrix")
sfbm <- as_SFBM(corr22, backingfile=file.path(tmp_dir, "corr22"), compact=TRUE)
cat(sprintf("  SFBM chr22 construido: %d x %d\n", nrow(corr22), ncol(corr22)))

# ── Diagnóstico del n_eff guardado ─────────────────────────────────
cat("\n══════════ n_eff GUARDADO (per-SNP calibrado) ══════════\n")
cat(sprintf("  median=%.0f, mean=%.0f, min=%.0f, max=%.0f\n",
    median(info_snp$n_eff), mean(info_snp$n_eff),
    min(info_snp$n_eff), max(info_snp$n_eff)))
cat(sprintf("  quantiles: %s\n",
    paste(sprintf("%.0f", quantile(info_snp$n_eff, c(.05,.25,.5,.75,.95))),
          collapse=", ")))

# ── Candidatos de n_eff ────────────────────────────────────────────
K_prev      <- 0.2419
neff_binary <- 78304
neff_kk     <- neff_binary * K_prev * (1 - K_prev)   # ≈ 14,359
neff_median <- median(info_snp$n_eff)                 # ~57,000

candidatos <- list(
  "A_perSNP_guardado"   = info_snp$n_eff,             # lo que ya falló
  "B_KK_constante_14359"= rep(neff_kk, nrow(info_snp)),
  "C_median_constante"  = rep(neff_median, nrow(info_snp)),
  "D_binary_78304"      = rep(neff_binary, nrow(info_snp))
)

# ── Probar cada candidato ──────────────────────────────────────────
probar_neff <- function(nombre, neff_vec) {
  cat(sprintf("\n══════════ [%s] ══════════\n", nombre))
  df_beta <- data.frame(beta=info_snp$beta, beta_se=info_snp$beta_se,
                        n_eff=neff_vec)
  # ratio sd_ss/sd_geno NO lo tenemos aquí sin sd_val; reportamos sd_ss
  sd_ss <- 2 / sqrt(df_beta$n_eff * df_beta$beta_se^2 + df_beta$beta^2)
  cat(sprintf("  n_eff median=%.0f, sd_ss median=%.4f\n",
              median(neff_vec), median(sd_ss)))

  t0 <- Sys.time()
  ma <- tryCatch(
    snp_ldpred2_auto(corr=sfbm, df_beta=df_beta,
      h2_init=0.1, vec_p_init=seq_log(1e-4, 0.5, length.out=10),
      burn_in=300L, num_iter=300L, report_step=100L,
      ncores=ncores, sparse=TRUE,
      allow_jump_sign=FALSE, shrink_corr=0.95),
    error=function(e) { cat(sprintf("  ERROR: %s\n", conditionMessage(e))); NULL })
  dt <- as.numeric(difftime(Sys.time(), t0, units="mins"))

  if (is.null(ma)) return(invisible(NULL))

  h2 <- sapply(ma, function(x) x$h2_est)
  p  <- sapply(ma, function(x) x$p_est)
  keep <- which(!is.na(h2) & h2 > 0.001 & h2 < 0.8 &
                !is.na(p)  & p  > 1e-6   & p  < 0.99)
  cat(sprintf("  Cadenas válidas: %d/10 (%.1f min)\n", length(keep), dt))
  if (length(keep) > 0) {
    cat(sprintf("  h2 median=%.5f, p median=%.2e\n",
                median(h2[keep]), median(p[keep])))
  }
  cat(sprintf("  h2_est: %s\n", paste(round(h2,5), collapse=", ")))
}

for (nm in names(candidatos)) probar_neff(nm, candidatos[[nm]])

cat("\n══════════ DIAG FIN ══════════\n")
cat("  Interpretación:\n")
cat("  - El n_eff con MÁS cadenas válidas es el correcto para genome-wide.\n")
cat("  - Si B (K(1-K)=14359) converge y A (per-SNP) no → usar K(1-K) en Step 2.\n")
REOF

apptainer exec \
  --bind /mnt/cephfs:/mnt/cephfs \
  --bind "$TMPDIR:/tmp,$TMPDIR:/var/tmp" \
  --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 --env MKL_NUM_THREADS=1 \
  --env BLIS_NUM_THREADS=1 --env NUMEXPR_NUM_THREADS=1 \
  --env BLOCKS_DIR="$BLOCKS_DIR" \
  "$SIF" /opt/miniforge3/envs/r-bio/bin/Rscript "$R_SCRIPT"

rm -rf "$TMPDIR"

echo ""
echo "============================================================"
echo "  DIAG n_eff FIN — $(date)"
echo "============================================================"
