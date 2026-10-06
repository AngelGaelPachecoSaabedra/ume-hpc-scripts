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
# PASO 2/2: Ensamblar SFBM genome-wide + LDpred2-auto + PRS
# Requiere: Step 1 completado (al menos algunos chr*_blocks.rds)
# v2 (2026-05-07): tolera chr fallidos, usa safe_as_SFBM
# ================================================================

SIF="/mnt/cephfs/biocontainers/images/jupyter-biotools-1.5.sif"
PLINK_DIR="/mnt/cephfs/hot_nvme/mcps/imputed-topmed/plink_files/maximally_unrelated"
OUT_DIR="/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes/ldpred2_gwide"
BLOCKS_DIR="${OUT_DIR}/blocks_resplit"

TMPDIR="/mnt/cephfs/m.2/scratch/$USER/ldp2_gw_run_tmp/$SLURM_JOB_ID"
mkdir -p "$TMPDIR"; chmod 700 "$TMPDIR"
export TMPDIR TMP="$TMPDIR" TEMP="$TMPDIR"

echo "============================================================"
echo "  LDpred2 gwide STEP 2 — Ensamblar + LDpred2-auto + PRS"
echo "  Job: $SLURM_JOB_ID   Nodo: $SLURMD_NODENAME"
echo "  $(date)"
echo "  Script version: v2 (2026-05-07)"
echo "============================================================"

# ──────────────────────────────────────────────────────────────────
# FIX #5: Solo verificar que existan archivos. El R script decide
# qué chr usar según status, así no se aborta por 1 chr fallido.
# ──────────────────────────────────────────────────────────────────
PRESENT=0
MISSING=0
FAILED=0
for CHR in $(seq 1 22); do
  if [[ -f "${BLOCKS_DIR}/chr${CHR}_blocks.rds" ]]; then
    PRESENT=$((PRESENT + 1))
  else
    echo "  FALTA: chr${CHR}_blocks.rds"
    MISSING=$((MISSING + 1))
  fi
done

if [[ $PRESENT -lt 18 ]]; then
  echo "ERROR: solo $PRESENT/22 chr presentes — Step 1 incompleto"
  exit 1
fi
echo "  $PRESENT/22 chr*_blocks.rds presentes"

R_SCRIPT="$TMPDIR/gwide_run.R"

cat << 'REOF' > "$R_SCRIPT"
Sys.setenv(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
           MKL_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1",
           BLIS_NUM_THREADS="1", NUMEXPR_NUM_THREADS="1")
if (requireNamespace("RhpcBLASctl", quietly=TRUE)) {
  RhpcBLASctl::blas_set_num_threads(1L)
  RhpcBLASctl::omp_set_num_threads(1L)
}
suppressMessages({
  library(bigsnpr); library(bigsparser); library(data.table); library(Matrix)
})
options(bigstatsr.check.parallel.blas=FALSE)
options(bigstatsr.ncores.max=48L)

bv <- packageVersion("bigsnpr")
cat(sprintf("  bigsnpr version: %s\n", as.character(bv)))
if (as.character(bv) == "1.12.21") {
  stop("bigsnpr 1.12.21 está roto para snp_ldpred2_auto. Usa 1.10.8.")
}

ncores     <- 24L
plink_dir  <- Sys.getenv("PLINK_DIR")
out_dir    <- Sys.getenv("OUT_DIR")
blocks_dir <- Sys.getenv("BLOCKS_DIR")
tmp_dir    <- Sys.getenv("TMPDIR")

betas_file <- file.path(out_dir, "ldpred2_gwide_results.rds")
if (file.exists(betas_file)) {
  cat("Ya existe resultado — saltando.\n")
  quit(save="no", status=0)
}

# ══════════════════════════════════════════════════════════════════
# FIX #3: safe_as_SFBM — workaround para columnas vacías en dsCMatrix.
# Inserta entradas microscópicas (1e-10) en la diagonal de columnas
# vacías para evitar que as_SFBM segfaulte.
# ══════════════════════════════════════════════════════════════════
safe_as_SFBM <- function(corr_mat, backingfile, compact=TRUE) {
  # Forzar dsCMatrix (column-compressed sparse symmetric)
  if (!is(corr_mat, "dsCMatrix")) {
    corr_mat <- as(corr_mat, "symmetricMatrix")
  }
  # Detectar columnas vacías mirando el slot @p (CSC pointers)
  col_nnz <- diff(corr_mat@p)
  empty_cols <- which(col_nnz == 0)

  if (length(empty_cols) > 0) {
    cat(sprintf("  safe_as_SFBM: %d columnas vacías detectadas — parchando\n",
                length(empty_cols)))

    # Manipular slots directamente (NO usar Matrix arithmetic, que
    # convertiría dsCMatrix → dgCMatrix y rompería simetría)
    n <- ncol(corr_mat)
    new_x_added  <- length(empty_cols)
    old_nnz      <- length(corr_mat@x)

    # Nuevos vectores @i, @x, @p
    new_i <- integer(old_nnz + new_x_added)
    new_x <- numeric(old_nnz + new_x_added)
    new_p <- integer(n + 1L)

    pos <- 0L
    for (j in seq_len(n)) {
      old_start <- corr_mat@p[j] + 1L  # 1-based
      old_end   <- corr_mat@p[j+1L]    # exclusive
      n_in_col  <- old_end - corr_mat@p[j]

      if (n_in_col > 0L) {
        idx_dst <- (pos+1L):(pos+n_in_col)
        new_i[idx_dst] <- corr_mat@i[old_start:old_end]
        new_x[idx_dst] <- corr_mat@x[old_start:old_end]
        pos <- pos + n_in_col
      } else {
        # Columna vacía: insertar (j, j, 1e-10) - diagonal microscópica
        pos <- pos + 1L
        new_i[pos] <- j - 1L  # 0-based en slot @i
        new_x[pos] <- 1e-10
      }
      new_p[j+1L] <- pos
    }

    # Reconstruir el dsCMatrix
    corr_mat@i <- new_i[seq_len(pos)]
    corr_mat@x <- new_x[seq_len(pos)]
    corr_mat@p <- new_p
  }

  as_SFBM(corr_mat, backingfile=backingfile, compact=compact)
}

# ══════════════════════════════════════════════════════════════════
# [1] Cargar bloques pre-calculados (saltar chr fallidos)
# ══════════════════════════════════════════════════════════════════
cat("\n[1] Cargando bloques pre-calculados...\n")

all_blocks   <- list()
all_info_snp <- list()
total_snps   <- 0L
total_blocks <- 0L
chr_ok       <- integer(0)
chr_failed   <- integer(0)

for (chr in 1:22) {
  f <- file.path(blocks_dir, sprintf("chr%d_blocks.rds", chr))
  if (!file.exists(f)) {
    cat(sprintf("  chr%d: archivo ausente — saltando\n", chr))
    chr_failed <- c(chr_failed, chr)
    next
  }
  d <- readRDS(f)

  if (is.null(d$status) || d$status != "ok") {
    cat(sprintf("  chr%d: status=%s — saltando\n", chr,
                ifelse(is.null(d$status), "NULL", d$status)))
    chr_failed <- c(chr_failed, chr)
    next
  }

  all_blocks <- c(all_blocks, d$blocks)
  all_info_snp[[length(all_info_snp) + 1]] <- as.data.frame(d$info_snp)
  total_snps   <- total_snps + d$n_snp
  total_blocks <- total_blocks + d$n_blocks
  chr_ok <- c(chr_ok, chr)
  cat(sprintf("  chr%d: %d bloques, %s SNPs\n", chr, d$n_blocks,
      format(d$n_snp, big.mark=",")))
}

if (length(chr_failed) > 0) {
  cat(sprintf("\n  ADVERTENCIA: %d chr fallidos (%s) — PRS será incompleto\n",
              length(chr_failed), paste(chr_failed, collapse=",")))
}
if (length(chr_ok) < 18L) {
  stop(sprintf("Solo %d chr OK — abortando. Reanudar Step 1 para chr fallidos.",
               length(chr_ok)))
}

info_snp_gw <- do.call(rbind, all_info_snp)
cat(sprintf("\n  Total: %s SNPs en %d bloques (%d chr OK)\n",
    format(total_snps, big.mark=","), total_blocks, length(chr_ok)))

# ══════════════════════════════════════════════════════════════════
# Diagnóstico: nnz total, distribución de thr_r2 y bloques
# ══════════════════════════════════════════════════════════════════
nnz_total <- sum(sapply(all_blocks, function(b) length(b@x)))
cat(sprintf("\n  nnz total: %s (~%.1f GB en disco SFBM)\n",
    format(nnz_total, big.mark=","), nnz_total * 12 / 1024^3))

cat("\n  Resumen por cromosoma:\n")
cat("  chr | thr_r2  | bloques | SNPs       | nnz\n")
cat("  ----|---------|---------|------------|------------\n")
for (i in seq_along(chr_ok)) {
  chr <- chr_ok[i]
  f <- file.path(blocks_dir, sprintf("chr%d_blocks.rds", chr))
  d <- readRDS(f)
  nnz_chr <- sum(sapply(d$blocks, function(b) length(b@x)))
  thr_str <- if (is.na(d$thr_r2_split)) "fallback" else sprintf("%.3f", d$thr_r2_split)
  cat(sprintf("  %3d | %-7s | %7d | %10s | %s\n",
    chr, thr_str, d$n_blocks,
    format(d$n_snp, big.mark=","), format(nnz_chr, big.mark=",")))
}

# Sanidad: verificar que nnz total < 2^31 (límite int32 de dsCMatrix)
if (nnz_total > 2^31 - 1) {
  stop(sprintf("nnz_total=%s excede 2^31-1 — overflow int32 en bdiag/SFBM. ",
       "Hay que reducir thr_r2 en snp_cor o filtrar más SNPs.",
       format(nnz_total, big.mark=",")))
}

# ══════════════════════════════════════════════════════════════════
# [2] Construir SFBM block-diagonal (con safe_as_SFBM)
# ══════════════════════════════════════════════════════════════════
cat("\n[2] Construyendo SFBM block-diagonal...\n")
t0 <- Sys.time()

corr_gw <- Matrix::bdiag(all_blocks)
corr_gw <- as(corr_gw, "symmetricMatrix")
rm(all_blocks); gc()

sfbm_file <- file.path(tmp_dir, "corr_gwide")
sfbm <- safe_as_SFBM(corr_gw, backingfile=sfbm_file, compact=TRUE)
t_sfbm <- as.numeric(difftime(Sys.time(), t0, units="mins"))
cat(sprintf("  SFBM: %.1f GB, %.1f min\n",
    file.info(paste0(sfbm_file, ".sbk"))$size / 1024^3, t_sfbm))
rm(corr_gw); gc()

# ══════════════════════════════════════════════════════════════════
# [3] LDpred2-auto genome-wide
# ══════════════════════════════════════════════════════════════════
cat("\n[3] LDpred2-auto GENOME-WIDE...\n")

# ══════════════════════════════════════════════════════════════════
# n_eff CORRECTO (doc oficial Privé): 4/(1/ncas + 1/ncon) por estrato.
# NO multiplicar por K(1-K) — eso NO es una corrección reconocida.
# El n_eff binario ya maneja la escala logística de REGENIE.
# La divergencia previa se debía a BLOQUES GIGANTES, no al n_eff;
# se resolvió con la re-partición (blocks_resplit/, ≤2,500 SNPs).
# La transformación a liability se hace DESPUÉS sobre h2 (Lee 2011).
# Consistente con Pipeline C (Custom) → comparables.
# ══════════════════════════════════════════════════════════════════
neff_f <- 4 * 17300 * 54196 / (17300 + 54196)   # 52,456
neff_m <- 4 *  8522 * 26733 / ( 8522 + 26733)   # 25,848
neff_ldpred2 <- neff_f + neff_m                  # 78,304

df_beta <- data.frame(
  beta    = info_snp_gw$beta,
  beta_se = info_snp_gw$beta_se,
  n_eff   = neff_ldpred2
)
cat(sprintf("    n_eff = %.0f (binario correcto, sin K(1-K))\n",
            neff_ldpred2))

n_chains <- 50L
burn_in  <- 1000L
num_iter <- 1000L

cat(sprintf("    %s variantes, %d bloques\n",
    format(nrow(df_beta), big.mark=","), total_blocks))
cat(sprintf("    burn_in=%d, num_iter=%d, chains=%d\n",
    burn_in, num_iter, n_chains))
cat("    allow_jump_sign=FALSE, shrink_corr=0.95\n")

t0 <- Sys.time()
multi_auto <- snp_ldpred2_auto(
  corr=sfbm, df_beta=df_beta,
  h2_init=0.1,
  vec_p_init=seq_log(1e-4, 0.5, length.out=n_chains),
  burn_in=burn_in, num_iter=num_iter, report_step=100L,
  ncores=ncores, sparse=TRUE,
  allow_jump_sign=FALSE,
  shrink_corr=0.95)
t_ldp <- as.numeric(difftime(Sys.time(), t0, units="mins"))
cat(sprintf("  LDpred2-auto: %.1f min\n", t_ldp))

h2_est <- sapply(multi_auto, function(x) x$h2_est)
p_est  <- sapply(multi_auto, function(x) x$p_est)
keep <- which(!is.na(h2_est) & h2_est > 0.001 & h2_est < 0.8 &
              !is.na(p_est)  & p_est  > 1e-6   & p_est  < 0.99)
cat(sprintf("  Cadenas válidas: %d/%d\n", length(keep), length(multi_auto)))

if (length(keep) > 0) {
  cat(sprintf("  h2 median: %.6f, p median: %.2e\n",
      median(h2_est[keep]), median(p_est[keep])))
}
cat("  h2_est:\n")
cat(sprintf("    %s\n", paste(round(h2_est, 6), collapse=", ")))

# RETRY
if (length(keep) < 3L) {
  cat("\n  [RETRY] shrink_corr=0.5, más iteraciones...\n")
  multi_auto <- snp_ldpred2_auto(
    corr=sfbm, df_beta=df_beta,
    h2_init=0.1,
    vec_p_init=seq_log(1e-5, 0.3, length.out=50),
    burn_in=2000L, num_iter=2000L, report_step=100L,
    ncores=ncores, sparse=TRUE,
    allow_jump_sign=FALSE,
    shrink_corr=0.5)
  h2_est <- sapply(multi_auto, function(x) x$h2_est)
  p_est  <- sapply(multi_auto, function(x) x$p_est)
  keep <- which(!is.na(h2_est) & h2_est > 0.001 & h2_est < 0.8 &
                !is.na(p_est)  & p_est  > 1e-6   & p_est  < 0.99)
  cat(sprintf("  RETRY cadenas: %d/%d\n", length(keep), length(multi_auto)))
}

# ══════════════════════════════════════════════════════════════════
# [4] Betas finales + PRS
# ══════════════════════════════════════════════════════════════════
cat("\n[4] Betas finales + PRS...\n")

if (length(keep) < 1L) {
  cat("  DIVERGENCIA TOTAL\n")
  saveRDS(list(status="diverged", h2_est=h2_est, p_est=p_est,
               n_snp=total_snps,
               chr_ok=chr_ok, chr_failed=chr_failed), betas_file)
} else {
  beta_mat <- sapply(multi_auto[keep], function(x) x$beta_est)
  if (is.null(dim(beta_mat))) beta_mat <- matrix(beta_mat, ncol=1)
  beta_final <- rowMeans(beta_mat, na.rm=TRUE)
  beta_final[!is.finite(beta_final)] <- 0.0

  h2_med <- median(h2_est[keep])
  p_med  <- median(p_est[keep])
  cat(sprintf("  h2=%.6f, p=%.2e, nonzero=%s\n",
      h2_med, p_med, format(sum(beta_final != 0), big.mark=",")))

  # PRS por cromosoma
  chr_list <- sort(unique(info_snp_gw$chr))
  prs_total <- NULL; fam <- NULL

  for (chr in chr_list) {
    idx_chr <- which(info_snp_gw$chr == chr)
    if (length(idx_chr) == 0) next

    prefix <- file.path(plink_dir, sprintf("mcps-freeze150k_qcd_chr%d_ivs", chr))
    bk_file <- file.path(tmp_dir, sprintf("prs_chr%d", chr))
    snp_readBed(paste0(prefix, ".bed"), backingfile=bk_file)
    obj <- snp_attach(paste0(bk_file, ".rds"))
    G <- obj$genotypes

    # Imputación defensiva también en cómputo de PRS
    n_na_g <- big_counts(G, ind.col=info_snp_gw$`_NUM_ID_`[idx_chr])[4, ] |> sum()
    if (n_na_g > 0) {
      G <- snp_fastImputeSimple(G, method="mean0", ncores=ncores)
    }

    prs_chr <- big_prodVec(X=G, y.col=beta_final[idx_chr],
      ind.col=info_snp_gw$`_NUM_ID_`[idx_chr], ncores=ncores)

    if (is.null(prs_total)) {
      prs_total <- prs_chr
      fam <- obj$fam
    } else {
      prs_total <- prs_total + prs_chr
    }
    cat(sprintf("    chr%d: %d SNPs, PRS added\n", chr, length(idx_chr)))
    rm(obj, G); gc()
  }

  saveRDS(list(
    status="ok",
    h2=h2_med, p=p_med,
    n_valid_chains=length(keep),
    n_snp_total=total_snps, n_blocks=total_blocks,
    h2_all=h2_est, p_all=p_est,
    neff_persnp_median=median(df_beta$n_eff),
    t_ldpred_min=t_ldp, t_sfbm_min=t_sfbm,
    beta_final=beta_final,
    info_snp=info_snp_gw[, c("chr","pos","a0","a1")],
    prs_total=prs_total, fam=fam,
    chr_ok=chr_ok, chr_failed=chr_failed,
    bigsnpr_version=as.character(bv)
  ), betas_file)

  cat(sprintf("\n ✓ LDpred2 genome-wide completado\n"))
  cat(sprintf("   h2=%.6f, p=%.2e, chains=%d/%d\n",
      h2_med, p_med, length(keep), length(multi_auto)))
  cat(sprintf("   PRS: %s individuos\n", format(length(prs_total), big.mark=",")))
  if (length(chr_failed) > 0) {
    cat(sprintf("   AVISO: PRS computado SIN chr %s\n",
                paste(chr_failed, collapse=",")))
  }
}
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
echo "  STEP 2 FIN — $(date)"
echo "============================================================"
