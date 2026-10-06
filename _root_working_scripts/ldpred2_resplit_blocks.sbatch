#!/bin/bash
#SBATCH --job-name=ldp2_resplit
#SBATCH --account=researchers
#SBATCH --qos=vip
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=90G
#SBATCH --output=/mnt/cephfs/orgs/home/angel.pacheco/logs/ldp2_resplit_%j.out
#SBATCH --error=/mnt/cephfs/orgs/home/angel.pacheco/logs/ldp2_resplit_%j.err

set -euo pipefail

# ================================================================
# RE-PARTICIÓN: sub-dividir bloques grandes (>MAX_BLOCK SNPs) en
# sub-bloques de ~TARGET SNPs. NO recalcula snp_cor (lo costoso);
# solo toma sub-matrices diagonales de los bloques ya guardados.
#
# Motivo: bloques de 16,496 SNPs (LD admixed extenso) hacen diverger
# el Gibbs sampler de LDpred2-auto. Bloques de ~2,500 convergen
# (confirmado con chr22).
#
# Genera blocks_resplit/chrN_blocks.rds (mismo formato que Step 1).
# ================================================================

SIF="/mnt/cephfs/biocontainers/images/jupyter-biotools-1.5.sif"
OUT_DIR="/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes/ldpred2_gwide"
BLOCKS_DIR="${OUT_DIR}/blocks"
RESPLIT_DIR="${OUT_DIR}/blocks_resplit"

mkdir -p "$RESPLIT_DIR"

TMPDIR="/mnt/cephfs/m.2/scratch/$USER/ldp2_resplit_$SLURM_JOB_ID"
mkdir -p "$TMPDIR"; chmod 700 "$TMPDIR"
export TMPDIR TMP="$TMPDIR" TEMP="$TMPDIR"

echo "============================================================"
echo "  RE-PARTICIÓN de bloques — Job: $SLURM_JOB_ID  $(date)"
echo "============================================================"

R_SCRIPT="$TMPDIR/resplit.R"
cat << 'REOF' > "$R_SCRIPT"
Sys.setenv(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
suppressMessages({ library(bigsnpr); library(Matrix) })

blocks_dir  <- Sys.getenv("BLOCKS_DIR")
resplit_dir <- Sys.getenv("RESPLIT_DIR")

# Parámetros de re-partición
MAX_BLOCK <- 4000L   # bloques mayores a esto se sub-dividen
TARGET    <- 2500L   # tamaño objetivo de sub-bloque

cat(sprintf("  MAX_BLOCK=%d, TARGET=%d\n\n", MAX_BLOCK, TARGET))

# Sub-dividir un bloque grande en piezas diagonales de ~TARGET
subdiv <- function(blk, target) {
  m <- nrow(blk)
  if (m <= target) return(list(blk))
  n_sub <- ceiling(m / target)
  sizes <- rep(floor(m / n_sub), n_sub)
  sizes[seq_len(m - sum(sizes))] <- sizes[seq_len(m - sum(sizes))] + 1L
  starts <- c(1L, cumsum(sizes[-n_sub]) + 1L)
  ends   <- cumsum(sizes)
  lapply(seq_len(n_sub), function(i) blk[starts[i]:ends[i], starts[i]:ends[i], drop=FALSE])
}

for (chr in 1:22) {
  f <- file.path(blocks_dir, sprintf("chr%d_blocks.rds", chr))
  if (!file.exists(f)) { cat(sprintf("  chr%d: ausente\n", chr)); next }
  d <- readRDS(f)
  if (is.null(d$status) || d$status != "ok") {
    cat(sprintf("  chr%d: status no-ok, copiando tal cual\n", chr))
    saveRDS(d, file.path(resplit_dir, sprintf("chr%d_blocks.rds", chr)))
    next
  }

  old_sizes <- sapply(d$blocks, nrow)
  new_blocks <- list()
  for (b in d$blocks) {
    subs <- subdiv(b, TARGET)
    new_blocks <- c(new_blocks, subs)
  }
  new_sizes <- sapply(new_blocks, nrow)

  # Verificación de integridad: mismo total de SNPs
  stopifnot(sum(new_sizes) == sum(old_sizes))
  stopifnot(sum(new_sizes) == d$n_snp)

  d$blocks      <- new_blocks
  d$block_sizes <- new_sizes
  d$n_blocks    <- length(new_blocks)
  d$resplit     <- TRUE
  d$resplit_max <- MAX_BLOCK
  d$resplit_target <- TARGET

  saveRDS(d, file.path(resplit_dir, sprintf("chr%d_blocks.rds", chr)))

  cat(sprintf("  chr%2d: %d bloques (max %d) → %d bloques (max %d)\n",
      chr, length(old_sizes), max(old_sizes),
      length(new_blocks), max(new_sizes)))
}

cat("\n ✓ Re-partición completada en blocks_resplit/\n")
REOF

apptainer exec \
  --bind /mnt/cephfs:/mnt/cephfs \
  --bind "$TMPDIR:/tmp,$TMPDIR:/var/tmp" \
  --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 \
  --env BLOCKS_DIR="$BLOCKS_DIR" --env RESPLIT_DIR="$RESPLIT_DIR" \
  "$SIF" /opt/miniforge3/envs/r-bio/bin/Rscript "$R_SCRIPT"

rm -rf "$TMPDIR"

echo ""
echo "============================================================"
echo "  RE-PARTICIÓN FIN — $(date)"
echo "============================================================"
