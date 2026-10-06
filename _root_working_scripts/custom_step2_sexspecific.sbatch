#!/bin/bash
#SBATCH --job-name=ldp2c_run_sex
#SBATCH --account=researchers
#SBATCH --qos=vip
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=150G
#SBATCH --output=/mnt/cephfs/orgs/home/angel.pacheco/logs/ldp2c_run_sex_%j.out
#SBATCH --error=/mnt/cephfs/orgs/home/angel.pacheco/logs/ldp2c_run_sex_%j.err
#SBATCH --mail-user=angelpachecosaavedra@gmail.com
#SBATCH --mail-type=END,FAIL

set -euo pipefail

# ================================================================
# PIPELINE C SEX — STEP 2: LDpred2-auto + PRS POR SEXO
# Lanzar DOS veces:
#   sbatch --export=SEX=females custom_step2_sexspecific.sbatch
#   sbatch --export=SEX=males   custom_step2_sexspecific.sbatch
# Requiere Step 1 sex (blocks_<SEX>/chrN_blocks.rds).
#
# Scorea a TODOS los individuos con las betas de este sexo. El Step 3
# (combine) selecciona despues el score del propio sexo de cada persona.
#
# FIXES: shrink_corr=0.95 (valido en 1.10.8, verificado), length.out=48,
#        num_iter=1000, mem=150G, operadores '<' correctos (no '<<').
# ================================================================

SEX=${SEX:-females}

SIF="/mnt/cephfs/biocontainers/images/jupyter-biotools-1.5.sif"
PLINK_DIR="/mnt/cephfs/hot_nvme/mcps/imputed-topmed/plink_files/maximally_unrelated"
OUT_DIR="/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes/ldpred2_custom"
BLOCKS_DIR="${OUT_DIR}/blocks_${SEX}"

TMPDIR="/mnt/cephfs/m.2/scratch/$USER/ldp2c_run_${SEX}_$SLURM_JOB_ID"
mkdir -p "$TMPDIR"; chmod 700 "$TMPDIR"
export TMPDIR TMP="$TMPDIR" TEMP="$TMPDIR"

echo "============================================================"
echo "  STEP 2 SEX — LDpred2 ${SEX} — Job $SLURM_JOB_ID"
echo "  Nodo: $SLURMD_NODENAME   $(date)"
echo "============================================================"

PRESENT=0
for CHR in $(seq 1 22) X; do
  [[ -f "${BLOCKS_DIR}/chr${CHR}_blocks.rds" ]] && PRESENT=$((PRESENT+1)) || echo "  FALTA: chr${CHR}"
done
echo "  $PRESENT/23 bloques presentes (${SEX})"
if [[ $PRESENT -lt 20 ]]; then echo "ERROR: faltan demasiados"; exit 1; fi

R_SCRIPT="$TMPDIR/run.R"
cat << 'REOF' > "$R_SCRIPT"
Sys.setenv(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
           MKL_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1",
           BLIS_NUM_THREADS="1", NUMEXPR_NUM_THREADS="1")
suppressMessages({ library(bigsnpr); library(bigsparser); library(data.table); library(Matrix) })
options(bigstatsr.check.parallel.blas=FALSE); options(bigstatsr.ncores.max=48L)

bv <- packageVersion("bigsnpr")
cat(sprintf("  bigsnpr: %s\n", as.character(bv)))
if (as.character(bv)=="1.12.21") stop("bigsnpr 1.12.21 roto.")

sex        <- Sys.getenv("SEX")
ncores     <- 24L
plink_dir  <- Sys.getenv("PLINK_DIR")
out_dir    <- Sys.getenv("OUT_DIR")
blocks_dir <- Sys.getenv("BLOCKS_DIR")
tmp_dir    <- Sys.getenv("TMPDIR")

betas_file <- file.path(out_dir, sprintf("ldpred2_custom_results_%s.rds", sex))
if (file.exists(betas_file)) { cat("Ya existe resultado, saltando.\n"); quit(save="no", status=0) }

safe_as_SFBM <- function(corr_mat, backingfile, compact=TRUE) {
  if (!is(corr_mat, "dsCMatrix")) corr_mat <- as(corr_mat, "symmetricMatrix")
  col_nnz <- diff(corr_mat@p); empty_cols <- which(col_nnz==0)
  if (length(empty_cols) > 0) {
    cat(sprintf("  safe_as_SFBM: %d columnas vacias, parchando\n", length(empty_cols)))
    n <- ncol(corr_mat); old_nnz <- length(corr_mat@x)
    new_i <- integer(old_nnz+length(empty_cols)); new_x <- numeric(old_nnz+length(empty_cols))
    new_p <- integer(n+1L); pos <- 0L
    for (j in seq_len(n)) {
      os <- corr_mat@p[j]+1L; oe <- corr_mat@p[j+1L]; nic <- oe-corr_mat@p[j]
      if (nic > 0L) { idx <- (pos+1L):(pos+nic); new_i[idx] <- corr_mat@i[os:oe]; new_x[idx] <- corr_mat@x[os:oe]; pos <- pos+nic }
      else { pos <- pos+1L; new_i[pos] <- j-1L; new_x[pos] <- 1e-10 }
      new_p[j+1L] <- pos
    }
    corr_mat@i <- new_i[seq_len(pos)]; corr_mat@x <- new_x[seq_len(pos)]; corr_mat@p <- new_p
  }
  as_SFBM(corr_mat, backingfile=backingfile, compact=compact)
}

cat("\n[1] Cargando bloques...\n")
all_blocks <- list(); all_info <- list(); total_snps <- 0L; total_blocks <- 0L
chr_ok <- c(); chr_failed <- c()
for (chr in c(1:22, 23)) {
  cl <- if (chr==23L) "X" else as.character(chr)
  f <- file.path(blocks_dir, sprintf("chr%s_blocks.rds", cl))
  if (!file.exists(f)) { chr_failed <- c(chr_failed, cl); next }
  d <- readRDS(f)
  if (is.null(d$status) || d$status!="ok") { cat(sprintf("  chr%s: no-ok\n", cl)); chr_failed <- c(chr_failed, cl); next }
  all_blocks <- c(all_blocks, d$blocks)
  all_info[[length(all_info)+1]] <- as.data.frame(d$info_snp)
  total_snps <- total_snps + d$n_snp; total_blocks <- total_blocks + d$n_blocks
  chr_ok <- c(chr_ok, cl)
  cat(sprintf("  chr%s: %d bloques, %s SNPs\n", cl, d$n_blocks, format(d$n_snp, big.mark=",")))
}
info_gw <- do.call(rbind, all_info)
cat(sprintf("\n  Total: %s SNPs, %d bloques, %d chr OK\n",
    format(total_snps, big.mark=","), total_blocks, length(chr_ok)))
if (length(chr_failed) > 0) cat(sprintf("  Fallidos: %s\n", paste(chr_failed, collapse=",")))

nnz_total <- sum(sapply(all_blocks, function(b) length(b@x)))
cat(sprintf("  nnz total: %s\n", format(nnz_total, big.mark=",")))
if (nnz_total > 2^31-1) stop("nnz excede int32")

cat("\n[2] SFBM block-diagonal...\n")
corr_gw <- as(Matrix::bdiag(all_blocks), "symmetricMatrix"); rm(all_blocks); gc()
sfbm <- safe_as_SFBM(corr_gw, backingfile=file.path(tmp_dir, "corr_gw"), compact=TRUE)
cat(sprintf("  SFBM: %.1f GB\n", file.info(paste0(file.path(tmp_dir,"corr_gw"),".sbk"))$size/1024^3))
rm(corr_gw); gc()

cat("\n[3] LDpred2-auto GENOME-WIDE...\n")
df_beta <- data.frame(beta=info_gw$beta, beta_se=info_gw$beta_se, n_eff=info_gw$n_eff)
cat(sprintf("    n_eff = %.0f, %s variantes, %d bloques\n",
    median(df_beta$n_eff), format(nrow(df_beta), big.mark=","), total_blocks))

t0 <- Sys.time()
ma <- snp_ldpred2_auto(corr=sfbm, df_beta=df_beta, h2_init=0.1,
    vec_p_init=seq_log(1e-4, 0.5, length.out=48),
    burn_in=1000L, num_iter=1000L, report_step=100L,
    ncores=ncores, sparse=TRUE, allow_jump_sign=FALSE, shrink_corr=0.95)
t_ldp <- as.numeric(difftime(Sys.time(), t0, units="mins"))
cat(sprintf("  LDpred2: %.1f min\n", t_ldp))

h2 <- sapply(ma, function(x) x$h2_est); p <- sapply(ma, function(x) x$p_est)
keep <- which(!is.na(h2) & h2 > 0.001 & h2 < 0.8 & !is.na(p) & p > 1e-6 & p < 0.99)
cat(sprintf("  Cadenas validas: %d/48\n", length(keep)))
if (length(keep)>0) cat(sprintf("  h2 median: %.6f, p median: %.2e\n", median(h2[keep]), median(p[keep])))

if (length(keep) < 3L) {
  cat("\n  [RETRY] shrink_corr=0.5...\n")
  ma <- snp_ldpred2_auto(corr=sfbm, df_beta=df_beta, h2_init=0.05,
      vec_p_init=seq_log(1e-5, 0.3, length.out=48),
      burn_in=2000L, num_iter=1000L, report_step=100L,
      ncores=ncores, sparse=TRUE, allow_jump_sign=FALSE, shrink_corr=0.5)
  h2 <- sapply(ma, function(x) x$h2_est); p <- sapply(ma, function(x) x$p_est)
  keep <- which(!is.na(h2) & h2 > 0.001 & h2 < 0.8 & !is.na(p) & p > 1e-6 & p < 0.99)
  cat(sprintf("  RETRY: %d/48\n", length(keep)))
}

cat("\n[4] Betas + PRS (todos los individuos)...\n")
if (length(keep) < 1L) {
  cat("  DIVERGENCIA TOTAL\n")
  saveRDS(list(status="diverged", sex=sex, h2_all=h2, p_all=p, chr_ok=chr_ok, chr_failed=chr_failed), betas_file)
} else {
  bm <- sapply(ma[keep], function(x) x$beta_est); if (is.null(dim(bm))) bm <- matrix(bm, ncol=1)
  beta_final <- rowMeans(bm, na.rm=TRUE); beta_final[!is.finite(beta_final)] <- 0
  h2_med <- median(h2[keep]); p_med <- median(p[keep])
  cat(sprintf("  h2=%.6f, p=%.2e, nonzero=%s\n", h2_med, p_med, format(sum(beta_final!=0), big.mark=",")))

  chr_list <- sort(unique(info_gw$chr)); prs_total <- NULL; fam <- NULL
  for (chr in chr_list) {
    cl <- if (chr==23L) "X" else as.character(chr)
    idx <- which(info_gw$chr == chr); if (length(idx)==0) next
    prefix <- file.path(plink_dir, sprintf("mcps-freeze150k_qcd_chr%s_ivs", cl))
    bk <- file.path(tmp_dir, sprintf("prs_chr%s", cl))
    snp_readBed(paste0(prefix, ".bed"), backingfile=bk)
    obj <- snp_attach(paste0(bk, ".rds")); G <- obj$genotypes
    n_na <- big_counts(G, ind.col=info_gw$`_NUM_ID_`[idx])[4,] |> sum()
    if (n_na > 0) G <- snp_fastImputeSimple(G, "mean0", ncores=ncores)
    prs_chr <- big_prodVec(X=G, y.col=beta_final[idx], ind.col=info_gw$`_NUM_ID_`[idx], ncores=ncores)
    if (is.null(prs_total)) { prs_total <- prs_chr; fam <- obj$fam } else prs_total <- prs_total + prs_chr
    cat(sprintf("    chr%s: %d SNPs\n", cl, length(idx)))
    rm(obj, G); gc()
  }

  saveRDS(list(status="ok", sex=sex, h2=h2_med, p=p_med, n_valid_chains=length(keep),
               n_snp_total=total_snps, n_blocks=total_blocks, h2_all=h2, p_all=p,
               t_ldpred_min=t_ldp, beta_final=beta_final,
               info_snp=info_gw[, c("chr","pos","a0","a1")],
               prs_total=prs_total, fam=fam, chr_ok=chr_ok, chr_failed=chr_failed,
               bigsnpr_version=as.character(bv)), betas_file)

  # TSV simple: IID + PRS de este sexo, para TODOS los individuos (sin filtrar)
  df_out <- data.frame(IID=as.character(fam$sample.ID), PRS=prs_total, stringsAsFactors=FALSE)
  colnames(df_out)[2] <- sprintf("PRS_%s", sex)
  write.table(df_out, file.path(out_dir, sprintf("prs_%s.tsv", sex)),
              sep="\t", row.names=FALSE, quote=FALSE)
  cat(sprintf("\n  OK %s: h2=%.6f, %d/48 cadenas, %s indiv -> prs_%s.tsv\n",
      sex, h2_med, length(keep), format(length(prs_total), big.mark=","), sex))
}
REOF

apptainer exec --bind /mnt/cephfs:/mnt/cephfs --bind "$TMPDIR:/tmp,$TMPDIR:/var/tmp" \
    --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 --env MKL_NUM_THREADS=1 \
    --env BLIS_NUM_THREADS=1 --env NUMEXPR_NUM_THREADS=1 \
    --env SEX="$SEX" --env PLINK_DIR="$PLINK_DIR" --env OUT_DIR="$OUT_DIR" --env BLOCKS_DIR="$BLOCKS_DIR" \
    "$SIF" /opt/miniforge3/envs/r-bio/bin/Rscript "$R_SCRIPT"

rm -rf "$TMPDIR"
echo "============================================================"
echo "  STEP 2 SEX ${SEX} FIN — $(date)"
echo "============================================================"
