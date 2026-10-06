#!/bin/bash
#SBATCH --job-name=ldp2c_blk_sex
#SBATCH --account=researchers
#SBATCH --qos=vip
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=150G
#SBATCH --output=/mnt/cephfs/orgs/home/angel.pacheco/logs/ldp2c_blk_sex_%A_%a.out
#SBATCH --error=/mnt/cephfs/orgs/home/angel.pacheco/logs/ldp2c_blk_sex_%A_%a.err
#SBATCH --mail-user=angelpachecosaavedra@gmail.com
#SBATCH --mail-type=FAIL

set -euo pipefail

# ================================================================
# PIPELINE C SEX — STEP 1: bloques LD por cromosoma y por sexo
# Lanzar DOS veces (una por sexo):
#   sbatch --array=1-23%4 --export=SEX=females custom_step1_sexspecific.sbatch
#   sbatch --array=1-23%4 --export=SEX=males   custom_step1_sexspecific.sbatch
# Requiere Step 0 sex (sumstats_qc_<SEX>.rds).
#
# Config VALIDADA: size=3/1000 (3 cM reales) + subdiv() a <=2,500.
# ================================================================

CHR=$SLURM_ARRAY_TASK_ID
SEX=${SEX:-females}

SIF="/mnt/cephfs/biocontainers/images/jupyter-biotools-1.5.sif"
PLINK_DIR="/mnt/cephfs/hot_nvme/mcps/imputed-topmed/plink_files/maximally_unrelated"
OUT_DIR="/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes/ldpred2_custom"
BLOCKS_DIR="${OUT_DIR}/blocks_${SEX}"
mkdir -p "$BLOCKS_DIR"

TMPDIR="/mnt/cephfs/scratch/$USER/ldp2c_blk_tmp/${SLURM_ARRAY_JOB_ID}_${CHR}_${SEX}"
mkdir -p "$TMPDIR"; chmod 700 "$TMPDIR"
export TMPDIR TMP="$TMPDIR" TEMP="$TMPDIR"

if [[ "$CHR" -eq 23 ]]; then CHR_LABEL="X"; else CHR_LABEL="$CHR"; fi

echo "============================================================"
echo "  STEP 1 SEX — chr${CHR_LABEL} (${SEX})  Job ${SLURM_ARRAY_JOB_ID}_${CHR}"
echo "  Nodo: $SLURMD_NODENAME   $(date)"
echo "============================================================"

R_SCRIPT="$TMPDIR/blocks.R"
cat << 'REOF' > "$R_SCRIPT"
Sys.setenv(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
           MKL_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1",
           BLIS_NUM_THREADS="1", NUMEXPR_NUM_THREADS="1")
suppressMessages({ library(bigsnpr); library(Matrix); library(data.table) })
options(bigstatsr.check.parallel.blas=FALSE); options(bigstatsr.ncores.max=24L)

bv <- packageVersion("bigsnpr")
cat(sprintf("  bigsnpr: %s\n", as.character(bv)))
if (as.character(bv)=="1.12.21") stop("bigsnpr 1.12.21 roto.")

chr        <- as.integer(Sys.getenv("CHR"))
sex        <- Sys.getenv("SEX")
chr_label  <- if (chr==23L) "X" else as.character(chr)
ncores     <- 24L
plink_dir  <- Sys.getenv("PLINK_DIR")
out_dir    <- Sys.getenv("OUT_DIR")
blocks_dir <- Sys.getenv("BLOCKS_DIR")
tmp_dir    <- Sys.getenv("TMPDIR")
TARGET_BLOCK <- 2500L

out_file <- file.path(blocks_dir, sprintf("chr%s_blocks.rds", chr_label))
if (file.exists(out_file)) { cat("Ya existe, saltando.\n"); quit(save="no", status=0) }

save_failed <- function(msg, stage) {
  cat(sprintf("  *** FALLO %s: %s ***\n", stage, msg))
  saveRDS(list(chr=chr, chr_label=chr_label, sex=sex, status="failed", stage=stage, msg=msg,
               n_snp=0L, n_blocks=0L, blocks=list(), info_snp=NULL), out_file)
  quit(save="no", status=0)
}

sumstats <- readRDS(file.path(out_dir, sprintf("sumstats_qc_%s.rds", sex)))
sumstats_chr <- sumstats[sumstats$chr == chr, ]
cat(sprintf("\n==== chr%s (%s): %s variantes pruned ====\n",
    chr_label, sex, format(nrow(sumstats_chr), big.mark=",")))
if (nrow(sumstats_chr) < 50) save_failed("muy pocas variantes", "input")

prefix <- file.path(plink_dir, sprintf("mcps-freeze150k_qcd_chr%s_ivs", chr_label))
if (!file.exists(paste0(prefix, ".bed"))) save_failed("no existe bed", "read_bed")

bk_file <- file.path(tmp_dir, sprintf("chr%s", chr_label))
snp_readBed(paste0(prefix, ".bed"), backingfile=bk_file)
obj <- snp_attach(paste0(bk_file, ".rds")); G <- obj$genotypes
if (chr == 23L) cat("  chrX: hombres pseudo-diploides 0/2\n")

n_na <- big_counts(G, ind.col=seq_len(ncol(G)))[4, ] |> sum()
if (n_na > 0) { cat(sprintf("  NAs=%s, imputando\n", format(n_na,big.mark=","))); G <- snp_fastImputeSimple(G, "mean0", ncores=ncores) }

map_df <- data.frame(chr=obj$map$chromosome, pos=obj$map$physical.pos,
                     a1=obj$map$allele1, a0=obj$map$allele2, stringsAsFactors=FALSE)
info_snp <- snp_match(sumstats_chr, map_df, strand_flip=FALSE, join_by_pos=TRUE, match.min.prop=0.01)
cat(sprintf("  Matched: %s\n", format(nrow(info_snp), big.mark=",")))

ind_col <- info_snp$`_NUM_ID_`
sd_vals <- big_scale()(G, ind.col=ind_col)$scale
keep_var <- which(sd_vals > 1e-6)
if (length(keep_var) < nrow(info_snp)) { info_snp <- info_snp[keep_var, ]; ind_col <- info_snp$`_NUM_ID_` }
if (nrow(info_snp) < 50) save_failed("muy pocas tras QC", "qc")
n_snp <- nrow(info_snp)

sd_geno <- sd_vals[keep_var]
sd_ss <- 2 / sqrt(info_snp$n_eff * info_snp$beta_se^2 + info_snp$beta^2)
cat(sprintf("  sd_ss/sd_geno ratio: %.3f\n", median(sd_ss/sd_geno)))

cat("\n  Interpolando mapa genetico b38 a cM...\n")
map_file <- sprintf("/mnt/cephfs/panel_referencia/maps/chr%s.b38.gmap.gz", chr_label)
map_ref <- fread(map_file, col.names=c("pos","chr","cM"), showProgress=FALSE)
map_ref <- map_ref[order(map_ref$pos), ]
cM_interpolated <- approx(x=map_ref$pos, y=map_ref$cM, xout=info_snp$pos, rule=2)$y
cat(sprintf("  cM min: %.2f | cM max: %.2f\n", min(cM_interpolated), max(cM_interpolated)))

cat(sprintf("  snp_cor ventana 3 cM (%s vars)...\n", format(length(ind_col), big.mark=",")))
t0 <- Sys.time()
# CRITICO: size = 3/1000 = 3 cM reales (size=3 cuelga = 3000 cM densa)
corr0 <- tryCatch(
  snp_cor(Gna=G, ind.col=ind_col, infos.pos=cM_interpolated,
          size=3/1000, ncores=ncores, thr_r2=0.005),
  error=function(e) save_failed(conditionMessage(e), "snp_cor"))
t_cor <- as.numeric(difftime(Sys.time(), t0, units="mins"))
m_chr <- ncol(corr0)
cat(sprintf("  snp_cor: %.1f min, dim=%d, nnz=%s\n", t_cor, m_chr, format(length(corr0@x), big.mark=",")))

max_size_seq <- round(seq_log(m_chr/30, m_chr/5, length.out=20))
max_size_seq <- unique(max_size_seq[max_size_seq>=100 & max_size_seq<=m_chr])
if (length(max_size_seq)==0) max_size_seq <- min(m_chr, 5000)
is_valid <- function(rs) !is.null(rs) && nrow(rs)>0 && is.numeric(rs$cost) && is.numeric(rs$n_block) && any(!is.na(rs$cost) & rs$n_block>0)
res_split <- NULL; thr_used <- NA
for (thr in c(0.02,0.05,0.1,0.2,0.3)) {
  rs <- tryCatch(snp_ldsplit(corr0, thr_r2=thr, min_size=50, max_size=max_size_seq, max_K=500), error=function(e) NULL)
  if (is_valid(rs)) { res_split <- rs; thr_used <- thr; break }
}
if (is.null(res_split)) {
  blocks0 <- list(corr0); cat("  ldsplit fallo -> 1 bloque (se re-particiona)\n")
} else {
  vi <- which(!is.na(res_split$cost) & !is.na(res_split$n_block) & res_split$n_block>0)
  res_split <- res_split[vi, ]
  bi <- which.min(res_split$cost + 0.5*log(res_split$n_block))
  bsz <- res_split$all_size[[bi]]; nb <- length(bsz)
  bs <- c(1L, cumsum(bsz[-nb])+1L); be <- cumsum(bsz)
  blocks0 <- lapply(seq_len(nb), function(b) corr0[bs[b]:be[b], bs[b]:be[b], drop=FALSE])
  cat(sprintf("  ldsplit: %d bloques (thr_r2=%.3f, max=%d)\n", nb, thr_used, max(bsz)))
}

# CRITICO: subdiv a <=2,500 (bloques mayores -> MCMC diverge a NaN)
subdiv <- function(blk, target) {
  m <- nrow(blk); if (m<=target) return(list(blk))
  n_sub <- ceiling(m/target); sizes <- rep(floor(m/n_sub), n_sub)
  extra <- m-sum(sizes); if (extra>0) sizes[seq_len(extra)] <- sizes[seq_len(extra)]+1L
  st <- c(1L, cumsum(sizes[-n_sub])+1L); en <- cumsum(sizes)
  lapply(seq_len(n_sub), function(i) blk[st[i]:en[i], st[i]:en[i], drop=FALSE])
}
blocks <- list(); for (b in blocks0) blocks <- c(blocks, subdiv(b, TARGET_BLOCK))
block_sizes <- sapply(blocks, nrow)
stopifnot(sum(block_sizes) == m_chr)
cat(sprintf("  Re-particionado: %d bloques (max=%d)\n", length(blocks), max(block_sizes)))

saveRDS(list(chr=chr, chr_label=chr_label, sex=sex, status="ok", n_snp=n_snp,
             n_blocks=length(blocks), block_sizes=block_sizes, blocks=blocks,
             info_snp=info_snp[, c("chr","pos","a0","a1","beta","beta_se","n_eff","_NUM_ID_")],
             thr_r2_split=thr_used, t_cor_min=t_cor, bigsnpr_version=as.character(bv)), out_file)
cat(sprintf("\n  OK chr%s (%s): %d bloques, %s SNPs, %.1f min\n",
    chr_label, sex, length(blocks), format(n_snp, big.mark=","), t_cor))
REOF

apptainer exec --bind /mnt/cephfs:/mnt/cephfs --bind "$TMPDIR:/tmp,$TMPDIR:/var/tmp" \
    --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 --env MKL_NUM_THREADS=1 \
    --env BLIS_NUM_THREADS=1 --env NUMEXPR_NUM_THREADS=1 \
    --env CHR="$CHR" --env SEX="$SEX" --env PLINK_DIR="$PLINK_DIR" \
    --env OUT_DIR="$OUT_DIR" --env BLOCKS_DIR="$BLOCKS_DIR" \
    "$SIF" /opt/miniforge3/envs/r-bio/bin/Rscript "$R_SCRIPT"

rm -rf "$TMPDIR"
echo "============================================================"
echo "  STEP 1 SEX chr${CHR_LABEL} (${SEX}) FIN — $(date)"
echo "============================================================"
