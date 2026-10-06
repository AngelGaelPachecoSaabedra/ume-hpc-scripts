#!/bin/bash
#SBATCH --job-name=ldp2c_blocks
#SBATCH --account=researchers
#SBATCH --qos=vip
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=150G
#SBATCH --output=/mnt/cephfs/orgs/home/angel.pacheco/logs/ldp2c_blocks_%A_%a.out
#SBATCH --error=/mnt/cephfs/orgs/home/angel.pacheco/logs/ldp2c_blocks_%A_%a.err
#SBATCH --mail-user=angelpachecosaavedra@gmail.com
#SBATCH --mail-type=FAIL

set -euo pipefail

# ================================================================
# PIPELINE C — STEP 1: bloques LD por cromosoma (array paralelo)
# Lanzar: sbatch --array=1-23%4 custom_step1_blocks.sbatch
# (23 = chrX).
# Requiere Step 0 completado (sumstats_qc.rds).
#
# Cada tarea: match -> interpolacion a mapa genetico (3 cM) ->
# snp_cor (ventana 3 cM REAL = size 3/1000) -> snp_ldsplit (retry)
# -> re-particion a <=2,500 SNPs (subdiv) -> guarda chrN_blocks.rds.
#
# CONFIG VALIDADA (produjo los 23 bloques del job 5693):
#   - size = 3 / 1000   (3 cM reales; size=3 cuelga = 3000 cM densa)
#   - subdiv() a <=2,500 (bloques >2,500 -> MCMC diverge a NaN)
# ================================================================

CHR=$SLURM_ARRAY_TASK_ID

SIF="/mnt/cephfs/biocontainers/images/jupyter-biotools-1.5.sif"
PLINK_DIR="/mnt/cephfs/hot_nvme/mcps/imputed-topmed/plink_files/maximally_unrelated"
OUT_DIR="/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes/ldpred2_custom"
BLOCKS_DIR="${OUT_DIR}/blocks"
mkdir -p "$BLOCKS_DIR"

TMPDIR="/mnt/cephfs/m.2/scratch/$USER/ldp2c_blocks_tmp/${SLURM_ARRAY_JOB_ID}_${CHR}"
mkdir -p "$TMPDIR"; chmod 700 "$TMPDIR"
export TMPDIR TMP="$TMPDIR" TEMP="$TMPDIR"

# Etiqueta de cromosoma (23 -> X)
if [[ "$CHR" -eq 23 ]]; then CHR_LABEL="X"; else CHR_LABEL="$CHR"; fi

echo "============================================================"
echo "  PIPELINE C STEP 1 — chr${CHR_LABEL}  Job ${SLURM_ARRAY_JOB_ID}_${CHR}"
echo "  Nodo: $SLURMD_NODENAME   $(date)"
echo "============================================================"

R_SCRIPT="$TMPDIR/blocks_chr${CHR}.R"
cat << 'REOF' > "$R_SCRIPT"
Sys.setenv(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
           MKL_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1",
           BLIS_NUM_THREADS="1", NUMEXPR_NUM_THREADS="1")
suppressMessages({ library(bigsnpr); library(Matrix); library(data.table) })
options(bigstatsr.check.parallel.blas=FALSE)
options(bigstatsr.ncores.max=24L)

bv <- packageVersion("bigsnpr")
cat(sprintf("  bigsnpr: %s\n", as.character(bv)))
if (as.character(bv)=="1.12.21") stop("bigsnpr 1.12.21 roto.")

chr        <- as.integer(Sys.getenv("CHR"))
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
  saveRDS(list(chr=chr, chr_label=chr_label, status="failed", stage=stage, msg=msg,
               n_snp=0L, n_blocks=0L, blocks=list(), info_snp=NULL), out_file)
  quit(save="no", status=0)
}

# Cargar sumstats pre-pruned de Step 0
sumstats <- readRDS(file.path(out_dir, "sumstats_qc.rds"))
sumstats_chr <- sumstats[sumstats$chr == chr, ]
cat(sprintf("\n==== chr%s: %s variantes GWAS pruned ====\n",
    chr_label, format(nrow(sumstats_chr), big.mark=",")))
if (nrow(sumstats_chr) < 50) save_failed("muy pocas variantes", "input")

prefix <- file.path(plink_dir, sprintf("mcps-freeze150k_qcd_chr%s_ivs", chr_label))
bed_file <- paste0(prefix, ".bed")
if (!file.exists(bed_file)) save_failed(paste("no existe", bed_file), "read_bed")

bk_file <- file.path(tmp_dir, sprintf("chr%s", chr_label))
snp_readBed(bed_file, backingfile=bk_file)
obj <- snp_attach(paste0(bk_file, ".rds"))
G <- obj$genotypes

# chrX: hombres ya pseudo-diploides (0/2) en bfile — no recodificar.
if (chr == 23L) cat("  chrX: hombres pseudo-diploides 0/2 (verificado HAP cts)\n")

# Imputacion defensiva
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

# sd_ss diagnostico (n_eff ya correcto en sumstats)
sd_geno <- sd_vals[keep_var]
sd_ss <- 2 / sqrt(info_snp$n_eff * info_snp$beta_se^2 + info_snp$beta^2)
cat(sprintf("  sd_ss/sd_geno ratio: %.3f\n", median(sd_ss/sd_geno)))

# ── Inyeccion de Mapa Genetico (cM) ────────────────────────────────
cat("\n  Interpolando mapa genetico b38 a cM...\n")
map_file <- sprintf("/mnt/cephfs/panel_referencia/maps/chr%s.b38.gmap.gz", chr_label)

# Cargar forzando los nombres de columna y ordenar por posicion fisica
map_ref <- fread(map_file, col.names = c("pos", "chr", "cM"), showProgress = FALSE)
map_ref <- map_ref[order(map_ref$pos), ]

# Interpolar posiciones fisicas del GWAS a distancias geneticas (cM)
# rule = 2: SNPs fuera del rango del mapa heredan el cM mas cercano
cM_interpolated <- approx(x = map_ref$pos,
                          y = map_ref$cM,
                          xout = info_snp$pos,
                          rule = 2)$y

cat(sprintf("  cM min: %.2f | cM max: %.2f\n", min(cM_interpolated), max(cM_interpolated)))
cat(sprintf("  snp_cor en ventana genetica de 3 cM (%s vars)...\n", format(length(ind_col), big.mark=",")))
t0 <- Sys.time()

# ── Calculo de correlacion ─────────────────────────────────────────
# CRITICO: size = 3 / 1000  -> ventana de 3 cM REALES.
# bigsnpr multiplica size por 1000 internamente cuando infos.pos esta en cM
# (ver doc de snp_cor). Usar size = 3 da 3000 cM = sin ventana = matriz
# densa que NO termina (cuelgues de 19-23 h en los jobs 5657 y 5722).
corr0 <- tryCatch(
  snp_cor(Gna = G,
          ind.col = ind_col,
          infos.pos = cM_interpolated,
          size = 3 / 1000,             # <-- 3 cM reales
          ncores = ncores,
          thr_r2 = 0.005),
  error=function(e) save_failed(conditionMessage(e), "snp_cor")
)

t_cor <- as.numeric(difftime(Sys.time(), t0, units="mins"))
m_chr <- ncol(corr0)
cat(sprintf("  snp_cor: %.1f min, dim=%d, nnz=%s\n", t_cor, m_chr, format(length(corr0@x), big.mark=",")))

# ── snp_ldsplit con retry escalado ─────────────────────────────────
max_size_seq <- round(seq_log(m_chr/30, m_chr/5, length.out=20))
max_size_seq <- unique(max_size_seq[max_size_seq>=100 & max_size_seq<=m_chr])
if (length(max_size_seq)==0) max_size_seq <- min(m_chr, 5000)
is_valid <- function(rs) !is.null(rs) && nrow(rs)>0 && is.numeric(rs$cost) &&
                         is.numeric(rs$n_block) && any(!is.na(rs$cost) & rs$n_block>0)
res_split <- NULL; thr_used <- NA
for (thr in c(0.02,0.05,0.1,0.2,0.3)) {
  rs <- tryCatch(snp_ldsplit(corr0, thr_r2=thr, min_size=50, max_size=max_size_seq, max_K=500),
                 error=function(e) NULL)
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

# ── Re-particion forzada a <=TARGET_BLOCK ──────────────────────────
# CRITICO: bloques >2,500 SNPs en cohorte admixta hacen que el MCMC de
# LDpred2-auto diverja a NaN (validado en Pipeline B). subdiv() acota.
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

saveRDS(list(chr=chr, chr_label=chr_label, status="ok", n_snp=n_snp,
             n_blocks=length(blocks), block_sizes=block_sizes, blocks=blocks,
             info_snp=info_snp[, c("chr","pos","a0","a1","beta","beta_se","n_eff","_NUM_ID_")],
             thr_r2_split=thr_used, t_cor_min=t_cor,
             bigsnpr_version=as.character(bv)), out_file)
cat(sprintf("\n  OK chr%s: %d bloques, %s SNPs, snp_cor %.1f min\n",
    chr_label, length(blocks), format(n_snp, big.mark=","), t_cor))
REOF

apptainer exec --bind /mnt/cephfs:/mnt/cephfs --bind "$TMPDIR:/tmp,$TMPDIR:/var/tmp" \
    --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 --env MKL_NUM_THREADS=1 \
    --env BLIS_NUM_THREADS=1 --env NUMEXPR_NUM_THREADS=1 \
    --env CHR="$CHR" --env PLINK_DIR="$PLINK_DIR" \
    --env OUT_DIR="$OUT_DIR" --env BLOCKS_DIR="$BLOCKS_DIR" \
    "$SIF" /opt/miniforge3/envs/r-bio/bin/Rscript "$R_SCRIPT"

rm -rf "$TMPDIR"
echo "============================================================"
echo "  STEP 1 chr${CHR_LABEL} FIN — $(date)"
echo "============================================================"
