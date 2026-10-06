#!/bin/bash
#SBATCH --job-name=ldpred2_gwide
#SBATCH --account=researchers
#SBATCH --qos=vip
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=170G
#SBATCH --output=/mnt/cephfs/orgs/home/angel.pacheco/logs/ldpred2_gwide_%j.out
#SBATCH --error=/mnt/cephfs/orgs/home/angel.pacheco/logs/ldpred2_gwide_%j.err
#SBATCH --mail-user=angelpachecosaavedra@gmail.com
#SBATCH --mail-type=END,FAIL

set -euo pipefail

# ================================================================
# LDpred2 GENOME-WIDE con snp_ldsplit (bloques LD independientes)
# bigsnpr 1.12.21 — Solución definitiva para chr19/chr22
#
# Arquitectura:
#   1. Calcula snp_cor por cromosoma (reutilizable)
#   2. snp_ldsplit → bloques LD independientes por chr
#   3. Zerifica correlaciones entre bloques
#   4. Concatena en una SFBM genome-wide
#   5. Corre UN solo snp_ldpred2_auto con todos los SNPs
#
# Ventajas:
#   - h2 y p estimados globalmente (no per-chr)
#   - Errores contenidos dentro de bloques
#   - shrink_corr + use_MLE=FALSE + allow_jump_sign=FALSE
#   - n_eff per-SNP calibrado con sd_val
# ================================================================

SIF="/mnt/cephfs/orgs/home/angel.pacheco/bigsnpr_nativo.sif"
GWAS_DIR="/mnt/cephfs/orgs/home/angel.pacheco/gwas_miami"
OUT_DIR="/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes/ldpred2_gwide"
PREV_TMPDIR="/mnt/cephfs/scratch/angel.pacheco/ldpred2_hm3_tmp/2931"

mkdir -p "$OUT_DIR"

TMPDIR="/mnt/cephfs/scratch/$USER/ldpred2_gwide_tmp/$SLURM_JOB_ID"
mkdir -p "$TMPDIR"; chmod 700 "$TMPDIR"
export TMPDIR TMP="$TMPDIR" TEMP="$TMPDIR"

echo "============================================================"
echo "  LDpred2 GENOME-WIDE (snp_ldsplit, v1.12.21)"
echo "  Job: $SLURM_JOB_ID   Nodo: $SLURMD_NODENAME"
echo "  $(date)"
echo "============================================================"

R_SCRIPT="$TMPDIR/ldpred2_gwide.R"

cat << 'REOF' > "$R_SCRIPT"
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
options(bigstatsr.ncores.max=48L)

cat(sprintf("  bigsnpr version: %s\n", as.character(packageVersion("bigsnpr"))))

ncores   <- 24L
prev_tmp <- Sys.getenv("PREV_TMPDIR", "/mnt/cephfs/scratch/angel.pacheco/ldpred2_hm3_tmp/2931")
out_dir  <- Sys.getenv("OUT_DIR", "/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes/ldpred2_gwide")
gwas_dir <- Sys.getenv("GWAS_DIR", "/mnt/cephfs/orgs/home/angel.pacheco/gwas_miami")
tmp_dir  <- Sys.getenv("TMPDIR")

betas_file <- file.path(out_dir, "ldpred2_gwide_results.rds")
if (file.exists(betas_file)) {
  cat("Ya existe resultado — saltando.\n")
  quit(save="no", status=0)
}

# ══════════════════════════════════════════════════════════════════
# [0] HapMap3 map
# ══════════════════════════════════════════════════════════════════
map_hm3 <- readRDS(file.path(
  "/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes/ldpred2",
  "map_hm3_ldpred2.rds"))
if ("pos_hg38" %in% names(map_hm3)) map_hm3$pos_use <- map_hm3$pos_hg38

# ══════════════════════════════════════════════════════════════════
# [1] GWAS meta-análisis IVW — todos los cromosomas
# ══════════════════════════════════════════════════════════════════
cat("\n[1] Cargando GWAS completo...\n")
gf <- fread(file.path(gwas_dir, "combined-gwas-results-all-ctrl-females.txt.gz"),
  select=c("CHROM","GENPOS","ID","ALLELE0","ALLELE1","A1FREQ","INFO","N","BETA","SE","P"),
  showProgress=FALSE)
gm <- fread(file.path(gwas_dir, "combined-gwas-results-all-ctrl-males.txt.gz"),
  select=c("CHROM","GENPOS","ID","ALLELE0","ALLELE1","A1FREQ","INFO","N","BETA","SE","P"),
  showProgress=FALSE)
gwas <- merge(gf, gm, by="ID", suffixes=c(".f",".m"))
rm(gf, gm); gc()

w_f <- 1/gwas$SE.f^2; w_m <- 1/gwas$SE.m^2
neff_f <- 4*17300*54196/(17300+54196)
neff_m <- 4*8522*26733/(8522+26733)
neff_ldpred2 <- neff_f + neff_m

sumstats_all <- data.frame(
  chr=as.integer(gwas$CHROM.f), pos=gwas$GENPOS.f,
  a1=gwas$ALLELE1.f, a0=gwas$ALLELE0.f,
  freq=(gwas$A1FREQ.f*gwas$N.f + gwas$A1FREQ.m*gwas$N.m)/(gwas$N.f+gwas$N.m),
  info=pmin(gwas$INFO.f, gwas$INFO.m),
  beta=(gwas$BETA.f*w_f + gwas$BETA.m*w_m)/(w_f+w_m),
  beta_se=sqrt(1/(w_f+w_m)), n_eff=neff_ldpred2, stringsAsFactors=FALSE)
sumstats_all$p <- 2*pnorm(-abs(sumstats_all$beta/sumstats_all$beta_se))
sumstats_all <- sumstats_all[sumstats_all$info>=0.8 & sumstats_all$freq>=0.01 &
  sumstats_all$freq<=0.99 & !is.na(sumstats_all$beta) &
  !is.na(sumstats_all$beta_se) & sumstats_all$beta_se>0, ]

key_hm3  <- paste(map_hm3$chr, map_hm3$pos_use, sep=":")
key_gwas <- paste(sumstats_all$chr, sumstats_all$pos, sep=":")
sumstats_all <- sumstats_all[key_gwas %in% key_hm3, ]
rm(gwas, w_f, w_m); gc()
cat(sprintf("  GWAS HM3-matched total: %s variantes\n",
    format(nrow(sumstats_all), big.mark=",")))

# ══════════════════════════════════════════════════════════════════
# [2] Por cromosoma: matching + snp_cor + snp_ldsplit
# ══════════════════════════════════════════════════════════════════
cat("\n[2] Procesando por cromosoma: matching + snp_cor + snp_ldsplit...\n")

all_corr_blocks <- list()   # lista de matrices dsCMatrix por bloque
all_info_snp    <- list()   # info_snp acumulado
all_sd_val      <- c()      # sd_val acumulado
global_offset   <- 0L       # offset para indexar bloques globalmente

for (chr in 1:22) {
  cat(sprintf("\n  --- chr%d ---\n", chr))

  sumstats_chr <- sumstats_all[sumstats_all$chr == chr, ]
  if (nrow(sumstats_chr) < 50) {
    cat(sprintf("    Saltando chr%d (%d variantes)\n", chr, nrow(sumstats_chr)))
    next
  }

  rds_chr <- file.path(prev_tmp, sprintf("chr%d.rds", chr))
  if (!file.exists(rds_chr)) {
    cat(sprintf("    No existe %s, saltando\n", basename(rds_chr)))
    next
  }

  obj <- snp_attach(rds_chr); G <- obj$genotypes
  map_df <- data.frame(chr=obj$map$chromosome, pos=obj$map$physical.pos,
    a1=obj$map$allele1, a0=obj$map$allele2, stringsAsFactors=FALSE)

  info_snp <- snp_match(sumstats_chr, map_df, strand_flip=FALSE,
    join_by_pos=TRUE, match.min.prop=0.01)
  ind_col <- info_snp$`_NUM_ID_`

  sd_vals <- big_scale()(G, ind.col=ind_col)$scale
  keep_var <- which(sd_vals > 1e-6)
  if (length(keep_var) < nrow(info_snp)) {
    info_snp <- info_snp[keep_var, ]
    ind_col  <- info_snp$`_NUM_ID_`
    sd_vals  <- sd_vals[keep_var]
  }
  cat(sprintf("    Matched: %s\n", format(nrow(info_snp), big.mark=",")))

  # snp_cor
  t0 <- Sys.time()
  corr0 <- snp_cor(Gna=G, ind.col=ind_col, infos.pos=info_snp$pos/1000,
    size=3000, ncores=ncores, thr_r2=0.005)
  t_cor <- as.numeric(difftime(Sys.time(), t0, units="mins"))
  cat(sprintf("    snp_cor: %.1f min, nnz=%s\n", t_cor,
      format(length(corr0@x), big.mark=",")))

  # snp_ldsplit — encontrar bloques independientes
  m_chr <- ncol(corr0)
  max_size_seq <- round(seq_log(m_chr/30, m_chr/5, length.out=8))
  max_size_seq <- max_size_seq[max_size_seq >= 100]
  if (length(max_size_seq) == 0) max_size_seq <- m_chr

  res_split <- snp_ldsplit(corr0, thr_r2=0.02, min_size=50,
    max_size=max_size_seq, max_K=500)

  # Elegir el split con mejor compromiso costo/bloques
  # Privé recomienda: cost bajo, pero no demasiados bloques
  best_idx <- which.min(res_split$cost + 0.5 * log(res_split$n_block))
  best_split <- res_split[best_idx, ]
  block_sizes <- best_split$all_size[[1]]
  n_blocks <- length(block_sizes)
  cat(sprintf("    snp_ldsplit: %d bloques, cost=%.1f\n",
      n_blocks, best_split$cost))

  # Zerificar correlaciones entre bloques
  block_starts <- c(1L, cumsum(block_sizes[-n_blocks]) + 1L)
  block_ends   <- cumsum(block_sizes)

  # Crear máscara: solo mantener entradas dentro de bloques
  corr_blocked <- corr0
  for (b in seq_len(n_blocks)) {
    bs <- block_starts[b]; be <- block_ends[b]
    # No necesitamos zerificar explícitamente — snp_ldsplit ya nos dice
    # los límites. Creamos sub-matrices y las añadimos a la lista.
  }

  # Añadir las sub-matrices de bloque a la lista global
  for (b in seq_len(n_blocks)) {
    bs <- block_starts[b]; be <- block_ends[b]
    block_corr <- corr0[bs:be, bs:be, drop=FALSE]
    all_corr_blocks[[length(all_corr_blocks) + 1]] <- block_corr
  }

  # Acumular info_snp y sd_val en el orden de los bloques
  all_info_snp[[length(all_info_snp) + 1]] <- info_snp
  all_sd_val <- c(all_sd_val, sd_vals)

  rm(corr0, obj, G); gc()
}

# ══════════════════════════════════════════════════════════════════
# [3] Construir SFBM block-diagonal genome-wide
# ══════════════════════════════════════════════════════════════════
cat("\n[3] Construyendo SFBM block-diagonal...\n")

# Concatenar info_snp
info_snp_gw <- do.call(rbind, all_info_snp)
sd_val_gw   <- all_sd_val
n_snp_total <- nrow(info_snp_gw)
cat(sprintf("  Total variantes genome-wide: %s en %d bloques\n",
    format(n_snp_total, big.mark=","), length(all_corr_blocks)))

# Construir la block-diagonal sparse matrix
corr_gw <- Matrix::bdiag(all_corr_blocks)
corr_gw <- as(corr_gw, "symmetricMatrix")
rm(all_corr_blocks); gc()

cat(sprintf("  Matrix block-diagonal: %dx%d, nnz=%s\n",
    nrow(corr_gw), ncol(corr_gw), format(length(corr_gw@x), big.mark=",")))

# Convertir a SFBM
sfbm_file <- file.path(tmp_dir, "corr_gwide")
sfbm <- as_SFBM(corr_gw, backingfile=sfbm_file, compact=TRUE)
cat(sprintf("  SFBM: %.1f GB\n", file.info(paste0(sfbm_file, ".sbk"))$size / 1024^3))
rm(corr_gw); gc()

# ══════════════════════════════════════════════════════════════════
# [4] Calibrar n_eff per-SNP + LDpred2-auto genome-wide
# ══════════════════════════════════════════════════════════════════
cat("\n[4] Preparando df_beta con n_eff per-SNP...\n")

df_beta <- info_snp_gw[, c("beta","beta_se","n_eff")]

# Calibrar n_eff per-SNP con sd_val
sd_ss_orig <- with(df_beta, 2 / sqrt(n_eff * beta_se^2 + beta^2))
ratio_orig <- median(sd_ss_orig / sd_val_gw)
cat(sprintf("  sd_ss/sd_val ratio ANTES: %.3f\n", ratio_orig))

df_beta$n_eff <- (2 / (sd_val_gw * df_beta$beta_se))^2 - (df_beta$beta / df_beta$beta_se)^2
df_beta$n_eff <- pmax(df_beta$n_eff, 1000)

sd_ss_new <- with(df_beta, 2 / sqrt(n_eff * beta_se^2 + beta^2))
ratio_new <- median(sd_ss_new / sd_val_gw)
cat(sprintf("  sd_ss/sd_val ratio DESPUÉS: %.3f\n", ratio_new))
cat(sprintf("  n_eff per-SNP: median=%.0f, range=[%.0f, %.0f]\n",
    median(df_beta$n_eff), min(df_beta$n_eff), max(df_beta$n_eff)))

# ── LDpred2-auto genome-wide ──────────────────────────────────────
cat("\n[5] LDpred2-auto GENOME-WIDE...\n")
cat(sprintf("    %s variantes, %d bloques\n",
    format(n_snp_total, big.mark=","), length(all_info_snp)))

n_chains <- 50L
burn_in  <- 1000L
num_iter <- 1000L

cat(sprintf("    burn_in=%d, num_iter=%d, chains=%d\n",
    burn_in, num_iter, n_chains))
cat("    allow_jump_sign=FALSE, shrink_corr=0.95, use_MLE=FALSE\n")

t0 <- Sys.time()
multi_auto <- snp_ldpred2_auto(
  corr=sfbm, df_beta=df_beta,
  h2_init=0.3,
  vec_p_init=seq_log(1e-4, 0.5, length.out=n_chains),
  burn_in=burn_in, num_iter=num_iter, report_step=100L,
  ncores=ncores, sparse=TRUE,
  allow_jump_sign=FALSE,
  shrink_corr=0.95,
  use_MLE=FALSE)
t_ldp <- as.numeric(difftime(Sys.time(), t0, units="mins"))
cat(sprintf("  LDpred2-auto: %.1f min\n", t_ldp))

h2_est <- sapply(multi_auto, function(x) x$h2_est)
p_est  <- sapply(multi_auto, function(x) x$p_est)
keep <- which(!is.na(h2_est) & h2_est > 0.001 & h2_est < 0.8 &
              !is.na(p_est)  & p_est  > 1e-6   & p_est  < 0.99)
cat(sprintf("  Cadenas válidas: %d/%d\n", length(keep), length(multi_auto)))
cat(sprintf("  h2_est range: [%.6f, %.6f]\n", min(h2_est[keep], na.rm=TRUE), max(h2_est[keep], na.rm=TRUE)))
cat(sprintf("  h2_est median: %.6f\n", median(h2_est[keep])))
cat(sprintf("  p_est median: %.2e\n", median(p_est[keep])))

# Diagnóstico completo
cat("  h2_est todas las cadenas:\n")
cat(sprintf("    %s\n", paste(round(h2_est, 6), collapse=", ")))

if (length(keep) < 3L) {
  cat("\n  ⚠ Pocas cadenas válidas. Intentando RETRY con shrink_corr=0.5...\n")
  multi_auto <- snp_ldpred2_auto(
    corr=sfbm, df_beta=df_beta,
    h2_init=0.1,
    vec_p_init=seq_log(1e-5, 0.3, length.out=50),
    burn_in=2000L, num_iter=2000L, report_step=100L,
    ncores=ncores, sparse=TRUE,
    allow_jump_sign=FALSE,
    shrink_corr=0.5,
    use_MLE=FALSE)
  h2_est <- sapply(multi_auto, function(x) x$h2_est)
  p_est  <- sapply(multi_auto, function(x) x$p_est)
  keep <- which(!is.na(h2_est) & h2_est > 0.001 & h2_est < 0.8 &
                !is.na(p_est)  & p_est  > 1e-6   & p_est  < 0.99)
  cat(sprintf("  RETRY cadenas válidas: %d/%d\n", length(keep), length(multi_auto)))
}

# ══════════════════════════════════════════════════════════════════
# [6] Extraer betas y calcular PRS
# ══════════════════════════════════════════════════════════════════
cat("\n[6] Extrayendo betas finales...\n")

if (length(keep) < 1L) {
  cat("  DIVERGENCIA TOTAL — sin cadenas válidas\n")
  saveRDS(list(status="diverged", h2_est=h2_est, p_est=p_est,
               n_snp=n_snp_total), betas_file)
} else {
  beta_mat <- sapply(multi_auto[keep], function(x) x$beta_est)
  if (is.null(dim(beta_mat))) beta_mat <- matrix(beta_mat, ncol=1)
  beta_final <- rowMeans(beta_mat, na.rm=TRUE)
  beta_final[!is.finite(beta_final)] <- 0.0

  h2_med <- median(h2_est[keep])
  p_med  <- median(p_est[keep])

  cat(sprintf("  h2=%.6f, p=%.2e\n", h2_med, p_med))
  cat(sprintf("  nonzero betas: %s / %s\n",
      format(sum(beta_final != 0), big.mark=","),
      format(length(beta_final), big.mark=",")))

  # Calcular PRS por cromosoma (necesitamos los genotipos)
  cat("\n[7] Calculando PRS por cromosoma...\n")

  # Mapear betas a cromosomas
  chr_list <- unique(info_snp_gw$chr)
  offset <- 0L
  prs_total <- NULL

  for (chr in sort(chr_list)) {
    idx_chr <- which(info_snp_gw$chr == chr)
    if (length(idx_chr) == 0) next

    rds_chr <- file.path(prev_tmp, sprintf("chr%d.rds", chr))
    obj <- snp_attach(rds_chr); G <- obj$genotypes

    beta_chr <- beta_final[idx_chr]
    ind_col_chr <- info_snp_gw$`_NUM_ID_`[idx_chr]

    prs_chr <- big_prodVec(X=G, y.col=beta_chr,
      ind.col=ind_col_chr, ncores=ncores)

    if (is.null(prs_total)) {
      prs_total <- prs_chr
      fam <- obj$fam
    } else {
      prs_total <- prs_total + prs_chr
    }
    cat(sprintf("    chr%d: %d SNPs, PRS added\n", chr, length(idx_chr)))
  }

  # Guardar
  saveRDS(list(
    status="ok",
    h2=h2_med, p=p_med,
    n_valid_chains=length(keep),
    n_snp_total=n_snp_total,
    h2_all=h2_est, p_all=p_est,
    neff_persnp_median=median(df_beta$n_eff),
    sd_ss_sd_val_ratio=ratio_new,
    t_ldpred_min=t_ldp,
    beta_final=beta_final,
    info_snp=info_snp_gw[, c("chr","pos","a0","a1")],
    prs_total=prs_total, fam=fam
  ), betas_file)

  cat(sprintf("\n ✓ LDpred2 genome-wide completado\n"))
  cat(sprintf("   h2=%.6f, p=%.2e, valid_chains=%d/%d\n",
      h2_med, p_med, length(keep), length(multi_auto)))
  cat(sprintf("   PRS calculado para %s individuos\n",
      format(length(prs_total), big.mark=",")))
}
REOF

apptainer exec \
  --bind /mnt/cephfs:/mnt/cephfs \
  --bind "$TMPDIR:/tmp,$TMPDIR:/var/tmp" \
  --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 --env MKL_NUM_THREADS=1 \
  --env PREV_TMPDIR="$PREV_TMPDIR" --env OUT_DIR="$OUT_DIR" --env GWAS_DIR="$GWAS_DIR" \
  "$SIF" /usr/bin/Rscript "$R_SCRIPT"

rm -rf "$TMPDIR"

echo ""
echo "============================================================"
echo "  LDpred2 genome-wide FIN — $(date)"
echo "============================================================"
