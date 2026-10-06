#!/bin/bash
#SBATCH --job-name=debug_ldp2
#SBATCH --account=researchers
#SBATCH --qos=vip
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=60G
#SBATCH --output=/mnt/cephfs/orgs/home/angel.pacheco/logs/debug_ldp2_%j.out
#SBATCH --error=/mnt/cephfs/orgs/home/angel.pacheco/logs/debug_ldp2_%j.err

set -euo pipefail

SIF="/mnt/cephfs/biocontainers/images/jupyter-biotools-1.4.sif"

TMPDIR="/mnt/cephfs/scratch/$USER/debug_ldp2/$SLURM_JOB_ID"
mkdir -p "$TMPDIR"
export TMPDIR TMP="$TMPDIR" TEMP="$TMPDIR"

R_SCRIPT="$TMPDIR/debug.R"
cat << 'REOF' > "$R_SCRIPT"
library(bigsnpr)
library(data.table)
library(Matrix)

# FIX: CERO paralelismo para evitar el assert_cores
ncores <- 1
options(bigstatsr.check.parallel.blas = FALSE)
options(bigstatsr.ncores.max = 1)

neff_f <- 4 * 17300 * 54196 / (17300 + 54196)
neff_m <- 4 * 8522 * 26733 / ( 8522 + 26733)
n_eff_val <- neff_f + neff_m

tmp_dir <- Sys.getenv("TMPDIR")

# Cargar solo chr22
bed <- "/mnt/cephfs/hot_nvme/mcps/imputed-topmed/plink_files/maximally_unrelated/mcps-freeze150k_qcd_chr22_ivs.bed"
bk <- file.path(tmp_dir, "debug_chr22")
snp_readBed(bed, backingfile = bk)
obj <- snp_attach(paste0(bk, ".rds"))
G <- obj$genotypes

# Cargar GWAS
gwas_dir <- "/mnt/cephfs/orgs/home/angel.pacheco/gwas_miami"
gf <- fread(file.path(gwas_dir, "combined-gwas-results-all-ctrl-females.txt.gz"),
    select = c("CHROM","GENPOS","ID","ALLELE0","ALLELE1","A1FREQ","INFO","N","BETA","SE","P"), showProgress = FALSE)
gm <- fread(file.path(gwas_dir, "combined-gwas-results-all-ctrl-males.txt.gz"),
    select = c("CHROM","GENPOS","ID","ALLELE0","ALLELE1","A1FREQ","INFO","N","BETA","SE","P"), showProgress = FALSE)
gwas <- merge(gf, gm, by = "ID", suffixes = c(".f", ".m"))
rm(gf, gm); gc()

w_f <- 1/gwas$SE.f^2; w_m <- 1/gwas$SE.m^2
sumstats <- data.frame(
    chr = gwas$CHROM.f, pos = gwas$GENPOS.f,
    a1 = gwas$ALLELE1.f, a0 = gwas$ALLELE0.f,
    freq = (gwas$A1FREQ.f * gwas$N.f + gwas$A1FREQ.m * gwas$N.m) / (gwas$N.f + gwas$N.m),
    info = pmin(gwas$INFO.f, gwas$INFO.m),
    beta = (gwas$BETA.f * w_f + gwas$BETA.m * w_m) / (w_f + w_m),
    beta_se = sqrt(1 / (w_f + w_m)),
    n_eff = n_eff_val, stringsAsFactors = FALSE
)
rm(gwas); gc()

sumstats <- sumstats[sumstats$info >= 0.8 & sumstats$freq >= 0.01 & sumstats$freq <= 0.99 &
    !is.na(sumstats$beta) & !is.na(sumstats$beta_se) & sumstats$beta_se > 0, ]
sumstats$chr <- as.integer(as.character(sumstats$chr))
sumstats <- sumstats[sumstats$chr == 22, ]

# Match
map_df <- data.frame(chr = obj$map$chromosome, pos = obj$map$physical.pos,
    a1 = obj$map$allele1, a0 = obj$map$allele2, stringsAsFactors = FALSE)
info_snp <- snp_match(sumstats, map_df, strand_flip = FALSE, join_by_pos = TRUE, match.min.prop = 0.01)

ind_col <- info_snp$`_NUM_ID_`

cat("\n========== LD MATRIX ==========\n")
corr0 <- snp_cor(Gna = G, ind.col = ind_col,
    infos.pos = info_snp$pos / 1000, size = 3000, ncores = ncores, thr_r2 = 0.01) 
corr0 <- corr0 + Matrix::Diagonal(nrow(corr0), 1e-10)

sfbm_file <- file.path(tmp_dir, "debug_corr_chr22")
corr <- as_SFBM(corr0, backingfile = sfbm_file, compact = TRUE)
rm(corr0); gc()

df_beta <- info_snp[, c("beta", "beta_se", "n_eff")]

# === TEST 1: h2_init = 0.01 (original) ===
cat("\n========== TEST 1: h2_init=0.01 (3 chains) ==========\n")
multi_auto <- snp_ldpred2_auto(
    corr = corr, df_beta = df_beta, h2_init = 0.01,
    vec_p_init = seq_log(1e-4, 0.5, length.out = 3),
    burn_in = 200, num_iter = 200, report_step = 50,
    ncores = 1, sparse = TRUE
)
h2_est <- sapply(multi_auto, function(x) x$h2_est)
cat(sprintf("h2_est: %s\n", paste(round(h2_est, 6), collapse=", ")))

# === TEST 2: h2_init = 0.2 ===
cat("\n========== TEST 2: h2_init=0.2 (3 chains) ==========\n")
multi_auto2 <- snp_ldpred2_auto(
    corr = corr, df_beta = df_beta, h2_init = 0.2,
    vec_p_init = seq_log(1e-4, 0.5, length.out = 3),
    burn_in = 200, num_iter = 200, report_step = 50,
    ncores = 1, sparse = TRUE
)
h2_est2 <- sapply(multi_auto2, function(x) x$h2_est)
cat(sprintf("h2_est: %s\n", paste(round(h2_est2, 6), collapse=", ")))

# === TEST 3: allow_jump_sign = FALSE ===
cat("\n========== TEST 3: h2_init=0.01, allow_jump_sign=FALSE ==========\n")
multi_auto3 <- snp_ldpred2_auto(
    corr = corr, df_beta = df_beta, h2_init = 0.01,
    vec_p_init = seq_log(1e-4, 0.5, length.out = 3),
    burn_in = 200, num_iter = 200, report_step = 50,
    ncores = 1, sparse = TRUE, allow_jump_sign = FALSE
)
h2_est3 <- sapply(multi_auto3, function(x) x$h2_est)
cat(sprintf("h2_est: %s\n", paste(round(h2_est3, 6), collapse=", ")))

cat("\n========== DIAGNÓSTICO COMPLETO ==========\n")
REOF

apptainer exec --bind /mnt/cephfs:/mnt/cephfs --bind "$TMPDIR:/tmp,$TMPDIR:/var/tmp" \
    "$SIF" /opt/miniforge3/envs/r-bio/bin/Rscript "$R_SCRIPT"

rm -rf "$TMPDIR"
