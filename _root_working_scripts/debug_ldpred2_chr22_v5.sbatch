#!/bin/bash
#SBATCH --job-name=debug_ldp2
#SBATCH --account=researchers
#SBATCH --qos=vip
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=120G
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
library(methods)

ncores <- 24
options(bigstatsr.check.parallel.blas = FALSE)
options(bigstatsr.ncores.max = 24)

# ================================================================
# safe_as_SFBM: El fix maestro
# ================================================================
safe_as_SFBM <- function(corr0, backingfile, compact = TRUE, eps = 1e-10) {
    stopifnot(is(corr0, "dsCMatrix"))
    col_nnz <- diff(corr0@p)
    empty <- which(col_nnz == 0L)

    if (length(empty) == 0L) {
        cat("      [SFBM] Sin columnas vacías\n")
        return(as_SFBM(corr0, backingfile = backingfile, compact = compact))
    }

    n      <- corr0@Dim[1]
    n_add  <- length(empty)
    old_p  <- corr0@p
    old_nnz <- length(corr0@i)
    new_nnz <- old_nnz + n_add

    cat(sprintf("      [SFBM] %d columnas vacías de %s → reparando slots\n", n_add, format(n, big.mark=",")))

    cum_add <- cumsum(col_nnz == 0L)
    new_p   <- old_p + c(0L, cum_add)

    new_i <- integer(new_nnz)
    new_x <- numeric(new_nnz)

    boundaries <- c(0L, empty)
    n_blocks   <- length(boundaries)

    for (b in seq_len(n_blocks)) {
        block_first <- boundaries[b] + 1L
        block_last  <- if (b < n_blocks) boundaries[b + 1L] - 1L else n
        shift       <- b - 1L

        if (block_last >= block_first) {
            old_start <- old_p[block_first] + 1L
            old_end   <- old_p[block_last + 1L]
            if (old_end >= old_start) {
                cnt <- old_end - old_start + 1L
                ns  <- old_start + shift
                new_i[ns:(ns + cnt - 1L)] <- corr0@i[old_start:old_end]
                new_x[ns:(ns + cnt - 1L)] <- corr0@x[old_start:old_end]
            }
        }

        if (b <= n_add) {
            ec <- empty[b]
            dp <- new_p[ec] + 1L
            new_i[dp] <- ec - 1L
            new_x[dp] <- eps
        }
    }
    cat(sprintf("      [SFBM] Reparado: %s → %s entries\n", format(old_nnz, big.mark=","), format(new_nnz, big.mark=",")))

    corr_fixed <- new("dsCMatrix", i = new_i, p = new_p, x = new_x, Dim = corr0@Dim, Dimnames = corr0@Dimnames, uplo = corr0@uplo, factors = list())
    as_SFBM(corr_fixed, backingfile = backingfile, compact = compact)
}

neff_f <- 4 * 17300 * 54196 / (17300 + 54196)
neff_m <- 4 * 8522 * 26733 / ( 8522 + 26733)
K_prev <- (17300 + 8522) / (17300 + 54196 + 8522 + 26733)
neff_ldpred2 <- (neff_f + neff_m) * K_prev * (1 - K_prev)
n_eff_val <- neff_ldpred2
cat(sprintf("n_eff combinado: %.1f\n", n_eff_val))

tmp_dir <- Sys.getenv("TMPDIR")

# Cargar solo chr22
bed <- "/mnt/cephfs/hot_nvme/mcps/imputed-topmed/plink_files/maximally_unrelated/mcps-freeze150k_qcd_chr22_ivs.bed"
bk <- file.path(tmp_dir, "debug_chr22")
snp_readBed(bed, backingfile = bk)
obj <- snp_attach(paste0(bk, ".rds"))
G <- obj$genotypes
cat(sprintf("Genotipos: %d individuos x %d variantes\n", nrow(G), ncol(G)))

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
# Filtrar a HapMap3
hm3_file <- "/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes/ldpred2/map_hm3_ldpred2.rds"
map_hm3 <- readRDS(hm3_file)
key_hm3 <- paste(map_hm3$chr, map_hm3$pos_hg38, sep = ":")
key_gwas <- paste(sumstats$chr, sumstats$pos, sep = ":")
sumstats <- sumstats[key_gwas %in% key_hm3, ]
cat(sprintf("Tras filtro HapMap3: %d variantes\n", nrow(sumstats)))
cat(sprintf("GWAS chr22 tras QC: %d variantes\n", nrow(sumstats)))

# Match
map_df <- data.frame(chr = obj$map$chromosome, pos = obj$map$physical.pos,
    a1 = obj$map$allele1, a0 = obj$map$allele2, stringsAsFactors = FALSE)
info_snp <- snp_match(sumstats, map_df, strand_flip = FALSE, join_by_pos = TRUE, match.min.prop = 0.01)
cat(sprintf("Matched: %d variantes\n", nrow(info_snp)))

# === DIAGNÓSTICOS ===
cat("\n========== DIAGNÓSTICOS ==========\n")
cat(sprintf("beta: median=%.6f, mean=%.6f, sd=%.6f, range=[%.4f, %.4f]\n",
    median(info_snp$beta), mean(info_snp$beta), sd(info_snp$beta), min(info_snp$beta), max(info_snp$beta)))
cat(sprintf("beta_se: median=%.6f, mean=%.6f, range=[%.4f, %.4f]\n",
    median(info_snp$beta_se), mean(info_snp$beta_se), min(info_snp$beta_se), max(info_snp$beta_se)))
cat(sprintf("n_eff: %.0f (constante)\n", n_eff_val))

sd_val <- 1 / sqrt(n_eff_val * info_snp$beta_se^2 + info_snp$beta^2)
cat(sprintf("\nsd_val (GWAS): median=%.6f, mean=%.6f, range=[%.4f, %.4f]\n",
    median(sd_val), mean(sd_val), min(sd_val), max(sd_val)))

ind_col <- info_snp$`_NUM_ID_`
sd_geno <- big_scale()(G, ind.col = ind_col)$scale / sqrt(nrow(G))
cat(sprintf("sd_geno: median=%.6f, mean=%.6f, range=[%.4f, %.4f]\n",
    median(sd_geno), mean(sd_geno), min(sd_geno), max(sd_geno)))

ratio <- sd_val / sd_geno
cat(sprintf("ratio sd_val/sd_geno: median=%.4f, mean=%.4f, range=[%.4f, %.4f]\n",
    median(ratio), mean(ratio), min(ratio), max(ratio)))
cat(sprintf("  ratio < 0.5: %d (%.1f%%)\n", sum(ratio < 0.5), 100*mean(ratio < 0.5)))
cat(sprintf("  ratio > 2.0: %d (%.1f%%)\n", sum(ratio > 2.0), 100*mean(ratio > 2.0)))

chisq <- (info_snp$beta / info_snp$beta_se)^2
cat(sprintf("\nChi-sq: median=%.4f, mean=%.4f, max=%.1f\n", median(chisq), mean(chisq), max(chisq)))

# === LD matrix ===
cat("\n========== LD MATRIX ==========\n")
cat(sprintf("Calculando LD matrix (%d vars)...\n", length(ind_col)))
corr0 <- snp_cor(Gna = G, ind.col = ind_col, infos.pos = info_snp$pos / 1000, size = 3000, ncores = ncores, thr_r2 = 0.02)
cat(sprintf("LD matrix: %d x %d, nnz=%s, clase=%s\n", nrow(corr0), ncol(corr0), format(length(corr0@x), big.mark=","), class(corr0)))

sfbm_file <- file.path(tmp_dir, "debug_corr_chr22")
corr <- safe_as_SFBM(corr0, backingfile = sfbm_file, compact = TRUE)
rm(corr0); gc()

df_beta <- info_snp[, c("beta", "beta_se", "n_eff")]

# === TEST 1: h2_init = 0.01 (original) ===
cat("\n========== TEST 1: h2_init=0.01 (3 chains) ==========\n")
multi_auto <- snp_ldpred2_auto(
    corr = corr, df_beta = df_beta, h2_init = 0.01,
    vec_p_init = seq_log(1e-4, 0.5, length.out = 3),
    burn_in = 200, num_iter = 200, report_step = 50,
    ncores = 3, sparse = TRUE
)
h2_est <- sapply(multi_auto, function(x) x$h2_est)
cat(sprintf("h2_est: %s\n", paste(round(h2_est, 6), collapse=", ")))

# === TEST 2: h2_init = 0.2 ===
cat("\n========== TEST 2: h2_init=0.2 (3 chains) ==========\n")
multi_auto2 <- snp_ldpred2_auto(
    corr = corr, df_beta = df_beta, h2_init = 0.2,
    vec_p_init = seq_log(1e-4, 0.5, length.out = 3),
    burn_in = 200, num_iter = 200, report_step = 50,
    ncores = 3, sparse = TRUE
)
h2_est2 <- sapply(multi_auto2, function(x) x$h2_est)
cat(sprintf("h2_est: %s\n", paste(round(h2_est2, 6), collapse=", ")))

# === TEST 3: allow_jump_sign = FALSE ===
cat("\n========== TEST 3: h2_init=0.01, allow_jump_sign=FALSE ==========\n")
multi_auto3 <- snp_ldpred2_auto(
    corr = corr, df_beta = df_beta, h2_init = 0.01,
    vec_p_init = seq_log(1e-4, 0.5, length.out = 3),
    burn_in = 200, num_iter = 200, report_step = 50,
    ncores = 3, sparse = TRUE, allow_jump_sign = FALSE
)
h2_est3 <- sapply(multi_auto3, function(x) x$h2_est)
cat(sprintf("h2_est: %s\n", paste(round(h2_est3, 6), collapse=", ")))

cat("\n========== DIAGNÓSTICO COMPLETO ==========\n")
REOF

apptainer exec --bind /mnt/cephfs:/mnt/cephfs --bind "$TMPDIR:/tmp,$TMPDIR:/var/tmp" \
    "$SIF" /opt/miniforge3/envs/r-bio/bin/Rscript "$R_SCRIPT"

rm -rf "$TMPDIR"
