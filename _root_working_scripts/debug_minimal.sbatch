#!/bin/bash
#SBATCH --job-name=debug_min
#SBATCH --account=researchers
#SBATCH --qos=vip
#SBATCH --partition=normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=60G
#SBATCH --time=01:00:00
#SBATCH --output=/mnt/cephfs/orgs/home/angel.pacheco/logs/debug_min_%j.out
#SBATCH --error=/mnt/cephfs/orgs/home/angel.pacheco/logs/debug_min_%j.err

SIF="/mnt/cephfs/biocontainers/images/bigsnpr-1.10.8.sif"
TMPDIR="/mnt/cephfs/scratch/$USER/debug_min/$SLURM_JOB_ID"
mkdir -p "$TMPDIR"
export TMPDIR TMP="$TMPDIR" TEMP="$TMPDIR"

R_SCRIPT="$TMPDIR/debug.R"
cat << 'RSCRIPT' > "$R_SCRIPT"
library(bigsnpr)

# Versión de bigsnpr
cat(sprintf("bigsnpr version: %s\n", packageVersion("bigsnpr")))
cat(sprintf("bigsparser version: %s\n", packageVersion("bigsparser")))

# Crear datos simulados que SABEMOS que funcionan
set.seed(1)
N <- 5000; M <- 1000
G_sim <- matrix(rbinom(N*M, 2, runif(M, 0.1, 0.5)), N, M)
h2_true <- 0.3
causal <- sample(M, 50)
beta_true <- rep(0, M)
beta_true[causal] <- rnorm(50, 0, sqrt(h2_true/50))
y <- G_sim %*% beta_true + rnorm(N, 0, sqrt(1-h2_true))

# GWAS simple
beta_hat <- sapply(1:M, function(j) coef(lm(y ~ G_sim[,j]))[2])
se_hat <- sapply(1:M, function(j) summary(lm(y ~ G_sim[,j]))$coefficients[2,2])

cat(sprintf("\n=== SIMULACIÓN (debería funcionar) ===\n"))
cat(sprintf("beta_hat: median=%.4f, sd=%.4f\n", median(beta_hat), sd(beta_hat)))

# LD con bigsnpr
tmp <- tempfile()
fake_bk <- snp_fake(N, M)
G_fake <- fake_bk$genotypes
G_fake[] <- G_sim
corr_sim <- snp_cor(G_fake, ncores = 1, thr_r2 = 0.02)
corr_sfbm <- as_SFBM(corr_sim, compact = TRUE)

df_beta_sim <- data.frame(beta = beta_hat, beta_se = se_hat, n_eff = N)

auto_sim <- snp_ldpred2_auto(corr_sfbm, df_beta_sim, h2_init = 0.3,
    vec_p_init = seq_log(1e-4, 0.5, 3), burn_in = 200, num_iter = 200,
    ncores = 1, sparse = TRUE)
h2_sim <- sapply(auto_sim, function(x) x$h2_est)
cat(sprintf("h2_est SIMULADO: %s\n", paste(round(h2_sim, 4), collapse=", ")))

# Ahora con datos REALES chr22 HapMap3
cat(sprintf("\n=== DATOS REALES CHR22 ===\n"))
library(data.table)
tmp_dir <- Sys.getenv("TMPDIR")

bed <- "/mnt/cephfs/hot_nvme/mcps/imputed-topmed/plink_files/maximally_unrelated/mcps-freeze150k_qcd_chr22_ivs.bed"
bk <- file.path(tmp_dir, "chr22")
snp_readBed(bed, backingfile = bk)
obj <- snp_attach(paste0(bk, ".rds"))
G <- obj$genotypes

gwas_dir <- "/mnt/cephfs/orgs/home/angel.pacheco/gwas_miami"
gf <- fread(file.path(gwas_dir, "combined-gwas-results-all-ctrl-females.txt.gz"),
    select=c("CHROM","GENPOS","ALLELE0","ALLELE1","A1FREQ","INFO","N","BETA","SE"), showProgress=FALSE)
gm <- fread(file.path(gwas_dir, "combined-gwas-results-all-ctrl-males.txt.gz"),
    select=c("CHROM","GENPOS","ALLELE0","ALLELE1","A1FREQ","INFO","N","BETA","SE"), showProgress=FALSE)
gf <- gf[CHROM == 22]; gm <- gm[CHROM == 22]
gwas <- merge(gf, gm, by=c("CHROM","GENPOS","ALLELE0","ALLELE1"), suffixes=c(".f",".m"))
rm(gf,gm); gc()

w_f <- 1/gwas$SE.f^2; w_m <- 1/gwas$SE.m^2
neff_f <- 4*17300*54196/(17300+54196)
neff_m <- 4*8522*26733/(8522+26733)
K <- (17300+8522)/(17300+54196+8522+26733)
n_eff <- (neff_f+neff_m)*K*(1-K)

sumstats <- data.frame(chr=22L, pos=gwas$GENPOS,
    a1=gwas$ALLELE1, a0=gwas$ALLELE0,
    beta=(gwas$BETA.f*w_f+gwas$BETA.m*w_m)/(w_f+w_m),
    beta_se=sqrt(1/(w_f+w_m)), n_eff=n_eff, stringsAsFactors=FALSE)
sumstats <- sumstats[!is.na(sumstats$beta) & sumstats$beta_se > 0, ]

# Filtrar HapMap3
hm3 <- readRDS("/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes/ldpred2/map_hm3_ldpred2.rds")
sumstats <- sumstats[paste(sumstats$chr,sumstats$pos,sep=":") %in% paste(hm3$chr,hm3$pos_hg38,sep=":"), ]
cat(sprintf("HapMap3 chr22: %d variantes\n", nrow(sumstats)))

map_df <- data.frame(chr=obj$map$chromosome, pos=obj$map$physical.pos,
    a1=obj$map$allele1, a0=obj$map$allele2, stringsAsFactors=FALSE)
info_snp <- snp_match(sumstats, map_df, strand_flip=FALSE, join_by_pos=TRUE, match.min.prop=0.01)
cat(sprintf("Matched: %d\n", nrow(info_snp)))

# Verificar que el QC interno de LDpred2 pasa
sd_ss <- with(info_snp, 1/sqrt(n_eff * beta_se^2 + beta^2))
sd_geno <- big_scale()(G, ind.col=info_snp$`_NUM_ID_`)$scale
ratio <- sd_ss / sd_geno
cat(sprintf("ratio sd_ss/sd_geno: median=%.4f\n", median(ratio)))

# LD
corr0 <- snp_cor(G, ind.col=info_snp$`_NUM_ID_`, infos.pos=info_snp$pos/1000,
    size=3000, ncores=4, thr_r2=0.02)
cat(sprintf("LD nnz: %s\n", format(length(corr0@x), big.mark=",")))
corr <- as_SFBM(corr0, compact=TRUE)
rm(corr0); gc()

df_beta <- info_snp[, c("beta","beta_se","n_eff")]

# Test con 1 sola chain y report_step=10 para ver progreso
cat("\n=== UNA CHAIN, report_step=10 ===\n")
auto1 <- snp_ldpred2_auto(corr, df_beta, h2_init=0.01,
    vec_p_init=0.1, burn_in=500, num_iter=500, report_step=10,
    ncores=1, sparse=TRUE)
cat(sprintf("h2_est: %s\n", auto1[[1]]$h2_est))
cat(sprintf("p_est: %s\n", auto1[[1]]$p_est))

# Ver path del h2
if (!is.na(auto1[[1]]$h2_est)) {
    cat("h2 path (últimos 20):\n")
    h2_path <- auto1[[1]]$path_p_est
    print(tail(h2_path, 20))
}
RSCRIPT

apptainer exec --bind /mnt/cephfs:/mnt/cephfs --bind "$TMPDIR:/tmp,$TMPDIR:/var/tmp" \
    "$SIF" /opt/miniforge3/envs/r-bio/bin/Rscript "$R_SCRIPT"
rm -rf "$TMPDIR"
