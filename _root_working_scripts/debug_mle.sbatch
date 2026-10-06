#!/bin/bash
#SBATCH --job-name=debug_mle
#SBATCH --account=researchers
#SBATCH --qos=vip
#SBATCH --partition=normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=30G
#SBATCH --time=00:30:00
#SBATCH --output=/mnt/cephfs/orgs/home/angel.pacheco/logs/debug_mle_%j.out
#SBATCH --error=/mnt/cephfs/orgs/home/angel.pacheco/logs/debug_mle_%j.err

SIF="/mnt/cephfs/biocontainers/images/jupyter-biotools-1.4.sif"
TMPDIR="/mnt/cephfs/scratch/$USER/debug_mle/$SLURM_JOB_ID"
mkdir -p "$TMPDIR"
export TMPDIR TMP="$TMPDIR" TEMP="$TMPDIR"

apptainer exec --bind /mnt/cephfs:/mnt/cephfs --bind "$TMPDIR:/tmp,$TMPDIR:/var/tmp" \
    "$SIF" /opt/miniforge3/envs/r-bio/bin/Rscript -e '
library(bigsnpr)
set.seed(1)
N <- 5000; M <- 1000
G_sim <- matrix(rbinom(N*M, 2, runif(M, 0.1, 0.5)), N, M)
h2_true <- 0.3
causal <- sample(M, 50)
beta_true <- rep(0, M)
beta_true[causal] <- rnorm(50, 0, sqrt(h2_true/50))
y <- G_sim %*% beta_true + rnorm(N, 0, sqrt(1-h2_true))
beta_hat <- sapply(1:M, function(j) coef(lm(y ~ G_sim[,j]))[2])
se_hat <- sapply(1:M, function(j) summary(lm(y ~ G_sim[,j]))$coefficients[2,2])

fake <- snp_fake(N, M)
fake$genotypes[] <- G_sim
corr <- snp_cor(fake$genotypes, ncores=1, thr_r2=0.02)
corr_sfbm <- as_SFBM(corr, compact=TRUE)
df_beta <- data.frame(beta=beta_hat, beta_se=se_hat, n_eff=N)

cat("=== TEST 1: use_MLE=TRUE (default) ===\n")
a1 <- snp_ldpred2_auto(corr_sfbm, df_beta, h2_init=0.3,
    vec_p_init=0.1, burn_in=500, num_iter=500, ncores=1, use_MLE=TRUE)
cat(sprintf("h2_est: %s\n", a1[[1]]$h2_est))

cat("\n=== TEST 2: use_MLE=FALSE ===\n")
a2 <- snp_ldpred2_auto(corr_sfbm, df_beta, h2_init=0.3,
    vec_p_init=0.1, burn_in=500, num_iter=500, ncores=1, use_MLE=FALSE)
cat(sprintf("h2_est: %s\n", a2[[1]]$h2_est))

cat("\n=== TEST 3: use_MLE=FALSE, shrink_corr=0.95 ===\n")
a3 <- snp_ldpred2_auto(corr_sfbm, df_beta, h2_init=0.3,
    vec_p_init=0.1, burn_in=500, num_iter=500, ncores=1, use_MLE=FALSE, shrink_corr=0.95)
cat(sprintf("h2_est: %s\n", a3[[1]]$h2_est))
'
rm -rf "$TMPDIR"
