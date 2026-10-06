#!/bin/bash
#SBATCH --job-name=sbrc_run
#SBATCH --account=researchers
#SBATCH --qos=vip
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=120G
#SBATCH --time=UNLIMITED
#SBATCH --output=/mnt/cephfs/orgs/home/angel.pacheco/logs/sbrc_run_%j.out
#SBATCH --error=/mnt/cephfs/orgs/home/angel.pacheco/logs/sbrc_run_%j.err
#SBATCH --mail-user=angelpachecosaavedra@gmail.com
#SBATCH --mail-type=BEGIN,END,FAIL

set -euo pipefail

SBRC_SIF="/mnt/cephfs/biocontainers/images/sbayesrc-0.2.6.sif"
SIF14="/mnt/cephfs/biocontainers/images/jupyter-biotools-1.4.sif"

PLINK_DIR="/mnt/cephfs/hot_nvme/mcps/imputed-topmed/plink_files/maximally_unrelated"
OUT_BASE="/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes/sbayesrc"
LD_DIR="${OUT_BASE}/ldm_79k"

# ¡NUEVO ARCHIVO CON LIFTOVER Y MATCH PERFECTO!
MA_FILE="${OUT_BASE}/T2DGGI_ALL_hg38_matched.ma" 

mkdir -p "${OUT_BASE}/results"

TMPDIR="/mnt/cephfs/scratch/$USER/sbrc_run_tmp/$SLURM_JOB_ID"
mkdir -p "$TMPDIR"
chmod 700 "$TMPDIR"
export TMPDIR TMP="$TMPDIR" TEMP="$TMPDIR"

echo "============================================================"
echo "  PIPELINE A — ETAPA 2: SBayesR (LD + MCMC)"
echo "  Nodo: $SLURMD_NODENAME"
echo "  RAM: $(free -g | awk 'NR==2{print $2}')G"
echo "  $(date)"
echo "============================================================"

# Verificar que el .ma existe
if [ ! -f "$MA_FILE" ]; then
    echo " No se encontró $MA_FILE"
    exit 1
fi

echo "  Summary stats: $MA_FILE ($(wc -l < "$MA_FILE") líneas)"

R_SCRIPT="$TMPDIR/run_sbayesr.R"
cat << 'REOF' > "$R_SCRIPT"
library(SBayesRC)

ma_file   <- "/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes/sbayesrc/T2DGGI_ALL_hg38_matched.ma"
geno_dir  <- "/mnt/cephfs/hot_nvme/mcps/imputed-topmed/plink_files/maximally_unrelated/mcps-freeze150k_qcd_chr{CHR}_ivs"
out_dir   <- "/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes/sbayesrc/results"
ld_dir    <- "/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes/sbayesrc/ldm_79k"
threads   <- 24

cat("\n=== PASO 1: Generando bloques LD ===\n")
tryCatch({
    SBayesRC::LDstep1(mafile = ma_file, genoPrefix = geno_dir, outDir = ld_dir, genoCHR = "1-22", log2file = TRUE)
    cat(" Paso 1 completado\n")
}, error = function(e) { cat(sprintf("❌ Error en LDstep1: %s\n", e$message)) })

cat("\n=== PASO 2: Generando matrices LD por bloque ===\n")
ld_script <- file.path(ld_dir, "ld.sh")
if (file.exists(ld_script)) {
    system(sprintf("bash %s", ld_script))
    cat(" Paso 2 completado\n")
} else {
    SBayesRC::LDstep2(outDir = ld_dir)
}

cat("\n=== PASO 3: Eigen-decomposición ===\n")
tryCatch({
    SBayesRC::LDstep3(outDir = ld_dir, threads = threads)
    cat(" Paso 3 completado\n")
}, error = function(e) { cat(sprintf("❌ Error en LDstep3: %s\n", e$message)) })

cat("\n=== PASO 4: SBayesR MCMC ===\n")
tryCatch({
    # Usamos annot = NULL para forzar el modo SBayesR clásico
    SBayesRC::sbayesrc(mafile = ma_file, LDdir = ld_dir, outDir = out_dir, annot = NULL, threads = threads)
    cat("\n SBayesR MCMC completado\n")
}, error = function(e) { cat(sprintf("\n❌ Error en SBayesR: %s\n", e$message)) })

cat("\n=== PASO 5: Exportando betas para scoring ===\n")
snp_res <- file.path(out_dir, "sbayesrc.snpRes")
if (!file.exists(snp_res)) {
    snp_files <- list.files(out_dir, pattern = "\\.snpRes$", full.names = TRUE)
    if (length(snp_files) > 0) snp_res <- snp_files[1]
}

if (file.exists(snp_res)) {
    df_snp <- read.table(snp_res, header = TRUE)
    df_nonzero <- df_snp[df_snp$A1Effect != 0, ]
    out_weights <- file.path(dirname(out_dir), "sbayesrc_weights_for_scoring.tsv")
    write.table(df_snp[, c("Id", "A1", "A2", "A1Effect", "A1Frq")], out_weights, sep = "\t", row.names = FALSE, quote = FALSE, col.names = c("SNP", "A1", "A2", "BETA", "FREQ"))
    cat(sprintf("\n Pesos para scoring: %s\n", out_weights))
}
cat(sprintf("\n=== SBayesR completado: %s ===\n", Sys.time()))
REOF

if [ -f "$SBRC_SIF" ]; then
    apptainer exec --bind /mnt/cephfs:/mnt/cephfs --bind "$TMPDIR:/tmp,$TMPDIR:/var/tmp" "$SBRC_SIF" Rscript "$R_SCRIPT"
else
    apptainer exec --bind /mnt/cephfs:/mnt/cephfs --bind "$TMPDIR:/tmp,$TMPDIR:/var/tmp" "$SIF14" /opt/miniforge3/envs/r-bio/bin/Rscript "$R_SCRIPT"
fi

rm -rf "$TMPDIR"
