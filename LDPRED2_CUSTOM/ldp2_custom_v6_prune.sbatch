#!/bin/bash
#SBATCH --job-name=ldp2_custom_prune
#SBATCH --account=researchers
#SBATCH --qos=vip
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=32G
#SBATCH --output=/mnt/cephfs/orgs/home/angel.pacheco/logs/ldp2_custom_prune_%j.out
#SBATCH --error=/mnt/cephfs/orgs/home/angel.pacheco/logs/ldp2_custom_prune_%j.err
#SBATCH --mail-user=angelpachecosaavedra@gmail.com
#SBATCH --mail-type=END,FAIL

set -euo pipefail

# ================================================================
# PIPELINE C v6 — PASO 1: Pruning MCPS (preparatorio)
#
# Genera:
#   - sumstats_pruned.rds  (sumstats ya filtrados por LD pruning en MCPS)
#   - pruned_ids.rds       (IDs que sobrevivieron al pruning)
#
# Estos outputs alimentan los workers paralelos por cromosoma.
# Este paso es rápido (~30 min) y corre una sola vez.
# ================================================================

SIF="/mnt/cephfs/biocontainers/images/jupyter-biotools-1.4.sif"
PLINK_DIR="/mnt/cephfs/hot_nvme/mcps/imputed-topmed/plink_files/maximally_unrelated"
GWAS_DIR="/mnt/cephfs/orgs/home/angel.pacheco/gwas_miami"
OUT_DIR="/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes/ldpred2_custom"

mkdir -p "$OUT_DIR"

TMPDIR="/mnt/cephfs/scratch/$USER/ldp2_custom_prune_tmp/$SLURM_JOB_ID"
mkdir -p "$TMPDIR"; chmod 700 "$TMPDIR"
export TMPDIR TMP="$TMPDIR" TEMP="$TMPDIR"

echo "============================================================"
echo "  PIPELINE C v6 — PASO 1: Pruning MCPS"
echo "  Nodo: $SLURMD_NODENAME"
echo "  $(date)"
echo "============================================================"

R_SCRIPT="$TMPDIR/prune_mcps.R"
cat << REOF > "$R_SCRIPT"
Sys.setenv(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
           MKL_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1")
if (requireNamespace("RhpcBLASctl", quietly=TRUE)) {
  RhpcBLASctl::blas_set_num_threads(1L)
  RhpcBLASctl::omp_set_num_threads(1L)
}
suppressMessages({
  library(data.table)
})

ncores    <- 24L
out_dir   <- "$OUT_DIR"
gwas_dir  <- "$GWAS_DIR"
plink_dir <- "$PLINK_DIR"
tmp_dir   <- Sys.getenv("TMPDIR")

# ================================================================
# [1] GWAS meta-análisis
# ================================================================
cat("\n[1] Cargando GWAS Meta-análisis...\n")

gf <- fread(file.path(gwas_dir, "combined-gwas-results-all-ctrl-females.txt.gz"),
            select=c("CHROM","GENPOS","ID","ALLELE0","ALLELE1","A1FREQ","INFO","N","BETA","SE","P"),
            showProgress=FALSE)
gm <- fread(file.path(gwas_dir, "combined-gwas-results-all-ctrl-males.txt.gz"),
            select=c("CHROM","GENPOS","ID","ALLELE0","ALLELE1","A1FREQ","INFO","N","BETA","SE","P"),
            showProgress=FALSE)
gwas <- merge(gf, gm, by="ID", suffixes=c(".f", ".m"))
rm(gf, gm); gc()

w_f <- 1/gwas\$SE.f^2; w_m <- 1/gwas\$SE.m^2

# n_eff binario con corrección log-odds → escala lineal
neff_f <- 4*17300*54196/(17300+54196)
neff_m <- 4*8522*26733/(8522+26733)
neff_ldpred2 <- neff_f + neff_m
#K_prev <- (17300+8522)/(17300+54196+8522+26733)
#neff_ldpred2 <- (neff_f + neff_m) * K_prev * (1 - K_prev)

sumstats <- data.frame(
  chr = gwas\$CHROM.f, pos = gwas\$GENPOS.f,
  a1 = gwas\$ALLELE1.f, a0 = gwas\$ALLELE0.f,
  freq = (gwas\$A1FREQ.f*gwas\$N.f + gwas\$A1FREQ.m*gwas\$N.m)/(gwas\$N.f+gwas\$N.m),
  info = pmin(gwas\$INFO.f, gwas\$INFO.m),
  beta = (gwas\$BETA.f*w_f + gwas\$BETA.m*w_m)/(w_f+w_m),
  beta_se = sqrt(1/(w_f+w_m)), n_eff = neff_ldpred2,
  stringsAsFactors=FALSE
)
sumstats\$p <- 2*pnorm(-abs(sumstats\$beta/sumstats\$beta_se))
rm(gwas); gc()

sumstats <- sumstats[sumstats\$info>=0.8 & sumstats\$freq>=0.01 & sumstats\$freq<=0.99 &
                     !is.na(sumstats\$beta) & !is.na(sumstats\$beta_se) & sumstats\$beta_se>0, ]
sumstats\$chr_orig <- as.character(sumstats\$chr)
sumstats\$chr <- as.character(sumstats\$chr)
sumstats\$chr[sumstats\$chr %in% c("X","chrX")] <- "23"
sumstats\$chr <- as.integer(sumstats\$chr)

cat(sprintf("  GWAS tras QC: %s variantes (chrX: %s)\n",
    format(nrow(sumstats), big.mark=","),
    format(sum(sumstats\$chr == 23), big.mark=",")))

# ================================================================
# [2] Generar IDs candidatos para extract
# ================================================================
cat("\n[2] Generando IDs candidatos...\n")
gwas_ids_file <- file.path(tmp_dir, "gwas_ids_for_pruning.txt")

is_autosome <- sumstats\$chr != 23
is_chrX     <- sumstats\$chr == 23

auto_ids_fwd <- paste0("chr", sumstats\$chr[is_autosome], ":",
                       sumstats\$pos[is_autosome], ":",
                       sumstats\$a1[is_autosome], ":",
                       sumstats\$a0[is_autosome])
auto_ids_rev <- paste0("chr", sumstats\$chr[is_autosome], ":",
                       sumstats\$pos[is_autosome], ":",
                       sumstats\$a0[is_autosome], ":",
                       sumstats\$a1[is_autosome])
chrx_ids_fwd <- paste0("chrX:", sumstats\$pos[is_chrX], ":",
                       sumstats\$a1[is_chrX], ":",
                       sumstats\$a0[is_chrX])
chrx_ids_rev <- paste0("chrX:", sumstats\$pos[is_chrX], ":",
                       sumstats\$a0[is_chrX], ":",
                       sumstats\$a1[is_chrX])

all_ids <- unique(c(auto_ids_fwd, auto_ids_rev, chrx_ids_fwd, chrx_ids_rev))
writeLines(all_ids, gwas_ids_file)
cat(sprintf("  IDs candidatos: %s\n", format(length(all_ids), big.mark=",")))
rm(auto_ids_fwd, auto_ids_rev, chrx_ids_fwd, chrx_ids_rev, all_ids); gc()

# ================================================================
# [3] LD-Pruning por cromosoma (plink2)
# ================================================================
cat("\n[3] LD-Pruning en MCPS (r²<0.9, ventana 500/50)...\n")

pruned_ids_all <- c()
for (chr in 1:23) {
  chr_label <- ifelse(chr == 23, "X", as.character(chr))
  prefix <- file.path(plink_dir, sprintf("mcps-freeze150k_qcd_chr%s_ivs", chr_label))
  if (!file.exists(paste0(prefix, ".bed"))) {
    cat(sprintf("    chr%s: NO EXISTE .bed, saltando\n", chr_label))
    next
  }
  prune_out <- file.path(tmp_dir, sprintf("prune_chr%s", chr_label))
  system(sprintf(
    "plink2 --bfile %s --extract %s --indep-pairwise 500 50 0.9 --out %s --threads %d --memory 8000 2>/dev/null",
    prefix, gwas_ids_file, prune_out, ncores))

  prune_in <- paste0(prune_out, ".prune.in")
  if (file.exists(prune_in)) {
    ids <- readLines(prune_in)
    pruned_ids_all <- c(pruned_ids_all, ids)
    cat(sprintf("    chr%s: %s variantes seleccionadas\n",
        chr_label, format(length(ids), big.mark=",")))
  } else {
    cat(sprintf("    chr%s: pruning falló\n", chr_label))
  }
}

cat(sprintf("\n  Total 'Filtro MCPS': %s variantes\n",
    format(length(pruned_ids_all), big.mark=",")))

# ================================================================
# [4] Filtrar sumstats a los pruned IDs
# ================================================================
pruned_parts <- strsplit(pruned_ids_all, ":")
pruned_chr_raw <- sapply(pruned_parts, \`[\`, 1)
pruned_chr <- as.integer(gsub("chr", "", gsub("chrX", "chr23", pruned_chr_raw)))
pruned_pos <- as.integer(sapply(pruned_parts, \`[\`, 2))

pruned_key <- paste(pruned_chr, pruned_pos, sep=":")
gwas_key <- paste(sumstats\$chr, sumstats\$pos, sep=":")
sumstats <- sumstats[gwas_key %in% pruned_key, ]

n_auto <- sum(sumstats\$chr != 23)
n_chrx <- sum(sumstats\$chr == 23)
cat(sprintf("\n  sumstats post-pruning: %s (%s autosomas + %s chrX)\n",
    format(nrow(sumstats), big.mark=","),
    format(n_auto, big.mark=","),
    format(n_chrx, big.mark=",")))

# ================================================================
# [5] Guardar outputs para workers
# ================================================================
saveRDS(sumstats, file.path(out_dir, "sumstats_pruned.rds"))
saveRDS(pruned_ids_all, file.path(out_dir, "pruned_ids.rds"))

# Resumen por cromosoma
cat("\n  Variantes por cromosoma post-pruning:\n")
for (chr in sort(unique(sumstats\$chr))) {
  n_chr <- sum(sumstats\$chr == chr)
  lbl <- ifelse(chr == 23, "X", as.character(chr))
  cat(sprintf("    chr%s: %s\n", lbl, format(n_chr, big.mark=",")))
}

cat(sprintf("\n✓ Outputs guardados en %s:\n", out_dir))
cat("    - sumstats_pruned.rds\n")
cat("    - pruned_ids.rds\n")
cat("\n  Siguiente paso: sbatch ldp2_custom_v6_launcher.sbatch\n")
REOF

apptainer exec --bind /mnt/cephfs:/mnt/cephfs \
    --bind "$TMPDIR:/tmp,$TMPDIR:/var/tmp" \
    --env OMP_NUM_THREADS=1 \
    --env OPENBLAS_NUM_THREADS=1 \
    --env MKL_NUM_THREADS=1 \
    "$SIF" /opt/miniforge3/envs/r-bio/bin/Rscript "$R_SCRIPT"

rm -rf "$TMPDIR"

echo ""
echo "============================================================"
echo "  Pruning completado — $(date)"
echo "  Siguiente: sbatch ldp2_custom_v6_launcher.sbatch"
echo "============================================================"
