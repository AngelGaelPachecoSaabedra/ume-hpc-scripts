#!/bin/bash
#SBATCH --job-name=ldp2c_sex0
#SBATCH --account=researchers
#SBATCH --qos=vip
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=150G
#SBATCH --output=/mnt/cephfs/orgs/home/angel.pacheco/logs/ldp2c_sex0_%j.out
#SBATCH --error=/mnt/cephfs/orgs/home/angel.pacheco/logs/ldp2c_sex0_%j.err
#SBATCH --mail-user=angelpachecosaavedra@gmail.com
#SBATCH --mail-type=END,FAIL

set -euo pipefail

# ================================================================
# PIPELINE C SEX — STEP 0: sumstats QC + pruning r2<0.85 POR SEXO
# Genera: sumstats_qc_females.rds y sumstats_qc_males.rds
# Sin IVW: cada sexo conserva sus betas crudos del GWAS REGENIE.
# n_eff por sexo = 4*cases*ctrls/(cases+ctrls), SIN K(1-K).
# ================================================================

SIF="/mnt/cephfs/biocontainers/images/jupyter-biotools-1.5.sif"
PLINK_DIR="/mnt/cephfs/hot_nvme/mcps/imputed-topmed/plink_files/maximally_unrelated"
GWAS_DIR="/mnt/cephfs/orgs/home/angel.pacheco/gwas_miami"
OUT_DIR="/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes/ldpred2_custom"
mkdir -p "$OUT_DIR"

TMPDIR="/mnt/cephfs/m.2/scratch/$USER/ldp2c_sex0_$SLURM_JOB_ID"
mkdir -p "$TMPDIR"; chmod 700 "$TMPDIR"
export TMPDIR TMP="$TMPDIR" TEMP="$TMPDIR"

echo "============================================================"
echo "  PIPELINE C SEX STEP 0 — sumstats por sexo — Job $SLURM_JOB_ID"
echo "  $(date)"
echo "============================================================"

R_SCRIPT="$TMPDIR/prune_sex.R"
cat << 'REOF' > "$R_SCRIPT"
Sys.setenv(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
suppressMessages({ library(data.table) })

ncores    <- 32L
gwas_dir  <- Sys.getenv("GWAS_DIR")
plink_dir <- Sys.getenv("PLINK_DIR")
out_dir   <- Sys.getenv("OUT_DIR")
tmp_dir   <- Sys.getenv("TMPDIR")

process_sex <- function(file_gz, sex_label, neff_const) {
  cat(sprintf("\n==== Procesando %s (n_eff=%.0f) ====\n", sex_label, neff_const))
  g <- fread(file.path(gwas_dir, file_gz),
             select=c("CHROM","GENPOS","ID","ALLELE0","ALLELE1","A1FREQ","INFO","N","BETA","SE","P"),
             showProgress=FALSE)
  sumstats <- data.frame(
    chr=g$CHROM, pos=g$GENPOS, a1=g$ALLELE1, a0=g$ALLELE0,
    freq=g$A1FREQ, info=g$INFO, beta=g$BETA, beta_se=g$SE,
    n_eff=neff_const, stringsAsFactors=FALSE)
  sumstats$p <- 2*pnorm(-abs(sumstats$beta/sumstats$beta_se))
  rm(g); gc()

  sumstats <- sumstats[sumstats$info>=0.8 & sumstats$freq>=0.01 & sumstats$freq<=0.99 &
                       !is.na(sumstats$beta) & !is.na(sumstats$beta_se) & sumstats$beta_se>0, ]
  chr_raw <- as.character(sumstats$chr); chr_raw[chr_raw %in% c("X","chrX","23")] <- "23"
  sumstats$chr <- suppressWarnings(as.integer(chr_raw))
  sumstats <- sumstats[!is.na(sumstats$chr) & sumstats$chr>=1 & sumstats$chr<=23, ]
  cat(sprintf("  Tras QC: %s variantes\n", format(nrow(sumstats), big.mark=",")))

  # Pruning r2<0.85 por cromosoma
  chr_lab <- ifelse(sumstats$chr==23L, "X", as.character(sumstats$chr))
  ids_fwd <- paste0("chr", chr_lab, ":", sumstats$pos, ":", sumstats$a1, ":", sumstats$a0)
  ids_rev <- paste0("chr", chr_lab, ":", sumstats$pos, ":", sumstats$a0, ":", sumstats$a1)
  ids_file <- file.path(tmp_dir, sprintf("gwas_ids_%s.txt", sex_label))
  writeLines(unique(c(ids_fwd, ids_rev)), ids_file)
  rm(ids_fwd, ids_rev, chr_lab); gc()

  pruned <- c()
  for (chr in 1:23) {
    cl <- if (chr==23L) "X" else as.character(chr)
    prefix <- file.path(plink_dir, sprintf("mcps-freeze150k_qcd_chr%s_ivs", cl))
    if (!file.exists(paste0(prefix, ".bed"))) next
    po <- file.path(tmp_dir, sprintf("prune_%s_chr%s", sex_label, cl))
    system(sprintf("plink2 --bfile %s --extract %s --indep-pairwise 500 50 0.85 --out %s --threads %d --memory 80000 2>/dev/null",
                   prefix, ids_file, po, ncores))
    pin <- paste0(po, ".prune.in")
    if (file.exists(pin)) { ids <- readLines(pin); pruned <- c(pruned, ids)
      cat(sprintf("    chr%s: %s\n", cl, format(length(ids), big.mark=","))) }
  }
  pp <- strsplit(pruned, ":")
  pcr <- gsub("chr","", sapply(pp, `[`, 1)); pcr[pcr=="X"] <- "23"
  pc <- suppressWarnings(as.integer(pcr)); ppos <- as.integer(sapply(pp, `[`, 2))
  ok <- !is.na(pc) & !is.na(ppos)
  pkey <- paste(pc[ok], ppos[ok], sep=":")
  gkey <- paste(sumstats$chr, sumstats$pos, sep=":")
  sumstats <- sumstats[gkey %in% pkey, ]
  cat(sprintf("  %s filtrado: %s (%s autosomas + %s chrX)\n", sex_label,
      format(nrow(sumstats), big.mark=","),
      format(sum(sumstats$chr<=22), big.mark=","),
      format(sum(sumstats$chr==23), big.mark=",")))

  saveRDS(sumstats, file.path(out_dir, sprintf("sumstats_qc_%s.rds", sex_label)))
  cat(sprintf("  Guardado: sumstats_qc_%s.rds\n", sex_label))
}

neff_f <- 4*17300*54196/(17300+54196)   # ~52,457
neff_m <- 4*8522*26733/(8522+26733)     # ~25,849
process_sex("combined-gwas-results-all-ctrl-females.txt.gz", "females", neff_f)
process_sex("combined-gwas-results-all-ctrl-males.txt.gz",   "males",   neff_m)
REOF

apptainer exec --bind /mnt/cephfs:/mnt/cephfs --bind "$TMPDIR:/tmp,$TMPDIR:/var/tmp" \
    --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 \
    --env GWAS_DIR="$GWAS_DIR" --env PLINK_DIR="$PLINK_DIR" --env OUT_DIR="$OUT_DIR" \
    "$SIF" /opt/miniforge3/envs/r-bio/bin/Rscript "$R_SCRIPT"

rm -rf "$TMPDIR"
echo "============================================================"
echo "  SEX STEP 0 FIN — $(date)"
echo "============================================================"
