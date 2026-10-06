#!/bin/bash
#SBATCH --job-name=ldp2c_prune
#SBATCH --account=researchers
#SBATCH --qos=vip
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=150G
#SBATCH --output=/mnt/cephfs/orgs/home/angel.pacheco/logs/ldp2c_prune_%j.out
#SBATCH --error=/mnt/cephfs/orgs/home/angel.pacheco/logs/ldp2c_prune_%j.err
#SBATCH --mail-user=angelpachecosaavedra@gmail.com
#SBATCH --mail-type=END,FAIL

set -euo pipefail

# ================================================================
# PIPELINE C — STEP 0: LD-Pruning r²<0.9 en MCPS (1 sola vez)
# Genera:
#   - pruned_keys.txt   : chr:pos de variantes que pasan pruning (para R)
#   - sumstats_qc.rds   : sumstats meta-IVW tras QC (para Step 1 y 2)
# Esto evita repetir el pruning en cada tarea del array de Step 1.
# ================================================================

SIF="/mnt/cephfs/biocontainers/images/jupyter-biotools-1.5.sif"
PLINK_DIR="/mnt/cephfs/hot_nvme/mcps/imputed-topmed/plink_files/maximally_unrelated"
GWAS_DIR="/mnt/cephfs/orgs/home/angel.pacheco/gwas_miami"
OUT_DIR="/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes/ldpred2_custom"
mkdir -p "$OUT_DIR"

TMPDIR="/mnt/cephfs/m.2/scratch/$USER/ldp2c_prune_$SLURM_JOB_ID"
mkdir -p "$TMPDIR"; chmod 700 "$TMPDIR"
export TMPDIR TMP="$TMPDIR" TEMP="$TMPDIR"

echo "============================================================"
echo "  PIPELINE C STEP 0 — Pruning r²<0.9 — Job $SLURM_JOB_ID"
echo "  $(date)"
echo "============================================================"

R_SCRIPT="$TMPDIR/prune.R"
cat << 'REOF' > "$R_SCRIPT"
Sys.setenv(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
suppressMessages({ library(data.table) })

ncores    <- 32L
gwas_dir  <- Sys.getenv("GWAS_DIR")
plink_dir <- Sys.getenv("PLINK_DIR")
out_dir   <- Sys.getenv("OUT_DIR")
tmp_dir   <- Sys.getenv("TMPDIR")

# ── GWAS meta-análisis IVW (idéntico a v6) ─────────────────────────
cat("[1] GWAS meta-análisis IVW...\n")
gf <- fread(file.path(gwas_dir, "combined-gwas-results-all-ctrl-females.txt.gz"),
            select=c("CHROM","GENPOS","ID","ALLELE0","ALLELE1","A1FREQ","INFO","N","BETA","SE","P"),
            showProgress=FALSE)
gm <- fread(file.path(gwas_dir, "combined-gwas-results-all-ctrl-males.txt.gz"),
            select=c("CHROM","GENPOS","ID","ALLELE0","ALLELE1","A1FREQ","INFO","N","BETA","SE","P"),
            showProgress=FALSE)
gwas <- merge(gf, gm, by="ID", suffixes=c(".f",".m"))
rm(gf, gm); gc()

w_f <- 1/gwas$SE.f^2; w_m <- 1/gwas$SE.m^2
# n_eff CORRECTO (sin K(1-K)) — consistente con HapMap3 re-corrido
neff_f <- 4*17300*54196/(17300+54196)
neff_m <- 4*8522*26733/(8522+26733)
neff_ldpred2 <- neff_f + neff_m   # 78,304
cat(sprintf("  n_eff = %.0f\n", neff_ldpred2))

sumstats <- data.frame(
    chr = gwas$CHROM.f, pos = gwas$GENPOS.f,
    a1 = gwas$ALLELE1.f, a0 = gwas$ALLELE0.f,
    freq = (gwas$A1FREQ.f*gwas$N.f + gwas$A1FREQ.m*gwas$N.m)/(gwas$N.f+gwas$N.m),
    info = pmin(gwas$INFO.f, gwas$INFO.m),
    beta = (gwas$BETA.f*w_f + gwas$BETA.m*w_m)/(w_f+w_m),
    beta_se = sqrt(1/(w_f+w_m)), n_eff = neff_ldpred2,
    stringsAsFactors = FALSE)
sumstats$p <- 2*pnorm(-abs(sumstats$beta/sumstats$beta_se))
rm(gwas); gc()

sumstats <- sumstats[sumstats$info>=0.8 & sumstats$freq>=0.01 & sumstats$freq<=0.99 &
                     !is.na(sumstats$beta) & !is.na(sumstats$beta_se) & sumstats$beta_se>0, ]
chr_raw <- as.character(sumstats$chr)
chr_raw[chr_raw %in% c("X","chrX","23")] <- "23"
sumstats$chr <- suppressWarnings(as.integer(chr_raw))
sumstats <- sumstats[!is.na(sumstats$chr) & sumstats$chr >= 1 & sumstats$chr <= 23, ]
cat(sprintf("  GWAS tras QC: %s (%s autosomas + %s chrX)\n",
    format(nrow(sumstats), big.mark=","),
    format(sum(sumstats$chr<=22), big.mark=","),
    format(sum(sumstats$chr==23), big.mark=",")))

# ── LD-pruning r²<0.9 por cromosoma ────────────────────────────────
cat("\n[2] LD-pruning r²<0.85...\n")
gwas_ids_file <- file.path(tmp_dir, "gwas_ids.txt")
chr_lab <- ifelse(sumstats$chr==23L, "X", as.character(sumstats$chr))
ids_fwd <- paste0("chr", chr_lab, ":", sumstats$pos, ":", sumstats$a1, ":", sumstats$a0)
ids_rev <- paste0("chr", chr_lab, ":", sumstats$pos, ":", sumstats$a0, ":", sumstats$a1)
writeLines(unique(c(ids_fwd, ids_rev)), gwas_ids_file)
rm(ids_fwd, ids_rev, chr_lab); gc()

pruned_ids_all <- c()
for (chr in 1:23) {
    chr_label <- if (chr==23L) "X" else as.character(chr)
    prefix <- file.path(plink_dir, sprintf("mcps-freeze150k_qcd_chr%s_ivs", chr_label))
    if (!file.exists(paste0(prefix, ".bed"))) next
    prune_out <- file.path(tmp_dir, sprintf("prune_chr%s", chr_label))
    system(sprintf(
        "plink2 --bfile %s --extract %s --indep-pairwise 500 50 0.85 --out %s --threads %d --memory 80000 2>/dev/null",
        prefix, gwas_ids_file, prune_out, ncores))
    pin <- paste0(prune_out, ".prune.in")
    if (file.exists(pin)) {
        ids <- readLines(pin); pruned_ids_all <- c(pruned_ids_all, ids)
        cat(sprintf("    chr%s: %s\n", chr_label, format(length(ids), big.mark=",")))
    }
}
cat(sprintf("  Total pruned: %s\n", format(length(pruned_ids_all), big.mark=",")))

# Filtrar sumstats a variantes pruned
pp <- strsplit(pruned_ids_all, ":")
pc_raw <- gsub("chr","", sapply(pp, `[`, 1)); pc_raw[pc_raw=="X"] <- "23"
pc <- suppressWarnings(as.integer(pc_raw))
pos <- as.integer(sapply(pp, `[`, 2))
ok <- !is.na(pc) & !is.na(pos)
pruned_key <- paste(pc[ok], pos[ok], sep=":")
gwas_key <- paste(sumstats$chr, sumstats$pos, sep=":")
sumstats <- sumstats[gwas_key %in% pruned_key, ]
cat(sprintf("  Sumstats filtrado: %s (%s autosomas + %s chrX)\n",
    format(nrow(sumstats), big.mark=","),
    format(sum(sumstats$chr<=22), big.mark=","),
    format(sum(sumstats$chr==23), big.mark=",")))

# ── Guardar para Step 1 y 2 ────────────────────────────────────────
saveRDS(sumstats, file.path(out_dir, "sumstats_qc.rds"))
cat(sprintf("\n Guardado: %s/sumstats_qc.rds (%s variantes)\n",
    out_dir, format(nrow(sumstats), big.mark=",")))
REOF

apptainer exec --bind /mnt/cephfs:/mnt/cephfs --bind "$TMPDIR:/tmp,$TMPDIR:/var/tmp" \
    --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 \
    --env GWAS_DIR="$GWAS_DIR" --env PLINK_DIR="$PLINK_DIR" --env OUT_DIR="$OUT_DIR" \
    "$SIF" /opt/miniforge3/envs/r-bio/bin/Rscript "$R_SCRIPT"

rm -rf "$TMPDIR"
echo "============================================================"
echo "  STEP 0 FIN — $(date)"
echo "============================================================"
