#!/bin/bash
#SBATCH --job-name=ldp2_custom_merge
#SBATCH --account=researchers
#SBATCH --qos=vip
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=/mnt/cephfs/orgs/home/angel.pacheco/logs/ldp2_custom_merge_%j.out
#SBATCH --error=/mnt/cephfs/orgs/home/angel.pacheco/logs/ldp2_custom_merge_%j.err
#SBATCH --mail-user=angelpachecosaavedra@gmail.com
#SBATCH --mail-type=END,FAIL

set -euo pipefail

# ================================================================
# PIPELINE C v6 — MERGE FINAL
# Solo correr DESPUÉS de que todos los workers chr{1..23} terminen
# ================================================================

SIF="/mnt/cephfs/biocontainers/images/jupyter-biotools-1.4.sif"
OUT_DIR="/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes/ldpred2_custom"

echo "============================================================"
echo "  PIPELINE C v6 MERGE FINAL"
echo "  Nodo: $SLURMD_NODENAME"
echo "  $(date)"
echo "============================================================"

TMPDIR="/mnt/cephfs/scratch/$USER/ldp2_custom_merge_tmp/$SLURM_JOB_ID"
mkdir -p "$TMPDIR"; chmod 700 "$TMPDIR"
R_SCRIPT="$TMPDIR/merge.R"

cat << REOF > "$R_SCRIPT"
suppressMessages({ library(data.table) })
out_dir <- "$OUT_DIR"

# ── 1) Cargar todos los chrs ────────────────────────────────────────
cat("[1] Cargando resultados por cromosoma...\n")
chr_data <- list()
labels <- c(as.character(1:22), "X")
for (lbl in labels) {
  f <- file.path(out_dir, sprintf("ldpred2_custom_betas_chr%s.rds", lbl))
  if (!file.exists(f)) {
    cat(sprintf("  chr%s: FALTA — excluido\n", lbl))
    next
  }
  chr_data[[lbl]] <- readRDS(f)
}
cat(sprintf("  Cromosomas cargados: %d/23\n", length(chr_data)))

# ── 2) Resumen ──────────────────────────────────────────────────────
cat("\n[2] Resumen:\n")
cat(sprintf("  %-6s %-10s %-8s %-10s %-10s %-10s %s\n",
            "chr", "n_matched", "chains", "h2", "p", "total_min", "status"))
for (lbl in names(chr_data)) {
  d <- chr_data[[lbl]]
  if (d\$status == "too_few_variants") {
    cat(sprintf("  %-6s %-10s %-8s %-10s %-10s %-10s %s\n",
                lbl, format(d\$n_matched, big.mark=","),
                "-", "-", "-", "-", d\$status))
  } else {
    cat(sprintf("  %-6s %-10s %-8d %-10s %-10s %-10.1f %s\n",
                lbl,
                format(d\$n_matched, big.mark=","),
                d\$n_valid_chains,
                if (is.na(d\$h2)) "NA" else sprintf("%.4f", d\$h2),
                if (is.na(d\$p)) "NA" else sprintf("%.2e", d\$p),
                d\$t_cor_min + d\$t_ldpred_min,
                d\$status))
  }
}

# ── 3) PRS total ────────────────────────────────────────────────────
cat("\n[3] Construyendo PRS total...\n")

# Usar fam del primer chr con datos válidos
fam_ref <- NULL
for (lbl in names(chr_data)) {
  if (!is.null(chr_data[[lbl]]\$fam)) {
    fam_ref <- chr_data[[lbl]]\$fam
    cat(sprintf("  fam tomado de chr%s (%d muestras)\n", lbl, nrow(fam_ref)))
    break
  }
}
if (is.null(fam_ref)) stop("Ningún chr tiene fam — ¿todos fallaron?")

# Validar que todos los chrs con prs_chr tengan el mismo fam
for (lbl in names(chr_data)) {
  d <- chr_data[[lbl]]
  if (!is.null(d\$fam) && !identical(d\$fam\$sample.ID, fam_ref\$sample.ID))
    stop(sprintf("chr%s tiene fam diferente", lbl))
}

prs_total  <- rep(0.0, nrow(fam_ref))
chrs_ok    <- c(); chrs_skip <- c(); chrs_few <- c()
matched_total <- 0L

for (lbl in names(chr_data)) {
  d <- chr_data[[lbl]]
  if (d\$status == "ok") {
    prs_total <- prs_total + d\$prs_chr
    chrs_ok   <- c(chrs_ok, lbl)
    matched_total <- matched_total + d\$n_matched
  } else if (d\$status == "diverged_zeroed") {
    chrs_skip <- c(chrs_skip, lbl)
  } else if (d\$status == "too_few_variants") {
    chrs_few <- c(chrs_few, lbl)
  }
}

cat(sprintf("  Chrs incluidos (OK)       : %s\n", paste(chrs_ok, collapse=",")))
if (length(chrs_skip) > 0)
  cat(sprintf("  Chrs excluidos (divergieron): %s\n", paste(chrs_skip, collapse=",")))
if (length(chrs_few) > 0)
  cat(sprintf("  Chrs excluidos (pocas vars) : %s\n", paste(chrs_few, collapse=",")))
cat(sprintf("  Variantes usadas: %s\n", format(matched_total, big.mark=",")))

# ── 4) Link IID → PATID ─────────────────────────────────────────────
cat("\n[4] Vinculando IID → PATID + sexo...\n")
covar_file <- "/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes/MCPS BASELINE.csv"
link_file  <- "/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes/RGN_LINK_IID.csv"

df_link      <- read.csv(link_file)
iid_to_patid <- setNames(as.character(df_link\$PATID), as.character(df_link\$IID))
df_baseline  <- read.csv(covar_file)
df_baseline  <- df_baseline[!is.na(df_baseline\$MALE), ]
patid_to_sex <- setNames(as.integer(df_baseline\$MALE), as.character(df_baseline\$PATID))

patids <- iid_to_patid[as.character(fam_ref\$sample.ID)]
patids[is.na(patids)] <- "UNKNOWN"
sexes  <- patid_to_sex[patids]
sexes[is.na(sexes)] <- -1L

df_out <- data.frame(
  PATID=patids, IID=fam_ref\$sample.ID,
  SEX=sexes, PRS_TOTAL=prs_total,
  stringsAsFactors=FALSE)
df_out <- df_out[df_out\$SEX >= 0L, ]

n_na <- sum(is.na(df_out\$PRS_TOTAL))
if (n_na > 0) cat(sprintf("  ⚠ %d muestras excluidas por PRS=NA\n", n_na))
df_out <- df_out[!is.na(df_out\$PRS_TOTAL), ]

# ── 5) Standardizar por sexo ────────────────────────────────────────
cat("\n[5] Standardizando por sexo...\n")
df_out\$PRS_Z      <- NA_real_
df_out\$PERCENTILE <- NA_real_
for (sv in c(0L, 1L)) {
  mask <- df_out\$SEX == sv
  vals <- df_out\$PRS_TOTAL[mask]
  sd_v <- sd(vals, na.rm=TRUE)
  if (length(vals) > 1L && !is.na(sd_v) && sd_v > 0) {
    df_out\$PRS_Z[mask]      <- (vals - mean(vals, na.rm=TRUE)) / sd_v
    df_out\$PERCENTILE[mask] <- rank(vals) / length(vals) * 100
    cat(sprintf("  SEX=%d: n=%d, mean=%.4f, sd=%.4f\n",
                sv, length(vals), mean(vals, na.rm=TRUE), sd_v))
  } else {
    cat(sprintf("  SEX=%d: no se pudo standardizar (sd=%.4f)\n", sv, sd_v))
  }
}

# ── 6) Export ──────────────────────────────────────────────────────
df_out <- df_out[order(-df_out\$PRS_TOTAL), ]
out_file <- file.path(out_dir, "prs_diabetes_ldpred2_custom_mcps_final.tsv")
write.table(df_out, out_file, sep="\t", row.names=FALSE, quote=FALSE)

# Guardar log agregado
saveRDS(list(
  chr_data      = lapply(chr_data, function(d) {
    d[!names(d) %in% c("beta_final", "prs_chr", "fam", "info_snp", "ind_col")]
  }),
  chrs_ok       = chrs_ok,
  chrs_diverged = chrs_skip,
  chrs_few      = chrs_few,
  matched_total = matched_total,
  n_samples     = nrow(df_out)
), file.path(out_dir, "ldpred2_custom_results_log.rds"))

cat(sprintf("\n✓ PIPELINE C v6 COMPLETADO\n"))
cat(sprintf("  Output    : %s\n", out_file))
cat(sprintf("  Muestras  : %s\n", format(nrow(df_out), big.mark=",")))
cat(sprintf("  Variantes : %s (de %d chrs)\n",
            format(matched_total, big.mark=","), length(chrs_ok)))
if (length(chrs_skip) > 0)
  cat(sprintf("  Excluidos : chrs %s (divergencia)\n", paste(chrs_skip, collapse=",")))
if (length(chrs_few) > 0)
  cat(sprintf("  Excluidos : chrs %s (pocas variantes)\n", paste(chrs_few, collapse=",")))
REOF

apptainer exec --bind /mnt/cephfs:/mnt/cephfs \
    "$SIF" /opt/miniforge3/envs/r-bio/bin/Rscript "$R_SCRIPT"

rm -rf "$TMPDIR"

echo ""
echo "============================================================"
echo "  MERGE COMPLETADO — $(date)"
echo "============================================================"
