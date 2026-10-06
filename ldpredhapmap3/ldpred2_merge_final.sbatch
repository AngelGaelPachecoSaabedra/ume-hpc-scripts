#!/bin/bash
#SBATCH --job-name=ldpred2_merge
#SBATCH --account=researchers
#SBATCH --qos=vip
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --output=/mnt/cephfs/orgs/home/angel.pacheco/logs/ldpred2_merge_%j.out
#SBATCH --error=/mnt/cephfs/orgs/home/angel.pacheco/logs/ldpred2_merge_%j.err
#SBATCH --mail-user=angelpachecosaavedra@gmail.com
#SBATCH --mail-type=END,FAIL

set -euo pipefail

# ================================================================
# MERGE FINAL — combina los 22 .rds por cromosoma en el PRS total
# Solo correr DESPUÉS de que todos los worker chr{1..22} terminen.
# ================================================================

SIF="/mnt/cephfs/biocontainers/images/jupyter-biotools-1.4.sif"
OUT_DIR="/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes/ldpred2"

echo "============================================================"
echo "  LDpred2 MERGE FINAL"
echo "  Nodo: $SLURMD_NODENAME"
echo "  $(date)"
echo "============================================================"

TMPDIR="/mnt/cephfs/scratch/$USER/ldpred2_merge_tmp/$SLURM_JOB_ID"
mkdir -p "$TMPDIR"; chmod 700 "$TMPDIR"
R_SCRIPT="$TMPDIR/merge.R"

cat << 'REOF' > "$R_SCRIPT"
suppressMessages({ library(data.table) })
out_dir <- "/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes/ldpred2"

# ── 1) Cargar todos los chrs ────────────────────────────────────────
cat("[1] Cargando resultados por cromosoma...\n")
chr_data <- list()
for (chr in 1:22) {
  f <- file.path(out_dir, sprintf("ldpred2_betas_chr%d.rds", chr))
  if (!file.exists(f)) {
    cat(sprintf("  ⚠ chr%d: FALTA %s — excluido del merge\n", chr, basename(f)))
    next
  }
  chr_data[[as.character(chr)]] <- readRDS(f)
}
cat(sprintf("  Cromosomas cargados: %d/22\n", length(chr_data)))

# ── 2) Resumen por cromosoma ────────────────────────────────────────
cat("\n[2] Resumen:\n")
cat(sprintf("  %-6s %-10s %-8s %-10s %-10s %-10s %s\n",
            "chr", "n_matched", "chains", "h2", "p", "total_min", "status"))
for (cl in names(chr_data)) {
  d <- chr_data[[cl]]
  cat(sprintf("  %-6s %-10s %-8d %-10s %-10s %-10.1f %s\n",
              cl,
              format(d$n_matched, big.mark=","),
              d$n_valid_chains,
              if (is.na(d$h2)) "NA" else sprintf("%.4f", d$h2),
              if (is.na(d$p))  "NA" else sprintf("%.2e", d$p),
              d$t_cor_min + d$t_ldpred_min,
              d$status))
}

# ── 3) Construir PRS_TOTAL ──────────────────────────────────────────
cat("\n[3] Construyendo PRS total...\n")
# Validar que todos tengan mismo fam
fam_ref <- chr_data[[1]]$fam
for (cl in names(chr_data)) {
  if (!identical(chr_data[[cl]]$fam$sample.ID, fam_ref$sample.ID))
    stop(sprintf("chr%s tiene fam diferente", cl))
}

prs_total <- rep(0.0, nrow(fam_ref))
chrs_ok   <- c()
chrs_skip <- c()
matched_total <- 0L

for (cl in names(chr_data)) {
  d <- chr_data[[cl]]
  if (d$status == "ok") {
    prs_total <- prs_total + d$prs_chr
    chrs_ok   <- c(chrs_ok, cl)
    matched_total <- matched_total + d$n_matched
  } else {
    chrs_skip <- c(chrs_skip, cl)
    # d$prs_chr es vector de ceros para diverged_zeroed — no aporta
  }
}
cat(sprintf("  Chrs incluidos (OK): %s\n", paste(chrs_ok, collapse=",")))
if (length(chrs_skip) > 0)
  cat(sprintf("  Chrs excluidos (divergieron): %s\n", paste(chrs_skip, collapse=",")))
cat(sprintf("  Variantes usadas: %s\n", format(matched_total, big.mark=",")))
cat(sprintf("  NAs en prs_total: %d\n", sum(is.na(prs_total))))

# ── 4) Link IID → PATID + sexo ──────────────────────────────────────
cat("\n[4] Vinculando IID → PATID + sexo...\n")
covar_file <- "/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes/MCPS BASELINE.csv"
link_file  <- "/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes/RGN_LINK_IID.csv"

df_link      <- read.csv(link_file)
iid_to_patid <- setNames(as.character(df_link$PATID), as.character(df_link$IID))
df_baseline  <- read.csv(covar_file)
df_baseline  <- df_baseline[!is.na(df_baseline$MALE), ]
patid_to_sex <- setNames(as.integer(df_baseline$MALE), as.character(df_baseline$PATID))

patids <- iid_to_patid[as.character(fam_ref$sample.ID)]
patids[is.na(patids)] <- "UNKNOWN"
sexes  <- patid_to_sex[patids]
sexes[is.na(sexes)] <- -1L

df_out <- data.frame(
  PATID     = patids,
  IID       = fam_ref$sample.ID,
  SEX       = sexes,
  PRS_TOTAL = prs_total,
  stringsAsFactors = FALSE
)
df_out <- df_out[df_out$SEX >= 0L, ]

# Filtrar NAs (por si acaso)
n_na <- sum(is.na(df_out$PRS_TOTAL))
if (n_na > 0) cat(sprintf("  ⚠ %d muestras con PRS=NA excluidas\n", n_na))
df_out <- df_out[!is.na(df_out$PRS_TOTAL), ]

# ── 5) Standardizar por sexo (con guardas robustas) ─────────────────
cat("\n[5] Standardizando por sexo...\n")
df_out$PRS_Z      <- NA_real_
df_out$PERCENTILE <- NA_real_
for (sv in c(0L, 1L)) {
  mask <- df_out$SEX == sv
  vals <- df_out$PRS_TOTAL[mask]
  sd_v <- sd(vals, na.rm=TRUE)
  if (length(vals) > 1L && !is.na(sd_v) && sd_v > 0) {
    df_out$PRS_Z[mask]      <- (vals - mean(vals, na.rm=TRUE)) / sd_v
    df_out$PERCENTILE[mask] <- rank(vals) / length(vals) * 100
    cat(sprintf("  SEX=%d: n=%d, mean=%.4f, sd=%.4f\n",
                sv, length(vals), mean(vals, na.rm=TRUE), sd_v))
  } else {
    cat(sprintf("  ⚠ SEX=%d: no se pudo standardizar (sd=%.4f)\n", sv, sd_v))
  }
}

# ── 6) Export ──────────────────────────────────────────────────────
df_out <- df_out[order(-df_out$PRS_TOTAL), ]
out_file <- file.path(out_dir, "prs_diabetes_ldpred2_final.tsv")
write.table(df_out, out_file, sep="\t", row.names=FALSE, quote=FALSE)

# Guardar también el log agregado
saveRDS(list(
  chr_data      = lapply(chr_data, function(d) d[c("chr","n_matched","h2","p",
                                                    "n_valid_chains","status",
                                                    "burn_in","num_iter","h2_init",
                                                    "t_cor_min","t_ldpred_min")]),
  chrs_ok       = chrs_ok,
  chrs_skip     = chrs_skip,
  matched_total = matched_total,
  n_samples     = nrow(df_out)
), file.path(out_dir, "ldpred2_hm3_results_log.rds"))

cat(sprintf("\n✓ PIPELINE COMPLETADO\n"))
cat(sprintf("  Output    : %s\n", out_file))
cat(sprintf("  Muestras  : %s\n", format(nrow(df_out), big.mark=",")))
cat(sprintf("  Variantes : %s (de %d chrs)\n", format(matched_total, big.mark=","), length(chrs_ok)))
if (length(chrs_skip) > 0)
  cat(sprintf("  Excluidos : chrs %s (divergencia)\n", paste(chrs_skip, collapse=",")))
REOF

apptainer exec --bind /mnt/cephfs:/mnt/cephfs \
    "$SIF" /opt/miniforge3/envs/r-bio/bin/Rscript "$R_SCRIPT"

rm -rf "$TMPDIR"

echo ""
echo "============================================================"
echo "  MERGE COMPLETADO — $(date)"
echo "============================================================"
