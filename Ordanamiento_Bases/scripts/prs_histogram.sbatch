#!/bin/bash
#SBATCH --job-name=prs_hist
#SBATCH --account=researchers
#SBATCH --qos=vip
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:15:00
#SBATCH --output=/mnt/cephfs/orgs/home/angel.pacheco/logs/prs_hist_%j.out
#SBATCH --error=/mnt/cephfs/orgs/home/angel.pacheco/logs/prs_hist_%j.err

set -euo pipefail

SIF="/mnt/cephfs/biocontainers/images/jupyter-biotools-1.2.sif"

TMPDIR="/mnt/cephfs/scratch/$USER/hist_tmp/$SLURM_JOB_ID"
mkdir -p "$TMPDIR"
chmod 700 "$TMPDIR"
export TMPDIR TMP="$TMPDIR" TEMP="$TMPDIR"

echo "==== HISTOGRAMA PRS EN R (ggplot2) ===="

R_SCRIPT="$TMPDIR/prs_histogram.R"
cat << 'EOF' > "$R_SCRIPT"
library(ggplot2)
library(scales)

out_dir <- "/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes"

# Cargar datos
cat("Cargando datos...\n")
df <- read.delim(file.path(out_dir, "prs_diabetes_t2d_final.tsv"))

# Estandarizar
prs_mean <- mean(df$PRS_T2D)
prs_sd   <- sd(df$PRS_T2D)
df$PRS_Z <- (df$PRS_T2D - prs_mean) / prs_sd

cat(sprintf("Muestras: %s\n", format(nrow(df), big.mark = ",")))
cat(sprintf("Media: %.4f, SD: %.4f\n", prs_mean, prs_sd))
cat(sprintf("Rango z-score: [%.2f, %.2f]\n", min(df$PRS_Z), max(df$PRS_Z)))

# Percentiles
p5  <- quantile(df$PRS_Z, 0.05)
p95 <- quantile(df$PRS_Z, 0.95)
n_high <- sum(df$PRS_Z >= 2)
n_low  <- sum(df$PRS_Z <= -2)

# Categoría de riesgo para colorear
df$risk <- ifelse(df$PRS_Z <= -2, "Bajo riesgo (\u2264 -2\u03c3)",
           ifelse(df$PRS_Z >= 2, "Alto riesgo (\u2265 2\u03c3)",
           "Riesgo intermedio"))

# Histograma
p <- ggplot(df, aes(x = PRS_Z, fill = risk)) +
  geom_histogram(bins = 80, color = "white", linewidth = 0.3, alpha = 0.85) +

  scale_fill_manual(
    values = c("Bajo riesgo (\u2264 -2\u03c3)" = "#4CAF50",
               "Riesgo intermedio"               = "#2196F3",
               "Alto riesgo (\u2265 2\u03c3)"   = "#F44336"),
    name = NULL
  ) +

  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +

  geom_vline(xintercept = p5, linetype = "dotted", color = "#4CAF50", linewidth = 0.9) +
  geom_vline(xintercept = p95, linetype = "dotted", color = "#F44336", linewidth = 0.9) +

  annotate("label", x = max(df$PRS_Z) - 0.3, y = Inf, vjust = 1.5,
    label = sprintf("Alto riesgo (\u22652\u03c3): %s (%.1f%%)\nBajo riesgo (\u2264-2\u03c3): %s (%.1f%%)",
      format(n_high, big.mark = ","), n_high / nrow(df) * 100,
      format(n_low, big.mark = ","),  n_low / nrow(df) * 100),
    hjust = 1, size = 3.5, fill = "white", alpha = 0.8, label.size = 0.3) +

  annotate("text", x = p5, y = 0, label = sprintf("P5 (%.2f\u03c3)", p5),
    vjust = -0.5, hjust = 1.1, size = 3, color = "#4CAF50") +
  annotate("text", x = p95, y = 0, label = sprintf("P95 (%.2f\u03c3)", p95),
    vjust = -0.5, hjust = -0.1, size = 3, color = "#F44336") +

  labs(
    title = "Distribuci\u00f3n del Polygenic Risk Score \u2014 Diabetes Tipo 2",
    subtitle = sprintf("Cohorte MCPS (n = %s) | PGS000014 (Khera et al. 2018)",
      format(nrow(df), big.mark = ",")),
    x = "PRS Estandarizado (z-score)",
    y = "N\u00famero de Pacientes"
  ) +

  theme_minimal(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(color = "gray40", size = 11),
    plot.margin   = margin(15, 20, 15, 15),
    legend.position = "bottom",
    legend.background = element_rect(fill = "white", color = "gray80", linewidth = 0.3),
    panel.grid.minor = element_blank()
  ) +

  scale_y_continuous(labels = comma)

# Guardar
out_png <- file.path(out_dir, "prs_diabetes_histogram.png")
out_pdf <- file.path(out_dir, "prs_diabetes_histogram.pdf")

ggsave(out_png, p, width = 12, height = 7, dpi = 200)
ggsave(out_pdf, p, width = 12, height = 7)

cat(sprintf("PNG: %s\n", out_png))
cat(sprintf("PDF: %s\n", out_pdf))
cat("Histograma generado exitosamente.\n")
EOF

apptainer exec \
    --bind /mnt/cephfs:/mnt/cephfs \
    --bind "$TMPDIR:/tmp,$TMPDIR:/var/tmp" \
    "$SIF" \
    /opt/miniforge3/envs/r-bio/bin/Rscript "$R_SCRIPT"

rm -rf "$TMPDIR"
echo "==== COMPLETADO ===="
