#!/bin/bash
#SBATCH --job-name=prs_hist2
#SBATCH --account=researchers
#SBATCH --qos=vip
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:15:00
#SBATCH --output=/mnt/cephfs/orgs/home/angel.pacheco/logs/prs_hist2_%j.out
#SBATCH --error=/mnt/cephfs/orgs/home/angel.pacheco/logs/prs_hist2_%j.err

set -euo pipefail

SIF="/mnt/cephfs/biocontainers/images/jupyter-biotools-1.2.sif"

TMPDIR="/mnt/cephfs/scratch/$USER/hist_tmp/$SLURM_JOB_ID"
mkdir -p "$TMPDIR"
chmod 700 "$TMPDIR"
export TMPDIR TMP="$TMPDIR" TEMP="$TMPDIR"

echo "==== HISTOGRAMA PRS MCPS GWAS (ggplot2) ===="

R_SCRIPT="$TMPDIR/prs_histogram_mcps.R"
cat << 'EOF' > "$R_SCRIPT"
library(ggplot2)
library(scales)

out_dir <- "/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes"

# Cargar datos
cat("Cargando datos...\n")
df <- read.delim(file.path(out_dir, "prs_diabetes_mcps_gwas_v4_final.tsv"))

# Etiquetar sexo
df$SEX_LABEL <- ifelse(df$SEX == 0, "Mujeres", "Hombres")

cat(sprintf("Muestras totales: %s\n", format(nrow(df), big.mark = ",")))
cat(sprintf("Mujeres: %s\n", format(sum(df$SEX == 0), big.mark = ",")))
cat(sprintf("Hombres: %s\n", format(sum(df$SEX == 1), big.mark = ",")))

for (s in c("Mujeres", "Hombres")) {
    sub <- df[df$SEX_LABEL == s, ]
    cat(sprintf("  %s: mean_Z=%.4f, sd_Z=%.4f, range=[%.2f, %.2f]\n",
                s, mean(sub$PRS_Z), sd(sub$PRS_Z), min(sub$PRS_Z), max(sub$PRS_Z)))
}

# Percentiles y conteos extremos
p5  <- quantile(df$PRS_Z, 0.05)
p95 <- quantile(df$PRS_Z, 0.95)
n_high <- sum(df$PRS_Z >= 2)
n_low  <- sum(df$PRS_Z <= -2)

# Categoría de riesgo
df$risk <- ifelse(df$PRS_Z <= -2, "Bajo riesgo (\u2264 -2\u03c3)",
           ifelse(df$PRS_Z >= 2, "Alto riesgo (\u2265 2\u03c3)",
           "Riesgo intermedio"))

# ================================================================
# GRÁFICA 1: Histograma global con zonas de riesgo
# ================================================================
p1 <- ggplot(df, aes(x = PRS_Z, fill = risk)) +
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
    subtitle = sprintf("Cohorte MCPS (n = %s) | GWAS propio estratificado por sexo | M\u00e9todo C+T",
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

ggsave(file.path(out_dir, "prs_mcps_gwas_histogram.png"), p1, width = 12, height = 7, dpi = 200)
ggsave(file.path(out_dir, "prs_mcps_gwas_histogram.pdf"), p1, width = 12, height = 7)

# ================================================================
# GRÁFICA 2: Facetado por sexo
# ================================================================
p2 <- ggplot(df, aes(x = PRS_Z, fill = risk)) +
  geom_histogram(bins = 80, color = "white", linewidth = 0.3, alpha = 0.85) +
  facet_wrap(~SEX_LABEL, ncol = 1, scales = "free_y") +

  scale_fill_manual(
    values = c("Bajo riesgo (\u2264 -2\u03c3)" = "#4CAF50",
               "Riesgo intermedio"               = "#2196F3",
               "Alto riesgo (\u2265 2\u03c3)"   = "#F44336"),
    name = NULL
  ) +

  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +

  labs(
    title = "PRS Diabetes T2 por Sexo \u2014 GWAS MCPS",
    subtitle = "Panel LD: maximally unrelated (79,612) | QC: INFO \u2265 0.8, MAF \u2265 0.01",
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
    strip.text = element_text(face = "bold", size = 12),
    panel.grid.minor = element_blank()
  ) +
  scale_y_continuous(labels = comma)

ggsave(file.path(out_dir, "prs_mcps_gwas_histogram_by_sex.png"), p2, width = 12, height = 9, dpi = 200)
ggsave(file.path(out_dir, "prs_mcps_gwas_histogram_by_sex.pdf"), p2, width = 12, height = 9)

# ================================================================
# GRÁFICA 3: Contribución por cromosoma
# ================================================================
chr_cols <- grep("^PRS_chr", names(df), value = TRUE)
if (length(chr_cols) > 0) {
    chr_means <- data.frame(
        CHR = gsub("PRS_chr", "", chr_cols),
        Mean_Abs = sapply(chr_cols, function(col) mean(abs(df[[col]])))
    )
    chr_means$CHR <- factor(chr_means$CHR, levels = c(as.character(1:22), "X"))

    p3 <- ggplot(chr_means, aes(x = CHR, y = Mean_Abs, fill = Mean_Abs)) +
      geom_col(width = 0.7) +
      scale_fill_gradient(low = "#4DBBD5", high = "#E64B35") +
      labs(
        title = "Contribuci\u00f3n promedio por cromosoma al PRS",
        subtitle = "Valor absoluto medio | GWAS MCPS estratificado por sexo",
        x = "Cromosoma",
        y = "Contribuci\u00f3n absoluta media"
      ) +
      theme_minimal(base_size = 13) +
      theme(
        legend.position = "none",
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(color = "gray40", size = 10),
        panel.grid.minor = element_blank()
      )

    ggsave(file.path(out_dir, "prs_mcps_gwas_chr_contribution.png"), p3, width = 12, height = 5, dpi = 200)
    ggsave(file.path(out_dir, "prs_mcps_gwas_chr_contribution.pdf"), p3, width = 12, height = 5)
    cat("Grafica de cromosomas generada.\n")
} else {
    cat("No se encontraron columnas PRS_chr*, omitiendo grafica de cromosomas.\n")
}

cat("\nGraficas guardadas en:\n")
cat(sprintf("  %s/prs_mcps_gwas_histogram.png\n", out_dir))
cat(sprintf("  %s/prs_mcps_gwas_histogram_by_sex.png\n", out_dir))
cat(sprintf("  %s/prs_mcps_gwas_chr_contribution.png\n", out_dir))
cat("Completado.\n")
EOF

apptainer exec \
    --bind /mnt/cephfs:/mnt/cephfs \
    --bind "$TMPDIR:/tmp,$TMPDIR:/var/tmp" \
    "$SIF" \
    /opt/miniforge3/envs/r-bio/bin/Rscript "$R_SCRIPT"

rm -rf "$TMPDIR"
echo "==== COMPLETADO ===="
