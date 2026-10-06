#!/bin/bash
#SBATCH --job-name=prs_amr
#SBATCH --account=researchers
#SBATCH --qos=vip
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:15:00
#SBATCH --output=/mnt/cephfs/orgs/home/angel.pacheco/logs/prs_amr_%j.out
#SBATCH --error=/mnt/cephfs/orgs/home/angel.pacheco/logs/prs_amr_%j.err

set -euo pipefail

SIF="/mnt/cephfs/biocontainers/images/jupyter-biotools-1.2.sif"

TMPDIR="/mnt/cephfs/scratch/$USER/hist_tmp/$SLURM_JOB_ID"
mkdir -p "$TMPDIR"
chmod 700 "$TMPDIR"
export TMPDIR TMP="$TMPDIR" TEMP="$TMPDIR"

echo "==== HISTOGRAMA PRS POR TERTILES DE ANCESTRÍA AMR ===="

R_SCRIPT="$TMPDIR/prs_amr_tertiles.R"
cat << 'EOF' > "$R_SCRIPT"
library(ggplot2)
library(scales)

out_dir <- "/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes"

# ================================================================
# Cargar PRS y ancestría
# ================================================================
cat("Cargando PRS...\n")
df_prs <- read.delim(file.path(out_dir, "prs_diabetes_mcps_gwas_final.tsv"))

cat("Cargando ancestría...\n")
df_anc <- read.delim(file.path(out_dir, "ancestry-proportions-k4.txt"))

# Filtrar solo MCPS
df_anc <- df_anc[df_anc$population == "MCPS", ]

cat(sprintf("PRS: %s pacientes\n", format(nrow(df_prs), big.mark = ",")))
cat(sprintf("Ancestr\u00eda MCPS: %s pacientes\n", format(nrow(df_anc), big.mark = ",")))

# ================================================================
# Cruzar por IID
# ================================================================
df <- merge(df_prs, df_anc[, c("IID", "Score_AMR", "Score_EUR", "Score_AFR", "Score_EAS")],
            by = "IID", all.x = FALSE)

cat(sprintf("Pacientes cruzados (PRS + ancestr\u00eda): %s\n", format(nrow(df), big.mark = ",")))

# ================================================================
# Crear tertiles de AMR
# ================================================================
df$AMR_tertile <- cut(df$Score_AMR,
                      breaks = quantile(df$Score_AMR, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE),
                      include.lowest = TRUE,
                      labels = c("T1: Baja AMR", "T2: Media AMR", "T3: Alta AMR"))

# Estadísticas por tertil
cat("\nEstad\u00edsticas por tertil de ancestr\u00eda AMR:\n")
for (t in levels(df$AMR_tertile)) {
    sub <- df[df$AMR_tertile == t, ]
    cat(sprintf("  %s (n=%s): AMR=[%.2f-%.2f], PRS_Z mean=%.4f, sd=%.4f\n",
                t, format(nrow(sub), big.mark = ","),
                min(sub$Score_AMR), max(sub$Score_AMR),
                mean(sub$PRS_Z), sd(sub$PRS_Z)))
}

# ================================================================
# GRÁFICA 1: Histograma por tertil de AMR (facetado)
# ================================================================
p1 <- ggplot(df, aes(x = PRS_Z, fill = AMR_tertile)) +
  geom_histogram(bins = 80, color = "white", linewidth = 0.3, alpha = 0.85) +
  facet_wrap(~AMR_tertile, ncol = 1, scales = "free_y") +

  scale_fill_manual(
    values = c("T1: Baja AMR" = "#4DBBD5",
               "T2: Media AMR" = "#E6AB02",
               "T3: Alta AMR"  = "#E64B35"),
    name = NULL
  ) +

  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +

  labs(
    title = "PRS Diabetes T2 estratificado por Ancestr\u00eda Amerindia (AMR)",
    subtitle = sprintf("Cohorte MCPS (n = %s) | GWAS propio | Tertiles de Score_AMR",
      format(nrow(df), big.mark = ",")),
    x = "PRS Estandarizado (z-score)",
    y = "N\u00famero de Pacientes",
    caption = "T1: menor proporci\u00f3n AMR | T3: mayor proporci\u00f3n AMR"
  ) +

  theme_minimal(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(color = "gray40", size = 11),
    plot.caption  = element_text(color = "gray50", size = 9),
    plot.margin   = margin(15, 20, 15, 15),
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 12),
    panel.grid.minor = element_blank()
  ) +
  scale_y_continuous(labels = comma)

ggsave(file.path(out_dir, "prs_mcps_gwas_amr_tertiles.png"), p1, width = 12, height = 9, dpi = 200)
ggsave(file.path(out_dir, "prs_mcps_gwas_amr_tertiles.pdf"), p1, width = 12, height = 9)

# ================================================================
# GRÁFICA 2: Overlay de densidades por tertil
# ================================================================
p2 <- ggplot(df, aes(x = PRS_Z, fill = AMR_tertile, color = AMR_tertile)) +
  geom_density(alpha = 0.3, linewidth = 0.8) +

  scale_fill_manual(
    values = c("T1: Baja AMR" = "#4DBBD5",
               "T2: Media AMR" = "#E6AB02",
               "T3: Alta AMR"  = "#E64B35"),
    name = "Tertil AMR"
  ) +
  scale_color_manual(
    values = c("T1: Baja AMR" = "#4DBBD5",
               "T2: Media AMR" = "#E6AB02",
               "T3: Alta AMR"  = "#E64B35"),
    name = "Tertil AMR"
  ) +

  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +

  labs(
    title = "Densidad del PRS por Tertil de Ancestr\u00eda Amerindia",
    subtitle = "Comparaci\u00f3n directa de distribuciones",
    x = "PRS Estandarizado (z-score)",
    y = "Densidad"
  ) +

  theme_minimal(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(color = "gray40", size = 11),
    plot.margin   = margin(15, 20, 15, 15),
    legend.position = c(0.85, 0.85),
    legend.background = element_rect(fill = "white", color = "gray80", linewidth = 0.3),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(out_dir, "prs_mcps_gwas_amr_density.png"), p2, width = 10, height = 6, dpi = 200)
ggsave(file.path(out_dir, "prs_mcps_gwas_amr_density.pdf"), p2, width = 10, height = 6)

# ================================================================
# GRÁFICA 3: Boxplot PRS por tertil AMR
# ================================================================
p3 <- ggplot(df, aes(x = AMR_tertile, y = PRS_Z, fill = AMR_tertile)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.3, outlier.alpha = 0.3) +

  scale_fill_manual(
    values = c("T1: Baja AMR" = "#4DBBD5",
               "T2: Media AMR" = "#E6AB02",
               "T3: Alta AMR"  = "#E64B35"),
    name = NULL
  ) +

  labs(
    title = "PRS Diabetes T2 por Tertil de Ancestr\u00eda AMR",
    subtitle = "Distribuci\u00f3n del z-score por grupo de ancestr\u00eda amerindia",
    x = "",
    y = "PRS Estandarizado (z-score)"
  ) +

  theme_minimal(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(color = "gray40", size = 11),
    plot.margin   = margin(15, 20, 15, 15),
    legend.position = "none",
    panel.grid.minor = element_blank()
  )

ggsave(file.path(out_dir, "prs_mcps_gwas_amr_boxplot.png"), p3, width = 8, height = 6, dpi = 200)
ggsave(file.path(out_dir, "prs_mcps_gwas_amr_boxplot.pdf"), p3, width = 8, height = 6)

# ================================================================
# GRÁFICA 4: Scatter PRS vs Score_AMR
# ================================================================
p4 <- ggplot(df, aes(x = Score_AMR, y = PRS_Z)) +
  geom_point(alpha = 0.05, size = 0.3, color = "#2196F3") +
  geom_smooth(method = "lm", color = "#E64B35", linewidth = 1, se = TRUE) +

  labs(
    title = "Relaci\u00f3n entre Ancestr\u00eda AMR y PRS Diabetes T2",
    subtitle = sprintf("r = %.4f | Cohorte MCPS (n = %s)",
      cor(df$Score_AMR, df$PRS_Z), format(nrow(df), big.mark = ",")),
    x = "Proporci\u00f3n de Ancestr\u00eda Amerindia (Score_AMR)",
    y = "PRS Estandarizado (z-score)"
  ) +

  theme_minimal(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(color = "gray40", size = 11),
    plot.margin   = margin(15, 20, 15, 15),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(out_dir, "prs_mcps_gwas_amr_scatter.png"), p4, width = 10, height = 7, dpi = 200)
ggsave(file.path(out_dir, "prs_mcps_gwas_amr_scatter.pdf"), p4, width = 10, height = 7)

cat("\nGr\u00e1ficas guardadas:\n")
cat(sprintf("  %s/prs_mcps_gwas_amr_tertiles.png\n", out_dir))
cat(sprintf("  %s/prs_mcps_gwas_amr_density.png\n", out_dir))
cat(sprintf("  %s/prs_mcps_gwas_amr_boxplot.png\n", out_dir))
cat(sprintf("  %s/prs_mcps_gwas_amr_scatter.png\n", out_dir))
cat("Completado.\n")
EOF

apptainer exec \
    --bind /mnt/cephfs:/mnt/cephfs \
    --bind "$TMPDIR:/tmp,$TMPDIR:/var/tmp" \
    "$SIF" \
    /opt/miniforge3/envs/r-bio/bin/Rscript "$R_SCRIPT"

rm -rf "$TMPDIR"
echo "==== COMPLETADO ===="
