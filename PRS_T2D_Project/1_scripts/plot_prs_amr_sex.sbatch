#!/bin/bash
#SBATCH --job-name=prs_amr_sex
#SBATCH --account=researchers
#SBATCH --qos=vip
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:15:00
#SBATCH --output=/mnt/cephfs/orgs/home/angel.pacheco/logs/prs_amr_sex_%j.out
#SBATCH --error=/mnt/cephfs/orgs/home/angel.pacheco/logs/prs_amr_sex_%j.err

set -euo pipefail

SIF="/mnt/cephfs/biocontainers/images/jupyter-biotools-1.2.sif"

TMPDIR="/mnt/cephfs/scratch/$USER/hist_tmp/$SLURM_JOB_ID"
mkdir -p "$TMPDIR"
chmod 700 "$TMPDIR"
export TMPDIR TMP="$TMPDIR" TEMP="$TMPDIR"

echo "==== HISTOGRAMA PRS POR TERTILES AMR × SEXO ===="

R_SCRIPT="$TMPDIR/prs_amr_sex.R"
cat << 'EOF' > "$R_SCRIPT"
library(ggplot2)
library(scales)

out_dir <- "/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes"

# ================================================================
# Cargar y cruzar datos
# ================================================================
cat("Cargando datos...\n")
df_prs <- read.delim(file.path(out_dir, "prs_diabetes_mcps_gwas_final.tsv"))
df_anc <- read.delim(file.path(out_dir, "ancestry-proportions-k4.txt"))
df_anc <- df_anc[df_anc$population == "MCPS", ]

df <- merge(df_prs, df_anc[, c("IID", "Score_AMR", "Score_EUR", "Score_AFR", "Score_EAS")],
            by = "IID", all.x = FALSE)

# Etiquetar sexo
df$SEX_LABEL <- ifelse(df$SEX == 0, "Mujeres", "Hombres")

# Tertiles de AMR
df$AMR_tertile <- cut(df$Score_AMR,
                      breaks = quantile(df$Score_AMR, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE),
                      include.lowest = TRUE,
                      labels = c("T1: Baja AMR", "T2: Media AMR", "T3: Alta AMR"))

cat(sprintf("Pacientes cruzados: %s\n", format(nrow(df), big.mark = ",")))
cat(sprintf("Mujeres: %s | Hombres: %s\n",
    format(sum(df$SEX == 0), big.mark = ","),
    format(sum(df$SEX == 1), big.mark = ",")))

# Estadísticas por sexo × tertil
cat("\nEstad\u00edsticas por Sexo x Tertil AMR:\n")
for (s in c("Mujeres", "Hombres")) {
    cat(sprintf("\n  === %s ===\n", s))
    for (t in levels(df$AMR_tertile)) {
        sub <- df[df$SEX_LABEL == s & df$AMR_tertile == t, ]
        cat(sprintf("    %s (n=%s): AMR=[%.2f-%.2f], PRS_Z mean=%.4f, sd=%.4f\n",
                    t, format(nrow(sub), big.mark = ","),
                    min(sub$Score_AMR), max(sub$Score_AMR),
                    mean(sub$PRS_Z), sd(sub$PRS_Z)))
    }
}

# Colores
amr_colors <- c("T1: Baja AMR" = "#4DBBD5",
                "T2: Media AMR" = "#E6AB02",
                "T3: Alta AMR"  = "#E64B35")

# ================================================================
# GRÁFICA 1: Histograma facetado Sexo (columnas) × Tertil AMR (filas)
# ================================================================
p1 <- ggplot(df, aes(x = PRS_Z, fill = AMR_tertile)) +
  geom_histogram(bins = 80, color = "white", linewidth = 0.3, alpha = 0.85) +
  facet_grid(AMR_tertile ~ SEX_LABEL, scales = "free_y") +

  scale_fill_manual(values = amr_colors, name = NULL) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +

  labs(
    title = "PRS Diabetes T2 por Sexo y Tertil de Ancestr\u00eda Amerindia",
    subtitle = sprintf("Cohorte MCPS (n = %s) | GWAS propio estratificado por sexo | M\u00e9todo C+T",
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
    strip.text = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank()
  ) +
  scale_y_continuous(labels = comma)

ggsave(file.path(out_dir, "prs_mcps_gwas_amr_sex_grid.png"), p1, width = 14, height = 10, dpi = 200)
ggsave(file.path(out_dir, "prs_mcps_gwas_amr_sex_grid.pdf"), p1, width = 14, height = 10)

# ================================================================
# GRÁFICA 2: Densidad overlay por tertil, facetado por sexo
# ================================================================
p2 <- ggplot(df, aes(x = PRS_Z, fill = AMR_tertile, color = AMR_tertile)) +
  geom_density(alpha = 0.25, linewidth = 0.8) +
  facet_wrap(~SEX_LABEL, ncol = 1) +

  scale_fill_manual(values = amr_colors, name = "Tertil AMR") +
  scale_color_manual(values = amr_colors, name = "Tertil AMR") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +

  labs(
    title = "Densidad PRS por Tertil AMR y Sexo",
    subtitle = "Comparaci\u00f3n directa de distribuciones por grupo",
    x = "PRS Estandarizado (z-score)",
    y = "Densidad"
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
  )

ggsave(file.path(out_dir, "prs_mcps_gwas_amr_sex_density.png"), p2, width = 10, height = 8, dpi = 200)
ggsave(file.path(out_dir, "prs_mcps_gwas_amr_sex_density.pdf"), p2, width = 10, height = 8)

# ================================================================
# GRÁFICA 3: Boxplot sexo × tertil AMR
# ================================================================
p3 <- ggplot(df, aes(x = AMR_tertile, y = PRS_Z, fill = AMR_tertile)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.3, outlier.alpha = 0.3) +
  facet_wrap(~SEX_LABEL) +

  scale_fill_manual(values = amr_colors, name = NULL) +

  labs(
    title = "PRS Diabetes T2 por Tertil AMR y Sexo",
    subtitle = "Distribuci\u00f3n del z-score por grupo",
    x = "",
    y = "PRS Estandarizado (z-score)"
  ) +

  theme_minimal(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(color = "gray40", size = 11),
    plot.margin   = margin(15, 20, 15, 15),
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 12),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 9)
  )

ggsave(file.path(out_dir, "prs_mcps_gwas_amr_sex_boxplot.png"), p3, width = 12, height = 6, dpi = 200)
ggsave(file.path(out_dir, "prs_mcps_gwas_amr_sex_boxplot.pdf"), p3, width = 12, height = 6)

# ================================================================
# GRÁFICA 4: Scatter PRS vs AMR, facetado por sexo
# ================================================================
p4 <- ggplot(df, aes(x = Score_AMR, y = PRS_Z)) +
  geom_point(alpha = 0.05, size = 0.3, color = "#2196F3") +
  geom_smooth(method = "lm", color = "#E64B35", linewidth = 1, se = TRUE) +
  facet_wrap(~SEX_LABEL) +

  labs(
    title = "Relaci\u00f3n PRS vs Ancestr\u00eda AMR por Sexo",
    subtitle = sprintf("Cohorte MCPS (n = %s) | r Mujeres = %.4f | r Hombres = %.4f",
      format(nrow(df), big.mark = ","),
      cor(df$PRS_Z[df$SEX == 0], df$Score_AMR[df$SEX == 0]),
      cor(df$PRS_Z[df$SEX == 1], df$Score_AMR[df$SEX == 1])),
    x = "Proporci\u00f3n de Ancestr\u00eda Amerindia (Score_AMR)",
    y = "PRS Estandarizado (z-score)"
  ) +

  theme_minimal(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(color = "gray40", size = 11),
    plot.margin   = margin(15, 20, 15, 15),
    strip.text = element_text(face = "bold", size = 12),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(out_dir, "prs_mcps_gwas_amr_sex_scatter.png"), p4, width = 12, height = 6, dpi = 200)
ggsave(file.path(out_dir, "prs_mcps_gwas_amr_sex_scatter.pdf"), p4, width = 12, height = 6)

cat("\nGr\u00e1ficas guardadas:\n")
cat(sprintf("  %s/prs_mcps_gwas_amr_sex_grid.png\n", out_dir))
cat(sprintf("  %s/prs_mcps_gwas_amr_sex_density.png\n", out_dir))
cat(sprintf("  %s/prs_mcps_gwas_amr_sex_boxplot.png\n", out_dir))
cat(sprintf("  %s/prs_mcps_gwas_amr_sex_scatter.png\n", out_dir))
cat("Completado.\n")
EOF

apptainer exec \
    --bind /mnt/cephfs:/mnt/cephfs \
    --bind "$TMPDIR:/tmp,$TMPDIR:/var/tmp" \
    "$SIF" \
    /opt/miniforge3/envs/r-bio/bin/Rscript "$R_SCRIPT"

rm -rf "$TMPDIR"
echo "==== COMPLETADO ===="
