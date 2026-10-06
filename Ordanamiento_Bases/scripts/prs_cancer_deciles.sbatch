#!/bin/bash
#SBATCH --job-name=prs_cancer_deciles
#SBATCH --account=researchers
#SBATCH --qos=vip
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=/mnt/cephfs/orgs/home/angel.pacheco/logs/prs_cancer_deciles_%j.out
#SBATCH --error=/mnt/cephfs/orgs/home/angel.pacheco/logs/prs_cancer_deciles_%j.err

set -euo pipefail

SIF="/mnt/cephfs/biocontainers/images/jupyter-biotools-1.5.sif"

TMPDIR="/mnt/cephfs/scratch/$USER/prs_deciles_tmp/$SLURM_JOB_ID"
mkdir -p "$TMPDIR"
chmod 700 "$TMPDIR"
export TMPDIR TMP="$TMPDIR" TEMP="$TMPDIR"

PRS_DIR="/mnt/cephfs/orgs/home/angel.pacheco/prs_cancer"
CC_FILE="/mnt/cephfs/orgs/home/angel.pacheco/Ordanamiento_Bases/MCPS_Cohorte_Casos_Controles_Oncologicos_F145K.csv"
OUT_DIR="$PRS_DIR"

echo "================================================================"
echo "FASE F+G: INTEGRACIÓN PRS + CASE-CONTROL Y GRÁFICAS POR DECIL"
echo "================================================================"

R_SCRIPT="$TMPDIR/prs_cancer_deciles.R"
cat << 'EOF' > "$R_SCRIPT"
library(ggplot2)
library(scales)
library(dplyr)

prs_dir  <- "/mnt/cephfs/orgs/home/angel.pacheco/prs_cancer"
cc_file  <- "/mnt/cephfs/orgs/home/angel.pacheco/Ordanamiento_Bases/MCPS_Cohorte_Casos_Controles_Oncologicos_F145K.csv"
out_dir  <- prs_dir

# ================================================================
# FASE F: INTEGRACIÓN PRS + CASO-CONTROL
# ================================================================
cat("Cargando PRS calculado...\n")
df_prs <- read.delim(file.path(prs_dir, "prs_cancer_mcps_all_samples.tsv"),
                     na.strings = "NA")

cat("Cargando cohorte caso-control oncológica...\n")
df_cc  <- read.csv(cc_file)

cat(sprintf("PRS: %s muestras\n", format(nrow(df_prs), big.mark = ",")))
cat(sprintf("Caso-Control: %s participantes\n", format(nrow(df_cc), big.mark = ",")))

# Join por IID_EXOME = IID
df <- merge(df_cc, df_prs[, c("IID",
                                "PRS_PGS003766_PROSTATE", "PRS_Z_PGS003766_PROSTATE",
                                "PRS_PGS000028_BREAST",   "PRS_Z_PGS000028_BREAST",
                                "PRS_PGS000078_LUNG",     "PRS_Z_PGS000078_LUNG",
                                "PRS_PGS000073_CERVICAL", "PRS_Z_PGS000073_CERVICAL")],
            by.x = "IID_EXOME", by.y = "IID", all.x = FALSE)

cat(sprintf("Participantes con PRS + fenotipo: %s\n", format(nrow(df), big.mark = ",")))

# Guardar tabla integrada
write.table(df, file.path(out_dir, "prs_cancer_integrated.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
cat("Tabla integrada guardada.\n\n")

# ================================================================
# FASE G: GRÁFICAS POR DECIL
# Configuración: 4 scores × contexto de evaluación
# Estrategia: para cada PGS, evaluar el score en:
#   (a) casos del sitio específico vs todos los controles
#   (b) casos del sitio específico vs controles del mismo sexo
# ================================================================

SCORE_CONFIGS <- list(
  list(
    pgs_id   = "PGS003766",
    col_raw  = "PRS_PGS003766_PROSTATE",
    col_z    = "PRS_Z_PGS003766_PROSTATE",
    site     = "PROSTATE",
    sex_val  = 1,          # 1=hombres
    sex_lab  = "Hombres",
    color    = "#1565C0",
    title    = "PRS Cáncer de Próstata (PGS003766)",
    subtitle = "Conti et al. 2021 | 451 SNPs | Cohorte MCPS Hombres"
  ),
  list(
    pgs_id   = "PGS000028",
    col_raw  = "PRS_PGS000028_BREAST",
    col_z    = "PRS_Z_PGS000028_BREAST",
    site     = "BREAST",
    sex_val  = 0,          # 0=mujeres
    sex_lab  = "Mujeres",
    color    = "#C62828",
    title    = "PRS Cáncer de Mama (PGS000028)",
    subtitle = "Mavaddat et al. 2019 | 83 SNPs | Cohorte MCPS Mujeres"
  ),
  list(
    pgs_id   = "PGS000078",
    col_raw  = "PRS_PGS000078_LUNG",
    col_z    = "PRS_Z_PGS000078_LUNG",
    site     = "LUNG",
    sex_val  = NA,         # ambos sexos
    sex_lab  = "Ambos Sexos",
    color    = "#2E7D32",
    title    = "PRS Cáncer de Pulmón (PGS000078)",
    subtitle = "McKay et al. 2017 | 109 SNPs | Cohorte MCPS"
  ),
  list(
    pgs_id   = "PGS000073",
    col_raw  = "PRS_PGS000073_CERVICAL",
    col_z    = "PRS_Z_PGS000073_CERVICAL",
    site     = "CERVICAL",
    sex_val  = 0,
    sex_lab  = "Mujeres",
    color    = "#6A1B9A",
    title    = "PRS Cáncer Cervical (PGS000073)",
    subtitle = "Chen et al. 2024 | 10 SNPs | Cohorte MCPS Mujeres"
  )
)

# Función: generar barplot de % casos por decil
plot_decile_barplot <- function(df_sub, col_z, site, pgs_id, color, title, subtitle, suffix = "") {
  # Asignar decil dentro de la subcohorte (1=menor riesgo, 10=mayor)
  df_sub <- df_sub[!is.na(df_sub[[col_z]]), ]
  df_sub$DECIL <- cut(df_sub[[col_z]],
                      breaks = quantile(df_sub[[col_z]], probs = seq(0, 1, 0.1), na.rm = TRUE),
                      include.lowest = TRUE,
                      labels = 1:10)

  # % casos en cada decil
  df_decil <- df_sub %>%
    group_by(DECIL) %>%
    summarise(
      n_total = n(),
      n_casos = sum(PHENOTYPE_CANCER_STATUS == "CASO", na.rm = TRUE),
      n_casos_site = sum(PHENOTYPE_CANCER_STATUS == "CASO" & CANCER_SITE == site, na.rm = TRUE),
      pct_casos      = n_casos      / n_total * 100,
      pct_casos_site = n_casos_site / n_total * 100,
      .groups = "drop"
    )

  # Referencia: % promedio del decil 5 (intermedio)
  ref_all  <- df_decil$pct_casos[df_decil$DECIL == 5]
  ref_site <- df_decil$pct_casos_site[df_decil$DECIL == 5]

  cat(sprintf("\n%s:\n", title))
  print(df_decil)

  # Guardar tabla de deciles
  write.table(df_decil,
    file.path(out_dir, sprintf("prs_deciles_%s_%s%s.tsv", pgs_id, site, suffix)),
    sep = "\t", row.names = FALSE, quote = FALSE)

  # ---- GRÁFICA: % CUALQUIER CÁNCER por decil ----
  p1 <- ggplot(df_decil, aes(x = DECIL, y = pct_casos)) +
    geom_col(fill = color, alpha = 0.8, width = 0.7) +
    geom_hline(yintercept = mean(df_decil$pct_casos), linetype = "dashed",
               color = "black", linewidth = 0.7) +
    geom_text(aes(label = sprintf("%.1f%%", pct_casos)),
              vjust = -0.4, size = 3.2, color = "gray20") +
    labs(
      title    = title,
      subtitle = subtitle,
      x        = "Decil de PRS (1 = menor riesgo, 10 = mayor riesgo)",
      y        = "% Cualquier Cáncer",
      caption  = sprintf("Línea punteada = promedio de la cohorte (%.1f%%)",
                          mean(df_decil$pct_casos))
    ) +
    scale_x_discrete(drop = FALSE) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title    = element_text(face = "bold", size = 15),
      plot.subtitle = element_text(color = "gray40", size = 11),
      plot.caption  = element_text(color = "gray50", size = 9),
      plot.margin   = margin(15, 20, 15, 15),
      panel.grid.minor = element_blank()
    )

  ggsave(file.path(out_dir, sprintf("prs_decil_%s_cualquier_cancer%s.png", pgs_id, suffix)),
         p1, width = 12, height = 7, dpi = 200)
  ggsave(file.path(out_dir, sprintf("prs_decil_%s_cualquier_cancer%s.pdf", pgs_id, suffix)),
         p1, width = 12, height = 7)

  # ---- GRÁFICA: % CÁNCER DEL SITIO ESPECÍFICO por decil ----
  p2 <- ggplot(df_decil, aes(x = DECIL, y = pct_casos_site)) +
    geom_col(fill = color, alpha = 0.85, width = 0.7) +
    geom_hline(yintercept = mean(df_decil$pct_casos_site), linetype = "dashed",
               color = "black", linewidth = 0.7) +
    geom_text(aes(label = sprintf("%.2f%%", pct_casos_site)),
              vjust = -0.4, size = 3.2, color = "gray20") +
    labs(
      title    = paste0(title, " — Casos de sitio específico"),
      subtitle = subtitle,
      x        = "Decil de PRS (1 = menor riesgo, 10 = mayor riesgo)",
      y        = sprintf("%% Casos (%s)", site),
      caption  = sprintf("Línea punteada = prevalencia promedio (%.2f%%)",
                          mean(df_decil$pct_casos_site))
    ) +
    scale_x_discrete(drop = FALSE) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title    = element_text(face = "bold", size = 15),
      plot.subtitle = element_text(color = "gray40", size = 11),
      plot.caption  = element_text(color = "gray50", size = 9),
      plot.margin   = margin(15, 20, 15, 15),
      panel.grid.minor = element_blank()
    )

  ggsave(file.path(out_dir, sprintf("prs_decil_%s_sitio_especifico%s.png", pgs_id, suffix)),
         p2, width = 12, height = 7, dpi = 200)
  ggsave(file.path(out_dir, sprintf("prs_decil_%s_sitio_especifico%s.pdf", pgs_id, suffix)),
         p2, width = 12, height = 7)

  invisible(df_decil)
}

# ================================================================
# GENERAR GRÁFICAS PARA CADA SCORE
# ================================================================
for (cfg in SCORE_CONFIGS) {
  cat(sprintf("\n>> Procesando %s (%s)...\n", cfg$pgs_id, cfg$site))

  # Filtrar subcohorte según sexo del score
  if (!is.na(cfg$sex_val)) {
    df_sub <- df[df$MALE == cfg$sex_val & !is.na(df$MALE), ]
  } else {
    df_sub <- df[!is.na(df$MALE), ]
  }

  # Excluir muestras sin PRS calculado (NA — sexo incompatible en scoring)
  df_sub <- df_sub[!is.na(df_sub[[cfg$col_z]]), ]

  cat(sprintf("   Subcohorte: %s participantes (%s)\n",
              format(nrow(df_sub), big.mark = ","), cfg$sex_lab))
  cat(sprintf("   Casos (%s): %s | Controles: %s\n",
              cfg$site,
              format(sum(df_sub$CANCER_SITE == cfg$site, na.rm = TRUE), big.mark = ","),
              format(sum(df_sub$PHENOTYPE_CANCER_STATUS == "CONTROL", na.rm = TRUE), big.mark = ",")))

  if (nrow(df_sub) < 100) {
    cat("   ADVERTENCIA: Subcohorte muy pequeña, omitiendo gráfica.\n")
    next
  }

  plot_decile_barplot(
    df_sub   = df_sub,
    col_z    = cfg$col_z,
    site     = cfg$site,
    pgs_id   = cfg$pgs_id,
    color    = cfg$color,
    title    = cfg$title,
    subtitle = cfg$subtitle
  )
}

# ================================================================
# GRÁFICA COMBINADA: 4 scores en una figura (decil vs OR relativo)
# ================================================================
cat("\n>> Generando figura combinada de ORs por decil...\n")

combined_rows <- list()
for (cfg in SCORE_CONFIGS) {
  if (!is.na(cfg$sex_val)) {
    df_sub <- df[df$MALE == cfg$sex_val & !is.na(df$MALE), ]
  } else {
    df_sub <- df[!is.na(df$MALE), ]
  }
  df_sub <- df_sub[!is.na(df_sub[[cfg$col_z]]), ]
  if (nrow(df_sub) < 100) next

  df_sub$DECIL <- cut(df_sub[[cfg$col_z]],
                      breaks = quantile(df_sub[[cfg$col_z]], probs = seq(0, 1, 0.1), na.rm = TRUE),
                      include.lowest = TRUE, labels = 1:10)

  df_d <- df_sub %>%
    group_by(DECIL) %>%
    summarise(
      n_casos_site = sum(CANCER_SITE == cfg$site, na.rm = TRUE),
      n_total = n(),
      .groups = "drop"
    ) %>%
    mutate(
      prev = n_casos_site / n_total,
      score = cfg$pgs_id,
      site  = cfg$site,
      color = cfg$color
    )

  # OR relativo al decil 5
  ref_prev <- df_d$prev[df_d$DECIL == 5]
  if (length(ref_prev) == 0 || ref_prev == 0) next
  df_d$OR_rel <- df_d$prev / ref_prev

  combined_rows[[length(combined_rows) + 1]] <- df_d
}

if (length(combined_rows) > 0) {
  df_combined <- bind_rows(combined_rows)
  df_combined$DECIL <- as.integer(as.character(df_combined$DECIL))
  df_combined$score_label <- factor(df_combined$score,
    levels = c("PGS003766","PGS000028","PGS000078","PGS000073"),
    labels = c("Próstata\n(PGS003766)","Mama\n(PGS000028)",
               "Pulmón\n(PGS000078)","Cervical\n(PGS000073)"))

  color_map <- c("PGS003766"="#1565C0","PGS000028"="#C62828",
                 "PGS000078"="#2E7D32","PGS000073"="#6A1B9A")

  p_comb <- ggplot(df_combined, aes(x = DECIL, y = OR_rel,
                                     color = score, group = score)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 3) +
    geom_hline(yintercept = 1.0, linetype = "dashed", color = "gray50", linewidth = 0.6) +
    scale_color_manual(
      values = color_map,
      labels = c("Próstata (PGS003766)","Mama (PGS000028)",
                 "Pulmón (PGS000078)","Cervical (PGS000073)"),
      name = "Cáncer / Score"
    ) +
    scale_x_continuous(breaks = 1:10) +
    labs(
      title    = "OR Relativo por Decil de PRS — 4 Cánceres (MCPS)",
      subtitle = "Referencia: decil 5 | PGS Catalog betas | Cohorte MCPS 140K",
      x        = "Decil de PRS (1 = menor riesgo, 10 = mayor riesgo)",
      y        = "OR relativo al decil 5",
      caption  = "PGS003766: Próstata (hombres) | PGS000028: Mama (mujeres) | PGS000078: Pulmón (todos) | PGS000073: Cervical (mujeres)"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title    = element_text(face = "bold", size = 15),
      plot.subtitle = element_text(color = "gray40", size = 11),
      plot.caption  = element_text(color = "gray50", size = 8.5),
      plot.margin   = margin(15, 20, 15, 15),
      legend.position = "right",
      legend.background = element_rect(fill = "white", color = "gray80", linewidth = 0.3),
      panel.grid.minor = element_blank()
    )

  ggsave(file.path(out_dir, "prs_cancer_combined_OR_deciles.png"),
         p_comb, width = 13, height = 7, dpi = 200)
  ggsave(file.path(out_dir, "prs_cancer_combined_OR_deciles.pdf"),
         p_comb, width = 13, height = 7)
  cat("Figura combinada guardada.\n")
}

cat("\n==== FASE F+G COMPLETADAS ====\n")
EOF

apptainer exec \
    --bind /mnt/cephfs:/mnt/cephfs \
    --bind "$TMPDIR:/tmp,$TMPDIR:/var/tmp" \
    "$SIF" \
    /opt/miniforge3/envs/r-bio/bin/Rscript "$R_SCRIPT"

rm -rf "$TMPDIR"
echo "[$(date +'%H:%M:%S')] ==== JOB COMPLETADO ===="
