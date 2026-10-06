#!/bin/bash
#SBATCH --job-name=sbrc_compare
#SBATCH --account=researchers
#SBATCH --qos=vip
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --exclude=mcps-hpc-01
#SBATCH --output=/mnt/cephfs/orgs/home/angel.pacheco/logs/sbrc_compare_%j.out
#SBATCH --error=/mnt/cephfs/orgs/home/angel.pacheco/logs/sbrc_compare_%j.err

set -euo pipefail

SIF="/mnt/cephfs/biocontainers/images/jupyter-biotools-1.4.sif"
TMPDIR="/mnt/cephfs/scratch/$USER/sbrc_cmp_tmp/$SLURM_JOB_ID"
mkdir -p "$TMPDIR"
chmod 700 "$TMPDIR"
export TMPDIR TMP="$TMPDIR" TEMP="$TMPDIR"

echo "============================================================"
echo "  COMPARACIÓN SBayesRC: T2DGGI (Pipeline A) vs MCPS (A2)"
echo "  Nodo: $SLURMD_NODENAME"
echo "  $(date)"
echo "============================================================"

R_SCRIPT="$TMPDIR/compare_sbrc.R"
cat << 'EOF' > "$R_SCRIPT"
suppressMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
  library(patchwork)
})

out_dir <- "/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes/sbayesrc"
fig_dir <- file.path(out_dir, "figures_comparison")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ================================================================
# 1. CARGAR BETAS DE AMBOS MODELOS
# ================================================================
cat(">> Cargando resultados SBayesRC...\n")

# Pipeline A — T2DGGI (consorcio global, N~1.4M)
f_t2dggi <- file.path(out_dir, "results/sbayesrc.txt")
# Pipeline A2 — MCPS (GWAS propio, N=78,304)
f_mcps   <- file.path(out_dir, "results_mcps/sbayesrc_mcps.txt")

stopifnot(file.exists(f_t2dggi), file.exists(f_mcps))

df_t2d <- fread(f_t2dggi, select = c("SNP", "A1", "BETA", "PIP"))
df_mcp <- fread(f_mcps,   select = c("SNP", "A1", "BETA", "PIP"))

setnames(df_t2d, c("BETA", "PIP"), c("BETA_T2DGGI", "PIP_T2DGGI"))
setnames(df_mcp, c("BETA", "PIP"), c("BETA_MCPS", "PIP_MCPS"))

cat(sprintf("  T2DGGI: %s variantes\n", format(nrow(df_t2d), big.mark = ",")))
cat(sprintf("  MCPS:   %s variantes\n", format(nrow(df_mcp), big.mark = ",")))

# Merge por SNP + A1 (mismo effect allele)
df <- merge(df_t2d, df_mcp[, .(SNP, A1, BETA_MCPS, PIP_MCPS)],
            by = c("SNP", "A1"), all = FALSE)
cat(sprintf("  Variantes en común (mismo SNP + A1): %s\n",
            format(nrow(df), big.mark = ",")))

# Clasificar variantes
df[, `:=`(
  nz_t2d = (BETA_T2DGGI != 0),
  nz_mcp = (BETA_MCPS   != 0)
)]

cat(sprintf("  BETA != 0 en T2DGGI: %s (%.1f%%)\n",
    format(sum(df$nz_t2d), big.mark = ","), sum(df$nz_t2d)/nrow(df)*100))
cat(sprintf("  BETA != 0 en MCPS:   %s (%.1f%%)\n",
    format(sum(df$nz_mcp), big.mark = ","), sum(df$nz_mcp)/nrow(df)*100))

# ================================================================
# 2. PARÁMETROS MCMC (hardcoded de los logs)
# ================================================================
mcmc_params <- data.table(
  Parameter = rep(c("h²_SNP", "nnz", "σ²", "Var_e", "Var_g",
                     "NumSnp (cat 3)", "NumSnp (cat 4)", "NumSnp (cat 5)",
                     "Vg (cat 3)", "Vg (cat 4)", "Vg (cat 5)"), each = 2),
  GWAS = rep(c("T2DGGI", "MCPS"), 11),
  Mean = c(
    # h²_SNP
    0.1249, 0.1690,
    # nnz
    9228.675, 1232.505,
    # σ²
    9.515e-04, 1.556e-03,
    # Var_e
    1.2405, 1.0000,
    # Var_g
    0.1249, 0.1690,
    # NumSnp3
    8484.71, 13.80,
    # NumSnp4
    329.06, 886.91,
    # NumSnp5
    12.43, 19.14,
    # Vg3
    0.666, 0.001,
    # Vg4
    0.261, 0.834,
    # Vg5
    0.070, 0.161
  ),
  SD = c(
    # h²_SNP
    0.00144, 0.00470,
    # nnz
    339.54, 243.71,
    # σ²
    5.00e-05, 1.14e-04,
    # Var_e
    0.00666, 3.89e-05,
    # Var_g
    0.00144, 0.00470,
    # NumSnp3
    207.60, 13.64,
    # NumSnp4
    35.34, 43.41,
    # NumSnp5
    5.56, 7.28,
    # Vg3
    0.0136, 0.0013,
    # Vg4
    0.0219, 0.0335,
    # Vg5
    0.0174, 0.0333
  )
)

mcmc_params[, GWAS := factor(GWAS, levels = c("T2DGGI", "MCPS"))]

# Colores consistentes con tu estilo
col_t2d <- "#4DBBD5"   # Azul — T2DGGI
col_mcp <- "#E64B35"   # Rojo — MCPS

# ================================================================
# GRÁFICA 1: Heritabilidad y Varianza por Componente (barras)
# ================================================================
cat("\n>> Gráfica 1: Arquitectura genética comparada...\n")

df_arch <- mcmc_params[Parameter %in% c("Vg (cat 3)", "Vg (cat 4)", "Vg (cat 5)")]
df_arch[, Component := factor(
  gsub("Vg \\(cat ([345])\\)", "Componente \\1", Parameter),
  levels = c("Componente 3", "Componente 4", "Componente 5")
)]
# Etiquetar componentes con su significado
df_arch[, Component_Label := factor(
  ifelse(Component == "Componente 3", "Comp 3\n(efecto pequeño)",
  ifelse(Component == "Componente 4", "Comp 4\n(efecto medio)",
                                      "Comp 5\n(efecto grande)")),
  levels = c("Comp 3\n(efecto pequeño)", "Comp 4\n(efecto medio)", "Comp 5\n(efecto grande)")
)]

p1 <- ggplot(df_arch, aes(x = Component_Label, y = Mean, fill = GWAS)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.75),
           width = 0.65, color = "white", linewidth = 0.4) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD),
                position = position_dodge(width = 0.75), width = 0.2, linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.3f", Mean), y = Mean + SD + 0.02),
            position = position_dodge(width = 0.75), size = 3.5, fontface = "bold") +

  scale_fill_manual(values = c("T2DGGI" = col_t2d, "MCPS" = col_mcp), name = "GWAS") +

  labs(
    title = "Arquitectura Genética de T2D — Proporción de Varianza por Componente",
    subtitle = "SBayesRC sin anotaciones | Distribución de Vg en 3 componentes de efecto",
    x = "",
    y = "Proporción de Varianza Genética (Vg)",
    caption = "Barras de error: ±1 SD de las cadenas MCMC post-burnin"
  ) +

  theme_minimal(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(color = "gray40", size = 11),
    plot.caption  = element_text(color = "gray50", size = 9),
    plot.margin   = margin(15, 20, 15, 15),
    legend.position = c(0.88, 0.88),
    legend.background = element_rect(fill = "white", color = "gray80", linewidth = 0.3),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(fig_dir, "01_vg_components_comparison.png"), p1, width = 11, height = 7, dpi = 300)
ggsave(file.path(fig_dir, "01_vg_components_comparison.pdf"), p1, width = 11, height = 7)

# ================================================================
# GRÁFICA 2: h² SNP y Var_e (paneles gemelos)
# ================================================================
cat(">> Gráfica 2: h² SNP y varianza residual...\n")

df_hsq <- mcmc_params[Parameter == "h²_SNP"]
df_ve  <- mcmc_params[Parameter == "Var_e"]

p2a <- ggplot(df_hsq, aes(x = GWAS, y = Mean, fill = GWAS)) +
  geom_bar(stat = "identity", width = 0.55, color = "white", linewidth = 0.4) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.15, linewidth = 0.7) +
  geom_text(aes(label = sprintf("%.4f ± %.4f", Mean, SD)), vjust = -1, size = 4, fontface = "bold") +
  scale_fill_manual(values = c("T2DGGI" = col_t2d, "MCPS" = col_mcp)) +
  labs(title = "Heredabilidad SNP (h²)", x = "", y = "h²_SNP") +
  ylim(0, max(df_hsq$Mean + df_hsq$SD) * 1.25) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

p2b <- ggplot(df_ve, aes(x = GWAS, y = Mean, fill = GWAS)) +
  geom_bar(stat = "identity", width = 0.55, color = "white", linewidth = 0.4) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.15, linewidth = 0.7) +
  geom_text(aes(label = sprintf("%.4f ± %.4f", Mean, SD)), vjust = -1, size = 4, fontface = "bold") +
  scale_fill_manual(values = c("T2DGGI" = col_t2d, "MCPS" = col_mcp)) +
  labs(title = "Varianza Residual (Var_e)", x = "", y = "Var_e") +
  ylim(0, max(df_ve$Mean + df_ve$SD) * 1.15) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

p2 <- p2a + p2b +
  plot_annotation(
    title = "Parámetros Globales SBayesRC — T2DGGI vs MCPS",
    subtitle = "Estimaciones posteriores MCMC (media ± SD, post-burnin)",
    theme = theme(
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(color = "gray40", size = 12)
    )
  )

ggsave(file.path(fig_dir, "02_hsq_vare_comparison.png"), p2, width = 12, height = 6, dpi = 300)
ggsave(file.path(fig_dir, "02_hsq_vare_comparison.pdf"), p2, width = 12, height = 6)

# ================================================================
# GRÁFICA 3: Distribución de Betas (densidad overlay)
# ================================================================
cat(">> Gráfica 3: Distribución de betas no-cero...\n")

# Solo variantes con beta != 0 en al menos un modelo
df_nz <- df[nz_t2d | nz_mcp]

# Para densidad: apilar los betas
df_dens <- rbind(
  data.table(BETA = df_nz[nz_t2d == TRUE, BETA_T2DGGI], GWAS = "T2DGGI"),
  data.table(BETA = df_nz[nz_mcp == TRUE,  BETA_MCPS],  GWAS = "MCPS")
)

# Winsorizar para visualización (quitar extremos >99.9th percentile)
q999_t <- quantile(abs(df_dens[GWAS == "T2DGGI", BETA]), 0.999)
q999_m <- quantile(abs(df_dens[GWAS == "MCPS",   BETA]), 0.999)
q_max  <- max(q999_t, q999_m)
df_dens_plot <- df_dens[abs(BETA) <= q_max]

p3 <- ggplot(df_dens_plot, aes(x = BETA, fill = GWAS, color = GWAS)) +
  geom_density(alpha = 0.3, linewidth = 0.8) +

  scale_fill_manual(values = c("T2DGGI" = col_t2d, "MCPS" = col_mcp), name = "GWAS") +
  scale_color_manual(values = c("T2DGGI" = col_t2d, "MCPS" = col_mcp), name = "GWAS") +

  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +

  labs(
    title = "Distribución de Efectos Posteriores (β ≠ 0)",
    subtitle = sprintf("T2DGGI: %s variantes | MCPS: %s variantes (winsorizado al 99.9° percentil)",
      format(sum(df$nz_t2d), big.mark = ","),
      format(sum(df$nz_mcp), big.mark = ",")),
    x = "Efecto Posterior (β)",
    y = "Densidad",
    caption = "SBayesRC sin anotaciones | Betas promediados post-burnin"
  ) +

  theme_minimal(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(color = "gray40", size = 11),
    plot.caption  = element_text(color = "gray50", size = 9),
    plot.margin   = margin(15, 20, 15, 15),
    legend.position = c(0.88, 0.85),
    legend.background = element_rect(fill = "white", color = "gray80", linewidth = 0.3),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(fig_dir, "03_beta_density_comparison.png"), p3, width = 11, height = 7, dpi = 300)
ggsave(file.path(fig_dir, "03_beta_density_comparison.pdf"), p3, width = 11, height = 7)

# ================================================================
# GRÁFICA 4: Scatter β_T2DGGI vs β_MCPS (variantes con señal)
# ================================================================
cat(">> Gráfica 4: Concordancia de betas (scatter)...\n")

# Usar variantes con beta != 0 en AMBOS modelos
df_both <- df[nz_t2d & nz_mcp]
cat(sprintf("  Variantes con β ≠ 0 en AMBOS: %s\n", format(nrow(df_both), big.mark = ",")))

# Pearson r
r_val <- cor(df_both$BETA_T2DGGI, df_both$BETA_MCPS, use = "complete.obs")
cat(sprintf("  Correlación Pearson: r = %.4f\n", r_val))

# Subsample para scatter (evitar overplotting)
set.seed(42)
n_plot <- min(nrow(df_both), 200000)
df_scatter <- df_both[sample(.N, n_plot)]

p4 <- ggplot(df_scatter, aes(x = BETA_T2DGGI, y = BETA_MCPS)) +
  geom_point(alpha = 0.03, size = 0.3, color = "#2196F3") +
  geom_smooth(method = "lm", color = "#E64B35", linewidth = 1, se = TRUE, fill = "#E64B3530") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.5) +
  geom_hline(yintercept = 0, color = "gray70", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "gray70", linewidth = 0.3) +

  labs(
    title = "Concordancia de Efectos Posteriores — T2DGGI vs MCPS",
    subtitle = sprintf("r = %.4f | %s variantes con β ≠ 0 en ambos modelos (muestreadas %s para visualización)",
      r_val, format(nrow(df_both), big.mark = ","), format(n_plot, big.mark = ",")),
    x = "β SBayesRC — T2DGGI (N ~ 1.4M)",
    y = "β SBayesRC — MCPS (N = 78,304)",
    caption = "Línea roja: regresión lineal | Línea gris punteada: identidad (y = x)"
  ) +

  theme_minimal(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(color = "gray40", size = 11),
    plot.caption  = element_text(color = "gray50", size = 9),
    plot.margin   = margin(15, 20, 15, 15),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(fig_dir, "04_beta_concordance_scatter.png"), p4, width = 10, height = 9, dpi = 300)
ggsave(file.path(fig_dir, "04_beta_concordance_scatter.pdf"), p4, width = 10, height = 9)

# ================================================================
# GRÁFICA 5: PIP Comparison (posterior inclusion probabilities)
# ================================================================
cat(">> Gráfica 5: Comparación de PIP...\n")

# Variantes con PIP > 0.01 en al menos un modelo
df_pip <- df[PIP_T2DGGI > 0.01 | PIP_MCPS > 0.01]
cat(sprintf("  Variantes con PIP > 0.01 (al menos un GWAS): %s\n",
            format(nrow(df_pip), big.mark = ",")))

# Subsample
n_pip <- min(nrow(df_pip), 100000)
df_pip_sample <- df_pip[sample(.N, n_pip)]

p5 <- ggplot(df_pip_sample, aes(x = PIP_T2DGGI, y = PIP_MCPS)) +
  geom_point(alpha = 0.05, size = 0.5, color = "#7B2D8E") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.5) +

  # Marcar umbrales
  geom_hline(yintercept = 0.95, linetype = "dotted", color = col_mcp, linewidth = 0.5) +
  geom_vline(xintercept = 0.95, linetype = "dotted", color = col_t2d, linewidth = 0.5) +

  annotate("text", x = 0.95, y = 0.02, label = "PIP = 0.95", hjust = -0.1,
           size = 3, color = col_t2d, fontface = "italic") +
  annotate("text", x = 0.02, y = 0.95, label = "PIP = 0.95", vjust = -0.5,
           size = 3, color = col_mcp, fontface = "italic") +

  labs(
    title = "Probabilidad de Inclusión Posterior (PIP) — T2DGGI vs MCPS",
    subtitle = sprintf("%s variantes con PIP > 0.01 en al menos un modelo",
      format(nrow(df_pip), big.mark = ",")),
    x = "PIP — T2DGGI",
    y = "PIP — MCPS",
    caption = "PIP > 0.95: evidencia fuerte de causalidad | Línea punteada: identidad"
  ) +

  theme_minimal(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(color = "gray40", size = 11),
    plot.caption  = element_text(color = "gray50", size = 9),
    plot.margin   = margin(15, 20, 15, 15),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(fig_dir, "05_pip_comparison.png"), p5, width = 9, height = 8, dpi = 300)
ggsave(file.path(fig_dir, "05_pip_comparison.pdf"), p5, width = 9, height = 8)

# ================================================================
# GRÁFICA 6: Número de variantes por categoría (stacked bars)
# ================================================================
cat(">> Gráfica 6: Número de variantes por componente de mixture...\n")

df_nsnp <- data.table(
  Component = rep(c("Cat 2\n(ruido)", "Cat 3\n(pequeño)", "Cat 4\n(medio)", "Cat 5\n(grande)"), each = 2),
  GWAS = rep(c("T2DGGI", "MCPS"), 4),
  Count = c(
    405.14, 316.54,      # NumSnp2
    8484.71, 13.80,      # NumSnp3
    329.06, 886.91,      # NumSnp4
    12.43, 19.14         # NumSnp5
  ),
  SD = c(
    308.22, 235.17,
    207.60, 13.64,
    35.34, 43.41,
    5.56, 7.28
  )
)

df_nsnp[, Component := factor(Component, levels = c("Cat 2\n(ruido)", "Cat 3\n(pequeño)",
                                                     "Cat 4\n(medio)", "Cat 5\n(grande)"))]
df_nsnp[, GWAS := factor(GWAS, levels = c("T2DGGI", "MCPS"))]

p6 <- ggplot(df_nsnp, aes(x = Component, y = Count, fill = GWAS)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.75),
           width = 0.65, color = "white", linewidth = 0.4) +
  geom_errorbar(aes(ymin = pmax(Count - SD, 0), ymax = Count + SD),
                position = position_dodge(width = 0.75), width = 0.2, linewidth = 0.6) +

  scale_fill_manual(values = c("T2DGGI" = col_t2d, "MCPS" = col_mcp), name = "GWAS") +

  # Escala log para poder ver las categorías grandes y chicas
  scale_y_log10(labels = comma, breaks = c(1, 10, 100, 1000, 10000)) +

  labs(
    title = "Variantes Seleccionadas por Componente de Mixture",
    subtitle = "Número promedio de SNPs con efecto ≠ 0 por categoría (escala log)",
    x = "",
    y = "Número de variantes (log₁₀)",
    caption = "Cat 2: σ²×0.0001 | Cat 3: σ²×0.001 | Cat 4: σ²×0.01 | Cat 5: σ²×0.1"
  ) +

  theme_minimal(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(color = "gray40", size = 11),
    plot.caption  = element_text(color = "gray50", size = 9),
    plot.margin   = margin(15, 20, 15, 15),
    legend.position = c(0.88, 0.88),
    legend.background = element_rect(fill = "white", color = "gray80", linewidth = 0.3),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(fig_dir, "06_nsnp_by_component.png"), p6, width = 11, height = 7, dpi = 300)
ggsave(file.path(fig_dir, "06_nsnp_by_component.pdf"), p6, width = 11, height = 7)

# ================================================================
# GRÁFICA 7: Manhattan-style de |β| por cromosoma (facetado)
# ================================================================
cat(">> Gráfica 7: Manhattan de betas por cromosoma...\n")

# Extraer cromosoma de los IDs tipo chr10:100000235:C:T
df[, CHR := as.integer(gsub("^chr([0-9]+):.*", "\\1", SNP))]
df[, POS := as.numeric(gsub("^chr[0-9]+:([0-9]+):.*", "\\1", SNP))]

# Top variantes (|β| > percentil 99.5 en alguno)
q995_t <- quantile(abs(df[nz_t2d == TRUE, BETA_T2DGGI]), 0.995)
q995_m <- quantile(abs(df[nz_mcp == TRUE, BETA_MCPS]), 0.995)

df_top_t <- df[nz_t2d & abs(BETA_T2DGGI) > q995_t, .(CHR, POS, BETA = abs(BETA_T2DGGI), GWAS = "T2DGGI")]
df_top_m <- df[nz_mcp & abs(BETA_MCPS) > q995_m, .(CHR, POS, BETA = abs(BETA_MCPS), GWAS = "MCPS")]
df_manh <- rbind(df_top_t, df_top_m)
df_manh[, GWAS := factor(GWAS, levels = c("T2DGGI", "MCPS"))]

# Posición genómica acumulada
chr_lengths <- df_manh[, .(max_pos = max(POS)), by = CHR][order(CHR)]
chr_lengths[, cumstart := cumsum(shift(max_pos, fill = 0))]
df_manh <- merge(df_manh, chr_lengths[, .(CHR, cumstart)], by = "CHR")
df_manh[, BP_CUM := POS + cumstart]

# Centros de cromosoma para labels
chr_centers <- df_manh[, .(center = mean(BP_CUM)), by = CHR]

p7 <- ggplot(df_manh, aes(x = BP_CUM, y = BETA, color = factor(CHR %% 2))) +
  geom_point(alpha = 0.4, size = 0.6) +

  facet_wrap(~GWAS, ncol = 1, scales = "free_y") +

  scale_color_manual(values = c("0" = "#2196F3", "1" = "#90CAF9"), guide = "none") +
  scale_x_continuous(
    breaks = chr_centers$center,
    labels = chr_centers$CHR,
    expand = expansion(mult = 0.01)
  ) +

  labs(
    title = "Señal Genómica — Variantes con Mayor Efecto Posterior",
    subtitle = sprintf("Top 0.5%% de |β| por GWAS | T2DGGI: %s variantes, MCPS: %s variantes",
      format(nrow(df_top_t), big.mark = ","), format(nrow(df_top_m), big.mark = ",")),
    x = "Cromosoma",
    y = "|β| Posterior",
    caption = "SBayesRC sin anotaciones | Top 0.5% del efecto absoluto"
  ) +

  theme_minimal(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(color = "gray40", size = 11),
    plot.caption  = element_text(color = "gray50", size = 9),
    plot.margin   = margin(15, 20, 15, 15),
    strip.text    = element_text(face = "bold", size = 13),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

ggsave(file.path(fig_dir, "07_manhattan_top_betas.png"), p7, width = 14, height = 8, dpi = 300)
ggsave(file.path(fig_dir, "07_manhattan_top_betas.pdf"), p7, width = 14, height = 8)

# ================================================================
# GRÁFICA 8: Venn-style (UpSet simplificado) de señales
# ================================================================
cat(">> Gráfica 8: Venn de variantes con señal...\n")

# Contar categorías de overlap
n_only_t2d <- sum(df$nz_t2d & !df$nz_mcp)
n_only_mcp <- sum(!df$nz_t2d & df$nz_mcp)
n_both     <- sum(df$nz_t2d & df$nz_mcp)
n_neither  <- sum(!df$nz_t2d & !df$nz_mcp)

df_venn <- data.table(
  Category = factor(c("Solo T2DGGI", "Ambos", "Solo MCPS"),
                    levels = c("Solo T2DGGI", "Ambos", "Solo MCPS")),
  Count = c(n_only_t2d, n_both, n_only_mcp)
)

p8 <- ggplot(df_venn, aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat = "identity", color = "white", linewidth = 0.4, width = 0.65) +
  geom_text(aes(label = sprintf("%s\n(%.1f%%)", format(Count, big.mark = ","),
                                 Count / sum(Count) * 100)),
            vjust = -0.3, size = 4.2, fontface = "bold") +

  scale_fill_manual(values = c("Solo T2DGGI" = col_t2d,
                                "Ambos"       = "#7B2D8E",
                                "Solo MCPS"   = col_mcp)) +

  labs(
    title = "Solapamiento de Variantes con Efecto (β ≠ 0)",
    subtitle = sprintf("Total variantes alineadas: %s | Sin señal en ninguno: %s",
      format(nrow(df), big.mark = ","), format(n_neither, big.mark = ",")),
    x = "",
    y = "Número de variantes",
    caption = "Variantes incluidas en el modelo con β posterior ≠ 0"
  ) +

  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.15))) +

  theme_minimal(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(color = "gray40", size = 11),
    plot.caption  = element_text(color = "gray50", size = 9),
    plot.margin   = margin(15, 20, 15, 15),
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(fig_dir, "08_variant_overlap.png"), p8, width = 9, height = 7, dpi = 300)
ggsave(file.path(fig_dir, "08_variant_overlap.pdf"), p8, width = 9, height = 7)

# ================================================================
# GRÁFICA 9: Resumen tipo "forest plot" de parámetros globales
# ================================================================
cat(">> Gráfica 9: Forest plot de parámetros MCMC...\n")

df_forest <- mcmc_params[Parameter %in% c("h²_SNP", "σ²")]

p9 <- ggplot(df_forest, aes(x = Mean, y = Parameter, color = GWAS)) +
  geom_pointrange(aes(xmin = Mean - 1.96*SD, xmax = Mean + 1.96*SD),
                  position = position_dodge(width = 0.5),
                  size = 0.8, linewidth = 0.8) +

  scale_color_manual(values = c("T2DGGI" = col_t2d, "MCPS" = col_mcp), name = "GWAS") +

  labs(
    title = "Estimaciones Posteriores MCMC — Parámetros Globales",
    subtitle = "Media ± 95% intervalo credibilidad (±1.96 SD)",
    x = "Estimación posterior",
    y = ""
  ) +

  theme_minimal(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(color = "gray40", size = 11),
    plot.margin   = margin(15, 20, 15, 15),
    legend.position = "bottom",
    legend.background = element_rect(fill = "white", color = "gray80", linewidth = 0.3),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(fig_dir, "09_mcmc_forest_plot.png"), p9, width = 10, height = 5, dpi = 300)
ggsave(file.path(fig_dir, "09_mcmc_forest_plot.pdf"), p9, width = 10, height = 5)

# ================================================================
# RESUMEN FINAL
# ================================================================
cat("\n============================================================\n")
cat("  RESUMEN DE COMPARACIÓN SBayesRC\n")
cat("============================================================\n")
cat(sprintf("  GWAS T2DGGI (global):  h²=%.4f, nnz=%s, Var_e=%.4f\n",
    0.1249, format(9229, big.mark=","), 1.2405))
cat(sprintf("  GWAS MCPS (propio):    h²=%.4f, nnz=%s, Var_e=%.4f\n",
    0.1690, format(1233, big.mark=","), 1.0000))
cat(sprintf("\n  Variantes alineadas: %s\n", format(nrow(df), big.mark = ",")))
cat(sprintf("  β ≠ 0 en ambos:     %s\n", format(n_both, big.mark = ",")))
cat(sprintf("  Correlación betas:   r = %.4f\n", r_val))
cat(sprintf("\n  9 gráficas exportadas en: %s/\n", fig_dir))
cat("============================================================\n")
EOF

apptainer exec \
    --bind /mnt/cephfs:/mnt/cephfs \
    --bind "$TMPDIR:/tmp,$TMPDIR:/var/tmp" \
    "$SIF" \
    /opt/miniforge3/envs/r-bio/bin/Rscript "$R_SCRIPT"

rm -rf "$TMPDIR"
echo "==== COMPARACIÓN COMPLETADA ===="
