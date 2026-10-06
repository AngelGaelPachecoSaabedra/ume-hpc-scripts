# ─────────────────────────────────────────────────────────────────────────────
# Validacion visual PRS - PGS000363
# Una grafica individual por cada aspecto de la comparacion
#
# Input : prs_comparativa_pgs000363.csv
# Output: 5 PNG + 5 PDF en OUT_DIR
# ─────────────────────────────────────────────────────────────────────────────

Sys.setlocale("LC_ALL", "en_US.UTF-8")

suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
})

OUT_DIR <- "/mnt/cephfs/orgs/home/angel.pacheco/Ordanamiento_Bases/scripts"
IN_CSV  <- file.path(OUT_DIR, "prs_comparativa_pgs000363.csv")

cat("Cargando datos...\n")
df <- read.csv(IN_CSV, stringsAsFactors = FALSE)
n  <- nrow(df)
cat(sprintf("Muestras: %s\n\n", format(n, big.mark = ",")))

# ─── Metricas globales ───────────────────────────────────────────────────────
r_pearson  <- cor(df$PRS_NEW, df$PRS_OLD, method = "pearson")
r_spearman <- cor(df$PRS_NEW, df$PRS_OLD, method = "spearman")
mae        <- mean(abs(df$PRS_NEW - df$PRS_OLD))
max_diff   <- max(abs(df$PRS_NEW - df$PRS_OLD))

pct_exact <- mean(df$pct_diff == 0) * 100
pct_1     <- mean(df$pct_diff <= 1) * 100
pct_2     <- mean(df$pct_diff <= 2) * 100
pct_5     <- mean(df$pct_diff <= 5) * 100

# Estandarizar (referencia = parametros del baseline)
old_mean <- mean(df$PRS_OLD)
old_sd   <- sd(df$PRS_OLD)
df$PRS_NEW_Z <- (df$PRS_NEW - old_mean) / old_sd
df$PRS_OLD_Z <- (df$PRS_OLD - old_mean) / old_sd
df$error     <- df$PRS_NEW - df$PRS_OLD

# ─── Helper: categorias de riesgo por z-score ────────────────────────────────
risk_cat <- function(z) {
  ifelse(z <= -2, "Bajo riesgo (<= -2s)",
  ifelse(z >=  2, "Alto riesgo (>= 2s)",
                  "Riesgo intermedio"))
}

risk_colors <- c(
  "Bajo riesgo (<= -2s)" = "#4CAF50",
  "Riesgo intermedio"    = "#2196F3",
  "Alto riesgo (>= 2s)"  = "#F44336"
)

# ─── Tema base del pipeline ──────────────────────────────────────────────────
tema_prs <- theme_minimal(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(color = "gray40", size = 11),
    plot.caption  = element_text(color = "gray50", size = 9),
    plot.margin   = margin(15, 20, 15, 15),
    legend.background = element_rect(fill = "white", color = "gray80", linewidth = 0.3),
    panel.grid.minor = element_blank()
  )

save_plot <- function(p, name, w, h) {
  ggsave(file.path(OUT_DIR, paste0(name, ".png")), p, width = w, height = h, dpi = 200)
  ggsave(file.path(OUT_DIR, paste0(name, ".pdf")), p, width = w, height = h)
  cat(sprintf("  Guardado: %s.png\n", name))
}

# ═════════════════════════════════════════════════════════════════════════════
# GRAFICA 1: Distribucion PRS_NEW (Spark+CuPy)
# ═════════════════════════════════════════════════════════════════════════════
cat("Grafica 1: Distribucion PRS_NEW...\n")

df$risk_new <- risk_cat(df$PRS_NEW_Z)

p5_new  <- quantile(df$PRS_NEW_Z, 0.05)
p95_new <- quantile(df$PRS_NEW_Z, 0.95)
n_high_new <- sum(df$PRS_NEW_Z >= 2)
n_low_new  <- sum(df$PRS_NEW_Z <= -2)

p1 <- ggplot(df, aes(x = PRS_NEW_Z, fill = risk_new)) +
  geom_histogram(bins = 80, color = "white", linewidth = 0.3, alpha = 0.85) +

  scale_fill_manual(values = risk_colors, name = NULL) +

  geom_vline(xintercept = 0,      linetype = "dashed", color = "black",   linewidth = 0.8) +
  geom_vline(xintercept = p5_new,  linetype = "dotted", color = "#4CAF50", linewidth = 0.9) +
  geom_vline(xintercept = p95_new, linetype = "dotted", color = "#F44336", linewidth = 0.9) +

  annotate("label",
    x = max(df$PRS_NEW_Z) - 0.3, y = Inf, vjust = 1.5,
    label = sprintf("Alto riesgo (>= 2s): %s (%.1f%%)\nBajo riesgo (<= -2s): %s (%.1f%%)",
      format(n_high_new, big.mark = ","), n_high_new / n * 100,
      format(n_low_new,  big.mark = ","), n_low_new  / n * 100),
    hjust = 1, size = 3.5, fill = "white", alpha = 0.8) +

  annotate("text", x = p5_new,  y = 0,
    label = sprintf("P5 (%.2fs)",  p5_new),
    vjust = -0.5, hjust = 1.1, size = 3, color = "#4CAF50") +
  annotate("text", x = p95_new, y = 0,
    label = sprintf("P95 (%.2fs)", p95_new),
    vjust = -0.5, hjust = -0.1, size = 3, color = "#F44336") +

  labs(
    title    = "Distribucion del PRS - Pipeline Spark+CuPy (nuevo)",
    subtitle = sprintf("PGS000363 | Cohorte MCPS (n = %s) | z-score estandarizado con parametros del baseline",
      format(n, big.mark = ",")),
    x       = "PRS Estandarizado (z-score)",
    y       = "Numero de Pacientes",
    caption = "Pipeline: Spark + CuPy (GPU)"
  ) +
  tema_prs +
  theme(legend.position = "bottom") +
  scale_y_continuous(labels = comma)

save_plot(p1, "prs_pgs000363_distribucion_NEW", 12, 7)

# ═════════════════════════════════════════════════════════════════════════════
# GRAFICA 2: Distribucion PRS_OLD (Baseline)
# ═════════════════════════════════════════════════════════════════════════════
cat("Grafica 2: Distribucion PRS_OLD...\n")

df$risk_old <- risk_cat(df$PRS_OLD_Z)

p5_old  <- quantile(df$PRS_OLD_Z, 0.05)
p95_old <- quantile(df$PRS_OLD_Z, 0.95)
n_high_old <- sum(df$PRS_OLD_Z >= 2)
n_low_old  <- sum(df$PRS_OLD_Z <= -2)

p2 <- ggplot(df, aes(x = PRS_OLD_Z, fill = risk_old)) +
  geom_histogram(bins = 80, color = "white", linewidth = 0.3, alpha = 0.85) +

  scale_fill_manual(values = risk_colors, name = NULL) +

  geom_vline(xintercept = 0,      linetype = "dashed", color = "black",   linewidth = 0.8) +
  geom_vline(xintercept = p5_old,  linetype = "dotted", color = "#4CAF50", linewidth = 0.9) +
  geom_vline(xintercept = p95_old, linetype = "dotted", color = "#F44336", linewidth = 0.9) +

  annotate("label",
    x = max(df$PRS_OLD_Z) - 0.3, y = Inf, vjust = 1.5,
    label = sprintf("Alto riesgo (>= 2s): %s (%.1f%%)\nBajo riesgo (<= -2s): %s (%.1f%%)",
      format(n_high_old, big.mark = ","), n_high_old / n * 100,
      format(n_low_old,  big.mark = ","), n_low_old  / n * 100),
    hjust = 1, size = 3.5, fill = "white", alpha = 0.8) +

  annotate("text", x = p5_old,  y = 0,
    label = sprintf("P5 (%.2fs)",  p5_old),
    vjust = -0.5, hjust = 1.1, size = 3, color = "#4CAF50") +
  annotate("text", x = p95_old, y = 0,
    label = sprintf("P95 (%.2fs)", p95_old),
    vjust = -0.5, hjust = -0.1, size = 3, color = "#F44336") +

  labs(
    title    = "Distribucion del PRS - Baseline (previo)",
    subtitle = sprintf("PGS000363 | Cohorte MCPS (n = %s) | z-score estandarizado",
      format(n, big.mark = ",")),
    x       = "PRS Estandarizado (z-score)",
    y       = "Numero de Pacientes",
    caption = "Pipeline: baseline de referencia (corrida previa)"
  ) +
  tema_prs +
  theme(legend.position = "bottom") +
  scale_y_continuous(labels = comma)

save_plot(p2, "prs_pgs000363_distribucion_OLD", 12, 7)

# ═════════════════════════════════════════════════════════════════════════════
# GRAFICA 3: Scatter PRS_NEW vs PRS_OLD
# ═════════════════════════════════════════════════════════════════════════════
cat("Grafica 3: Scatter PRS_NEW vs PRS_OLD...\n")

lim <- range(c(df$PRS_OLD_Z, df$PRS_NEW_Z))

p3 <- ggplot(df, aes(x = PRS_OLD_Z, y = PRS_NEW_Z)) +
  geom_point(alpha = 0.05, size = 0.3, color = "#2196F3") +

  geom_abline(slope = 1, intercept = 0,
              color = "black", linetype = "dashed", linewidth = 0.8) +
  geom_smooth(method = "lm", color = "#E64B35", linewidth = 1, se = TRUE) +

  coord_equal(xlim = lim, ylim = lim) +

  annotate("label",
    x = lim[2] - diff(lim) * 0.02,
    y = lim[1] + diff(lim) * 0.02,
    vjust = 0,
    label = sprintf(
      "Pearson r = %.6f\nSpearman rho = %.6f\nMAE = %.2e\nMax diff = %.2e",
      r_pearson, r_spearman, mae, max_diff),
    hjust = 1, size = 3.5, fill = "white", alpha = 0.85) +

  labs(
    title    = "Comparacion PRS: Spark+CuPy vs Baseline",
    subtitle = sprintf("PGS000363 | n = %s | Pearson r = %.6f | Spearman rho = %.6f",
      format(n, big.mark = ","), r_pearson, r_spearman),
    x       = "PRS Baseline (z-score)",
    y       = "PRS Spark+CuPy (z-score)",
    caption = "Linea negra: identidad (y = x) | Linea roja: regresion lineal con IC 95%"
  ) +
  tema_prs +
  theme(legend.position = "none")

save_plot(p3, "prs_pgs000363_scatter_prs", 10, 9)

# ═════════════════════════════════════════════════════════════════════════════
# GRAFICA 4: Histograma del error (PRS_NEW - PRS_OLD)
# ═════════════════════════════════════════════════════════════════════════════
cat("Grafica 4: Histograma del error...\n")

p25_err <- quantile(abs(df$error), 0.25)
p75_err <- quantile(abs(df$error), 0.75)
p99_err <- quantile(abs(df$error), 0.99)

df$err_cat <- ifelse(abs(df$error) <= p25_err, "Error bajo (<= P25)",
             ifelse(abs(df$error) <= p75_err, "Error medio (P25-P75)",
                                              "Error alto (> P75)"))
df$err_cat <- factor(df$err_cat,
                     levels = c("Error bajo (<= P25)",
                                "Error medio (P25-P75)",
                                "Error alto (> P75)"))

p4 <- ggplot(df, aes(x = error, fill = err_cat)) +
  geom_histogram(bins = 80, color = "white", linewidth = 0.3, alpha = 0.85) +

  scale_fill_manual(
    values = c("Error bajo (<= P25)"   = "#4CAF50",
               "Error medio (P25-P75)" = "#2196F3",
               "Error alto (> P75)"    = "#F44336"),
    name = NULL
  ) +

  geom_vline(xintercept =  0,       linetype = "dashed", color = "black",   linewidth = 0.8) +
  geom_vline(xintercept =  p99_err, linetype = "dotted", color = "#F44336", linewidth = 0.9) +
  geom_vline(xintercept = -p99_err, linetype = "dotted", color = "#F44336", linewidth = 0.9) +

  annotate("label",
    x = max(df$error, na.rm = TRUE) - diff(range(df$error)) * 0.01,
    y = Inf, vjust = 1.5,
    label = sprintf(
      "MAE = %.2e\nMax diff = %.2e\nMedia error = %.2e\nP99 error = %.2e",
      mae, max_diff, mean(df$error), p99_err),
    hjust = 1, size = 3.5, fill = "white", alpha = 0.85) +

  annotate("text", x =  p99_err, y = 0,
    label = sprintf("P99 (%.2e)", p99_err),
    vjust = -0.5, hjust = -0.05, size = 3, color = "#F44336") +
  annotate("text", x = -p99_err, y = 0,
    label = sprintf("P99 (%.2e)", -p99_err),
    vjust = -0.5, hjust = 1.05, size = 3, color = "#F44336") +

  labs(
    title    = "Distribucion del Error - PRS Nuevo vs Baseline",
    subtitle = sprintf("PGS000363 | n = %s | Error = PRS_nuevo - PRS_baseline",
      format(n, big.mark = ",")),
    x       = "PRS_nuevo - PRS_baseline (escala original)",
    y       = "Numero de muestras",
    caption = "Verde: <= P25 | Azul: P25-P75 | Rojo: > P75 del error absoluto"
  ) +
  tema_prs +
  theme(legend.position = "bottom") +
  scale_y_continuous(labels = comma)

save_plot(p4, "prs_pgs000363_histograma_error", 12, 7)

# ═════════════════════════════════════════════════════════════════════════════
# GRAFICA 5: Scatter de percentiles pct_NEW vs pct_OLD
# ═════════════════════════════════════════════════════════════════════════════
cat("Grafica 5: Scatter de percentiles...\n")

df$pct_cat <- ifelse(df$pct_diff == 0,  "Concordancia exacta",
             ifelse(df$pct_diff <= 2,   "Discordancia +-1-2 pctil",
                                        "Discordancia >=3 pctil"))
df$pct_cat <- factor(df$pct_cat,
                     levels = c("Concordancia exacta",
                                "Discordancia +-1-2 pctil",
                                "Discordancia >=3 pctil"))

p5 <- ggplot(df, aes(x = pct_OLD, y = pct_NEW, color = pct_cat)) +
  geom_jitter(alpha = 0.06, size = 0.35, width = 0.3, height = 0.3) +

  geom_abline(slope = 1, intercept = 0,
              color = "black", linetype = "dashed", linewidth = 0.8) +

  scale_color_manual(
    values = c("Concordancia exacta"      = "#4CAF50",
               "Discordancia +-1-2 pctil" = "#2196F3",
               "Discordancia >=3 pctil"   = "#F44336"),
    name = "Concordancia"
  ) +

  coord_equal(xlim = c(1, 100), ylim = c(1, 100)) +
  scale_x_continuous(breaks = seq(0, 100, 20)) +
  scale_y_continuous(breaks = seq(0, 100, 20)) +

  annotate("label",
    x = 2, y = 98, vjust = 1,
    label = sprintf(
      "Exacto: %.1f%%\n+-1 pctil: %.1f%%\n+-2 pctil: %.1f%%\n+-5 pctil: %.1f%%",
      pct_exact, pct_1, pct_2, pct_5),
    hjust = 0, size = 3.5, fill = "white", alpha = 0.85) +

  labs(
    title    = "Concordancia de Percentiles - PRS Nuevo vs Baseline",
    subtitle = sprintf("PGS000363 | n = %s | Percentil asignado por ranking en cada pipeline",
      format(n, big.mark = ",")),
    x       = "Percentil Baseline (1-100)",
    y       = "Percentil Spark+CuPy (1-100)",
    caption = "Verde: match exacto | Azul: +-1-2 pctil | Rojo: discordancia >=3 pctil"
  ) +
  tema_prs +
  theme(
    legend.position = c(0.88, 0.15),
    legend.text  = element_text(size = 9),
    legend.title = element_text(size = 9, face = "bold")
  )

save_plot(p5, "prs_pgs000363_scatter_percentiles", 10, 9)

# ─────────────────────────────────────────────────────────────────────────────
cat(sprintf("\nArchivos guardados en: %s\n", OUT_DIR))
cat("  [1] prs_pgs000363_distribucion_NEW.png/.pdf\n")
cat("  [2] prs_pgs000363_distribucion_OLD.png/.pdf\n")
cat("  [3] prs_pgs000363_scatter_prs.png/.pdf\n")
cat("  [4] prs_pgs000363_histograma_error.png/.pdf\n")
cat("  [5] prs_pgs000363_scatter_percentiles.png/.pdf\n")
cat("COMPLETADO.\n")
