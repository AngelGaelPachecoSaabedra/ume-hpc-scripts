Sys.setlocale("LC_ALL", "en_US.UTF-8")

suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
})

OUT_DIR <- "/mnt/cephfs/orgs/home/angel.pacheco/Ordanamiento_Bases/scripts"
IN_CSV  <- file.path(OUT_DIR, "prs_comparativa_pgs000363.csv")

df <- read.csv(IN_CSV, stringsAsFactors = FALSE)
n  <- nrow(df)

df$diff <- df$PRS_NEW - df$PRS_OLD

# ─── Estadisticas ────────────────────────────────────────────────────────────
mu  <- mean(df$diff)
sig <- sd(df$diff)

cat("=== Diferencias PRS_NEW - PRS_OLD ===\n")
cat(sprintf("  n        = %s\n",   format(n, big.mark = ",")))
cat(sprintf("  Media    = %+.6e\n", mu))
cat(sprintf("  SD       = %.6e\n",  sig))
cat(sprintf("  Media+1SD= %+.6e\n", mu + sig))
cat(sprintf("  Media-1SD= %+.6e\n", mu - sig))
cat(sprintf("  Media+2SD= %+.6e\n", mu + 2*sig))
cat(sprintf("  Media-2SD= %+.6e\n", mu - 2*sig))
cat(sprintf("  Dentro +-1SD: %.2f%%\n", mean(abs(df$diff - mu) <= sig)   * 100))
cat(sprintf("  Dentro +-2SD: %.2f%%\n", mean(abs(df$diff - mu) <= 2*sig) * 100))

# ─── Grafica ─────────────────────────────────────────────────────────────────
p <- ggplot(df, aes(x = diff)) +
  geom_histogram(bins = 80, fill = "#2196F3", color = "white",
                 linewidth = 0.3, alpha = 0.85) +

  # media
  geom_vline(xintercept = mu,
             color = "black", linetype = "solid", linewidth = 1.0) +
  # media +- 1 SD
  geom_vline(xintercept = mu + sig,
             color = "#E64B35", linetype = "dashed", linewidth = 0.9) +
  geom_vline(xintercept = mu - sig,
             color = "#E64B35", linetype = "dashed", linewidth = 0.9) +
  # media +- 2 SD
  geom_vline(xintercept = mu + 2*sig,
             color = "#F44336", linetype = "dotted", linewidth = 0.9) +
  geom_vline(xintercept = mu - 2*sig,
             color = "#F44336", linetype = "dotted", linewidth = 0.9) +

  annotate("text", x = mu, y = Inf,
    label = sprintf("Media\n%+.2e", mu),
    vjust = 1.4, hjust = -0.1, size = 3.2, color = "black", fontface = "bold") +

  annotate("text", x = mu + sig, y = Inf,
    label = sprintf("+1SD\n%+.2e", mu + sig),
    vjust = 1.4, hjust = -0.1, size = 3, color = "#E64B35") +
  annotate("text", x = mu - sig, y = Inf,
    label = sprintf("-1SD\n%+.2e", mu - sig),
    vjust = 1.4, hjust = 1.1, size = 3, color = "#E64B35") +

  annotate("text", x = mu + 2*sig, y = Inf,
    label = sprintf("+2SD\n%+.2e", mu + 2*sig),
    vjust = 3.2, hjust = -0.1, size = 3, color = "#F44336") +
  annotate("text", x = mu - 2*sig, y = Inf,
    label = sprintf("-2SD\n%+.2e", mu - 2*sig),
    vjust = 3.2, hjust = 1.1, size = 3, color = "#F44336") +

  annotate("label",
    x = max(df$diff) - diff(range(df$diff)) * 0.01, y = Inf, vjust = 1.5,
    label = sprintf(
      "Media  = %+.4e\nSD     =  %.4e\n+-1SD: %.1f%% muestras\n+-2SD: %.1f%% muestras",
      mu, sig,
      mean(abs(df$diff - mu) <= sig)   * 100,
      mean(abs(df$diff - mu) <= 2*sig) * 100),
    hjust = 1, size = 3.5, fill = "white", alpha = 0.85) +

  labs(
    title    = "Diferencias PRS: Media y Desviacion Estandar",
    subtitle = sprintf("PGS000363 | n = %s | Diff = PRS_nuevo - PRS_baseline",
      format(n, big.mark = ",")),
    x       = "PRS_nuevo - PRS_baseline",
    y       = "Numero de muestras",
    caption = "Negro: media | Rojo discontinuo: +-1SD | Rojo punteado: +-2SD"
  ) +

  theme_minimal(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(color = "gray40", size = 11),
    plot.caption  = element_text(color = "gray50", size = 9),
    plot.margin   = margin(15, 20, 15, 15),
    panel.grid.minor = element_blank()
  ) +
  scale_y_continuous(labels = comma)

ggsave(file.path(OUT_DIR, "prs_pgs000363_diff_mean_sd.png"), p, width = 12, height = 7, dpi = 200)
ggsave(file.path(OUT_DIR, "prs_pgs000363_diff_mean_sd.pdf"), p, width = 12, height = 7)
cat("\nGuardado: prs_pgs000363_diff_mean_sd.png\n")
cat("COMPLETADO.\n")
