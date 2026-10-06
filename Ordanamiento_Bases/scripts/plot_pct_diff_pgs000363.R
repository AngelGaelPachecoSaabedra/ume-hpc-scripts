Sys.setlocale("LC_ALL", "en_US.UTF-8")

suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
})

OUT_DIR <- "/mnt/cephfs/orgs/home/angel.pacheco/Ordanamiento_Bases/scripts"
df <- read.csv(file.path(OUT_DIR, "prs_comparativa_pgs000363.csv"),
               stringsAsFactors = FALSE)
n <- nrow(df)

mu  <- mean(df$pct_diff)
med <- median(df$pct_diff)
sig <- sd(df$pct_diff)

cat(sprintf("n       = %s\n",   format(n, big.mark = ",")))
cat(sprintf("Media   = %+.6f\n", mu))
cat(sprintf("Mediana = %+.1f\n", med))
cat(sprintf("SD      =  %.6f\n", sig))

# Categoria por magnitud de diferencia
df$cat <- ifelse(df$pct_diff == 0, "Sin diferencia",
          ifelse(abs(df$pct_diff) <= 1, "+-1 percentil",
          ifelse(abs(df$pct_diff) <= 2, "+-2 percentiles",
                                        ">+-2 percentiles")))
df$cat <- factor(df$cat,
                 levels = c("Sin diferencia", "+-1 percentil",
                            "+-2 percentiles", ">+-2 percentiles"))

p <- ggplot(df, aes(x = pct_diff, fill = cat)) +
  geom_bar(color = "white", linewidth = 0.3, alpha = 0.85, width = 0.8) +

  scale_fill_manual(
    values = c("Sin diferencia"   = "#4CAF50",
               "+-1 percentil"    = "#2196F3",
               "+-2 percentiles"  = "#E6AB02",
               ">+-2 percentiles" = "#F44336"),
    name = NULL
  ) +

  # media
  geom_vline(xintercept = mu,
             color = "black", linetype = "solid", linewidth = 1.0) +
  # +- 1 SD
  geom_vline(xintercept = mu + sig,
             color = "#E64B35", linetype = "dashed", linewidth = 0.9) +
  geom_vline(xintercept = mu - sig,
             color = "#E64B35", linetype = "dashed", linewidth = 0.9) +
  # +- 2 SD
  geom_vline(xintercept = mu + 2*sig,
             color = "#F44336", linetype = "dotted", linewidth = 0.9) +
  geom_vline(xintercept = mu - 2*sig,
             color = "#F44336", linetype = "dotted", linewidth = 0.9) +

  annotate("text", x = mu, y = Inf,
    label = sprintf("Media\n%+.3f", mu),
    vjust = 1.4, hjust = -0.1, size = 3.2, color = "black", fontface = "bold") +

  annotate("text", x = mu + sig, y = Inf,
    label = sprintf("+1SD\n%+.2f", mu + sig),
    vjust = 1.4, hjust = -0.1, size = 3, color = "#E64B35") +
  annotate("text", x = mu - sig, y = Inf,
    label = sprintf("-1SD\n%+.2f", mu - sig),
    vjust = 1.4, hjust = 1.1, size = 3, color = "#E64B35") +

  annotate("text", x = mu + 2*sig, y = Inf,
    label = sprintf("+2SD\n%+.2f", mu + 2*sig),
    vjust = 3.2, hjust = -0.1, size = 3, color = "#F44336") +
  annotate("text", x = mu - 2*sig, y = Inf,
    label = sprintf("-2SD\n%+.2f", mu - 2*sig),
    vjust = 3.2, hjust = 1.1, size = 3, color = "#F44336") +

  annotate("label",
    x = max(df$pct_diff) + 0.3, y = Inf, vjust = 1.5,
    label = sprintf(
      "Media   = %+.4f\nMediana = %+.1f\nSD      =  %.4f\n\nSin dif: %.1f%%\n+-1 pctil: %.1f%%\n+-2 pctil: %.1f%%",
      mu, med, sig,
      mean(df$pct_diff == 0) * 100,
      mean(abs(df$pct_diff) <= 1) * 100,
      mean(abs(df$pct_diff) <= 2) * 100),
    hjust = 1, size = 3.5, fill = "white", alpha = 0.85) +

  scale_x_continuous(breaks = seq(min(df$pct_diff), max(df$pct_diff), by = 1)) +
  scale_y_continuous(labels = comma) +

  labs(
    title    = "Diferencia de Percentiles: Spark+CuPy vs Baseline",
    subtitle = sprintf("PGS000363 | n = %s | pct_diff = pct_NEW - pct_OLD (entero 0-100)",
      format(n, big.mark = ",")),
    x       = "Diferencia de percentil (pct_NEW - pct_OLD)",
    y       = "Numero de muestras",
    caption = "Negro: media | Rojo discontinuo: +-1SD | Rojo punteado: +-2SD"
  ) +

  theme_minimal(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(color = "gray40", size = 11),
    plot.caption  = element_text(color = "gray50", size = 9),
    plot.margin   = margin(15, 20, 15, 15),
    legend.position = "bottom",
    legend.background = element_rect(fill = "white", color = "gray80", linewidth = 0.3),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(OUT_DIR, "prs_pgs000363_pct_diff.png"), p, width = 14, height = 7, dpi = 200)
ggsave(file.path(OUT_DIR, "prs_pgs000363_pct_diff.pdf"), p, width = 14, height = 7)
cat("Guardado: prs_pgs000363_pct_diff.png\nCOMPLETADO.\n")
