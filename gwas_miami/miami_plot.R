#!/usr/bin/env Rscript
# miami_plot_nature.R — Miami plot completo conforme a Nature
# ============================================================
# Pipeline completo en R (sin GPU):
#   1. Carga de GWAS (males + females, formato REGENIE .gz)
#   2. Cálculo de coordenadas genómicas acumulativas
#   3. Carga y anotación GTF (gen más cercano por SNP significativo)
#   4. Selección de top SNPs por gen (1 por gen/sexo, ≤500 kb)
#   5. Miami plot con patchwork (males arriba, females abajo espejado)
#
# Especificaciones Nature aplicadas (ver nature_specs.md):
#   - 183 mm ancho × 150 mm alto (≤170 mm máx.)
#   - Sans-serif Arial/Helvetica (resolución automática)
#   - Texto 5–7 pt; panel label 8 pt bold minúscula
#   - Paleta Wong 2011 (azul #0072b2 / azul cielo #56b4e9 alternados)
#   - Sin gridlines, ejes con líneas y tick marks
#   - Nube rasterizada 450 dpi vía ggrastr::rasterise()
#   - Etiquetas con ggrepel (sin superposición)
#   - PDF vectorial vía ragg::agg_pdf (fuentes embebidas TrueType 42)
#
# Dependencias:
#   install.packages(c("ggplot2","dplyr","data.table","ggrastr",
#                       "ggrepel","patchwork","ragg","systemfonts"))
#
# Uso:
#   Rscript miami_plot_nature.R
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(ggrastr)
  library(ggrepel)
})

# ============================================================
# CONSTANTES
# ============================================================
MM_PER_IN <- 25.4

# Paleta Wong (2011) — color-blind safe
WONG <- c(
  black          = "#000000",
  orange         = "#e69f00",
  sky_blue       = "#56b4e9",
  bluish_green   = "#009e73",
  yellow         = "#f0e442",
  blue           = "#0072b2",
  vermillion     = "#d55e00",
  reddish_purple = "#cc79a7"
)

# Colores alternos por cromosoma (azul / azul cielo)
CHR_COLORS <- unname(WONG[c("blue", "sky_blue")])

# Etiquetas de cromosomas en orden
CHR_LABELS <- c(as.character(1:22), "X")

# ── Archivos de entrada ──
MALE_FILE   <- "combined-gwas-results-all-ctrl-males.txt.gz"
FEMALE_FILE <- "combined-gwas-results-all-ctrl-females.txt.gz"
GTF_FILE    <- "Homo_sapiens.GRCh38.105.gtf.gz"

# ── Parámetros ──
SIGNIFICANCE   <- 5e-8
SUGGESTIVE     <- 1e-5
MAX_GENE_DIST  <- 500000L   # bp
MAX_ANNOT_MALE <- 20L       # etiquetas por panel (legibilidad)
MAX_ANNOT_FEM  <- 20L
OUTPUT_PDF     <- "miami_t2d_nature.pdf"
OUTPUT_CSV     <- "top_snps_annotated.csv"

# ── Dimensiones Nature ──
WIDTH_MM  <- 183
HEIGHT_MM <- 150   # ≤170 mm
DPI       <- 450


# ============================================================
# 0. RESOLUCIÓN DE FUENTE
# ============================================================
resolve_font <- function() {
  candidates <- c("Arial", "Helvetica", "Liberation Sans", "Arimo",
                   "Nimbus Sans", "Nimbus Sans L", "DejaVu Sans")
  installed <- tryCatch(
    unique(systemfonts::system_fonts()$family),
    error = function(e) character(0)
  )
  for (cand in candidates) {
    if (cand %in% installed) {
      message("Fuente seleccionada: ", cand)
      return(cand)
    }
  }
  warning("No se encontro Arial/Helvetica/Liberation Sans/Arimo. ",
          "Usando 'sans' por defecto — verifica el PDF resultante.")
  "sans"
}


# ============================================================
# 1. CARGA DE DATOS GWAS
# ============================================================
load_gwas <- function(filepath) {
  message("Cargando ", filepath, "...")
  t0 <- proc.time()

  dt <- fread(filepath, select = c("CHROM", "GENPOS", "P"),
              colClasses = list(character = "CHROM"),
              fill = TRUE, showProgress = FALSE)

  # Normalizar cromosoma
  dt[, CHROM := sub("^chr", "", CHROM)]
  dt[CHROM == "23", CHROM := "X"]

  # Forzar numérico (P puede tener "NA" como texto en líneas malformadas)
  dt[, GENPOS := as.numeric(GENPOS)]
  dt[, P := as.numeric(P)]

  # Filtrar cromosomas válidos y p > 0
  valid <- c(as.character(1:22), "X")
  dt <- dt[CHROM %in% valid & !is.na(P) & P > 0 & !is.na(GENPOS)]

  # Renombrar
  setnames(dt, c("CHROM", "GENPOS", "P"), c("CHR", "POS", "pvalue"))

  # Orden natural de cromosomas (suppressWarnings para as.integer("X") → NA)
  chr_num <- suppressWarnings(as.integer(dt$CHR))
  chr_num[dt$CHR == "X"]  <- 23L
  chr_num[dt$CHR == "Y"]  <- 24L
  chr_num[dt$CHR == "MT"] <- 25L
  chr_num[dt$CHR == "M"]  <- 25L
  dt[, chr_num := chr_num]
  setorder(dt, chr_num, POS)
  dt[, chr_num := NULL]

  elapsed <- (proc.time() - t0)[3]
  message(sprintf("  %s SNPs cargados en %.1fs", format(nrow(dt), big.mark = ","), elapsed))
  dt
}


# ============================================================
# 2. COORDENADAS GENÓMICAS ACUMULATIVAS
# ============================================================
compute_coordinates <- function(dt, chr_labels = CHR_LABELS) {
  message("Calculando coordenadas acumulativas...")

  gap <- 10e6
  cumulative_offset <- 0
  chr_centers <- numeric(length(chr_labels))
  names(chr_centers) <- chr_labels

  dt[, pos_cumulative := 0.0]

  for (i in seq_along(chr_labels)) {
    chr_name <- chr_labels[i]
    idx <- which(dt$CHR == chr_name)
    if (length(idx) == 0L) {
      chr_centers[chr_name] <- cumulative_offset
      next
    }
    cur_pos <- dt$POS[idx]
    min_p <- min(cur_pos)
    max_p <- max(cur_pos)
    chr_centers[chr_name] <- cumulative_offset + (min_p + max_p) / 2
    set(dt, i = idx, j = "pos_cumulative", value = cur_pos + cumulative_offset)
    cumulative_offset <- cumulative_offset + max_p + gap
  }

  dt[, log10p := -log10(pvalue)]

  # Color alterno por cromosoma
  chr_idx <- match(dt$CHR, chr_labels)
  dt[, color := CHR_COLORS[(chr_idx %% 2L) + 1L]]

  list(dt = dt, chr_centers = chr_centers)
}


# ============================================================
# 3. CARGA Y ANOTACIÓN GTF
# ============================================================
load_gtf_genes <- function(gtf_path) {
  if (!file.exists(gtf_path)) {
    warning("GTF no encontrado: ", gtf_path)
    return(NULL)
  }
  message("Cargando GTF: ", gtf_path, " ...")
  t0 <- proc.time()

  con <- gzfile(gtf_path, open = "rt")
  on.exit(close(con))

  genes <- list()
  n <- 0L

  while (length(line <- readLines(con, n = 1L)) > 0) {
    if (startsWith(line, "#")) next
    fields <- strsplit(line, "\t", fixed = TRUE)[[1]]
    if (length(fields) < 9 || fields[3] != "gene") next

    chrom <- sub("^chr", "", fields[1])
    start <- as.integer(fields[4])
    end   <- as.integer(fields[5])

    # Extraer gene_name
    attr_str <- fields[9]
    m <- regmatches(attr_str, regexpr('gene_name "([^"]+)"', attr_str))
    if (length(m) == 0 || nchar(m) == 0) next
    gene_name <- sub('gene_name "', "", sub('"$', "", m))

    n <- n + 1L
    genes[[n]] <- data.table(chr = chrom, start = start, end = end,
                              gene_name = gene_name)
  }

  if (n == 0) return(NULL)
  genes_dt <- rbindlist(genes)
  setorder(genes_dt, chr, start)

  elapsed <- (proc.time() - t0)[3]
  message(sprintf("  %s genes cargados en %.1fs",
                  format(nrow(genes_dt), big.mark = ","), elapsed))
  genes_dt
}


annotate_snps <- function(snps_dt, genes_dt) {
  if (is.null(genes_dt) || nrow(snps_dt) == 0) {
    snps_dt[, `:=`(nearest_gene = "UNKNOWN", distance_to_gene = -1L)]
    return(snps_dt)
  }
  message(sprintf("Anotando %s SNPs significativos...",
                  format(nrow(snps_dt), big.mark = ",")))
  t0 <- proc.time()

  res_gene <- character(nrow(snps_dt))
  res_dist <- integer(nrow(snps_dt))

  for (chrom in unique(snps_dt$CHR)) {
    snp_idx <- which(snps_dt$CHR == chrom)
    g_idx   <- which(genes_dt$chr == chrom)
    if (length(g_idx) == 0) {
      res_gene[snp_idx] <- paste0("UNKNOWN_chr", chrom)
      res_dist[snp_idx] <- -1L
      next
    }

    snp_pos  <- snps_dt$POS[snp_idx]
    g_starts <- genes_dt$start[g_idx]
    g_ends   <- genes_dt$end[g_idx]
    g_names  <- genes_dt$gene_name[g_idx]

    # Para cada SNP: distancia vectorizada contra todos los genes del chr
    for (j in seq_along(snp_idx)) {
      p <- snp_pos[j]
      # Distancia: 0 si dentro, positiva si fuera
      d <- ifelse(p >= g_starts & p <= g_ends, 0L,
            ifelse(p < g_starts, g_starts - p, p - g_ends))
      best <- which.min(d)
      res_gene[snp_idx[j]] <- g_names[best]
      res_dist[snp_idx[j]] <- d[best]
    }
  }

  snps_dt[, `:=`(nearest_gene = res_gene, distance_to_gene = res_dist)]

  elapsed <- (proc.time() - t0)[3]
  within <- sum(res_dist == 0)
  message(sprintf("  Anotacion completada en %.1fs — %d dentro de genes",
                  elapsed, within))
  snps_dt
}


# ============================================================
# 4. SELECCIÓN DE TOP SNPs
# ============================================================
select_top_snps <- function(sig_dt, max_distance = MAX_GENE_DIST) {
  if (nrow(sig_dt) == 0) return(sig_dt[0])

  # Filtrar por distancia
  in_range <- sig_dt[distance_to_gene >= 0 & distance_to_gene <= max_distance]
  if (nrow(in_range) == 0) return(sig_dt[0])

  # Top SNP por gen y dataset (sexo)
  top <- in_range[, .SD[which.min(pvalue)], by = .(nearest_gene, dataset)]
  message(sprintf("  Top SNPs unicos por gen: %d (Males: %d, Females: %d)",
                  nrow(top),
                  sum(top$dataset == "Males"),
                  sum(top$dataset == "Females")))
  top
}


select_best_for_annotation <- function(top_dt, max_males = MAX_ANNOT_MALE,
                                        max_females = MAX_ANNOT_FEM) {
  if (nrow(top_dt) == 0) return(top_dt[0])
  males   <- head(top_dt[dataset == "Males"][order(pvalue)], max_males)
  females <- head(top_dt[dataset == "Females"][order(pvalue)], max_females)
  rbind(males, females)[order(pvalue)]
}


# ============================================================
# 5. MIAMI PLOT — panel único, eje Y compartido, chr labels al centro
# ============================================================
build_miami_plot <- function(males_res, females_res, sig_annot,
                              font_family, output_pdf = OUTPUT_PDF,
                              width_mm = WIDTH_MM, height_mm = HEIGHT_MM,
                              dpi = DPI) {

  males_dt   <- copy(males_res$dt)
  females_dt <- copy(females_res$dt)
  chr_centers <- males_res$chr_centers

  sig_y <- -log10(SIGNIFICANCE)
  sug_y <- -log10(SUGGESTIVE)

  # ── Preparar datos: males positivos, females negativos, con offset ──
  # El offset crea un gap alrededor de y=0 para las etiquetas de cromosoma
  y_offset <- 3
  males_dt[, y := log10p + y_offset]
  females_dt[, y := -(log10p + y_offset)]

  # Límites asimétricos: cada mitad se ajusta a su propio máximo
  y_top    <- max(males_dt$y, na.rm = TRUE) * 1.12
  y_bottom <- min(females_dt$y, na.rm = TRUE) * 1.12   # negativo

  # Breaks del eje Y: independientes arriba y abajo, sin ceros
  y_step <- if (max(y_top, abs(y_bottom)) > 50) 10 else
            if (max(y_top, abs(y_bottom)) > 25) 5 else 2

  # Males (arriba): breaks con offset, labels sin offset
  y_real_top <- seq(y_step, floor(y_top - y_offset), by = y_step)
  y_brk_top  <- y_real_top + y_offset

  # Females (abajo): breaks con offset, labels sin offset
  y_real_bot <- seq(y_step, floor(abs(y_bottom) - y_offset), by = y_step)
  y_brk_bot  <- -(y_real_bot + y_offset)

  y_breaks <- c(rev(y_brk_bot), y_brk_top)
  y_labels <- c(rev(as.character(y_real_bot)), as.character(y_real_top))

  # Data.frame para etiquetas de cromosoma en y = 0
  chr_label_df <- data.frame(
    x     = as.numeric(chr_centers),
    label = names(chr_centers),
    stringsAsFactors = FALSE
  )

  # ── Tema Nature ──
  nature_theme <- theme_minimal(base_family = font_family) +
    theme(
      panel.grid       = element_blank(),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background  = element_rect(fill = "white", colour = NA),
      axis.line.y      = element_line(colour = "black", linewidth = 0.6 / .pt),
      axis.line.x      = element_blank(),
      axis.ticks.y     = element_line(colour = "black", linewidth = 0.6 / .pt),
      axis.ticks.x     = element_blank(),
      axis.ticks.length = unit(1.2, "pt"),
      axis.text.y      = element_text(size = 6, colour = "black",
                                       family = font_family),
      axis.text.x      = element_blank(),
      axis.title       = element_text(size = 7, colour = "black",
                                       family = font_family),
      legend.position  = "none",
      plot.margin      = margin(6, 6, 6, 6, unit = "pt")
    )

  # ── Construir plot ──
  p <- ggplot() +

    # Nube de puntos males (arriba, y > 0) — rasterizada
    rasterise(
      geom_point(data = males_dt,
                 aes(x = pos_cumulative, y = y, colour = color),
                 size = 0.15, shape = 16, stroke = 0),
      dpi = dpi
    ) +

    # Nube de puntos females (abajo, y < 0) — rasterizada
    rasterise(
      geom_point(data = females_dt,
                 aes(x = pos_cumulative, y = y, colour = color),
                 size = 0.15, shape = 16, stroke = 0),
      dpi = dpi
    ) +

    scale_colour_identity() +

    # Líneas de significancia — males (arriba, con offset)
    geom_hline(yintercept = sig_y + y_offset, linetype = "dashed",
               linewidth = 0.6 / .pt, colour = "black") +
    geom_hline(yintercept = sug_y + y_offset, linetype = "dotted",
               linewidth = 0.5 / .pt, colour = "grey50") +

    # Líneas de significancia — females (abajo, espejadas con offset)
    geom_hline(yintercept = -(sig_y + y_offset), linetype = "dashed",
               linewidth = 0.6 / .pt, colour = "black") +
    geom_hline(yintercept = -(sug_y + y_offset), linetype = "dotted",
               linewidth = 0.5 / .pt, colour = "grey50") +

    # Etiquetas de cromosoma centradas en y = 0
    geom_text(data = chr_label_df,
              aes(x = x, y = 0, label = label),
              size = 5 / .pt, family = font_family, colour = "black",
              vjust = 0.5, hjust = 0.5) +

    # Labels de grupo
    annotate("text", x = -Inf, y = Inf,
             label = "Males (N=35,255)",
             hjust = -0.02, vjust = 1.4,
             size = 6 / .pt, family = font_family, colour = "black") +
    annotate("text", x = -Inf, y = -Inf,
             label = "Females (N=71,496)",
             hjust = -0.02, vjust = -0.6,
             size = 6 / .pt, family = font_family, colour = "black") +

    # Ejes
    scale_x_continuous(expand = expansion(mult = 0.01)) +
    scale_y_continuous(
      breaks = y_breaks,
      labels = y_labels,
      limits = c(y_bottom, y_top),
      expand = expansion(mult = 0.01)
    ) +
    labs(x = NULL, y = "-log10(P)") +
    nature_theme

  # ── Anotaciones de genes (ggrepel) ──
  if (!is.null(sig_annot) && nrow(sig_annot) > 0) {

    # Añadir pos_cumulative y log10p desde cada dataset
    sig_annot <- merge(
      sig_annot,
      males_dt[, .(CHR, POS, pos_cumulative_m = pos_cumulative, log10p_m = log10p)],
      by = c("CHR", "POS"), all.x = TRUE
    )
    sig_annot <- merge(
      sig_annot,
      females_dt[, .(CHR, POS, pos_cumulative_f = pos_cumulative, log10p_f = log10p)],
      by = c("CHR", "POS"), all.x = TRUE
    )

    # Males hits (y positivo con offset)
    m_hits <- sig_annot[dataset == "Males" & !is.na(pos_cumulative_m)]
    if (nrow(m_hits) > 0) {
      m_hits[, `:=`(x_plot = pos_cumulative_m, y_plot = log10p_m + y_offset)]
      p <- p +
        geom_text_repel(
          data = m_hits,
          aes(x = x_plot, y = y_plot, label = nearest_gene),
          size = 5 / .pt, fontface = "italic",
          family = font_family, colour = "black",
          segment.size = 0.3 / .pt, segment.colour = "grey40",
          max.overlaps = Inf, min.segment.length = 0,
          box.padding = 0.3, point.padding = 0.2,
          force = 2, force_pull = 0.5,
          nudge_y = 2, direction = "both",
          seed = 42
        )
    }

    # Females hits (y negativo con offset)
    f_hits <- sig_annot[dataset == "Females" & !is.na(pos_cumulative_f)]
    if (nrow(f_hits) > 0) {
      f_hits[, `:=`(x_plot = pos_cumulative_f, y_plot = -(log10p_f + y_offset))]
      p <- p +
        geom_text_repel(
          data = f_hits,
          aes(x = x_plot, y = y_plot, label = nearest_gene),
          size = 5 / .pt, fontface = "italic",
          family = font_family, colour = "black",
          segment.size = 0.3 / .pt, segment.colour = "grey40",
          max.overlaps = Inf, min.segment.length = 0,
          box.padding = 0.3, point.padding = 0.2,
          force = 2, force_pull = 0.5,
          nudge_y = -2, direction = "both",
          seed = 42
        )
    }
  }

  # ── Guardar PDF ──
  width_in  <- width_mm / MM_PER_IN
  height_in <- height_mm / MM_PER_IN

  device <- grDevices::cairo_pdf
  message("Usando cairo_pdf (embebido TrueType)")

  dir.create(dirname(output_pdf), recursive = TRUE, showWarnings = FALSE)
  ggsave(
    filename = output_pdf,
    plot     = p,
    width    = width_in,
    height   = height_in,
    units    = "in",
    device   = device,
    dpi      = dpi
  )
  message("PDF guardado: ", output_pdf)

  # Verificación con pdfinfo/pdffonts si están disponibles
  if (Sys.which("pdfinfo") != "") {
    message("\n--- pdfinfo ---")
    system2("pdfinfo", output_pdf)
  }
  if (Sys.which("pdffonts") != "") {
    message("\n--- pdffonts ---")
    system2("pdffonts", output_pdf)
  }

  invisible(p)
}


# ============================================================
# MAIN
# ============================================================
main <- function() {
  t_total <- proc.time()

  cat("============================================================\n")
  cat("  Miami Plot — R / Nature specs\n")
  cat("============================================================\n")

  # Verificar archivos
  for (f in c(MALE_FILE, FEMALE_FILE, GTF_FILE)) {
    if (!file.exists(f)) {
      stop("Archivo no encontrado: ", f, "\n  Directorio actual: ", getwd())
    }
  }

  # 0. Fuente
  font_family <- resolve_font()

  # 1. Cargar GWAS
  cat("\n--- PASO 1: Carga de datos GWAS ---\n")
  males_dt   <- load_gwas(MALE_FILE)
  females_dt <- load_gwas(FEMALE_FILE)

  # 2. Coordenadas
  cat("\n--- PASO 2: Coordenadas acumulativas ---\n")
  males_res   <- compute_coordinates(males_dt)
  females_res <- compute_coordinates(females_dt)

  # 3. SNPs significativos + anotación GTF
  cat("\n--- PASO 3: SNPs significativos + anotacion GTF ---\n")
  m_sig <- males_res$dt[pvalue < SIGNIFICANCE][, dataset := "Males"]
  f_sig <- females_res$dt[pvalue < SIGNIFICANCE][, dataset := "Females"]
  all_sig <- rbind(m_sig, f_sig)
  message(sprintf("  Significativos: Males=%s, Females=%s, Total=%s",
                  format(nrow(m_sig), big.mark = ","),
                  format(nrow(f_sig), big.mark = ","),
                  format(nrow(all_sig), big.mark = ",")))

  sig_annot <- NULL
  if (nrow(all_sig) > 0) {
    genes_dt <- load_gtf_genes(GTF_FILE)
    all_sig  <- annotate_snps(all_sig, genes_dt)

    # 4. Selección de top SNPs
    cat("\n--- PASO 4: Seleccion de top SNPs ---\n")
    top_snps <- select_top_snps(all_sig)
    if (nrow(top_snps) > 0) {
      sig_annot <- select_best_for_annotation(top_snps)
      message(sprintf("  SNPs para anotacion: %d", nrow(sig_annot)))

      # Guardar CSV
      fwrite(top_snps, OUTPUT_CSV)
      message("CSV guardado: ", OUTPUT_CSV)
    }
  }

  # 5. Miami plot
  cat("\n--- PASO 5: Generacion del Miami Plot ---\n")
  build_miami_plot(males_res, females_res, sig_annot, font_family)

  elapsed <- (proc.time() - t_total)[3]
  cat("\n============================================================\n")
  cat(sprintf("  Males:          %s SNPs\n", format(nrow(males_res$dt), big.mark = ",")))
  cat(sprintf("  Females:        %s SNPs\n", format(nrow(females_res$dt), big.mark = ",")))
  cat(sprintf("  Significativos: %s\n", format(nrow(all_sig), big.mark = ",")))
  cat(sprintf("  Genes anotados: %d\n", if (!is.null(sig_annot)) nrow(sig_annot) else 0))
  cat(sprintf("  Salida PDF:     %s\n", OUTPUT_PDF))
  cat(sprintf("  Salida CSV:     %s\n", OUTPUT_CSV))
  cat(sprintf("  Tiempo total:   %.1fs (%.1f min)\n", elapsed, elapsed / 60))
  cat("============================================================\n")
}

main()
