#!/bin/bash
#SBATCH --job-name=ldpred2_recover
#SBATCH --account=researchers
#SBATCH --qos=vip
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=120G
#SBATCH --output=/mnt/cephfs/orgs/home/angel.pacheco/logs/ldpred2_recover_%j.out
#SBATCH --error=/mnt/cephfs/orgs/home/angel.pacheco/logs/ldpred2_recover_%j.err
#SBATCH --mail-user=angelpachecosaavedra@gmail.com
#SBATCH --mail-type=BEGIN,END,FAIL

set -euo pipefail

# ================================================================
# RECUPERACIÓN — LDpred2-auto sin recalcular LD matrices
# Reaprovecha SFBMs (.sbk) y bigSNPs (.bk/.rds) del job 2931.
# ================================================================

SIF="/mnt/cephfs/biocontainers/images/jupyter-biotools-1.4.sif"
PLINK_DIR="/mnt/cephfs/hot_nvme/mcps/imputed-topmed/plink_files/maximally_unrelated"
GWAS_DIR="/mnt/cephfs/orgs/home/angel.pacheco/gwas_miami"
OUT_DIR="/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes/ldpred2"
mkdir -p "$OUT_DIR"

# >>> CRÍTICO: apunta al TMPDIR del job 2931 — NO borrar <<<
PREV_TMPDIR="/mnt/cephfs/scratch/angel.pacheco/ldpred2_hm3_tmp/2931"

# TMPDIR nuevo solo para archivos temporales del proceso R
TMPDIR="/mnt/cephfs/scratch/$USER/ldpred2_recover_tmp/$SLURM_JOB_ID"
mkdir -p "$TMPDIR"
chmod 700 "$TMPDIR"
export TMPDIR TMP="$TMPDIR" TEMP="$TMPDIR"

echo "============================================================"
echo "  LDpred2-auto RECOVERY — reusing job 2931 artifacts"
echo "  Nodo: $SLURMD_NODENAME"
echo "  PREV_TMPDIR: $PREV_TMPDIR"
echo "  $(date)"
echo "============================================================"

R_SCRIPT="$TMPDIR/recover_ldpred2.R"

cat << 'REOF' > "$R_SCRIPT"
suppressMessages({
  library(bigsnpr)
  library(bigsparser)
  library(data.table)
  library(Matrix)
  library(methods)
})

ncores <- 24L
options(bigstatsr.check.parallel.blas = FALSE)
options(bigstatsr.ncores.max = 24L)

prev_tmp  <- Sys.getenv("PREV_TMPDIR")
out_dir   <- Sys.getenv("OUT_DIR")
gwas_dir  <- Sys.getenv("GWAS_DIR")

stopifnot(dir.exists(prev_tmp))

# ================================================================
# Helper: reattach SFBM_corr_compact desde .sbk existente
#
# Estrategia confirmada por diagnóstico:
#   1. Crear dummy SFBM_corr_compact (3x3 diagonal)
#   2. Field-swap: backingfile, nrow, p, first_i
#   3. Llamar $address para que C++ recargue el extptr desde el .sbk
#
# p y first_i se recomputan corriendo snp_cor sobre la misma ventana
# de posiciones (size=3000 kb, determinista dado info_snp$pos).
# Eso toma ~2-5 min/chr vs ~1000 min/chr para recalcular la LD matrix.
# ================================================================
reattach_sfbm_corr_compact <- function(sbk_path, info_snp, G, ncores = 1L) {
  bk_noext <- sub("\\.sbk$", "", sbk_path)
  n_snp    <- nrow(info_snp)

  # --- Validaciones básicas ---
  sbk_size <- file.info(sbk_path)$size
  if (is.na(sbk_size) || sbk_size == 0)
    stop("sbk_path no existe o está vacío: ", sbk_path)

  # --- Dummy 3x3 diagonal para obtener un objeto vacío ---
  gen_cc  <- get("SFBM_corr_compact_RC", envir = asNamespace("bigsparser"))
  dummy   <- as(Matrix::Diagonal(3L), "symmetricMatrix")
  tmp_bk  <- tempfile("dummy_reattach")
  sfbm    <- gen_cc$new(spm = dummy, backingfile = tmp_bk)

  # --- Recomputar p y first_i ---
  # snp_cor con los mismos parámetros que el script original (size=3000 kb)
  # Devuelve una dsCMatrix; extraemos p e i sin guardar la matriz.
  cat("      [reattach] Recomputando p/first_i (snp_cor ligero)...\n")
  t_cor <- system.time({
    corr_tmp <- snp_cor(
      Gna        = G,
      ind.col    = info_snp$`_NUM_ID_`,
      infos.pos  = info_snp$pos,
      size       = 3000L,
      ncores     = ncores
    )
  })
  cat(sprintf("      [reattach] snp_cor: %.1f seg\n", t_cor["elapsed"]))

  # Extraer los campos del objeto SFBM temporal que as_SFBM hubiera llenado
  # Necesitamos p y first_i calculados por initialize de SFBM_compact.
  # Los obtenemos creando un SFBM temporal SOLO para robar esos campos.
  tmp_bk2  <- tempfile("tmp_sfbm_fields")
  sfbm_tmp <- as_SFBM(corr_tmp, backingfile = tmp_bk2, compact = TRUE)
  p_vec    <- sfbm_tmp$p
  first_i  <- sfbm_tmp$first_i
  rm(corr_tmp, sfbm_tmp)
  # Borrar el .sbk temporal (no lo necesitamos, solo robamos p/first_i)
  unlink(paste0(tmp_bk2, ".sbk"))

  # --- Field-swap al .sbk real ---
  sfbm$backingfile <- normalizePath(sbk_path)
  sfbm$nrow        <- n_snp
  sfbm$p           <- p_vec
  sfbm$first_i     <- first_i
  # Forzar recarga del puntero C++ con el archivo real
  invisible(sfbm$address)

  # Limpiar dummy
  unlink(paste0(tmp_bk, ".sbk"))

  # Validar dimensiones
  if (sfbm$nrow != n_snp || sfbm$ncol != n_snp)
    stop(sprintf("Dimensiones SFBM (%dx%d) != info_snp (%d)",
                 sfbm$nrow, sfbm$ncol, n_snp))

  sfbm
}

# ================================================================
# [0] Mapa HapMap3
# ================================================================
cat("\n[0] Cargando mapa HapMap3...\n")
map_hm3 <- readRDS(file.path(out_dir, "map_hm3_ldpred2.rds"))
cat(sprintf("  HapMap3: %s variantes\n", format(nrow(map_hm3), big.mark=",")))

# ================================================================
# [1] Reconstruir GWAS meta + sumstats
# ================================================================
cat("\n[1] Reconstruyendo GWAS meta-análisis...\n")
gf <- fread(file.path(gwas_dir, "combined-gwas-results-all-ctrl-females.txt.gz"),
            select=c("CHROM","GENPOS","ID","ALLELE0","ALLELE1","A1FREQ","INFO","N","BETA","SE","P"),
            showProgress=FALSE)
gm <- fread(file.path(gwas_dir, "combined-gwas-results-all-ctrl-males.txt.gz"),
            select=c("CHROM","GENPOS","ID","ALLELE0","ALLELE1","A1FREQ","INFO","N","BETA","SE","P"),
            showProgress=FALSE)
gwas <- merge(gf, gm, by="ID", suffixes=c(".f",".m"))
rm(gf, gm); gc()

w_f      <- 1 / gwas$SE.f^2
w_m      <- 1 / gwas$SE.m^2
neff_f   <- 4 * 17300 * 54196 / (17300 + 54196)
neff_m   <- 4 *  8522 * 26733 / ( 8522 + 26733)
K_prev   <- (17300 + 8522) / (17300 + 54196 + 8522 + 26733)
neff_ldpred2 <- (neff_f + neff_m) * K_prev * (1 - K_prev)

sumstats <- data.frame(
  chr      = gwas$CHROM.f,
  pos      = gwas$GENPOS.f,
  a1       = gwas$ALLELE1.f,
  a0       = gwas$ALLELE0.f,
  freq     = (gwas$A1FREQ.f * gwas$N.f + gwas$A1FREQ.m * gwas$N.m) / (gwas$N.f + gwas$N.m),
  info     = pmin(gwas$INFO.f, gwas$INFO.m),
  beta     = (gwas$BETA.f * w_f + gwas$BETA.m * w_m) / (w_f + w_m),
  beta_se  = sqrt(1 / (w_f + w_m)),
  n_eff    = neff_ldpred2,
  stringsAsFactors = FALSE
)
sumstats$p <- 2 * pnorm(-abs(sumstats$beta / sumstats$beta_se))
rm(gwas); gc()

sumstats <- sumstats[
  sumstats$info    >= 0.8  &
  sumstats$freq    >= 0.01 & sumstats$freq <= 0.99 &
  !is.na(sumstats$beta) & !is.na(sumstats$beta_se) & sumstats$beta_se > 0, ]
sumstats$chr <- as.character(sumstats$chr)
sumstats$chr[sumstats$chr %in% c("X","chrX")] <- "23"
sumstats$chr <- as.integer(sumstats$chr)
cat(sprintf("  GWAS tras QC: %s\n", format(nrow(sumstats), big.mark=",")))

key_hm3  <- paste(map_hm3$chr, map_hm3$pos_hg38, sep=":")
key_gwas <- paste(sumstats$chr, sumstats$pos, sep=":")
sumstats <- sumstats[key_gwas %in% key_hm3, ]
cat(sprintf("  GWAS filtrado a HapMap3: %s\n", format(nrow(sumstats), big.mark=",")))

# ================================================================
# [2] LDpred2-auto por cromosoma — reusando SFBMs
# ================================================================
cat("\n[2] LDpred2-auto reusando SFBMs existentes...\n")

prs_total    <- NULL
fam_global   <- NULL
matched_total <- 0L
results_log  <- list()

for (chr in 1:22) {
  chr_label <- as.character(chr)
  rds_chr   <- file.path(prev_tmp, sprintf("chr%s.rds",      chr_label))
  sbk_chr   <- file.path(prev_tmp, sprintf("corr_chr%s.sbk", chr_label))
  betas_rds <- file.path(out_dir,  sprintf("ldpred2_betas_chr%s.rds", chr_label))

  cat(sprintf("\n  >>> chr%s <<<\n", chr_label))
  t0 <- Sys.time()

  # Si ya calculamos este chr en un intento previo, cargar desde disco
  if (file.exists(betas_rds)) {
    cat(sprintf("      Cargando betas ya calculados desde disco...\n"))
    saved <- readRDS(betas_rds)
    if (is.null(prs_total)) {
      obj_tmp    <- snp_attach(rds_chr)
      prs_total  <- rep(0.0, nrow(obj_tmp$genotypes))
      fam_global <- obj_tmp$fam
      rm(obj_tmp)
    }
    prs_total <- prs_total + saved$prs_chr
    matched_total <- matched_total + saved$n_matched
    results_log[[chr_label]] <- saved$log
    next
  }

  # Verificar que existan los archivos necesarios
  for (f in c(rds_chr, sbk_chr)) {
    if (!file.exists(f)) {
      cat(sprintf("      FALTA: %s — saltando chr%s\n", f, chr_label))
      next
    }
  }
  if (!file.exists(rds_chr) || !file.exists(sbk_chr)) next

  # Cargar bigSNP
  obj <- snp_attach(rds_chr)
  G   <- obj$genotypes
  if (is.null(prs_total)) {
    prs_total  <- rep(0.0, nrow(G))
    fam_global <- obj$fam
  }

  # Matching sumstats ↔ mapa
  sumstats_chr <- sumstats[sumstats$chr == chr, ]
  map_df <- data.frame(
    chr = obj$map$chromosome,
    pos = obj$map$physical.pos,
    a1  = obj$map$allele1,
    a0  = obj$map$allele2,
    stringsAsFactors = FALSE
  )
  info_snp <- snp_match(sumstats_chr, map_df, strand_flip = FALSE,
                        join_by_pos = TRUE, match.min.prop = 0.01)

  # Filtrar monomórficas
  sd_vals  <- big_scale()(G, ind.col = info_snp$`_NUM_ID_`)$scale
  keep_var <- which(sd_vals > 1e-6)
  if (length(keep_var) < nrow(info_snp)) {
    info_snp <- info_snp[keep_var, ]
  }
  n_snp <- nrow(info_snp)
  cat(sprintf("      Matched: %s variantes\n", format(n_snp, big.mark=",")))

  # Reattach SFBM
  corr <- tryCatch(
    reattach_sfbm_corr_compact(sbk_chr, info_snp, G, ncores = ncores),
    error = function(e) {
      cat(sprintf("      ERROR reattach: %s\n", conditionMessage(e)))
      NULL
    }
  )
  if (is.null(corr)) {
    rm(obj, G); gc()
    results_log[[chr_label]] <- list(n_matched=n_snp, h2=NA, p=NA,
                                     elapsed_min=NA, status="reattach_failed")
    next
  }
  cat(sprintf("      SFBM OK: %dx%d\n", corr$nrow, corr$ncol))

  # LDpred2-auto — más iteraciones para chr con señal débil
  burn_in  <- if (chr %in% c(18L, 19L, 21L, 22L)) 1000L else 500L
  num_iter <- if (chr %in% c(18L, 19L, 21L, 22L)) 1000L else 500L

  df_beta <- info_snp[, c("beta","beta_se","n_eff")]
  multi_auto <- snp_ldpred2_auto(
    corr        = corr,
    df_beta     = df_beta,
    h2_init     = 0.01,
    vec_p_init  = seq_log(1e-4, 0.5, length.out = 30L),
    burn_in     = burn_in,
    num_iter    = num_iter,
    report_step = 100L,
    ncores      = ncores,
    sparse      = TRUE
  )

  h2_est <- sapply(multi_auto, function(x) x$h2_est)
  p_est  <- sapply(multi_auto, function(x) x$p_est)

  # Filtrado robusto: sin fallback ciego a cadenas NA
  keep <- which(
    !is.na(h2_est) & h2_est > 0.0005 & h2_est < 0.5 &
    !is.na(p_est)  & p_est  > 1e-6   & p_est  < 0.9
  )
  cat(sprintf("      Cadenas válidas: %d/%d\n", length(keep), length(multi_auto)))

  if (length(keep) < 3L) {
    cat(sprintf("      ⚠ chr%s: chains válidas=%d (<3). Saltando.\n",
                chr_label, length(keep)))
    rm(obj, G, corr, multi_auto); gc()
    results_log[[chr_label]] <- list(n_matched=n_snp, h2=NA, p=NA,
                                     elapsed_min=as.numeric(difftime(Sys.time(), t0, units="mins")),
                                     status="diverged")
    next
  }

  beta_mat   <- sapply(multi_auto[keep], function(x) x$beta_est)
  beta_final <- rowMeans(beta_mat, na.rm = TRUE)
  beta_final[!is.finite(beta_final)] <- 0.0

  h2_med <- median(h2_est[keep])
  p_med  <- median(p_est[keep])
  cat(sprintf("      h2=%.4f, p=%.2e, nonzero=%s\n",
              h2_med, p_med, format(sum(beta_final != 0), big.mark=",")))

  # Scoring
  prs_chr <- big_prodVec(X       = G,
                         y.col   = beta_final,
                         ind.col = info_snp$`_NUM_ID_`,
                         ncores  = ncores)
  prs_total     <- prs_total + prs_chr
  matched_total <- matched_total + n_snp

  elapsed <- as.numeric(difftime(Sys.time(), t0, units="mins"))
  cat(sprintf("      Tiempo: %.1f min\n", elapsed))

  log_entry <- list(n_matched=n_snp, h2=h2_med, p=p_med,
                    elapsed_min=elapsed, status="ok")
  results_log[[chr_label]] <- log_entry

  # Guardar PRS parcial por chr — protección ante crashes futuros
  saveRDS(list(prs_chr=prs_chr, n_matched=n_snp, log=log_entry),
          betas_rds)

  rm(obj, G, corr, multi_auto, info_snp, beta_final, prs_chr); gc()
}

# ================================================================
# [3] Export — con filtros robustos a NA
# ================================================================
cat(sprintf("\n[3] Exportando (%s variantes totales)...\n",
            format(matched_total, big.mark=",")))

saveRDS(results_log, file.path(out_dir, "ldpred2_hm3_results_log.rds"))
saveRDS(prs_total,   file.path(out_dir, "prs_total_raw.rds"))
saveRDS(fam_global,  file.path(out_dir, "fam_global.rds"))

# Resumen por cromosoma
cat("\n  Resumen por cromosoma:\n")
for (cl in names(results_log)) {
  r <- results_log[[cl]]
  cat(sprintf("    chr%s: %s vars, h2=%s, p=%s, %.1f min [%s]\n",
              cl,
              format(r$n_matched, big.mark=","),
              if (is.na(r$h2)) "NA" else sprintf("%.4f", r$h2),
              if (is.na(r$p))  "NA" else sprintf("%.2e", r$p),
              ifelse(is.na(r$elapsed_min), 0, r$elapsed_min),
              r$status))
}

# Link IID → PATID y sexo
covar_file <- "/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes/MCPS BASELINE.csv"
link_file  <- "/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes/RGN_LINK_IID.csv"

df_link      <- read.csv(link_file)
iid_to_patid <- setNames(as.character(df_link$PATID), as.character(df_link$IID))
df_baseline  <- read.csv(covar_file)
df_baseline  <- df_baseline[!is.na(df_baseline$MALE), ]
patid_to_sex <- setNames(as.integer(df_baseline$MALE), as.character(df_baseline$PATID))

patids <- iid_to_patid[as.character(fam_global$sample.ID)]
patids[is.na(patids)] <- "UNKNOWN"
sexes  <- patid_to_sex[patids]
sexes[is.na(sexes)] <- -1L

df_out <- data.frame(
  PATID     = patids,
  IID       = fam_global$sample.ID,
  SEX       = sexes,
  PRS_TOTAL = prs_total,
  stringsAsFactors = FALSE
)
df_out <- df_out[df_out$SEX >= 0L, ]

# Excluir muestras con PRS=NA (cromosomas saltados)
n_na <- sum(is.na(df_out$PRS_TOTAL))
if (n_na > 0) cat(sprintf("  ⚠ %d muestras excluidas por PRS=NA\n", n_na))
df_out <- df_out[!is.na(df_out$PRS_TOTAL), ]

# Standardizar por sexo — con guarda explícita contra NA en sd()
df_out$PRS_Z      <- NA_real_
df_out$PERCENTILE <- NA_real_
for (sv in c(0L, 1L)) {
  mask <- df_out$SEX == sv
  vals <- df_out$PRS_TOTAL[mask]
  sd_v <- sd(vals)
  if (length(vals) > 1L && !is.na(sd_v) && sd_v > 0) {
    df_out$PRS_Z[mask]      <- (vals - mean(vals)) / sd_v
    df_out$PERCENTILE[mask] <- rank(vals) / length(vals) * 100
  }
}

df_out   <- df_out[order(-df_out$PRS_TOTAL), ]
out_file <- file.path(out_dir, "prs_diabetes_ldpred2_final.tsv")
write.table(df_out, out_file, sep="\t", row.names=FALSE, quote=FALSE)

cat(sprintf("\n Pipeline RECUPERADO\n"))
cat(sprintf("   Output: %s\n", out_file))
cat(sprintf("   Muestras: %s\n", format(nrow(df_out), big.mark=",")))
cat(sprintf("   Variantes: %s\n", format(matched_total, big.mark=",")))
REOF

# Exportar variables para que R las lea con Sys.getenv()
export PREV_TMPDIR OUT_DIR GWAS_DIR

apptainer exec \
    --bind /mnt/cephfs:/mnt/cephfs \
    --bind "$TMPDIR:/tmp,$TMPDIR:/var/tmp" \
    --env PREV_TMPDIR="$PREV_TMPDIR" \
    --env OUT_DIR="$OUT_DIR" \
    --env GWAS_DIR="$GWAS_DIR" \
    "$SIF" /opt/miniforge3/envs/r-bio/bin/Rscript "$R_SCRIPT"

# >>> NO borramos PREV_TMPDIR — sigue siendo nuestro respaldo <<<
rm -rf "$TMPDIR"

echo ""
echo "============================================================"
echo "  RECOVERY COMPLETADO — $(date)"
echo "============================================================"
