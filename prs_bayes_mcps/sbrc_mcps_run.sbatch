#!/bin/bash
#SBATCH --job-name=sbrc_mcps
#SBATCH --account=researchers
#SBATCH --qos=vip
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=150G
#SBATCH --output=/mnt/cephfs/orgs/home/angel.pacheco/logs/sbrc_mcps_run_%j.out
#SBATCH --error=/mnt/cephfs/orgs/home/angel.pacheco/logs/sbrc_mcps_run_%j.err
#SBATCH --mail-user=angelpachecosaavedra@gmail.com
#SBATCH --mail-type=BEGIN,END,FAIL

set -euo pipefail

# ================================================================
# SBayesRC con GWAS MCPS — pipeline coherente
#
#   Betas (β)   : meta-análisis MCPS sex-stratified (REGENIE F+M)
#   Panel LD    : MCPS in-house (ldm_79k, N=79,612 maximally unrelated)
#   Variantes   : intersección GWAS × panel LD
#
# Output: results_mcps/sbayesrc_mcps.{txt,par,rds,mcmcsamples}
#
# Diferencia vs sbrc_run_2932 (T2DGGI): mismo panel LD, distinto GWAS.
# ================================================================

SIF="/mnt/cephfs/biocontainers/images/jupyter-biotools-1.4.sif"
SBRC_SIF="/mnt/cephfs/biocontainers/images/sbayesrc-0.2.6.sif"

OUT_BASE="/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes/sbayesrc"
LD_DIR="${OUT_BASE}/ldm_79k"
MA_FILE="${OUT_BASE}/MCPS_meta_matched.ma"

RESULTS_DIR="${OUT_BASE}/results_mcps"
mkdir -p "$RESULTS_DIR"

OUT_PREFIX="${RESULTS_DIR}/sbayesrc_mcps"

TMPDIR="/mnt/cephfs/scratch/$USER/sbrc_mcps_tmp/$SLURM_JOB_ID"
mkdir -p "$TMPDIR"; chmod 700 "$TMPDIR"
export TMPDIR TMP="$TMPDIR" TEMP="$TMPDIR"

echo "============================================================"
echo "  PIPELINE A2 — SBayesRC con GWAS MCPS"
echo "  Nodo: $SLURMD_NODENAME"
echo "  RAM: $(free -g | awk 'NR==2{print $2}')G"
echo "  $(date)"
echo "============================================================"

# ================================================================
# Pre-flight checks
# ================================================================
echo ""
echo "--- Pre-flight ---"

if [[ ! -f "$MA_FILE" ]]; then
  echo "ERROR: No existe $MA_FILE"
  echo "Primero corre: sbatch sbrc_mcps_generate_ma.sbatch"
  exit 1
fi

if [[ ! -d "$LD_DIR" ]]; then
  echo "ERROR: No existe $LD_DIR"
  exit 1
fi

if [[ ! -f "$LD_DIR/ldm.info" ]] || [[ ! -f "$LD_DIR/snp.info" ]]; then
  echo "ERROR: Falta ldm.info o snp.info en $LD_DIR"
  exit 1
fi

N_BLOCKS=$(tail -n +2 "$LD_DIR/ldm.info" | wc -l)
N_SNPS_LD=$(tail -n +2 "$LD_DIR/snp.info" | wc -l)
N_VARS_MA=$(tail -n +2 "$MA_FILE" | wc -l)

echo "  ✓ MA file       : $MA_FILE"
echo "    variantes     : $(printf "%'d" $N_VARS_MA)"
echo "  ✓ LD panel      : $LD_DIR"
echo "    bloques       : $(printf "%'d" $N_BLOCKS)"
echo "    variantes     : $(printf "%'d" $N_SNPS_LD)"
echo "  ✓ Output prefix : $OUT_PREFIX"

# ================================================================
# PASO 1: SBayesRC MCMC con GWAS MCPS
# ================================================================
echo ""
echo "============================================================"
echo "  SBayesRC MCMC (sin anotaciones, modo SBayesR)"
echo "============================================================"

R_SCRIPT="$TMPDIR/run_sbrc_mcps.R"

cat << REOF > "$R_SCRIPT"
suppressMessages({ library(SBayesRC); library(data.table) })

ma_file   <- "$MA_FILE"
ld_dir    <- "$LD_DIR"
out_prefix <- "$OUT_PREFIX"

cat("\n[1] Pre-validación del .ma...\n")
ma_head <- fread(ma_file, nrows=10)
print(ma_head)
cat(sprintf("\n  Columnas: %s\n", paste(names(ma_head), collapse=", ")))

# Verificar columnas esperadas
required_cols <- c("SNP","A1","A2","freq","b","se","p","N")
missing <- setdiff(required_cols, names(ma_head))
if (length(missing) > 0) {
  stop("Faltan columnas en .ma: ", paste(missing, collapse=", "))
}

# Estadísticas rápidas
ma_full <- fread(ma_file)
cat(sprintf("\n  Total variantes en .ma: %s\n", format(nrow(ma_full), big.mark=",")))
cat(sprintf("  N (sample size): %s\n", format(unique(ma_full\$N), big.mark=",")))
cat(sprintf("  beta range: [%.3e, %.3e]\n", min(ma_full\$b), max(ma_full\$b)))
cat(sprintf("  beta median |b|: %.3e\n", median(abs(ma_full\$b))))

# Distribución por chr
# --- NUEVO DIAGNÓSTICO DE SOLAPAMIENTO ---

ld_snps <- fread(file.path(ld_dir, "snp.info"), select = "ID")
overlap <- sum(ma_full\$SNP %in% ld_snps\$ID)



cat(sprintf("\n[1.1] Chequeo de solapamiento GWAS vs LD:\n"))
cat(sprintf("  SNPs en común: %s (%.2f%% del GWAS)\n", 
            format(overlap, big.mark=","), (overlap/nrow(ma_full))*100))

if (overlap < 100) {
  cat("\n!!! ERROR: Solapamiento casi nulo !!!\n")
  cat("  Ejemplo GWAS : ", ma_full\$SNP[1], "\n")
  cat("  Ejemplo LD   : ", ld_snps\$ID[1], "\n")
  stop("El formato de los SNPs no coincide (ej. prefijo 'chr' vs numérico puro).")
}
# ----------------------------------------



rm(ma_full); gc()
ma_file_imp <- paste0(out_prefix, "_imp.ma")
if (file.exists(ma_file_imp)) {
  cat("[1.5] .ma imputado ya existe — reutilizando\n")
  ma_file <- ma_file_imp
} else {
  file.copy(file.path(ld_dir, "snp.info"), file.path(ld_dir, "snp.info.16col"), overwrite=TRUE)
  file.copy(file.path(ld_dir, "snp.info.10col"), file.path(ld_dir, "snp.info"), overwrite=TRUE)

  # Imputar variantes faltantes del LD panel
  cat("[1.5] Imputando variantes faltantes...\n")
  SBayesRC::impute(
      mafile = ma_file, 
      LDdir  = ld_dir,
      output = paste0(out_prefix, "_imp.ma")  # <-- Cambio clave: 'output' en lugar de 'outPrefix'
  )
  file.copy(file.path(ld_dir, "snp.info.16col"), file.path(ld_dir, "snp.info"), overwrite=TRUE)
}
# Usar el .ma imputado en lugar del original
ma_file <- paste0(out_prefix, "_imp.ma")
cat(sprintf("  Variantes tras imputación: %s\n",
    format(length(readLines(ma_file)) - 1L, big.mark=",")))

cat("\n[2] Ejecutando SBayesRC::sbayesrc()...\n")
cat(sprintf("  Inicio: %s\n", Sys.time()))

# Parámetros: igual que el job 2932 que SÍ convergió
result <- SBayesRC::sbayesrc(
  mafile     = ma_file,
  LDdir      = ld_dir,
  outPrefix  = out_prefix,
  niter      = 3000,           # mismo que job 2932
  burn       = 1000,           # mismo que job 2932
  thresh     = 0.9,            # cutoff varianza eigen, igual al estándar
  bTune      = FALSE,
  log2file   = FALSE
)

cat(sprintf("  Fin: %s\n", Sys.time()))
cat("\n[3] Output esperado:\n")
output_files <- list.files(dirname(out_prefix), pattern=basename(out_prefix), full.names=TRUE)
for (f in output_files) {
  size_mb <- file.info(f)\$size / 1024^2
  cat(sprintf("  %s (%.1f MB)\n", basename(f), size_mb))
}

# Cargar pesos generados
snp_file <- paste0(out_prefix, ".txt")
if (file.exists(snp_file)) {
  snp_res <- fread(snp_file)
  n_total <- nrow(snp_res)
  effect_col <- intersect(c("A1Effect","BETA"), names(snp_res))[1]
  if (!is.na(effect_col)) {
    n_nz <- sum(snp_res[[effect_col]] != 0, na.rm=TRUE)
    cat(sprintf("\n[4] Resumen final:\n"))
    cat(sprintf("  Variantes total: %s\n", format(n_total, big.mark=",")))
    cat(sprintf("  Con beta != 0  : %s (%.2f%%)\n",
                format(n_nz, big.mark=","), n_nz/n_total*100))
    nz_beta <- snp_res[[effect_col]][snp_res[[effect_col]] != 0]
    cat(sprintf("  |beta| range   : [%.3e, %.3e]\n",
                min(abs(nz_beta), na.rm=TRUE),
                max(abs(nz_beta), na.rm=TRUE)))
  }
}

# Mostrar parámetros estimados
par_file <- paste0(out_prefix, ".par")
if (file.exists(par_file)) {
  cat("\n[5] Parámetros MCMC:\n")
  cat(paste(readLines(par_file), collapse="\n"), "\n")
}

cat("\n✓ SBayesRC con GWAS MCPS completado\n")
REOF

apptainer exec --bind /mnt/cephfs:/mnt/cephfs \
    --bind "$TMPDIR:/tmp,$TMPDIR:/var/tmp" \
    --env OMP_NUM_THREADS=24 \
    --env OPENBLAS_NUM_THREADS=24 \
    --env MKL_NUM_THREADS=24 \
    "$SBRC_SIF" Rscript "$R_SCRIPT"

# ================================================================
# PASO 2: Exportar weights para scoring (formato plink2 --score)
# ================================================================
echo ""
echo "============================================================"
echo "  Exportando weights para scoring"
echo "============================================================"

EXPORT_R="$TMPDIR/export_weights.R"
cat << REOF > "$EXPORT_R"
suppressMessages(library(data.table))

snp_file <- paste0("$OUT_PREFIX", ".txt")
if (!file.exists(snp_file)) stop("No existe ", snp_file)

df <- fread(snp_file)
cat(sprintf("Columnas detectadas: %s\n", paste(names(df), collapse=", ")))

# Detectar nombres de columnas (varían entre versiones)
id_col     <- intersect(c("Id","Name","SNP","ID"), names(df))[1]
effect_col <- intersect(c("A1Effect","BETA","Beta","beta"), names(df))[1]
freq_col   <- intersect(c("A1Frq","A1Freq","MAF","Freq"), names(df))[1]

stopifnot(!is.na(id_col), !is.na(effect_col))
cat(sprintf("ID col: %s, Effect col: %s\n", id_col, effect_col))

out_w <- "$RESULTS_DIR/sbayesrc_mcps_weights_for_scoring.tsv"
out_dt <- data.table(
  SNP  = df[[id_col]],
  A1   = df\$A1,
  A2   = df\$A2,
  BETA = df[[effect_col]]
)
if (!is.na(freq_col)) out_dt\$FREQ <- df[[freq_col]]
fwrite(out_dt, out_w, sep="\t")
cat(sprintf("✓ Weights exportados: %s\n", out_w))
cat(sprintf("  %s variantes total, %s con BETA!=0\n",
    format(nrow(out_dt), big.mark=","),
    format(sum(out_dt\$BETA != 0), big.mark=",")))
REOF

apptainer exec --bind /mnt/cephfs:/mnt/cephfs \
    "$SIF" /opt/miniforge3/envs/r-bio/bin/Rscript "$EXPORT_R"

rm -rf "$TMPDIR"

echo ""
echo "============================================================"
echo "  ✓ PIPELINE A2 COMPLETADO — $(date)"
echo ""
echo "  Outputs en $RESULTS_DIR/:"
echo "    sbayesrc_mcps.txt                       (betas + metadata)"
echo "    sbayesrc_mcps.par                       (params MCMC)"
echo "    sbayesrc_mcps.rds                       (objeto R)"
echo "    sbayesrc_mcps.mcmcsamples               (cadenas MCMC)"
echo "    sbayesrc_mcps_weights_for_scoring.tsv   (formato plink2)"
echo ""
echo "  Para scoring: usar plink2 --score con sbayesrc_mcps_weights_for_scoring.tsv"
echo "============================================================"
