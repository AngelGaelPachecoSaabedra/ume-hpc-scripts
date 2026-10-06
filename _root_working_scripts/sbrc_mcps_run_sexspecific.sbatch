#!/bin/bash
#SBATCH --job-name=sbrc_run_sex
#SBATCH --account=researchers
#SBATCH --qos=vip
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=150G
#SBATCH --output=/mnt/cephfs/orgs/home/angel.pacheco/logs/sbrc_run_sex_%j.out
#SBATCH --error=/mnt/cephfs/orgs/home/angel.pacheco/logs/sbrc_run_sex_%j.err
#SBATCH --mail-user=angelpachecosaavedra@gmail.com
#SBATCH --mail-type=END,FAIL

set -euo pipefail

# ================================================================
# SBayesRC POR SEXO (Pipeline D sex-specific)
# Lanzar DOS veces (pueden ir EN PARALELO gracias al panel aislado):
#   sbatch --export=SEX=females sbrc_mcps_run_sexspecific.sbatch
#   sbatch --export=SEX=males   sbrc_mcps_run_sexspecific.sbatch
#
# Requiere MCPS_<SEX>_matched.ma (generate_ma sex).
#
# >>> BLINDAJE DE CONCURRENCIA <<<
# El swap in-place de snp.info (16col<->10col) que pide impute() corromperia
# el panel LD COMPARTIDO si dos jobs lo tocan a la vez. Aqui cada job opera
# sobre una COPIA LOCAL aislada (symlinks a los binarios, copia de los .info),
# asi el panel compartido en ldm_79k NUNCA se modifica.
# ================================================================

SEX=${SEX:-females}

SIF="/mnt/cephfs/biocontainers/images/jupyter-biotools-1.4.sif"
SBRC_SIF="/mnt/cephfs/biocontainers/images/sbayesrc-0.2.6.sif"

OUT_BASE="/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes/sbayesrc"
LD_DIR="${OUT_BASE}/ldm_79k"                       # panel COMPARTIDO (read-only de facto)
MA_FILE="${OUT_BASE}/MCPS_${SEX}_matched.ma"
RESULTS_DIR="${OUT_BASE}/results_${SEX}"
OUT_PREFIX="${RESULTS_DIR}/sbayesrc_${SEX}"
mkdir -p "$RESULTS_DIR"

TMPDIR="/mnt/cephfs/scratch/$USER/sbrc_run_sex_tmp/${SLURM_JOB_ID}_${SEX}"
mkdir -p "$TMPDIR"; chmod 700 "$TMPDIR"
export TMPDIR TMP="$TMPDIR" TEMP="$TMPDIR"
trap 'rm -rf "$TMPDIR"' EXIT

echo "============================================================"
echo "  PIPELINE D sex — SBayesRC ${SEX} — Job $SLURM_JOB_ID"
echo "  Nodo: $SLURMD_NODENAME   $(date)"
echo "============================================================"

# ── Pre-flight ─────────────────────────────────────────────────────
echo ""
echo "--- Pre-flight ---"
[[ -f "$MA_FILE" ]] || { echo "ERROR: falta $MA_FILE (corre generate_ma sex)"; exit 1; }
[[ -d "$LD_DIR" ]]  || { echo "ERROR: falta $LD_DIR"; exit 1; }
[[ -f "$LD_DIR/ldm.info" && -f "$LD_DIR/snp.info" ]] || { echo "ERROR: falta ldm.info/snp.info"; exit 1; }
echo "  MA       : $MA_FILE  ($(printf "%'d" "$(tail -n +2 "$MA_FILE" | wc -l)") variantes)"
echo "  LD panel : $LD_DIR   ($(printf "%'d" "$(tail -n +2 "$LD_DIR/ldm.info" | wc -l)") bloques)"

# ── Panel LD AISLADO (symlinks a binarios, copia de .info) ─────────
echo ""
echo "--- Construyendo copia local aislada del panel LD ---"
LD_LOCAL="$TMPDIR/ld_local"; mkdir -p "$LD_LOCAL"
for f in "$LD_DIR"/*; do
  bn=$(basename "$f")
  case "$bn" in
    snp.info*|ldm.info) cp -f "$f" "$LD_LOCAL/" ;;   # copias EDITABLES (swap local)
    *)                  ln -sf "$f" "$LD_LOCAL/" ;;  # symlink a binarios grandes
  esac
done
echo "  Copia local: $LD_LOCAL (los .info se editan aqui; el panel compartido queda intacto)"

# ── SBayesRC MCMC ──────────────────────────────────────────────────
echo ""
echo "============================================================"
echo "  SBayesRC MCMC (${SEX}, panel aislado)"
echo "============================================================"

R_SCRIPT="$TMPDIR/run_sbrc.R"
cat << 'REOF' > "$R_SCRIPT"
suppressMessages({ library(SBayesRC); library(data.table) })

sex        <- Sys.getenv("SEX")
ma_file    <- Sys.getenv("MA_FILE")
ld_dir     <- Sys.getenv("LD_DIR_LOCAL")     # <- copia local aislada
out_prefix <- Sys.getenv("OUT_PREFIX")

cat(sprintf("\n[1] Pre-validacion del .ma (%s)...\n", sex))
ma_head <- fread(ma_file, nrows=10); print(ma_head)
required <- c("SNP","A1","A2","freq","b","se","p","N")
missing <- setdiff(required, names(ma_head))
if (length(missing) > 0) stop("Faltan columnas en .ma: ", paste(missing, collapse=", "))

ma_full <- fread(ma_file)
cat(sprintf("  Variantes: %s | N: %s | beta median |b|: %.3e\n",
    format(nrow(ma_full), big.mark=","), format(unique(ma_full$N), big.mark=","),
    median(abs(ma_full$b))))

# Chequeo de solapamiento GWAS vs LD (snp.info local, 16col)
ld_snps <- fread(file.path(ld_dir, "snp.info"), select="ID")
overlap <- sum(ma_full$SNP %in% ld_snps$ID)
cat(sprintf("\n[1.1] Solapamiento GWAS vs LD: %s (%.2f%% del GWAS)\n",
    format(overlap, big.mark=","), 100*overlap/nrow(ma_full)))
if (overlap < 100) {
  cat("  Ejemplo GWAS:", ma_full$SNP[1], "| Ejemplo LD:", ld_snps$ID[1], "\n")
  stop("Solapamiento casi nulo: formato de SNP no coincide.")
}
rm(ma_full); gc()

# ── Imputacion (swap snp.info SOLO en la copia local) ──────────────
ma_file_imp <- paste0(out_prefix, "_imp.ma")
if (file.exists(ma_file_imp)) {
  cat("[1.5] .ma imputado ya existe — reutilizando\n")
} else {
  snp_info    <- file.path(ld_dir, "snp.info")
  snp_info_16 <- file.path(ld_dir, "snp.info.16col")
  snp_info_10 <- file.path(ld_dir, "snp.info.10col")
  did_swap <- FALSE
  if (file.exists(snp_info_10)) {
    cat("[1.5] Swap LOCAL snp.info -> 10col para impute()\n")
    file.copy(snp_info, snp_info_16, overwrite=TRUE)   # backup local (16col actual)
    file.copy(snp_info_10, snp_info, overwrite=TRUE)   # 10col para impute
    did_swap <- TRUE
  } else {
    cat("[1.5] No hay snp.info.10col local; impute() usa snp.info tal cual\n")
  }
  SBayesRC::impute(mafile=ma_file, LDdir=ld_dir, output=ma_file_imp)
  if (did_swap) {
    file.copy(snp_info_16, snp_info, overwrite=TRUE)   # restaurar 16col (local)
    cat("[1.5] snp.info local restaurado a 16col\n")
  }
}
ma_file <- ma_file_imp
cat(sprintf("  Variantes tras imputacion: %s\n",
    format(length(readLines(ma_file)) - 1L, big.mark=",")))

cat("\n[2] SBayesRC::sbayesrc()...\n")
cat(sprintf("  Inicio: %s\n", Sys.time()))
result <- SBayesRC::sbayesrc(
  mafile    = ma_file,
  LDdir     = ld_dir,
  outPrefix = out_prefix,
  niter     = 3000,
  burn      = 1000,
  thresh    = 0.9,
  bTune     = FALSE,
  log2file  = FALSE
)
cat(sprintf("  Fin: %s\n", Sys.time()))

snp_file <- paste0(out_prefix, ".txt")
if (file.exists(snp_file)) {
  sr <- fread(snp_file)
  ec <- intersect(c("A1Effect","BETA"), names(sr))[1]
  if (!is.na(ec)) {
    nz <- sum(sr[[ec]] != 0, na.rm=TRUE)
    cat(sprintf("\n[3] %s: %s variantes, %s con beta!=0 (%.2f%%)\n",
        sex, format(nrow(sr), big.mark=","), format(nz, big.mark=","), 100*nz/nrow(sr)))
  }
}
par_file <- paste0(out_prefix, ".par")
if (file.exists(par_file)) { cat("\n[4] Params MCMC:\n"); cat(paste(readLines(par_file), collapse="\n"), "\n") }
cat(sprintf("\n OK SBayesRC %s completado\n", sex))
REOF

apptainer exec --bind /mnt/cephfs:/mnt/cephfs --bind "$TMPDIR:/tmp,$TMPDIR:/var/tmp" \
    --env OMP_NUM_THREADS=24 --env OPENBLAS_NUM_THREADS=24 --env MKL_NUM_THREADS=24 \
    --env SEX="$SEX" --env MA_FILE="$MA_FILE" --env LD_DIR_LOCAL="$LD_LOCAL" --env OUT_PREFIX="$OUT_PREFIX" \
    "$SBRC_SIF" Rscript "$R_SCRIPT"

# ── Exportar weights para scoring ──────────────────────────────────
echo ""
echo "============================================================"
echo "  Exportando weights (${SEX})"
echo "============================================================"
EXPORT_R="$TMPDIR/export_weights.R"
cat << 'REOF' > "$EXPORT_R"
suppressMessages(library(data.table))
out_prefix  <- Sys.getenv("OUT_PREFIX")
results_dir <- Sys.getenv("RESULTS_DIR")
sex         <- Sys.getenv("SEX")
snp_file <- paste0(out_prefix, ".txt")
if (!file.exists(snp_file)) stop("No existe ", snp_file)
df <- fread(snp_file)
id_col     <- intersect(c("Id","Name","SNP","ID"), names(df))[1]
effect_col <- intersect(c("A1Effect","BETA","Beta","beta"), names(df))[1]
freq_col   <- intersect(c("A1Frq","A1Freq","MAF","Freq"), names(df))[1]
stopifnot(!is.na(id_col), !is.na(effect_col))
out_w <- file.path(results_dir, sprintf("sbayesrc_%s_weights_for_scoring.tsv", sex))
out_dt <- data.table(SNP=df[[id_col]], A1=df$A1, A2=df$A2, BETA=df[[effect_col]])
if (!is.na(freq_col)) out_dt$FREQ <- df[[freq_col]]
fwrite(out_dt, out_w, sep="\t")
cat(sprintf("OK weights: %s (%s variantes, %s con BETA!=0)\n", out_w,
    format(nrow(out_dt), big.mark=","), format(sum(out_dt$BETA != 0), big.mark=",")))
REOF

apptainer exec --bind /mnt/cephfs:/mnt/cephfs \
    --env OUT_PREFIX="$OUT_PREFIX" --env RESULTS_DIR="$RESULTS_DIR" --env SEX="$SEX" \
    "$SIF" /opt/miniforge3/envs/r-bio/bin/Rscript "$EXPORT_R"

echo ""
echo "============================================================"
echo "  PIPELINE D sex ${SEX} COMPLETADO — $(date)"
echo "  Outputs en $RESULTS_DIR/:"
echo "    sbayesrc_${SEX}.txt / .par / weights_for_scoring.tsv"
echo "  Panel LD compartido en ldm_79k: INTACTO (se uso copia local)"
echo "============================================================"
