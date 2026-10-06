#!/bin/bash
#SBATCH --job-name=sbrc_mcps_ma
#SBATCH --account=researchers
#SBATCH --qos=vip
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=120G
#SBATCH --time=02:00:00
#SBATCH --output=/mnt/cephfs/orgs/home/angel.pacheco/logs/sbrc_mcps_ma_%j.out
#SBATCH --error=/mnt/cephfs/orgs/home/angel.pacheco/logs/sbrc_mcps_ma_%j.err
#SBATCH --mail-user=angelpachecosaavedra@gmail.com
#SBATCH --mail-type=END,FAIL

set -euo pipefail

# ================================================================
# v2: usa data.table merge en lugar de hash environment
# (el hash de 43M elementos con assign() crashea R con segfault)
# ================================================================

SIF="/mnt/cephfs/biocontainers/images/jupyter-biotools-1.4.sif"
GWAS_DIR="/mnt/cephfs/orgs/home/angel.pacheco/gwas_miami"
PLINK_DIR="/mnt/cephfs/hot_nvme/mcps/imputed-topmed/plink_files/maximally_unrelated"
OUT_BASE="/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes/sbayesrc"
MA_OUTPUT="${OUT_BASE}/MCPS_meta_matched.ma"

mkdir -p "$OUT_BASE"

TMPDIR="/mnt/cephfs/scratch/$USER/sbrc_mcps_ma_tmp/$SLURM_JOB_ID"
mkdir -p "$TMPDIR"; chmod 700 "$TMPDIR"
export TMPDIR TMP="$TMPDIR" TEMP="$TMPDIR"

echo "============================================================"
echo "  Generando .ma desde GWAS MCPS para SBayesRC (v2)"
echo "  Nodo: $SLURMD_NODENAME"
echo "  $(date)"
echo "============================================================"

# ================================================================
# Paso 1: Extraer todos los SNP IDs del .bim
# ================================================================
echo ""
echo "--- 1. Construyendo set de IDs del .bim ---"

BIM_IDS="$TMPDIR/all_bim_ids.txt"
> "$BIM_IDS"

for CHR in $(seq 1 22) X; do
  BIM="${PLINK_DIR}/mcps-freeze150k_qcd_chr${CHR}_ivs.bim"
  if [[ -f "$BIM" ]]; then
    awk '{print $2}' "$BIM" >> "$BIM_IDS"
  fi
done

N_BIM=$(wc -l < "$BIM_IDS")
echo "  IDs en .bim totales: $(printf "%'d" $N_BIM)"

# ================================================================
# Paso 2: Generar .ma con data.table merge
# ================================================================
R_SCRIPT="$TMPDIR/generate_ma.R"

cat << REOF > "$R_SCRIPT"
suppressMessages({ library(data.table) })
setDTthreads(8L)

gwas_dir <- "$GWAS_DIR"
ma_out   <- "$MA_OUTPUT"
bim_file <- "$BIM_IDS"

cat("\n[1] Cargando GWAS sex-stratified...\n")
gf <- fread(file.path(gwas_dir, "combined-gwas-results-all-ctrl-females.txt.gz"),
            select=c("CHROM","GENPOS","ID","ALLELE0","ALLELE1","A1FREQ","INFO","N","BETA","SE","P"),
            showProgress=FALSE)
gm <- fread(file.path(gwas_dir, "combined-gwas-results-all-ctrl-males.txt.gz"),
            select=c("CHROM","GENPOS","ID","ALLELE0","ALLELE1","A1FREQ","INFO","N","BETA","SE","P"),
            showProgress=FALSE)
cat(sprintf("  Females: %s, Males: %s\n",
    format(nrow(gf), big.mark=","), format(nrow(gm), big.mark=",")))

cat("\n[2] Meta-análisis IVW...\n")
gwas <- merge(gf, gm, by="ID", suffixes=c(".f", ".m"))
rm(gf, gm); gc()
cat(sprintf("  Tras merge: %s variantes\n", format(nrow(gwas), big.mark=",")))

w_f <- 1 / gwas\$SE.f^2
w_m <- 1 / gwas\$SE.m^2

# n_eff con corrección log-odds → escala lineal
neff_f <- 4*17300*54196/(17300+54196)
neff_m <- 4*8522*26733/(8522+26733)
#K_prev <- (17300+8522)/(17300+54196+8522+26733)
#neff_lineal <- (neff_f + neff_m) * K_prev * (1 - K_prev)
neff_lineal <- neff_f + neff_m                               # 78,304

#cat(sprintf("  K (prevalencia)    : %.4f\n", K_prev))
cat(sprintf("  N_eff combinado    : %.0f\n", neff_f + neff_m))
cat(sprintf("  N_eff escala lineal: %.0f\n", neff_lineal))

# Construir tabla con formato .ma
ma_df <- data.table(
  CHR  = gwas\$CHROM.f,
  POS  = gwas\$GENPOS.f,
  A1   = gwas\$ALLELE1.f,   # Effect allele
  A2   = gwas\$ALLELE0.f,   # Non-effect
  freq = (gwas\$A1FREQ.f * gwas\$N.f + gwas\$A1FREQ.m * gwas\$N.m) / (gwas\$N.f + gwas\$N.m),
  info = pmin(gwas\$INFO.f, gwas\$INFO.m),
  b    = (gwas\$BETA.f * w_f + gwas\$BETA.m * w_m) / (w_f + w_m),
  se   = sqrt(1 / (w_f + w_m)),
  N    = round(neff_lineal)
)
ma_df[, p := 2 * pnorm(-abs(b / se))]
rm(gwas); gc()

cat("\n[3] QC...\n")
n0 <- nrow(ma_df)
ma_df <- ma_df[info >= 0.8 &
               freq >= 0.01 & freq <= 0.99 &
               !is.na(b) & !is.na(se) & se > 0]
cat(sprintf("  %s → %s variantes (%.1f%% retenidas)\n",
    format(n0, big.mark=","), format(nrow(ma_df), big.mark=","),
    100*nrow(ma_df)/n0))

# ================================================================
# [4] Matching .ma vs .bim usando data.table merge (NO hash env)
# ================================================================
cat("\n[4] Construyendo SNP IDs y matching contra .bim...\n")

cat("  Leyendo .bim IDs como data.table...\n")
bim_dt <- fread(bim_file, header=FALSE, col.names="bim_id")
cat(sprintf("  %s IDs cargados\n", format(nrow(bim_dt), big.mark=",")))
setkey(bim_dt, bim_id)

# Construir IDs candidatos
ma_df[, chr_str := ifelse(CHR == 23 | as.character(CHR) %in% c("X","23"),
                          "chrX", paste0("chr", CHR))]
ma_df[, id_fwd := paste0(chr_str, ":", POS, ":", A1, ":", A2)]
ma_df[, id_rev := paste0(chr_str, ":", POS, ":", A2, ":", A1)]

cat("  Matching forward via data.table merge...\n")
ma_df[, fwd_match := id_fwd %in% bim_dt\$bim_id]

cat("  Matching reverse para los no-matched...\n")
ma_df[, rev_match := FALSE]
ma_df[fwd_match == FALSE, rev_match := id_rev %in% bim_dt\$bim_id]

# Asignar id_match
ma_df[, id_match := NA_character_]
ma_df[fwd_match == TRUE, id_match := id_fwd]
ma_df[fwd_match == FALSE & rev_match == TRUE, id_match := id_rev]

n_fwd <- sum(ma_df\$fwd_match, na.rm=TRUE)
n_rev <- sum(ma_df\$rev_match & !ma_df\$fwd_match, na.rm=TRUE)
n_unmatched <- sum(is.na(ma_df\$id_match))

cat(sprintf("\n  Match forward (A1:A2): %s\n", format(n_fwd, big.mark=",")))
cat(sprintf("  Match reverse (A2:A1): %s\n", format(n_rev, big.mark=",")))
cat(sprintf("  Total matched        : %s (%.1f%%)\n",
    format(n_fwd + n_rev, big.mark=","), 100*(n_fwd+n_rev)/nrow(ma_df)))
cat(sprintf("  Sin match            : %s (%.1f%%)\n",
    format(n_unmatched, big.mark=","), 100*n_unmatched/nrow(ma_df)))

# Filtrar a matched
ma_df <- ma_df[!is.na(id_match)]
ma_df[, SNP := id_match]

cat("\n[5] Escribiendo .ma final...\n")
ma_out_df <- ma_df[, .(SNP, A1, A2, freq, b, se, p, N)]

# Quitar NA/Inf en b o se
n_bad <- sum(!is.finite(ma_out_df\$b) | !is.finite(ma_out_df\$se) | ma_out_df\$se <= 0)
if (n_bad > 0) {
  cat(sprintf("  %d filas con b/se inválidos — removiendo\n", n_bad))
  ma_out_df <- ma_out_df[is.finite(b) & is.finite(se) & se > 0]
}

fwrite(ma_out_df, ma_out, sep=" ", quote=FALSE)
cat(sprintf("  ✓ Escritas %s variantes\n", format(nrow(ma_out_df), big.mark=",")))

cat("\n[6] Validación final:\n")
cat(sprintf("  Archivo : %s\n", ma_out))
cat(sprintf("  Tamaño  : %.1f MB\n", file.info(ma_out)\$size / 1024^2))
cat(sprintf("  N const : %s\n", format(unique(ma_out_df\$N), big.mark=",")))
cat(sprintf("  beta median |b| : %.4e\n", median(abs(ma_out_df\$b))))
cat(sprintf("  beta range      : [%.3e, %.3e]\n",
    min(ma_out_df\$b), max(ma_out_df\$b)))

# Distribución por chr
cat("\n  Variantes por cromosoma:\n")
chr_summary <- ma_out_df[, .(n=.N), by=.(chr=sub(":.*", "", SNP))]
chr_summary <- chr_summary[order(chr)]
for (i in 1:nrow(chr_summary)) {
  cat(sprintf("    %s: %s\n", chr_summary\$chr[i],
      format(chr_summary\$n[i], big.mark=",")))
}
REOF

apptainer exec --bind /mnt/cephfs:/mnt/cephfs \
    --bind "$TMPDIR:/tmp,$TMPDIR:/var/tmp" \
    --env OMP_NUM_THREADS=8 \
    --env OPENBLAS_NUM_THREADS=8 \
    --env MKL_NUM_THREADS=8 \
    "$SIF" /opt/miniforge3/envs/r-bio/bin/Rscript "$R_SCRIPT"

echo ""
echo "--- Validación rápida ---"
echo "  Header:"
head -1 "$MA_OUTPUT"
echo ""
echo "  Primeras 3 variantes:"
head -4 "$MA_OUTPUT" | tail -3 | column -t
echo ""
N_FINAL=$(tail -n +2 "$MA_OUTPUT" | wc -l)
echo "  Variantes en .ma final: $(printf "%'d" $N_FINAL)"

echo ""
echo "  Quick check chr20:"
awk '{print $2}' "${PLINK_DIR}/mcps-freeze150k_qcd_chr20_ivs.bim" | sort > "$TMPDIR/bim20.txt"
awk 'NR>1 && $1 ~ /^chr20:/' "$MA_OUTPUT" | awk '{print $1}' | sort > "$TMPDIR/ma20.txt"
N20=$(comm -12 "$TMPDIR/bim20.txt" "$TMPDIR/ma20.txt" | wc -l)
echo "    En común con .bim chr20: $(printf "%'d" $N20) (de $(wc -l < "$TMPDIR/ma20.txt") en .ma)"

rm -rf "$TMPDIR"

echo ""
echo "============================================================"
echo "   COMPLETADO"
echo ""
echo "  Output: $MA_OUTPUT"
echo ""
echo "  Siguiente: sbatch sbrc_mcps_run.sbatch"
echo "============================================================"
