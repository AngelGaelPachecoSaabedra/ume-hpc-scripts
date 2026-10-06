#!/bin/bash
#SBATCH --job-name=sbrc_ma_sex
#SBATCH --account=researchers
#SBATCH --qos=vip
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=120G
#SBATCH --output=/mnt/cephfs/orgs/home/angel.pacheco/logs/sbrc_ma_sex_%j.out
#SBATCH --error=/mnt/cephfs/orgs/home/angel.pacheco/logs/sbrc_ma_sex_%j.err
#SBATCH --mail-user=angelpachecosaavedra@gmail.com
#SBATCH --mail-type=END,FAIL

set -euo pipefail

# ================================================================
# SBayesRC — generar .ma POR SEXO (Pipeline D sex-specific)
# Lanzar DOS veces:
#   sbatch --export=SEX=females sbrc_mcps_generate_ma_sexspecific.sbatch
#   sbatch --export=SEX=males   sbrc_mcps_generate_ma_sexspecific.sbatch
#
# Diferencia vs version meta: SIN IVW. Cada sexo conserva sus betas
# crudos REGENIE y su propio n_eff. Output: MCPS_<sex>_matched.ma
# ================================================================

SEX=${SEX:-females}

SIF="/mnt/cephfs/biocontainers/images/jupyter-biotools-1.4.sif"
GWAS_DIR="/mnt/cephfs/orgs/home/angel.pacheco/gwas_miami"
PLINK_DIR="/mnt/cephfs/hot_nvme/mcps/imputed-topmed/plink_files/maximally_unrelated"
OUT_BASE="/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes/sbayesrc"
MA_OUTPUT="${OUT_BASE}/MCPS_${SEX}_matched.ma"
mkdir -p "$OUT_BASE"

TMPDIR="/mnt/cephfs/scratch/$USER/sbrc_ma_sex_tmp/${SLURM_JOB_ID}_${SEX}"
mkdir -p "$TMPDIR"; chmod 700 "$TMPDIR"
export TMPDIR TMP="$TMPDIR" TEMP="$TMPDIR"
trap 'rm -rf "$TMPDIR"' EXIT

echo "============================================================"
echo "  Generando .ma SBayesRC por sexo: ${SEX} — Job $SLURM_JOB_ID"
echo "  Nodo: $SLURMD_NODENAME   $(date)"
echo "============================================================"

# ── 1. IDs del .bim ────────────────────────────────────────────────
echo ""
echo "--- 1. Set de IDs del .bim ---"
BIM_IDS="$TMPDIR/all_bim_ids.txt"; > "$BIM_IDS"
for CHR in $(seq 1 22) X; do
  BIM="${PLINK_DIR}/mcps-freeze150k_qcd_chr${CHR}_ivs.bim"
  [[ -f "$BIM" ]] && awk '{print $2}' "$BIM" >> "$BIM_IDS"
done
echo "  IDs en .bim: $(printf "%'d" "$(wc -l < "$BIM_IDS")")"

# ── 2. Generar .ma (por sexo, sin IVW) ─────────────────────────────
R_SCRIPT="$TMPDIR/generate_ma.R"
cat << 'REOF' > "$R_SCRIPT"
suppressMessages({ library(data.table) })
setDTthreads(8L)

sex      <- Sys.getenv("SEX")
gwas_dir <- Sys.getenv("GWAS_DIR")
ma_out   <- Sys.getenv("MA_OUTPUT")
bim_file <- Sys.getenv("BIM_IDS")

# GWAS y n_eff de este sexo
if (sex == "females") {
  gwas_file <- "combined-gwas-results-all-ctrl-females.txt.gz"
  neff <- 4*17300*54196/(17300+54196)   # ~52,457
} else if (sex == "males") {
  gwas_file <- "combined-gwas-results-all-ctrl-males.txt.gz"
  neff <- 4*8522*26733/(8522+26733)     # ~25,849
} else stop("SEX debe ser 'females' o 'males'")

cat(sprintf("\n[1] Cargando GWAS %s (n_eff=%.0f)...\n", sex, neff))
g <- fread(file.path(gwas_dir, gwas_file),
           select=c("CHROM","GENPOS","ID","ALLELE0","ALLELE1","A1FREQ","INFO","N","BETA","SE","P"),
           showProgress=FALSE)
cat(sprintf("  Variantes crudas: %s\n", format(nrow(g), big.mark=",")))

cat("\n[2] Construyendo tabla .ma (SIN IVW, betas crudos del sexo)...\n")
ma_df <- data.table(
  CHR  = g$CHROM,
  POS  = g$GENPOS,
  A1   = g$ALLELE1,    # effect allele
  A2   = g$ALLELE0,    # non-effect
  freq = g$A1FREQ,
  info = g$INFO,
  b    = g$BETA,
  se   = g$SE,
  N    = round(neff)
)
ma_df[, p := 2*pnorm(-abs(b/se))]
rm(g); gc()

cat("\n[3] QC...\n")
n0 <- nrow(ma_df)
ma_df <- ma_df[info>=0.8 & freq>=0.01 & freq<=0.99 & !is.na(b) & !is.na(se) & se>0]
cat(sprintf("  %s -> %s variantes (%.1f%% retenidas)\n",
    format(n0, big.mark=","), format(nrow(ma_df), big.mark=","), 100*nrow(ma_df)/n0))

cat("\n[4] Matching contra .bim (data.table)...\n")
bim_dt <- fread(bim_file, header=FALSE, col.names="bim_id"); setkey(bim_dt, bim_id)
cat(sprintf("  %s IDs en .bim\n", format(nrow(bim_dt), big.mark=",")))

ma_df[, chr_str := ifelse(CHR==23 | as.character(CHR) %in% c("X","23"),
                          "chrX", paste0("chr", CHR))]
ma_df[, id_fwd := paste0(chr_str, ":", POS, ":", A1, ":", A2)]
ma_df[, id_rev := paste0(chr_str, ":", POS, ":", A2, ":", A1)]

ma_df[, fwd_match := id_fwd %in% bim_dt$bim_id]
ma_df[, rev_match := FALSE]
ma_df[fwd_match==FALSE, rev_match := id_rev %in% bim_dt$bim_id]
ma_df[, id_match := NA_character_]
ma_df[fwd_match==TRUE, id_match := id_fwd]
ma_df[fwd_match==FALSE & rev_match==TRUE, id_match := id_rev]

n_fwd <- sum(ma_df$fwd_match, na.rm=TRUE)
n_rev <- sum(ma_df$rev_match & !ma_df$fwd_match, na.rm=TRUE)
cat(sprintf("  Forward: %s | Reverse: %s | Total matched: %s (%.1f%%)\n",
    format(n_fwd, big.mark=","), format(n_rev, big.mark=","),
    format(n_fwd+n_rev, big.mark=","), 100*(n_fwd+n_rev)/nrow(ma_df)))

ma_df <- ma_df[!is.na(id_match)]; ma_df[, SNP := id_match]

cat("\n[5] Escribiendo .ma final...\n")
ma_out_df <- ma_df[, .(SNP, A1, A2, freq, b, se, p, N)]
ma_out_df <- ma_out_df[is.finite(b) & is.finite(se) & se>0]
fwrite(ma_out_df, ma_out, sep=" ", quote=FALSE)
cat(sprintf("  OK: %s variantes -> %s\n", format(nrow(ma_out_df), big.mark=","), ma_out))
cat(sprintf("  N const: %s | beta median |b|: %.4e\n",
    format(unique(ma_out_df$N), big.mark=","), median(abs(ma_out_df$b))))
REOF

apptainer exec --bind /mnt/cephfs:/mnt/cephfs --bind "$TMPDIR:/tmp,$TMPDIR:/var/tmp" \
    --env OMP_NUM_THREADS=8 --env OPENBLAS_NUM_THREADS=8 --env MKL_NUM_THREADS=8 \
    --env SEX="$SEX" --env GWAS_DIR="$GWAS_DIR" --env MA_OUTPUT="$MA_OUTPUT" --env BIM_IDS="$BIM_IDS" \
    "$SIF" /opt/miniforge3/envs/r-bio/bin/Rscript "$R_SCRIPT"

echo ""
echo "--- Validacion ---"
head -1 "$MA_OUTPUT"
head -4 "$MA_OUTPUT" | tail -3 | column -t
echo "  Variantes en .ma (${SEX}): $(printf "%'d" "$(tail -n +2 "$MA_OUTPUT" | wc -l)")"
echo ""
echo "============================================================"
echo "  .ma ${SEX} COMPLETADO — siguiente: sbrc_mcps_run_sexspecific.sbatch"
echo "============================================================"
