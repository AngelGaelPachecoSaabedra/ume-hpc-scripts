#!/bin/bash
#SBATCH --job-name=ldp2_gw_prep
#SBATCH --account=researchers
#SBATCH --qos=vip
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=150G
#SBATCH --output=/mnt/cephfs/orgs/home/angel.pacheco/logs/ldp2_gw_prep_%A_%a.out
#SBATCH --error=/mnt/cephfs/orgs/home/angel.pacheco/logs/ldp2_gw_prep_%A_%a.err
#SBATCH --mail-user=angelpachecosaavedra@gmail.com
#SBATCH --mail-type=FAIL

set -euo pipefail

# ================================================================
# PASO 1/2: Preparar bloques LD por cromosoma (paralelo)
# Cada worker: .bed → QC → snp_cor → snp_ldsplit → guarda bloques
# Lanzar: sbatch --array=1-22%6 ldpred2_gwide_step1.sbatch
# ================================================================

CHR=$SLURM_ARRAY_TASK_ID

# ──────────────────────────────────────────────────────────────────
# IMPORTANTE: usar jupyter-biotools-1.5.sif (bigsnpr 1.10.8 en r-bio).
# bigsnpr_nativo.sif tiene 1.12.21 que está roto.
# El R con bigsnpr 1.10.8 vive en el conda env r-bio:
#   /opt/miniforge3/envs/r-bio/bin/Rscript
# Verifica con:
#   apptainer exec $SIF /opt/miniforge3/envs/r-bio/bin/Rscript \
#     -e 'packageVersion("bigsnpr")'
# ──────────────────────────────────────────────────────────────────
SIF="/mnt/cephfs/biocontainers/images/jupyter-biotools-1.5.sif"

GWAS_DIR="/mnt/cephfs/orgs/home/angel.pacheco/gwas_miami"
PLINK_DIR="/mnt/cephfs/hot_nvme/mcps/imputed-topmed/plink_files/maximally_unrelated"
BLOCKS_DIR="/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes/ldpred2_gwide/blocks"

mkdir -p "$BLOCKS_DIR"

TMPDIR="/mnt/cephfs/scratch/$USER/ldp2_gw_prep_tmp/${SLURM_ARRAY_JOB_ID}_${CHR}"
mkdir -p "$TMPDIR"; chmod 700 "$TMPDIR"
export TMPDIR TMP="$TMPDIR" TEMP="$TMPDIR"

echo "============================================================"
echo "  LDpred2 gwide STEP 1 — chr${CHR} (QC + snp_cor + ldsplit)"
echo "  Job: ${SLURM_ARRAY_JOB_ID}_${CHR}   Nodo: ${SLURMD_NODENAME}"
echo "  $(date)"
echo "  Script version: v2 (2026-05-07) — FIX bigsnpr 1.10.8"
echo "============================================================"

R_SCRIPT="$TMPDIR/prep_chr${CHR}.R"

cat << 'REOF' > "$R_SCRIPT"
# ── Threading defensivo ────────────────────────────────────────────
Sys.setenv(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
           MKL_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1",
           BLIS_NUM_THREADS="1", NUMEXPR_NUM_THREADS="1")
if (requireNamespace("RhpcBLASctl", quietly=TRUE)) {
  RhpcBLASctl::blas_set_num_threads(1L)
  RhpcBLASctl::omp_set_num_threads(1L)
}
suppressMessages({
  library(bigsnpr); library(bigsparser); library(data.table); library(Matrix)
})
options(bigstatsr.check.parallel.blas=FALSE)
options(bigstatsr.ncores.max=24L)

# ══════════════════════════════════════════════════════════════════
# FIX #1: ABORTAR INMEDIATAMENTE si la versión de bigsnpr es la rota.
# 1.12.21 tiene bugs documentados en corMat (segfault chr6) y
# snp_ldsplit (chr4: n_block no-numeric). Pipeline funciona con 1.10.8.
# ══════════════════════════════════════════════════════════════════
bv <- packageVersion("bigsnpr")
cat(sprintf("  bigsnpr version: %s\n", as.character(bv)))
if (as.character(bv) == "1.12.21") {
  stop("bigsnpr 1.12.21 está roto para snp_cor/snp_ldsplit. ",
       "Usa el contenedor con 1.10.8.")
}
if (bv >= "1.11.0" && bv < "1.10.8") {
  warning("Versión no probada de bigsnpr: ", bv,
          ". Procediendo con cautela.")
}

chr        <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
ncores     <- 12L  # Coincide con --cpus-per-task. nb_cores() retornaría 11 (deja 1 para OS), pero snp_cor escala bien con todos.
plink_dir  <- Sys.getenv("PLINK_DIR")
gwas_dir   <- Sys.getenv("GWAS_DIR")
blocks_dir <- Sys.getenv("BLOCKS_DIR")
tmp_dir    <- Sys.getenv("TMPDIR")

chr_label <- as.character(chr)
out_file  <- file.path(blocks_dir, sprintf("chr%s_blocks.rds", chr_label))

if (file.exists(out_file)) {
  cat(sprintf("Ya existe %s — saltando.\n", basename(out_file)))
  quit(save="no", status=0)
}

cat(sprintf("\n==== chr%s (ncores=%d) ====\n", chr_label, ncores))

# ══════════════════════════════════════════════════════════════════
# FIX #5: Wrapper que captura cualquier fallo y guarda placeholder
# para que Step 2 pueda continuar con los chr exitosos.
# ══════════════════════════════════════════════════════════════════
save_failed <- function(msg, stage) {
  cat(sprintf("\n  *** FALLO en %s: %s ***\n", stage, msg))
  saveRDS(list(chr=chr, status="failed", stage=stage, msg=msg,
               n_snp=0L, n_blocks=0L, blocks=list(),
               info_snp=NULL, sd_val=NULL), out_file)
  quit(save="no", status=0)  # status=0 para no marcar el array job como FAIL
}

# ── [0] HapMap3 ────────────────────────────────────────────────────
map_hm3 <- tryCatch(
  readRDS(file.path(
    "/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes/ldpred2",
    "map_hm3_ldpred2.rds")),
  error=function(e) save_failed(conditionMessage(e), "load_hm3")
)
if ("pos_hg38" %in% names(map_hm3)) map_hm3$pos_use <- map_hm3$pos_hg38

# ── [1] GWAS meta chr ─────────────────────────────────────────────
cat(sprintf("\n[1] GWAS meta chr%s...\n", chr_label))
gf <- fread(file.path(gwas_dir, "combined-gwas-results-all-ctrl-females.txt.gz"),
  select=c("CHROM","GENPOS","ID","ALLELE0","ALLELE1","A1FREQ","INFO","N","BETA","SE","P"),
  showProgress=FALSE)
gm <- fread(file.path(gwas_dir, "combined-gwas-results-all-ctrl-males.txt.gz"),
  select=c("CHROM","GENPOS","ID","ALLELE0","ALLELE1","A1FREQ","INFO","N","BETA","SE","P"),
  showProgress=FALSE)
gf <- gf[gf$CHROM == chr, ]; gm <- gm[gm$CHROM == chr, ]
gwas <- merge(gf, gm, by="ID", suffixes=c(".f",".m"))
rm(gf, gm); gc()

w_f <- 1/gwas$SE.f^2; w_m <- 1/gwas$SE.m^2
neff_f <- 4*17300*54196/(17300+54196)
neff_m <- 4*8522*26733/(8522+26733)

# ══════════════════════════════════════════════════════════════════
# Approach Pipeline C (NO Pipeline B v8):
# Aquí n_eff inicial = 78,304 (SIN corrección K(1-K) global).
# La calibración per-SNP en sección [2c] recupera el n_eff "efectivo"
# por SNP usando (2/(sd*se))² - (beta/se)², que es el método de Privé.
# Esto es MUTUAMENTE EXCLUYENTE con × K(1-K) global; no mezclar.
# Ratio sd_ss_orig/sd_val esperado ≈ 0.86 antes de calibrar — OK.
# ══════════════════════════════════════════════════════════════════
neff_ldpred2 <- neff_f + neff_m   # ≈ 78,304 (semilla; se recalibra per-SNP)

sumstats <- data.frame(
  chr=gwas$CHROM.f, pos=gwas$GENPOS.f,
  a1=gwas$ALLELE1.f, a0=gwas$ALLELE0.f,
  freq=(gwas$A1FREQ.f*gwas$N.f + gwas$A1FREQ.m*gwas$N.m)/(gwas$N.f+gwas$N.m),
  info=pmin(gwas$INFO.f, gwas$INFO.m),
  beta=(gwas$BETA.f*w_f + gwas$BETA.m*w_m)/(w_f+w_m),
  beta_se=sqrt(1/(w_f+w_m)), n_eff=neff_ldpred2, stringsAsFactors=FALSE)
sumstats$p <- 2*pnorm(-abs(sumstats$beta/sumstats$beta_se))
sumstats <- sumstats[sumstats$info>=0.8 & sumstats$freq>=0.01 & sumstats$freq<=0.99 &
  !is.na(sumstats$beta) & !is.na(sumstats$beta_se) & sumstats$beta_se>0, ]
sumstats$chr <- as.integer(sumstats$chr)

key_hm3  <- paste(map_hm3$chr, map_hm3$pos_use, sep=":")
key_gwas <- paste(sumstats$chr, sumstats$pos, sep=":")
sumstats <- sumstats[key_gwas %in% key_hm3, ]
rm(gwas); gc()
cat(sprintf("  GWAS HM3: %s\n", format(nrow(sumstats), big.mark=",")))

# ── [2] Lectura .bed + Matching + QC ──────────────────────────────
cat("\n[2] Leyendo .bed, matching, QC...\n")
prefix <- file.path(plink_dir, sprintf("mcps-freeze150k_qcd_chr%s_ivs", chr_label))
bed_file <- paste0(prefix, ".bed")
if (!file.exists(bed_file)) save_failed(paste("No existe:", bed_file), "read_bed")

bk_file <- file.path(tmp_dir, sprintf("chr%s", chr_label))
snp_readBed(bed_file, backingfile=bk_file)
obj <- snp_attach(paste0(bk_file, ".rds"))
G   <- obj$genotypes

# ══════════════════════════════════════════════════════════════════
# FIX #2: Verificar y manejar NAs antes de snp_cor.
# snp_cor con NAs en G puede dar correlaciones espurias o (en algunas
# versiones) crashear. Imputamos con la media (modo 'mean0') que es
# lo más conservador para LD.
# ══════════════════════════════════════════════════════════════════
n_na_total <- big_counts(G, ind.col=seq_len(ncol(G)))[4, ] |> sum()
if (n_na_total > 0) {
  pct_na <- 100 * n_na_total / (as.numeric(nrow(G)) * ncol(G))
  cat(sprintf("  NAs en G: %s (%.4f%%) — imputando con snp_fastImputeSimple\n",
      format(n_na_total, big.mark=","), pct_na))
  G <- snp_fastImputeSimple(G, method="mean0", ncores=ncores)
} else {
  cat("  NAs en G: 0\n")
}

map_df <- data.frame(chr=obj$map$chromosome, pos=obj$map$physical.pos,
  a1=obj$map$allele1, a0=obj$map$allele2, stringsAsFactors=FALSE)
info_snp <- snp_match(sumstats, map_df, strand_flip=FALSE,
  join_by_pos=TRUE, match.min.prop=0.01)
ind_col <- info_snp$`_NUM_ID_`

# 2a. sd > 0
sd_vals <- big_scale()(G, ind.col=ind_col)$scale
keep_var <- which(sd_vals > 1e-6)
if (length(keep_var) < nrow(info_snp)) {
  info_snp <- info_snp[keep_var, ]
  ind_col  <- info_snp$`_NUM_ID_`
  sd_vals  <- sd_vals[keep_var]
}

# 2b. Filtro QC con sd_ss original
df_tmp <- info_snp[, c("beta","beta_se","n_eff")]
sd_ss_orig <- with(df_tmp, 2 / sqrt(n_eff * beta_se^2 + beta^2))
ratio_orig <- median(sd_ss_orig / sd_vals)
cat(sprintf("  sd_ss_orig/sd_val ratio: %.3f\n", ratio_orig))

keep_qc <- which(sd_ss_orig > 0.5 * sd_vals &
                 sd_ss_orig < sd_vals + 0.1 &
                 sd_ss_orig > 0.1 &
                 sd_vals > 0.05)
n_removed_qc <- length(sd_ss_orig) - length(keep_qc)
if (n_removed_qc > 0) {
  cat(sprintf("  Filtro QC: removidas %d (%.1f%%)\n",
      n_removed_qc, 100*n_removed_qc/length(sd_ss_orig)))
  info_snp <- info_snp[keep_qc, ]
  ind_col  <- info_snp$`_NUM_ID_`
  sd_vals  <- sd_vals[keep_qc]
  df_tmp   <- df_tmp[keep_qc, ]
} else {
  cat("  Filtro QC: 0 removidas\n")
}

# 2c. Calibrar n_eff per-SNP (Privé)
df_tmp$n_eff <- (2 / (sd_vals * df_tmp$beta_se))^2 - (df_tmp$beta / df_tmp$beta_se)^2
info_snp$n_eff <- pmax(df_tmp$n_eff, 1000)
rm(df_tmp, sd_ss_orig)

n_snp <- nrow(info_snp)
cat(sprintf("  Final: %s variantes, n_eff median=%.0f\n",
    format(n_snp, big.mark=","), median(info_snp$n_eff)))

if (n_snp < 50) {
  saveRDS(list(chr=chr, status="too_few", n_snp=n_snp,
               blocks=list(), info_snp=NULL, sd_val=NULL), out_file)
  quit(save="no", status=0)
}

# ── [3] snp_cor ────────────────────────────────────────────────────
cat("\n[3] snp_cor (size=3000 kb, thr_r2=0.005)...\n")
t0 <- Sys.time()
corr0 <- tryCatch(
  snp_cor(Gna=G, ind.col=ind_col, infos.pos=info_snp$pos/1000,
    size=3000, ncores=ncores, thr_r2=0.005),
  error=function(e) save_failed(conditionMessage(e), "snp_cor")
)
t_cor <- as.numeric(difftime(Sys.time(), t0, units="mins"))
cat(sprintf("  snp_cor: %.1f min, dim=%dx%d, nnz=%s\n",
    t_cor, nrow(corr0), ncol(corr0), format(length(corr0@x), big.mark=",")))

# ══════════════════════════════════════════════════════════════════
# FIX #4: Verificar integridad de corr0 antes de snp_ldsplit.
# Detectar columnas vacías (que romperían bdiag/as_SFBM en Step 2).
# ══════════════════════════════════════════════════════════════════
col_nnz <- diff(corr0@p)
n_empty_cols <- sum(col_nnz == 0)
if (n_empty_cols > 0) {
  cat(sprintf("  ADVERTENCIA: %d columnas vacías en corr0 — pueden romper SFBM\n",
              n_empty_cols))
  # No abortamos: el wrapper safe_as_SFBM en Step 2 las maneja.
}

# ── [4] snp_ldsplit → bloques independientes ───────────────────────
cat("\n[4] snp_ldsplit...\n")
m_chr <- ncol(corr0)

# length.out=20 (Privé vignette) en lugar de 8 — más opciones de partición
max_size_seq <- round(seq_log(m_chr/30, m_chr/5, length.out=20))
max_size_seq <- unique(max_size_seq[max_size_seq >= 100 & max_size_seq <= m_chr])
if (length(max_size_seq) == 0) max_size_seq <- min(m_chr, 5000)

cat(sprintf("  m_chr=%d, max_size candidates (n=%d): %s...%s\n",
            m_chr, length(max_size_seq),
            paste(head(max_size_seq, 3), collapse=","),
            paste(tail(max_size_seq, 2), collapse=",")))

# ══════════════════════════════════════════════════════════════════
# RETRY ESCALADO: MCPS es cohorte admixed con LD extenso. thr_r2=0.02
# (default Privé EUR) suele ser demasiado estricto y retorna res_split
# vacío. Probamos thr_r2 escalando; en EUR convergería en 0.02, en
# admixed típicamente entre 0.05–0.1.
# ══════════════════════════════════════════════════════════════════
thr_r2_seq <- c(0.02, 0.05, 0.1, 0.2, 0.3)
res_split  <- NULL
thr_used   <- NA_real_

is_valid_split <- function(rs) {
  if (is.null(rs) || nrow(rs) == 0) return(FALSE)
  if (!is.numeric(rs$cost) || !is.numeric(rs$n_block)) return(FALSE)
  any(!is.na(rs$cost) & !is.na(rs$n_block) & rs$n_block > 0)
}

for (try_thr in thr_r2_seq) {
  t0_split <- Sys.time()
  res_try <- tryCatch(
    snp_ldsplit(corr0, thr_r2=try_thr, min_size=50,
                max_size=max_size_seq, max_K=500),
    error=function(e) {
      cat(sprintf("    thr_r2=%.3f ERROR: %s\n", try_thr, conditionMessage(e)))
      NULL
    }
  )
  dt <- as.numeric(difftime(Sys.time(), t0_split, units="secs"))

  if (is_valid_split(res_try)) {
    n_valid_rows <- sum(!is.na(res_try$cost) & res_try$n_block > 0)
    cat(sprintf("  ✓ thr_r2=%.3f → %d candidatos válidos (%.1fs)\n",
                try_thr, n_valid_rows, dt))
    res_split <- res_try
    thr_used  <- try_thr
    break
  } else {
    cat(sprintf("    thr_r2=%.3f vacío/malformado (%.1fs) — subiendo umbral\n",
                try_thr, dt))
  }
}

# ══════════════════════════════════════════════════════════════════
# FALLBACK: si TODOS los thr_r2 fallaron, usar un único bloque
# (cromosoma completo). LDpred2 funciona así, solo es menos eficiente
# en memoria. Mejor que abortar.
# ══════════════════════════════════════════════════════════════════
if (is.null(res_split)) {
  cat("  TODOS los thr_r2 fallaron — fallback: 1 bloque = cromosoma completo\n")
  block_sizes <- m_chr
  n_blocks    <- 1L
  blocks_list <- list(corr0)
  best_cost   <- NA_real_
  thr_used    <- NA_real_
} else {
  valid_idx <- which(!is.na(res_split$cost) &
                     !is.na(res_split$n_block) &
                     res_split$n_block > 0)
  res_split <- res_split[valid_idx, ]

  best_idx   <- which.min(res_split$cost + 0.5 * log(res_split$n_block))
  best_split <- res_split[best_idx, ]
  block_sizes <- best_split$all_size[[1]]
  n_blocks    <- length(block_sizes)
  best_cost   <- best_split$cost
  cat(sprintf("  %d bloques (thr_r2=%.3f, cost=%.1f, max=%d, min=%d, mediana=%.0f)\n",
      n_blocks, thr_used, best_cost,
      max(block_sizes), min(block_sizes), median(block_sizes)))

  if (sum(block_sizes) != m_chr) {
    save_failed(sprintf("sum(block_sizes)=%d != m_chr=%d",
                        sum(block_sizes), m_chr), "block_sizes")
  }

  block_starts <- c(1L, cumsum(block_sizes[-n_blocks]) + 1L)
  block_ends   <- cumsum(block_sizes)

  blocks_list <- lapply(seq_len(n_blocks), function(b) {
    bs <- block_starts[b]; be <- block_ends[b]
    corr0[bs:be, bs:be, drop=FALSE]
  })
}

# ── [5] Guardar ────────────────────────────────────────────────────
saveRDS(list(
  chr=chr, status="ok", n_snp=n_snp,
  n_blocks=n_blocks, block_sizes=block_sizes,
  blocks=blocks_list,
  info_snp=info_snp[, c("chr","pos","a0","a1","beta","beta_se","n_eff","_NUM_ID_")],
  sd_val=sd_vals,
  t_cor_min=t_cor,
  thr_r2_split=thr_used,
  bigsnpr_version=as.character(bv)
), out_file)

cat(sprintf("\n ✓ chr%s: %d bloques (thr_r2=%s), %s SNPs, snp_cor %.1f min\n",
    chr_label, n_blocks,
    ifelse(is.na(thr_used), "fallback", sprintf("%.3f", thr_used)),
    format(n_snp, big.mark=","), t_cor))
REOF

apptainer exec \
  --bind /mnt/cephfs:/mnt/cephfs \
  --bind "$TMPDIR:/tmp,$TMPDIR:/var/tmp" \
  --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 --env MKL_NUM_THREADS=1 \
  --env BLIS_NUM_THREADS=1 --env NUMEXPR_NUM_THREADS=1 \
  --env SLURM_ARRAY_TASK_ID="$CHR" \
  --env PLINK_DIR="$PLINK_DIR" --env GWAS_DIR="$GWAS_DIR" --env BLOCKS_DIR="$BLOCKS_DIR" \
  "$SIF" /opt/miniforge3/envs/r-bio/bin/Rscript "$R_SCRIPT"
# NOTA: NO pasar SLURM_NODELIST ni SLURM_JOB_NUM_NODES al contenedor.
# Causan segfaults en snp_cor (bigparallelr intenta parsear hostnames
# y falla porque /etc/slurm/slurm.conf no está en el contenedor).
# El warning de scontrol que aparece en .err es cosmético — ignorarlo.

rm -rf "$TMPDIR"

echo ""
echo "============================================================"
echo "  STEP 1 chr${CHR} FIN — $(date)"
echo "============================================================"
