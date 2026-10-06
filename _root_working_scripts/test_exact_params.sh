#!/bin/bash
# Test DEFINITIVO v2 — fix paralelismo BLAS
# Recompute snp_cor con los parámetros EXACTOS del job 2931.
# Compara byte a byte el .sbk resultante con el .sbk del job original.

set -euo pipefail

SIF="/mnt/cephfs/biocontainers/images/jupyter-biotools-1.4.sif"
PREV_TMPDIR="/mnt/cephfs/scratch/angel.pacheco/ldpred2_hm3_tmp/2931"
GWAS_DIR="/mnt/cephfs/orgs/home/angel.pacheco/gwas_miami"
OUT_DIR="/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes/ldpred2"

R_SCRIPT=$(mktemp /tmp/test_exact_XXXX.R)
cat << 'REOF' > "$R_SCRIPT"
# ── FIX paralelismo: forzar BLAS a 1 thread ANTES de cargar bigstatsr ────
Sys.setenv(OMP_NUM_THREADS = "1")
Sys.setenv(OPENBLAS_NUM_THREADS = "1")
Sys.setenv(MKL_NUM_THREADS = "1")
Sys.setenv(VECLIB_MAXIMUM_THREADS = "1")
if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
  RhpcBLASctl::blas_set_num_threads(1L)
  RhpcBLASctl::omp_set_num_threads(1L)
}

suppressMessages({
  library(bigsnpr); library(bigsparser); library(data.table); library(Matrix)
})
options(bigstatsr.check.parallel.blas = FALSE)

prev_tmp <- Sys.getenv("PREV_TMPDIR")
out_dir  <- Sys.getenv("OUT_DIR")
gwas_dir <- Sys.getenv("GWAS_DIR")

# ============================================================
# safe_as_SFBM — COPIA EXACTA del job 2931
# ============================================================
safe_as_SFBM <- function(corr0, backingfile, compact = TRUE, eps = 1e-10) {
  stopifnot(is(corr0, "dsCMatrix"))
  col_nnz <- diff(corr0@p)
  empty <- which(col_nnz == 0L)
  if (length(empty) == 0L) {
    cat("      [SFBM] Sin columnas vacías\n")
    return(as_SFBM(corr0, backingfile = backingfile, compact = compact))
  }
  n <- corr0@Dim[1]; n_add <- length(empty); old_p <- corr0@p
  old_nnz <- length(corr0@i); new_nnz <- old_nnz + n_add
  cat(sprintf("      [SFBM] %d columnas vacías → reparando\n", n_add))
  cum_add <- cumsum(col_nnz == 0L)
  new_p   <- old_p + c(0L, cum_add)
  new_i   <- integer(new_nnz); new_x <- numeric(new_nnz)
  boundaries <- c(0L, empty); n_blocks <- length(boundaries)
  for (b in seq_len(n_blocks)) {
    block_first <- boundaries[b] + 1L
    block_last  <- if (b < n_blocks) boundaries[b + 1L] - 1L else n
    shift <- b - 1L
    if (block_last >= block_first) {
      old_start <- old_p[block_first] + 1L
      old_end   <- old_p[block_last + 1L]
      if (old_end >= old_start) {
        cnt <- old_end - old_start + 1L
        ns  <- old_start + shift
        new_i[ns:(ns+cnt-1L)] <- corr0@i[old_start:old_end]
        new_x[ns:(ns+cnt-1L)] <- corr0@x[old_start:old_end]
      }
    }
    if (b <= n_add) {
      ec <- empty[b]; dp <- new_p[ec] + 1L
      new_i[dp] <- ec - 1L; new_x[dp] <- eps
    }
  }
  corr_fixed <- new("dsCMatrix",
    i=new_i, p=new_p, x=new_x, Dim=corr0@Dim,
    Dimnames=corr0@Dimnames, uplo=corr0@uplo, factors=list())
  as_SFBM(corr_fixed, backingfile=backingfile, compact=compact)
}

# ============================================================
# [1] Reconstruir sumstats
# ============================================================
cat("[1] Sumstats chr22 (idéntico al job 2931)...\n")
map_hm3 <- readRDS(file.path(out_dir, "map_hm3_ldpred2.rds"))
if ("pos_hg38" %in% names(map_hm3)) map_hm3$pos_use <- map_hm3$pos_hg38

gf <- fread(file.path(gwas_dir, "combined-gwas-results-all-ctrl-females.txt.gz"),
            select=c("CHROM","GENPOS","ID","ALLELE0","ALLELE1","A1FREQ","INFO","N","BETA","SE","P"),
            showProgress=FALSE)
gm <- fread(file.path(gwas_dir, "combined-gwas-results-all-ctrl-males.txt.gz"),
            select=c("CHROM","GENPOS","ID","ALLELE0","ALLELE1","A1FREQ","INFO","N","BETA","SE","P"),
            showProgress=FALSE)
gwas <- merge(gf, gm, by="ID", suffixes=c(".f",".m"))
rm(gf, gm); gc()

w_f <- 1/gwas$SE.f^2; w_m <- 1/gwas$SE.m^2
neff_f <- 4*17300*54196/(17300+54196); neff_m <- 4*8522*26733/(8522+26733)
K_prev <- (17300+8522)/(17300+54196+8522+26733)
neff_ldpred2 <- (neff_f + neff_m) * K_prev * (1 - K_prev)

sumstats <- data.frame(
  chr=gwas$CHROM.f, pos=gwas$GENPOS.f,
  a1=gwas$ALLELE1.f, a0=gwas$ALLELE0.f,
  freq=(gwas$A1FREQ.f*gwas$N.f + gwas$A1FREQ.m*gwas$N.m)/(gwas$N.f+gwas$N.m),
  info=pmin(gwas$INFO.f, gwas$INFO.m),
  beta=(gwas$BETA.f*w_f + gwas$BETA.m*w_m)/(w_f+w_m),
  beta_se=sqrt(1/(w_f+w_m)),
  n_eff=neff_ldpred2, stringsAsFactors=FALSE)
sumstats$p <- 2*pnorm(-abs(sumstats$beta/sumstats$beta_se))
sumstats <- sumstats[sumstats$info>=0.8 & sumstats$freq>=0.01 & sumstats$freq<=0.99 &
                     !is.na(sumstats$beta) & !is.na(sumstats$beta_se) & sumstats$beta_se>0, ]
sumstats$chr <- as.character(sumstats$chr)
sumstats$chr[sumstats$chr %in% c("X","chrX")] <- "23"
sumstats$chr <- as.integer(sumstats$chr)

key_hm3 <- paste(map_hm3$chr, map_hm3$pos_use, sep=":")
sumstats <- sumstats[paste(sumstats$chr, sumstats$pos, sep=":") %in% key_hm3, ]
sumstats_chr22 <- sumstats[sumstats$chr == 22, ]
cat(sprintf("  sumstats chr22: %d\n", nrow(sumstats_chr22)))

# ============================================================
# [2] Matching + filtro monomórficas
# ============================================================
cat("\n[2] Matching chr22...\n")
obj <- snp_attach(file.path(prev_tmp, "chr22.rds"))
G <- obj$genotypes
map_df <- data.frame(chr=obj$map$chromosome, pos=obj$map$physical.pos,
                     a1=obj$map$allele1, a0=obj$map$allele2, stringsAsFactors=FALSE)
info_snp <- snp_match(sumstats_chr22, map_df, strand_flip=FALSE,
                      join_by_pos=TRUE, match.min.prop=0.01)
ind_col <- info_snp$`_NUM_ID_`
sd_vals <- big_scale()(G, ind.col=ind_col)$scale
keep_var <- which(sd_vals > 1e-6)
if (length(keep_var) < nrow(info_snp)) {
  info_snp <- info_snp[keep_var, ]
  ind_col  <- info_snp$`_NUM_ID_`
}
cat(sprintf("  Matched: %d (debería ser 14,360)\n", nrow(info_snp)))

# ============================================================
# [3] snp_cor con parámetros EXACTOS + BLAS=1
# ============================================================
cat("\n[3] snp_cor con parámetros EXACTOS (pos/1000, size=3000, thr_r2=0.005)...\n")
cat(sprintf("  BLAS threads: %s\n", Sys.getenv("OPENBLAS_NUM_THREADS")))

t0 <- Sys.time()
corr0 <- snp_cor(
  Gna       = G,
  ind.col   = ind_col,
  infos.pos = info_snp$pos / 1000,
  size      = 3000,
  ncores    = 24L,
  thr_r2    = 0.005
)
t_elapsed <- as.numeric(difftime(Sys.time(), t0, units="mins"))
cat(sprintf("  Tiempo snp_cor: %.1f min (original: 111.4 min para chr22)\n", t_elapsed))
cat(sprintf("  dim: %dx%d, nnz: %d (original: 9,933,050)\n",
            nrow(corr0), ncol(corr0), length(corr0@x)))

# ============================================================
# [4] safe_as_SFBM
# ============================================================
cat("\n[4] safe_as_SFBM...\n")
tmp_bk <- tempfile("test_exact")
sfbm_new <- safe_as_SFBM(corr0, backingfile=tmp_bk, compact=TRUE)

sbk_new_path <- paste0(tmp_bk, ".sbk")
sbk_new_size <- file.info(sbk_new_path)$size
sbk_real      <- file.path(prev_tmp, "corr_chr22.sbk")
sbk_real_size <- file.info(sbk_real)$size

cat(sprintf("\n[5] COMPARACIÓN:\n"))
cat(sprintf("  .sbk original : %d bytes\n", sbk_real_size))
cat(sprintf("  .sbk recomput : %d bytes\n", sbk_new_size))

if (sbk_new_size == sbk_real_size) {
  cat("  ✓ TAMAÑOS IDÉNTICOS\n\n")
  con1 <- file(sbk_real,     "rb"); a <- readBin(con1, "double", n=sbk_real_size/8, size=8); close(con1)
  con2 <- file(sbk_new_path, "rb"); b <- readBin(con2, "double", n=sbk_new_size/8, size=8); close(con2)
  max_diff  <- max(abs(a - b))
  mean_diff <- mean(abs(a - b))
  cat(sprintf("  max|a-b|  : %.6e\n", max_diff))
  cat(sprintf("  mean|a-b| : %.6e\n", mean_diff))
  if (max_diff < 1e-10) {
    cat("\n  ✓✓ .sbk BIT-IDÉNTICO — recomputable determinísticamente\n")
    cat("  → LOS .sbk DEL 2931 SON REUTILIZABLES vía field-swap\n")
  } else if (max_diff < 1e-6) {
    cat("\n  ≈ diferencias numéricas pequeñas (paralelismo no-determinista)\n")
    cat("  → .sbk del 2931 son equivalentes estadísticamente, reutilizables\n")
  } else {
    cat(sprintf("\n  ⚠ diferencias grandes (max=%.2e) — investigar\n", max_diff))
  }
  cat(sprintf("\n  p[length(p)] (nval total): %d\n", sfbm_new$p[length(sfbm_new$p)]))
  cat(sprintf("  first_i[1:10]: %s\n", paste(head(sfbm_new$first_i, 10), collapse=",")))
} else {
  cat(sprintf("  ✗ TAMAÑOS DIFERENTES: ratio = %.4f\n", sbk_real_size/sbk_new_size))
  cat("  → Algo cambió entre el job 2931 y ahora\n")
}

unlink(c(sbk_new_path, paste0(tmp_bk, ".rds")))
cat("\nTest definitivo v2 completado.\n")
REOF

export PREV_TMPDIR OUT_DIR GWAS_DIR

apptainer exec --bind /mnt/cephfs:/mnt/cephfs \
    --env PREV_TMPDIR="$PREV_TMPDIR" --env OUT_DIR="$OUT_DIR" --env GWAS_DIR="$GWAS_DIR" \
    --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 --env MKL_NUM_THREADS=1 \
    "$SIF" /opt/miniforge3/envs/r-bio/bin/Rscript "$R_SCRIPT"

rm -f "$R_SCRIPT"
