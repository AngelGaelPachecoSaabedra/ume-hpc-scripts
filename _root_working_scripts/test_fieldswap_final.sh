#!/bin/bash
# Test FINAL — confirma que field-swap con .sbk del 2931 + p/first_i recomputados
# produce LDpred2-auto convergente. 5-10 min para chr22.

set -euo pipefail

SIF="/mnt/cephfs/biocontainers/images/jupyter-biotools-1.4.sif"
PREV_TMPDIR="/mnt/cephfs/scratch/angel.pacheco/ldpred2_hm3_tmp/2931"
GWAS_DIR="/mnt/cephfs/orgs/home/angel.pacheco/gwas_miami"
OUT_DIR="/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes/ldpred2"

R_SCRIPT=$(mktemp /tmp/test_final_XXXX.R)
cat << 'REOF' > "$R_SCRIPT"
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
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

# safe_as_SFBM del job 2931
safe_as_SFBM <- function(corr0, backingfile, compact = TRUE, eps = 1e-10) {
  stopifnot(is(corr0, "dsCMatrix"))
  col_nnz <- diff(corr0@p)
  empty <- which(col_nnz == 0L)
  if (length(empty) == 0L) return(as_SFBM(corr0, backingfile = backingfile, compact = compact))
  n <- corr0@Dim[1]; n_add <- length(empty); old_p <- corr0@p
  old_nnz <- length(corr0@i); new_nnz <- old_nnz + n_add
  cum_add <- cumsum(col_nnz == 0L)
  new_p   <- old_p + c(0L, cum_add)
  new_i   <- integer(new_nnz); new_x <- numeric(new_nnz)
  boundaries <- c(0L, empty); n_blocks <- length(boundaries)
  for (b in seq_len(n_blocks)) {
    block_first <- boundaries[b] + 1L
    block_last  <- if (b < n_blocks) boundaries[b + 1L] - 1L else n
    shift <- b - 1L
    if (block_last >= block_first) {
      old_start <- old_p[block_first] + 1L; old_end <- old_p[block_last + 1L]
      if (old_end >= old_start) {
        cnt <- old_end - old_start + 1L; ns <- old_start + shift
        new_i[ns:(ns+cnt-1L)] <- corr0@i[old_start:old_end]
        new_x[ns:(ns+cnt-1L)] <- corr0@x[old_start:old_end]
      }
    }
    if (b <= n_add) { ec <- empty[b]; dp <- new_p[ec] + 1L; new_i[dp] <- ec - 1L; new_x[dp] <- eps }
  }
  corr_fixed <- new("dsCMatrix", i=new_i, p=new_p, x=new_x, Dim=corr0@Dim,
                    Dimnames=corr0@Dimnames, uplo=corr0@uplo, factors=list())
  as_SFBM(corr_fixed, backingfile=backingfile, compact=compact)
}

# ─────────────────────────────────────────────────────────────────
# Sumstats chr22 (breve)
# ─────────────────────────────────────────────────────────────────
cat("[1] Sumstats chr22...\n")
map_hm3 <- readRDS(file.path(out_dir, "map_hm3_ldpred2.rds"))
if ("pos_hg38" %in% names(map_hm3)) map_hm3$pos_use <- map_hm3$pos_hg38

gf <- fread(file.path(gwas_dir, "combined-gwas-results-all-ctrl-females.txt.gz"),
            select=c("CHROM","GENPOS","ID","ALLELE0","ALLELE1","A1FREQ","INFO","N","BETA","SE","P"),
            showProgress=FALSE)
gm <- fread(file.path(gwas_dir, "combined-gwas-results-all-ctrl-males.txt.gz"),
            select=c("CHROM","GENPOS","ID","ALLELE0","ALLELE1","A1FREQ","INFO","N","BETA","SE","P"),
            showProgress=FALSE)
gwas <- merge(gf, gm, by="ID", suffixes=c(".f",".m")); rm(gf, gm); gc()

w_f <- 1/gwas$SE.f^2; w_m <- 1/gwas$SE.m^2
neff_f <- 4*17300*54196/(17300+54196); neff_m <- 4*8522*26733/(8522+26733)
K_prev <- (17300+8522)/(17300+54196+8522+26733)
neff_ldpred2 <- (neff_f + neff_m) * K_prev * (1 - K_prev)
sumstats <- data.frame(
  chr=gwas$CHROM.f, pos=gwas$GENPOS.f, a1=gwas$ALLELE1.f, a0=gwas$ALLELE0.f,
  freq=(gwas$A1FREQ.f*gwas$N.f + gwas$A1FREQ.m*gwas$N.m)/(gwas$N.f+gwas$N.m),
  info=pmin(gwas$INFO.f, gwas$INFO.m),
  beta=(gwas$BETA.f*w_f + gwas$BETA.m*w_m)/(w_f+w_m),
  beta_se=sqrt(1/(w_f+w_m)), n_eff=neff_ldpred2, stringsAsFactors=FALSE)
sumstats$p <- 2*pnorm(-abs(sumstats$beta/sumstats$beta_se))
sumstats <- sumstats[sumstats$info>=0.8 & sumstats$freq>=0.01 & sumstats$freq<=0.99 &
                     !is.na(sumstats$beta) & !is.na(sumstats$beta_se) & sumstats$beta_se>0, ]
sumstats$chr <- as.integer(sumstats$chr)
key_hm3 <- paste(map_hm3$chr, map_hm3$pos_use, sep=":")
sumstats <- sumstats[paste(sumstats$chr, sumstats$pos, sep=":") %in% key_hm3, ]
sumstats_chr22 <- sumstats[sumstats$chr == 22, ]

obj <- snp_attach(file.path(prev_tmp, "chr22.rds"))
G <- obj$genotypes
map_df <- data.frame(chr=obj$map$chromosome, pos=obj$map$physical.pos,
                     a1=obj$map$allele1, a0=obj$map$allele2, stringsAsFactors=FALSE)
info_snp <- snp_match(sumstats_chr22, map_df, strand_flip=FALSE,
                      join_by_pos=TRUE, match.min.prop=0.01)
ind_col <- info_snp$`_NUM_ID_`
sd_vals <- big_scale()(G, ind.col=ind_col)$scale
keep_var <- which(sd_vals > 1e-6)
if (length(keep_var) < nrow(info_snp)) { info_snp <- info_snp[keep_var, ]; ind_col <- info_snp$`_NUM_ID_` }
cat(sprintf("  Matched: %d\n", nrow(info_snp)))

# ─────────────────────────────────────────────────────────────────
# [2] Recomputar corr0 + safe_as_SFBM (idéntico al 2931)
# ─────────────────────────────────────────────────────────────────
cat("\n[2] Recomputando corr0 + safe_as_SFBM (produce .sbk idéntico al 2931)...\n")
t0 <- Sys.time()
corr0 <- snp_cor(Gna=G, ind.col=ind_col, infos.pos=info_snp$pos/1000,
                 size=3000, ncores=24L, thr_r2=0.005)
cat(sprintf("  snp_cor: %.1f min\n", as.numeric(difftime(Sys.time(), t0, units="mins"))))

tmp_bk <- tempfile("corr_fresh")
corr_fresh <- safe_as_SFBM(corr0, backingfile=tmp_bk, compact=TRUE)
rm(corr0); gc()

# ─────────────────────────────────────────────────────────────────
# [3A] LDpred2-auto con SFBM fresco
# ─────────────────────────────────────────────────────────────────
cat("\n[3A] LDpred2-auto con SFBM fresco (del recompute)...\n")
df_beta <- info_snp[, c("beta","beta_se","n_eff")]
t0 <- Sys.time()
multi_fresh <- snp_ldpred2_auto(
  corr=corr_fresh, df_beta=df_beta, h2_init=0.01,
  vec_p_init=seq_log(1e-4, 0.5, length.out=30),
  burn_in=500L, num_iter=500L, ncores=24L, sparse=TRUE)
cat(sprintf("  tiempo: %.1f min\n", as.numeric(difftime(Sys.time(), t0, units="mins"))))

h2_f <- sapply(multi_fresh, function(x) x$h2_est)
p_f  <- sapply(multi_fresh, function(x) x$p_est)
keep_f <- which(!is.na(h2_f) & h2_f>0.0005 & h2_f<0.5 & !is.na(p_f) & p_f>1e-6 & p_f<0.9)
cat(sprintf("  FRESH — chains válidas: %d/30, h2_med=%s, p_med=%s\n",
            length(keep_f),
            if (length(keep_f)>0) sprintf("%.4f", median(h2_f[keep_f])) else "NA",
            if (length(keep_f)>0) sprintf("%.2e", median(p_f[keep_f])) else "NA"))

# ─────────────────────────────────────────────────────────────────
# [3B] LDpred2-auto con field-swap al .sbk del 2931
# ─────────────────────────────────────────────────────────────────
cat("\n[3B] LDpred2-auto con FIELD-SWAP al .sbk original del 2931...\n")

# Dummy + swap de campos con los valores correctos del recompute
gen_cc <- get("SFBM_corr_compact_RC", envir=asNamespace("bigsparser"))
dummy  <- as(Matrix::Diagonal(3L), "symmetricMatrix")
tmp_swap <- tempfile("dummy_swap")
sfbm_swap <- gen_cc$new(spm=dummy, backingfile=tmp_swap)

# Tomar p y first_i del SFBM fresco (calculados con los MISMOS datos/params)
sbk_2931 <- file.path(prev_tmp, "corr_chr22.sbk")
sfbm_swap$backingfile <- normalizePath(sbk_2931)
sfbm_swap$nrow        <- corr_fresh$nrow
sfbm_swap$p           <- corr_fresh$p
sfbm_swap$first_i     <- corr_fresh$first_i
invisible(sfbm_swap$address)

cat(sprintf("  swap OK: %dx%d, nval=%d (apunta a .sbk del 2931)\n",
            sfbm_swap$nrow, sfbm_swap$ncol, sfbm_swap$nval))

t0 <- Sys.time()
multi_swap <- snp_ldpred2_auto(
  corr=sfbm_swap, df_beta=df_beta, h2_init=0.01,
  vec_p_init=seq_log(1e-4, 0.5, length.out=30),
  burn_in=500L, num_iter=500L, ncores=24L, sparse=TRUE)
cat(sprintf("  tiempo: %.1f min\n", as.numeric(difftime(Sys.time(), t0, units="mins"))))

h2_s <- sapply(multi_swap, function(x) x$h2_est)
p_s  <- sapply(multi_swap, function(x) x$p_est)
keep_s <- which(!is.na(h2_s) & h2_s>0.0005 & h2_s<0.5 & !is.na(p_s) & p_s>1e-6 & p_s<0.9)
cat(sprintf("  SWAP — chains válidas: %d/30, h2_med=%s, p_med=%s\n",
            length(keep_s),
            if (length(keep_s)>0) sprintf("%.4f", median(h2_s[keep_s])) else "NA",
            if (length(keep_s)>0) sprintf("%.2e", median(p_s[keep_s])) else "NA"))

# ─────────────────────────────────────────────────────────────────
# VEREDICTO
# ─────────────────────────────────────────────────────────────────
cat("\n[4] VEREDICTO:\n")
if (length(keep_s) >= 3L && length(keep_f) >= 3L) {
  cat("  ✓✓ Ambos enfoques producen LDpred2-auto convergente.\n")
  cat("  → Field-swap FUNCIONA. Los .sbk del 2931 SON reutilizables.\n")
  cat("  → PLAN: recompute ligero de corr0→safe_as_SFBM (solo para p/first_i),\n")
  cat("    swap al .sbk del 2931, LDpred2-auto. Total ~2-3 días vs 2 semanas.\n")
} else if (length(keep_f) >= 3L && length(keep_s) < 3L) {
  cat("  ⚠ FRESH converge pero SWAP NO. Hay un bug sutil en el swap.\n")
  cat("  → Revisar nval o alineación del .sbk.\n")
} else {
  cat("  ✗ Ninguno converge. Revisar parámetros LDpred2 o datos.\n")
}

unlink(c(paste0(tmp_bk, ".sbk"), paste0(tmp_swap, ".sbk")))
cat("\nTest FINAL completado.\n")
REOF

export PREV_TMPDIR OUT_DIR GWAS_DIR
apptainer exec --bind /mnt/cephfs:/mnt/cephfs \
    --env PREV_TMPDIR="$PREV_TMPDIR" --env OUT_DIR="$OUT_DIR" --env GWAS_DIR="$GWAS_DIR" \
    --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 --env MKL_NUM_THREADS=1 \
    "$SIF" /opt/miniforge3/envs/r-bio/bin/Rscript "$R_SCRIPT"

rm -f "$R_SCRIPT"
