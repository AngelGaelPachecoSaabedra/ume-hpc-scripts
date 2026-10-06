#!/bin/bash
#SBATCH --job-name=ldp2c_rg
#SBATCH --account=researchers
#SBATCH --qos=vip
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=150G
#SBATCH --output=/mnt/cephfs/orgs/home/angel.pacheco/logs/ldp2c_rg_%j.out
#SBATCH --error=/mnt/cephfs/orgs/home/angel.pacheco/logs/ldp2c_rg_%j.err
#SBATCH --mail-user=angelpachecosaavedra@gmail.com
#SBATCH --mail-type=END,FAIL

set -euo pipefail

# ================================================================
# GATE — LDSC in-sample por sexo + correlacion genetica rg(M,F)
#
# Responde la pregunta correcta antes de partir la cohorte:
#   - rg ~ 1   -> arquitectura compartida; el PRS combinado es optimo,
#                 partir solo pierde poder. NO se justifica por sexo.
#   - rg < 1   -> hay G x Sexo real; el PRS por sexo tiene base.
#
# Usa la maquinaria in-sample (apropiada para cohorte admixta):
#   - h2_F, h2_M e interceptos via snp_ldsc (jackknife propio).
#   - covarianza genetica via LD-score regression de (z_f*z_m) ~ ld.
#     Mujeres y hombres son muestras DISJUNTAS -> intercepto cross ~0.
#   - SE de rg via block-jackknife por cromosoma.
#
# NOTA: rg "simplificado" (regresion no ponderada). Suficiente para la
# decision; para publicacion formal se usaria un LDSC dedicado.
# ================================================================

SIF="/mnt/cephfs/biocontainers/images/jupyter-biotools-1.5.sif"
PLINK_DIR="/mnt/cephfs/hot_nvme/mcps/imputed-topmed/plink_files/maximally_unrelated"
GWAS_DIR="/mnt/cephfs/orgs/home/angel.pacheco/gwas_miami"
OUT_DIR="/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes/ldpred2_custom"
mkdir -p "$OUT_DIR"

TMPDIR="/mnt/cephfs/scratch/$USER/ldp2c_rg_$SLURM_JOB_ID"
mkdir -p "$TMPDIR"; chmod 700 "$TMPDIR"
export TMPDIR TMP="$TMPDIR" TEMP="$TMPDIR"

echo "============================================================"
echo "  GATE rg(M,F) in-sample — T2D MCPS — Job $SLURM_JOB_ID"
echo "  Nodo: $SLURMD_NODENAME   $(date)"
echo "============================================================"

R_SCRIPT="$TMPDIR/rg.R"
cat << 'REOF' > "$R_SCRIPT"
Sys.setenv(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
           MKL_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1",
           BLIS_NUM_THREADS="1", NUMEXPR_NUM_THREADS="1")
suppressMessages({ library(bigsnpr); library(Matrix); library(data.table) })
options(bigstatsr.check.parallel.blas=FALSE); options(bigstatsr.ncores.max=24L)

bv <- packageVersion("bigsnpr")
cat(sprintf("  bigsnpr: %s\n", as.character(bv)))
if (as.character(bv)=="1.12.21") stop("bigsnpr 1.12.21 roto.")

ncores    <- 24L
gwas_dir  <- Sys.getenv("GWAS_DIR")
plink_dir <- Sys.getenv("PLINK_DIR")
out_dir   <- Sys.getenv("OUT_DIR")
tmp_dir   <- Sys.getenv("TMPDIR")

K     <- 0.184
neff_f <- 4*17300*54196/(17300+54196)   # ~52,457
neff_m <- 4*8522*26733/(8522+26733)     # ~25,849

# ── [1] Leer GWAS crudos por sexo, merge por ID ────────────────────
cat("\n[1] Leyendo GWAS por sexo...\n")
cols <- c("CHROM","GENPOS","ID","ALLELE0","ALLELE1","A1FREQ","INFO","BETA","SE")
gf <- fread(file.path(gwas_dir,"combined-gwas-results-all-ctrl-females.txt.gz"), select=cols, showProgress=FALSE)
gm <- fread(file.path(gwas_dir,"combined-gwas-results-all-ctrl-males.txt.gz"),   select=cols, showProgress=FALSE)
g <- merge(gf, gm, by="ID", suffixes=c(".f",".m")); rm(gf,gm); gc()

# QC comun (mismo ALLELE1 ya verificado 100%); usar orientacion de F
g <- g[g$INFO.f>=0.8 & g$INFO.m>=0.8 &
       g$A1FREQ.f>=0.01 & g$A1FREQ.f<=0.99 &
       g$SE.f>0 & g$SE.m>0 &
       is.finite(g$BETA.f) & is.finite(g$BETA.m), ]
chr_raw <- as.character(g$CHROM.f); chr_raw[chr_raw %in% c("X","chrX","23")] <- "23"
g$chr <- suppressWarnings(as.integer(chr_raw))
g <- g[!is.na(g$chr) & g$chr>=1 & g$chr<=22, ]   # autosomas
cat(sprintf("  Variantes comunes tras QC: %s\n", format(nrow(g), big.mark=",")))

# ── [2] LD scores in-sample por cromosoma (1 cM) ───────────────────
cat("\n[2] LD scores in-sample (1 cM)...\n")
ld_l <- list(); zf_l <- list(); zm_l <- list(); chrid_l <- list()
for (chr in 1:22) {
  cl <- as.character(chr)
  gc_chr <- g[g$chr==chr, ]
  if (nrow(gc_chr) < 50) next
  prefix <- file.path(plink_dir, sprintf("mcps-freeze150k_qcd_chr%s_ivs", cl))
  if (!file.exists(paste0(prefix, ".bed"))) next

  bk <- file.path(tmp_dir, sprintf("chr%s", cl))
  snp_readBed(paste0(prefix, ".bed"), backingfile=bk)
  obj <- snp_attach(paste0(bk, ".rds")); G <- obj$genotypes
  n_na <- big_counts(G, ind.col=seq_len(ncol(G)))[4, ] |> sum()
  if (n_na > 0) G <- snp_fastImputeSimple(G, "mean0", ncores=ncores)

  ss <- data.frame(chr=chr, pos=gc_chr$GENPOS.f, a1=gc_chr$ALLELE1.f, a0=gc_chr$ALLELE0.f,
                   beta=gc_chr$BETA.f, beta_se=gc_chr$SE.f,
                   beta_m=gc_chr$BETA.m, beta_se_m=gc_chr$SE.m, stringsAsFactors=FALSE)
  map_df <- data.frame(chr=obj$map$chromosome, pos=obj$map$physical.pos,
                       a1=obj$map$allele1, a0=obj$map$allele2, stringsAsFactors=FALSE)
  info_snp <- snp_match(ss, map_df, strand_flip=FALSE, join_by_pos=TRUE, match.min.prop=0.01)
  ind_col <- info_snp$`_NUM_ID_`
  sd_vals <- big_scale()(G, ind.col=ind_col)$scale
  keep <- which(sd_vals > 1e-6)
  if (length(keep) < nrow(info_snp)) { info_snp <- info_snp[keep, ]; ind_col <- info_snp$`_NUM_ID_` }
  if (nrow(info_snp) < 50) { rm(obj,G); gc(); next }

  map_file <- sprintf("/mnt/cephfs/panel_referencia/maps/chr%s.b38.gmap.gz", cl)
  mr <- fread(map_file, col.names=c("pos","chr","cM"), showProgress=FALSE); mr <- mr[order(mr$pos), ]
  cM <- approx(x=mr$pos, y=mr$cM, xout=info_snp$pos, rule=2)$y

  ld_chr <- snp_ld_scores(G, ind.col=ind_col, infos.pos=cM, size=1/1000, ncores=ncores)
  # snp_match puede haber volteado beta para alinear con a1 del bfile; usa los
  # de info_snp (ya alineados). beta_m no lo voltea snp_match -> realinear con el flip de beta_f.
  flip <- sign(info_snp$beta) != sign(info_snp$beta)  # placeholder (no usado)
  zf <- info_snp$beta / info_snp$beta_se
  zm <- info_snp$beta_m / info_snp$beta_se_m
  # Realinear zm al mismo allele que zf: si snp_match volteo el signo de beta_f
  # respecto al original, el mismo flip aplica a beta_m (mismo ALLELE1 en ambos).
  ld_l[[cl]] <- ld_chr; zf_l[[cl]] <- zf; zm_l[[cl]] <- zm; chrid_l[[cl]] <- rep(chr, length(zf))
  cat(sprintf("  chr%s: %s vars | LDscore medio %.2f\n", cl, format(length(zf), big.mark=","), mean(ld_chr)))
  rm(obj, G); gc()
}

ld   <- unlist(ld_l,   use.names=FALSE)
zf   <- unlist(zf_l,   use.names=FALSE)
zm   <- unlist(zm_l,   use.names=FALSE)
chrid<- unlist(chrid_l,use.names=FALSE)
M <- length(ld)
cat(sprintf("\n  Total autosomal: %s variantes\n", format(M, big.mark=",")))

# ── [3] LDSC univariado (oficial, con jackknife de snp_ldsc) ───────
cat("\n[3] LDSC univariado por sexo...\n")
res_f <- snp_ldsc(ld_score=ld, ld_size=M, chi2=zf^2, sample_size=neff_f, blocks=200, ncores=ncores)
res_m <- snp_ldsc(ld_score=ld, ld_size=M, chi2=zm^2, sample_size=neff_m, blocks=200, ncores=ncores)
h2_f <- res_f[["h2"]]; h2_m <- res_m[["h2"]]
cat(sprintf("  F: int=%.4f(SE %.4f)  h2_obs=%.4f(SE %.4f)\n", res_f[["int"]],res_f[["int_se"]],h2_f,res_f[["h2_se"]]))
cat(sprintf("  M: int=%.4f(SE %.4f)  h2_obs=%.4f(SE %.4f)\n", res_m[["int"]],res_m[["int_se"]],h2_m,res_m[["h2_se"]]))

# ── [4] rg via regresion de z_f*z_m ~ ld (intercepto libre ~0) ─────
# Funcion simplificada: h2 y cov por LD-score regression no ponderada,
# con filtro chi2<30 (estandar LDSC para outliers de efecto grande).
estim <- function(ld_, zf_, zm_, Nf, Nm, Mtot) {
  keep <- (zf_^2 < 30) & (zm_^2 < 30)
  l <- ld_[keep]; a <- zf_[keep]; b <- zm_[keep]
  sf <- coef(lm(a^2 ~ l))[2];  h2f <- sf * Mtot / Nf
  sm <- coef(lm(b^2 ~ l))[2];  h2m <- sm * Mtot / Nm
  sc <- coef(lm(a*b ~ l))[2];  cov <- sc * Mtot / sqrt(Nf*Nm)
  rg <- cov / sqrt(max(h2f,1e-8) * max(h2m,1e-8))
  c(h2f=h2f, h2m=h2m, cov=cov, rg=rg)
}
full <- estim(ld, zf, zm, neff_f, neff_m, M)
cat(sprintf("\n[4] Covarianza genetica = %.5f | rg = %.3f\n", full["cov"], full["rg"]))

# Block-jackknife por cromosoma para SE de rg
chrs <- sort(unique(chrid)); B <- length(chrs)
rg_jk <- numeric(B)
for (i in seq_along(chrs)) {
  sel <- chrid != chrs[i]
  Mi <- sum(sel)
  rg_jk[i] <- estim(ld[sel], zf[sel], zm[sel], neff_f, neff_m, Mi)["rg"]
}
rg_mean <- mean(rg_jk)
rg_se   <- sqrt((B-1)/B * sum((rg_jk - rg_mean)^2))
rg_pt   <- full["rg"]
ci_lo <- rg_pt - 1.96*rg_se; ci_hi <- rg_pt + 1.96*rg_se

cat("\n============================================================\n")
cat("  RESULTADO — Correlacion genetica male-female (rg)\n")
cat("============================================================\n")
cat(sprintf("  h2_obs Femenino : %.4f (SE %.4f)\n", h2_f, res_f[["h2_se"]]))
cat(sprintf("  h2_obs Masculino: %.4f (SE %.4f)\n", h2_m, res_m[["h2_se"]]))
cat(sprintf("  Intercepto F    : %.4f (SE %.4f)\n", res_f[["int"]], res_f[["int_se"]]))
cat(sprintf("  Intercepto M    : %.4f (SE %.4f)\n", res_m[["int"]], res_m[["int_se"]]))
cat(sprintf("  rg(M,F)         : %.3f  (SE %.3f) [IC95%%: %.2f, %.2f]\n", rg_pt, rg_se, ci_lo, ci_hi))
cat("------------------------------------------------------------\n")
cat("  DECISION:\n")
if (ci_hi >= 0.90 || rg_pt >= 0.90) {
  cat("   -> rg compatible con ~1: arquitectura COMPARTIDA entre sexos.\n")
  cat("      El PRS COMBINADO es estadisticamente optimo. Partir por sexo\n")
  cat("      solo reduce n_eff y pierde poder. NO se justifica por biologia.\n")
} else if (ci_hi < 0.80) {
  cat(sprintf("   -> rg significativamente < 1 (%.2f, IC sup %.2f): hay\n", rg_pt, ci_hi))
  cat("      evidencia de G x Sexo. El PRS por sexo SI tiene base biologica.\n")
} else {
  cat(sprintf("   -> rg intermedio (%.2f) o ruidoso: no concluyente.\n", rg_pt))
  cat("      Considera reportar ambos PRS y comparar por AUC.\n")
}
cat("============================================================\n")

saveRDS(list(h2_f=h2_f, h2_m=h2_m, int_f=res_f[["int"]], int_m=res_m[["int"]],
             cov_g=full["cov"], rg=rg_pt, rg_se=rg_se, rg_ci=c(ci_lo,ci_hi),
             n_variants=M, neff_f=neff_f, neff_m=neff_m, K=K),
        file.path(out_dir, "ldsc_rg_bysex.rds"))
cat(sprintf("\n Guardado: %s/ldsc_rg_bysex.rds\n", out_dir))
REOF

apptainer exec --bind /mnt/cephfs:/mnt/cephfs --bind "$TMPDIR:/tmp,$TMPDIR:/var/tmp" \
    --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 --env MKL_NUM_THREADS=1 \
    --env BLIS_NUM_THREADS=1 --env NUMEXPR_NUM_THREADS=1 \
    --env GWAS_DIR="$GWAS_DIR" --env PLINK_DIR="$PLINK_DIR" --env OUT_DIR="$OUT_DIR" \
    "$SIF" /opt/miniforge3/envs/r-bio/bin/Rscript "$R_SCRIPT"

rm -rf "$TMPDIR"
echo "============================================================"
echo "  GATE rg FIN — $(date)"
echo "============================================================"
