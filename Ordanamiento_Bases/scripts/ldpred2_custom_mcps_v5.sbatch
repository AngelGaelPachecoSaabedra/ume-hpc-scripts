#!/bin/bash
#SBATCH --job-name=ldp2_custom
#SBATCH --account=researchers
#SBATCH --qos=vip
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=120G
#SBATCH --output=/mnt/cephfs/orgs/home/angel.pacheco/logs/ldp2_custom_%j.out
#SBATCH --error=/mnt/cephfs/orgs/home/angel.pacheco/logs/ldp2_custom_%j.err
#SBATCH --mail-user=angelpachecosaavedra@gmail.com
#SBATCH --mail-type=BEGIN,END,FAIL

set -euo pipefail

# ================================================================
# PIPELINE C — LDpred2 CUSTOM MCPS PRUNED v5 (Chr 1-23)
#
# v4 fixes (vs v2/v3):
#   1. thr_r2=0.005 en snp_cor → previene integer overflow en cumsum
#   2. safe_as_SFBM: manipulación directa de slots dsCMatrix
#      para columnas vacías (sin aritmética de Matrix → sin dgCMatrix)
#   3. chrX IDs: chrX:pos:... además de chr23:pos:... (ya en v2)
#   4. r² 0.9, ventana 500/50 (ya en v2)
# ================================================================

SIF="/mnt/cephfs/biocontainers/images/jupyter-biotools-1.4.sif"

PLINK_DIR="/mnt/cephfs/hot_nvme/mcps/imputed-topmed/plink_files/maximally_unrelated"
GWAS_DIR="/mnt/cephfs/orgs/home/angel.pacheco/gwas_miami"
OUT_DIR="/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes/ldpred2"
mkdir -p "$OUT_DIR"

TMPDIR="/mnt/cephfs/scratch/$USER/ldp2_custom_tmp/$SLURM_JOB_ID"
mkdir -p "$TMPDIR"
chmod 700 "$TMPDIR"
export TMPDIR TMP="$TMPDIR" TEMP="$TMPDIR"

echo "============================================================"
echo "  PIPELINE C — LDpred2 CUSTOM MCPS PRUNED v4"
echo "  Nodo: $SLURMD_NODENAME"
echo "  RAM: $(free -g | awk 'NR==2{print $2}')G"
echo "  $(date)"
echo "============================================================"

R_SCRIPT="$TMPDIR/run_ldpred2_custom_v2.R"
cat << 'REOF' > "$R_SCRIPT"
library(bigsnpr)
library(data.table)
library(methods)

ncores <- 24
options(bigstatsr.check.parallel.blas = FALSE)
options(bigstatsr.ncores.max = 24)

# ================================================================
# safe_as_SFBM: convierte dsCMatrix → SFBM sin aritmética de Matrix
# (ver comentarios en ldpred2_v7_hapmap3.sbatch para explicación)
# ================================================================
safe_as_SFBM <- function(corr0, backingfile, compact = TRUE, eps = 1e-10) {
    stopifnot(is(corr0, "dsCMatrix"))
    col_nnz <- diff(corr0@p)
    empty <- which(col_nnz == 0L)

    if (length(empty) == 0L) {
        cat("      [SFBM] Sin columnas vacías\n")
        return(as_SFBM(corr0, backingfile = backingfile, compact = compact))
    }

    n      <- corr0@Dim[1]
    n_add  <- length(empty)
    old_p  <- corr0@p
    old_nnz <- length(corr0@i)
    new_nnz <- old_nnz + n_add

    cat(sprintf("      [SFBM] %d columnas vacías de %s → reparando slots\n",
        n_add, format(n, big.mark=",")))

    if (as.double(new_nnz) > .Machine$integer.max)
        stop(sprintf("Total nnz (%s) excedería límite 32-bit",
             format(as.double(new_nnz), big.mark=",")))

    cum_add <- cumsum(col_nnz == 0L)
    new_p   <- old_p + c(0L, cum_add)

    new_i <- integer(new_nnz)
    new_x <- numeric(new_nnz)

    boundaries <- c(0L, empty)
    n_blocks   <- length(boundaries)

    for (b in seq_len(n_blocks)) {
        block_first <- boundaries[b] + 1L
        block_last  <- if (b < n_blocks) boundaries[b + 1L] - 1L else n
        shift       <- b - 1L

        if (block_last >= block_first) {
            old_start <- old_p[block_first] + 1L
            old_end   <- old_p[block_last + 1L]
            if (old_end >= old_start) {
                cnt <- old_end - old_start + 1L
                ns  <- old_start + shift
                new_i[ns:(ns + cnt - 1L)] <- corr0@i[old_start:old_end]
                new_x[ns:(ns + cnt - 1L)] <- corr0@x[old_start:old_end]
            }
        }

        if (b <= n_add) {
            ec <- empty[b]
            dp <- new_p[ec] + 1L
            new_i[dp] <- ec - 1L
            new_x[dp] <- eps
        }
    }

    cat(sprintf("      [SFBM] Reparado: %s → %s entries\n",
        format(old_nnz, big.mark=","), format(new_nnz, big.mark=",")))

    corr_fixed <- new("dsCMatrix",
        i = new_i, p = new_p, x = new_x,
        Dim = corr0@Dim, Dimnames = corr0@Dimnames,
        uplo = corr0@uplo, factors = list())

    as_SFBM(corr_fixed, backingfile = backingfile, compact = compact)
}

out_dir   <- "/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes/ldpred2"
gwas_dir  <- "/mnt/cephfs/orgs/home/angel.pacheco/gwas_miami"
plink_dir <- "/mnt/cephfs/hot_nvme/mcps/imputed-topmed/plink_files/maximally_unrelated"
tmp_dir   <- Sys.getenv("TMPDIR")

# ================================================================
# [1] Cargar GWAS y meta-análisis IVW
# ================================================================
cat("\n[1] Cargando GWAS Meta-análisis...\n")

gf <- fread(file.path(gwas_dir, "combined-gwas-results-all-ctrl-females.txt.gz"),
            select=c("CHROM","GENPOS","ID","ALLELE0","ALLELE1","A1FREQ","INFO","N","BETA","SE","P"),
            showProgress=FALSE)
gm <- fread(file.path(gwas_dir, "combined-gwas-results-all-ctrl-males.txt.gz"),
            select=c("CHROM","GENPOS","ID","ALLELE0","ALLELE1","A1FREQ","INFO","N","BETA","SE","P"),
            showProgress=FALSE)
gwas <- merge(gf, gm, by="ID", suffixes=c(".f", ".m"))
rm(gf, gm); gc()

w_f <- 1/gwas$SE.f^2; w_m <- 1/gwas$SE.m^2

# n_eff para trait binario: 4 * cases * controls / (cases + controls)
neff_f <- 4 * 17300 * 54196 / (17300 + 54196)  # = 52,456
neff_m <- 4 *  8522 * 26733 / ( 8522 + 26733)  # = 25,848

# FIX: REGENIE log-odds → LDpred2 escala lineal
# n_eff_ldpred2 = n_eff_binario × K × (1-K), K = prevalencia
K_prev <- (17300 + 8522) / (17300 + 54196 + 8522 + 26733)  # = 0.2419
neff_ldpred2 <- (neff_f + neff_m) * K_prev * (1 - K_prev)   # ≈ 14,359

sumstats <- data.frame(
    chr = gwas$CHROM.f, pos = gwas$GENPOS.f,
    a1 = gwas$ALLELE1.f, a0 = gwas$ALLELE0.f,
    freq = (gwas$A1FREQ.f*gwas$N.f + gwas$A1FREQ.m*gwas$N.m)/(gwas$N.f+gwas$N.m),
    info = pmin(gwas$INFO.f, gwas$INFO.m),
    beta = (gwas$BETA.f*w_f + gwas$BETA.m*w_m)/(w_f+w_m),
    beta_se = sqrt(1/(w_f+w_m)), n_eff = neff_ldpred2,
    stringsAsFactors = FALSE
)
sumstats$p <- 2 * pnorm(-abs(sumstats$beta / sumstats$beta_se))
rm(gwas); gc()

# QC
sumstats <- sumstats[sumstats$info >= 0.8 & sumstats$freq >= 0.01 & sumstats$freq <= 0.99 &
                     !is.na(sumstats$beta) & !is.na(sumstats$beta_se) & sumstats$beta_se > 0, ]

# Guardar chr original ANTES de convertir X→23
# (necesitamos ambas formas para los IDs de PLINK)
sumstats$chr_orig <- as.character(sumstats$chr)

# Estandarizar chrX a 23 para uso interno
sumstats$chr <- as.character(sumstats$chr)
sumstats$chr[sumstats$chr %in% c("X", "chrX")] <- "23"
sumstats$chr <- as.integer(sumstats$chr)

cat(sprintf("  GWAS tras QC: %s variantes (incluyendo %s en chrX)\n",
    format(nrow(sumstats), big.mark=","),
    format(sum(sumstats$chr == 23), big.mark=",")))

# ================================================================
# [2] Creación del Filtro Propio (LD-Pruning en MCPS)
# ================================================================
cat("\n[2] Creando filtro propio (LD-Pruning r²<0.9 en datos MCPS)...\n")

gwas_ids_file <- file.path(tmp_dir, "gwas_ids_for_pruning.txt")

# FIX v2: Para autosomas, IDs son chr1:pos:a1:a0
# Para chrX, IDs son chrX:pos:a1:a0 (NO chr23)
is_autosome <- sumstats$chr != 23
is_chrX     <- sumstats$chr == 23

# Autosomas: chr{N}:pos:a1:a0 y chr{N}:pos:a0:a1
auto_ids_fwd <- paste0("chr", sumstats$chr[is_autosome], ":",
                        sumstats$pos[is_autosome], ":",
                        sumstats$a1[is_autosome], ":",
                        sumstats$a0[is_autosome])
auto_ids_rev <- paste0("chr", sumstats$chr[is_autosome], ":",
                        sumstats$pos[is_autosome], ":",
                        sumstats$a0[is_autosome], ":",
                        sumstats$a1[is_autosome])

# chrX: chrX:pos:a1:a0 y chrX:pos:a0:a1
chrx_ids_fwd <- paste0("chrX:", sumstats$pos[is_chrX], ":",
                        sumstats$a1[is_chrX], ":",
                        sumstats$a0[is_chrX])
chrx_ids_rev <- paste0("chrX:", sumstats$pos[is_chrX], ":",
                        sumstats$a0[is_chrX], ":",
                        sumstats$a1[is_chrX])

all_ids <- unique(c(auto_ids_fwd, auto_ids_rev, chrx_ids_fwd, chrx_ids_rev))
writeLines(all_ids, gwas_ids_file)
cat(sprintf("  IDs candidatos para extract: %s\n", format(length(all_ids), big.mark=",")))

rm(auto_ids_fwd, auto_ids_rev, chrx_ids_fwd, chrx_ids_rev, all_ids); gc()

# LD-Pruning por cromosoma
pruned_ids_all <- c()

for (chr in 1:23) {
    chr_label <- ifelse(chr == 23, "X", as.character(chr))
    prefix <- file.path(plink_dir, sprintf("mcps-freeze150k_qcd_chr%s_ivs", chr_label))

    if (!file.exists(paste0(prefix, ".bed"))) next

    prune_out <- file.path(tmp_dir, sprintf("prune_chr%s", chr_label))

    # FIX v2: r² = 0.9 en vez de 0.2
    # LDpred2 necesita variantes en LD para el modelo Bayesiano
    # Un pruning suave (r²<0.9) remueve solo redundancia extrema
    # y deja ~1-1.5M variantes (ideal para LDpred2)
    system(sprintf(
        "plink2 --bfile %s --extract %s --indep-pairwise 500 50 0.9 --out %s --threads %d --memory 8000 2>/dev/null",
        prefix, gwas_ids_file, prune_out, ncores
    ))

    prune_in <- paste0(prune_out, ".prune.in")
    if (file.exists(prune_in)) {
        ids <- readLines(prune_in)
        pruned_ids_all <- c(pruned_ids_all, ids)
        cat(sprintf("    chr%s: %s variantes seleccionadas\n",
            chr_label, format(length(ids), big.mark=",")))
    } else {
        cat(sprintf("    chr%s: pruning falló o 0 variantes\n", chr_label))
    }
}

cat(sprintf("\n  Total variantes en 'Filtro MCPS': %s\n",
    format(length(pruned_ids_all), big.mark=",")))

# Parsear IDs pruneados para filtrar sumstats
# IDs son chr1:12345:A:G o chrX:12345:A:G
pruned_parts <- strsplit(pruned_ids_all, ":")
pruned_chr_raw <- sapply(pruned_parts, `[`, 1)  # "chr1", "chrX", etc
pruned_chr <- as.integer(gsub("chr", "", gsub("chrX", "chr23", pruned_chr_raw)))
pruned_pos <- as.integer(sapply(pruned_parts, `[`, 2))

pruned_key <- paste(pruned_chr, pruned_pos, sep=":")
gwas_key <- paste(sumstats$chr, sumstats$pos, sep=":")

sumstats <- sumstats[gwas_key %in% pruned_key, ]

n_auto <- sum(sumstats$chr != 23)
n_chrx <- sum(sumstats$chr == 23)
cat(sprintf("  GWAS filtrado: %s variantes (%s autosomas + %s chrX)\n",
    format(nrow(sumstats), big.mark=","),
    format(n_auto, big.mark=","),
    format(n_chrx, big.mark=",")))

# ================================================================
# [3] LDpred2-auto por cromosoma (1-23)
# ================================================================
cat("\n[3] MCMC y Scoring por Cromosoma...\n")

prs_total <- NULL
fam_global <- NULL
matched_total <- 0
results_log <- list()

for (chr in 1:23) {
    chr_label <- ifelse(chr == 23, "X", as.character(chr))
    prefix <- file.path(plink_dir, sprintf("mcps-freeze150k_qcd_chr%s_ivs", chr_label))
    bed_file <- paste0(prefix, ".bed")

    if (!file.exists(bed_file)) next

    sumstats_chr <- sumstats[sumstats$chr == chr, ]
    if (nrow(sumstats_chr) < 50) {
        cat(sprintf("\n  chr%s: solo %d variantes, saltando\n", chr_label, nrow(sumstats_chr)))
        next
    }

    cat(sprintf("\n  >>> chr%s: %s variantes GWAS %s<<<\n",
        chr_label, format(nrow(sumstats_chr), big.mark=","),
        ifelse(chr == 23, "(chrX pseudo-diploide) ", "")))
    t0 <- Sys.time()

    # Leer .bed original por cromosoma
    bk_file <- file.path(tmp_dir, sprintf("chr%s", chr_label))
    snp_readBed(bed_file, backingfile = bk_file)
    obj <- snp_attach(paste0(bk_file, ".rds"))
    G <- obj$genotypes

    if (is.null(prs_total)) {
        prs_total <- rep(0, nrow(G))
        fam_global <- obj$fam
    }

    # Match variantes
    map_df <- data.frame(
        chr = obj$map$chromosome,
        pos = obj$map$physical.pos,
        a1  = obj$map$allele1,
        a0  = obj$map$allele2,
        stringsAsFactors = FALSE
    )
    info_snp <- snp_match(sumstats_chr, map_df, strand_flip=FALSE,
                          join_by_pos=TRUE, match.min.prop=0.01)

    cat(sprintf("      Matched: %s variantes\n", format(nrow(info_snp), big.mark=",")))
    matched_total <- matched_total + nrow(info_snp)

    if (nrow(info_snp) < 50) {
        cat("      Muy pocas, saltando\n")
        rm(obj, G); gc()
        next
    }

    # FIX: Filtrar variantes monomorphic/zero-variance ANTES de snp_cor
    # Estas causan columnas vacías en la sparse matrix → segfault en col_count_sym
    ind_col <- info_snp$`_NUM_ID_`
    sd_vals <- big_scale()(G, ind.col = ind_col)$scale
    keep_var <- which(sd_vals > 1e-6)
    if (length(keep_var) < nrow(info_snp)) {
        cat(sprintf("      Removidas %d variantes con varianza ~0\n",
            nrow(info_snp) - length(keep_var)))
        info_snp <- info_snp[keep_var, ]
        ind_col <- info_snp$`_NUM_ID_`
    }

    if (nrow(info_snp) < 50) {
        cat("      Muy pocas variantes tras filtro varianza, saltando\n")
        rm(obj, G); gc()
        next
    }

    # LD matrix sparse (dsCMatrix)
    cat(sprintf("      Calculando LD matrix (%s variantes)...\n",
        format(length(ind_col), big.mark=",")))
    corr0 <- snp_cor(
        Gna = G, ind.col = ind_col,
        infos.pos = info_snp$pos / 1000,
        size = 3000, ncores = ncores,
        thr_r2 = 0.005
    )

    cat(sprintf("      LD matrix: %s x %s, nnz=%s, clase=%s\n",
        format(nrow(corr0), big.mark=","), format(ncol(corr0), big.mark=","),
        format(length(corr0@x), big.mark=","), class(corr0)))

    # Convertir a SFBM con fix de columnas vacías
    sfbm_file <- file.path(tmp_dir, sprintf("corr_chr%s", chr_label))
    corr <- safe_as_SFBM(corr0, backingfile = sfbm_file, compact = TRUE)
    rm(corr0); gc()

    # LDpred2-auto
    df_beta <- info_snp[, c("beta", "beta_se", "n_eff")]
    multi_auto <- snp_ldpred2_auto(
        corr = corr, df_beta = df_beta, h2_init = 0.01,
        vec_p_init = seq_log(1e-4, 0.5, length.out = 30),
        burn_in = 500, num_iter = 500, report_step = 100,
        ncores = ncores, sparse = TRUE
    )

    # Filtrar chains convergentes
    h2_est <- sapply(multi_auto, function(x) x$h2_est)
    p_est  <- sapply(multi_auto, function(x) x$p_est)
    keep <- which(h2_est > 0.0005 & h2_est < 0.5 & p_est > 1e-6 & p_est < 0.9)
    if (length(keep) < 3) keep <- seq_along(multi_auto)
    beta_final <- rowMeans(sapply(multi_auto[keep], function(x) x$beta_est))

    # Scoring
    prs_chr <- big_prodVec(X=G, y.col=beta_final,
                           ind.col=info_snp$`_NUM_ID_`, ncores=ncores)
    prs_total <- prs_total + prs_chr

    elapsed <- as.numeric(difftime(Sys.time(), t0, units="mins"))
    cat(sprintf("      h2=%.4f, p=%.2e, chains=%d/%d, nonzero=%s, Tiempo=%.1f min\n",
        median(h2_est[keep]), median(p_est[keep]),
        length(keep), length(multi_auto),
        format(sum(beta_final != 0), big.mark=","), elapsed))

    results_log[[chr_label]] <- list(
        n_matched = nrow(info_snp),
        h2 = median(h2_est[keep]),
        p = median(p_est[keep]),
        n_nonzero = sum(beta_final != 0),
        elapsed_min = elapsed
    )

    # Limpieza agresiva
    rm(obj, G, corr, multi_auto, info_snp, beta_final, prs_chr); gc()
}

# ================================================================
# [4] Exportar Resultados
# ================================================================
cat(sprintf("\n[4] Exportando (%s variantes totales)...\n",
    format(matched_total, big.mark=",")))

# Resumen por cromosoma
cat("\n  Resumen por cromosoma:\n")
total_time <- 0
for (cn in names(results_log)) {
    r <- results_log[[cn]]
    cat(sprintf("    chr%s: %s vars, h2=%.4f, p=%.2e, %s nonzero, %.1f min\n",
        cn, format(r$n_matched, big.mark=","), r$h2, r$p,
        format(r$n_nonzero, big.mark=","), r$elapsed_min))
    total_time <- total_time + r$elapsed_min
}
cat(sprintf("  Tiempo total MCMC: %.1f horas\n", total_time / 60))

# Conectar IDs
covar_file <- "/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes/MCPS BASELINE.csv"
link_file  <- "/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes/RGN_LINK_IID.csv"

df_link <- read.csv(link_file)
iid_to_patid <- setNames(as.character(df_link$PATID), as.character(df_link$IID))
df_baseline <- read.csv(covar_file)
df_baseline <- df_baseline[!is.na(df_baseline$MALE), ]
patid_to_sex <- setNames(as.integer(df_baseline$MALE), as.character(df_baseline$PATID))

patids <- iid_to_patid[as.character(fam_global$sample.ID)]
patids[is.na(patids)] <- "UNKNOWN"
sexes <- patid_to_sex[patids]
sexes[is.na(sexes)] <- -1L

df_out <- data.frame(PATID=patids, IID=fam_global$sample.ID, SEX=sexes,
                     PRS_TOTAL=prs_total, stringsAsFactors=FALSE)
df_out <- df_out[df_out$SEX >= 0, ]

# Estandarizar Z-score por sexo
df_out$PRS_Z <- NA_real_
df_out$PERCENTILE <- NA_real_
for (sv in c(0, 1)) {
    mask <- df_out$SEX == sv
    vals <- df_out$PRS_TOTAL[mask]
    if (length(vals) > 0 && sd(vals) > 0) {
        df_out$PRS_Z[mask] <- (vals - mean(vals)) / sd(vals)
        df_out$PERCENTILE[mask] <- rank(vals) / length(vals) * 100
    }
}

df_out <- df_out[order(-df_out$PRS_TOTAL), ]
out_file <- file.path(out_dir, "prs_diabetes_ldpred2_custom_mcps_final.tsv")
write.table(df_out, out_file, sep="\t", row.names=FALSE, quote=FALSE)
saveRDS(results_log, file.path(out_dir, "ldpred2_custom_results_log.rds"))

cat(sprintf("\n Pipeline Custom MCPS completado\n"))
cat(sprintf("   Output: %s\n", out_file))
cat(sprintf("   Muestras: %s\n", format(nrow(df_out), big.mark=",")))
cat(sprintf("   Variantes totales: %s\n", format(matched_total, big.mark=",")))
cat(sprintf("   Cromosomas: 1-22 + X\n"))
REOF

apptainer exec --bind /mnt/cephfs:/mnt/cephfs --bind "$TMPDIR:/tmp,$TMPDIR:/var/tmp" \
    "$SIF" /opt/miniforge3/envs/r-bio/bin/Rscript "$R_SCRIPT"

rm -rf "$TMPDIR"
echo "============================================================"
echo "  PIPELINE C v4 COMPLETADO"
echo "============================================================"
