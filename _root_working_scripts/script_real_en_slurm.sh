#!/bin/bash
#SBATCH --job-name=ldpred2_hm3
#SBATCH --account=researchers
#SBATCH --qos=vip
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=120G
#SBATCH --output=/mnt/cephfs/orgs/home/angel.pacheco/logs/ldpred2_hm3_%j.out
#SBATCH --error=/mnt/cephfs/orgs/home/angel.pacheco/logs/ldpred2_hm3_%j.err
#SBATCH --mail-user=angelpachecosaavedra@gmail.com
#SBATCH --mail-type=BEGIN,END,FAIL

set -euo pipefail

# ================================================================
# PIPELINE B — LDpred2-auto — v8 HapMap3
#
# v7 fixes (vs v5/v6):
#   1. thr_r2=0.005 en snp_cor → previene integer overflow en cumsum
#      interno (la cohorte admixta MCPS genera LD más extenso que EUR,
#      produciendo >2^31 entries sin threshold).
#   2. safe_as_SFBM: manipulación directa de slots @p/@i/@x de la
#      dsCMatrix para rellenar columnas vacías. NO usa aritmética de
#      Matrix (que internamente convierte a dgCMatrix → overflow).
#   3. Filtrado de variantes sd≈0 previo a snp_cor (ya en v5).
#   4. as_SFBM(compact=TRUE) para guardado ultraligero (ya en v5).
# ================================================================

SIF="/mnt/cephfs/biocontainers/images/jupyter-biotools-1.4.sif"

PLINK_DIR="/mnt/cephfs/hot_nvme/mcps/imputed-topmed/plink_files/maximally_unrelated"
GWAS_DIR="/mnt/cephfs/orgs/home/angel.pacheco/gwas_miami"
OUT_DIR="/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes/ldpred2"
mkdir -p "$OUT_DIR"

TMPDIR="/mnt/cephfs/scratch/$USER/ldpred2_hm3_tmp/$SLURM_JOB_ID"
mkdir -p "$TMPDIR"
chmod 700 "$TMPDIR"
export TMPDIR TMP="$TMPDIR" TEMP="$TMPDIR"

echo "============================================================"
echo "  PIPELINE B — LDpred2-auto — v7 HapMap3"
echo "  Nodo: $SLURMD_NODENAME"
echo "  RAM: $(free -g | awk 'NR==2{print $2}')G"
echo "  $(date)"
echo "============================================================"

# ================================================================
# Todo en R — incluyendo descarga HapMap3 y filtrado
# ================================================================
echo ""
echo "--- Ejecutando LDpred2 (HapMap3 filtered) ---"

R_SCRIPT="$TMPDIR/run_ldpred2_hm3.R"
cat << 'REOF' > "$R_SCRIPT"
library(bigsnpr)
library(data.table)
library(methods)

ncores <- 24
options(bigstatsr.check.parallel.blas = FALSE)
options(bigstatsr.ncores.max = 24)

# ================================================================
# safe_as_SFBM: convierte dsCMatrix → SFBM sin aritmética de Matrix
#
# Problema: as_SFBM() segfaults con columnas vacías en dsCMatrix.
#   - Convertir a dgCMatrix causa integer overflow (>2^31 entries).
#   - Sumar Diagonal() internamente pasa por dgCMatrix (mismo overflow).
# Solución: manipular directamente los slots @p/@i/@x de la dsCMatrix
#   para insertar entries diagonales microscópicos (eps) en columnas
#   vacías, sin ninguna operación aritmética del paquete Matrix.
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

    # --- Nuevos column pointers ---
    cum_add <- cumsum(col_nnz == 0L)
    new_p   <- old_p + c(0L, cum_add)

    # --- Nuevos vectores i/x ---
    new_i <- integer(new_nnz)
    new_x <- numeric(new_nnz)

    # Copiar entries en bloques entre columnas vacías consecutivas.
    # Dentro de un bloque, todas las entries se desplazan la misma cantidad.
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

        # Insertar diagonal para la columna vacía
        if (b <= n_add) {
            ec <- empty[b]
            dp <- new_p[ec] + 1L
            new_i[dp] <- ec - 1L   # 0-based row index
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
# [0] Descargar mapa HapMap3
# ================================================================
cat("\n[0] Cargando mapa HapMap3...\n")

hm3_file <- file.path(out_dir, "map_hm3_ldpred2.rds")
if (!file.exists(hm3_file)) {
    hm3_url <- "https://ndownloader.figshare.com/files/25503788"
    download.file(hm3_url, hm3_file, mode = "wb")
}
map_hm3 <- readRDS(hm3_file)
cat(sprintf("  HapMap3 variantes: %s\n", format(nrow(map_hm3), big.mark = ",")))
cat(sprintf("  Columnas: %s\n", paste(names(map_hm3), collapse = ", ")))
cat(sprintf("  Cromosomas: %s\n", paste(sort(unique(map_hm3$chr)), collapse = ", ")))

if ("pos_hg38" %in% names(map_hm3)) {
    cat("  Posiciones hg38 disponibles\n")
    map_hm3$pos_use <- map_hm3$pos_hg38
} else if ("pos" %in% names(map_hm3)) {
    cat(sprintf("  Columnas disponibles: %s\n", paste(names(map_hm3), collapse = ", ")))
    cat(" Usando 'pos' como posición (verificar build)\n")
    map_hm3$pos_use <- map_hm3$pos
}

# ================================================================
# [1] Cargar GWAS y meta-análisis IVW
# ================================================================
cat("\n[1] GWAS Meta-análisis...\n")

gf <- fread(file.path(gwas_dir, "combined-gwas-results-all-ctrl-females.txt.gz"),
            select = c("CHROM","GENPOS","ID","ALLELE0","ALLELE1","A1FREQ","INFO","N","BETA","SE","P"),
            showProgress = FALSE)
gm <- fread(file.path(gwas_dir, "combined-gwas-results-all-ctrl-males.txt.gz"),
            select = c("CHROM","GENPOS","ID","ALLELE0","ALLELE1","A1FREQ","INFO","N","BETA","SE","P"),
            showProgress = FALSE)
gwas <- merge(gf, gm, by = "ID", suffixes = c(".f", ".m"))
rm(gf, gm); gc()

w_f <- 1/gwas$SE.f^2; w_m <- 1/gwas$SE.m^2

# n_eff para trait binario: 4 * cases * controls / (cases + controls)
# Mujeres: 17,300 casos / 54,196 controles
# Hombres:  8,522 casos / 26,733 controles
neff_f <- 4 * 17300 * 54196 / (17300 + 54196)  # = 52,456
neff_m <- 4 *  8522 * 26733 / ( 8522 + 26733)  # = 25,848

# FIX CRÍTICO: REGENIE reporta betas en log-odds (logistic regression).
# LDpred2 asume escala lineal: sd_ss = 1/sqrt(n_eff*SE² + β²) ≈ sd_geno.
# Sin ajuste, sd_ss/sd_geno ≈ sqrt(K*(1-K)) ≈ 0.43 → LDpred2 descarta
# todas las variantes y MCMC diverge (h2=NA).
# Solución: n_eff_ldpred2 = n_eff_binario × K × (1-K)
# donde K = prevalencia = 25822/106751 = 0.2419
K_prev <- (17300 + 8522) / (17300 + 54196 + 8522 + 26733)  # = 0.2419
neff_ldpred2 <- (neff_f + neff_m) * K_prev * (1 - K_prev)   # ≈ 14,359

sumstats <- data.frame(
    chr = gwas$CHROM.f,
    pos = gwas$GENPOS.f,
    a1 = gwas$ALLELE1.f, a0 = gwas$ALLELE0.f,
    freq = (gwas$A1FREQ.f * gwas$N.f + gwas$A1FREQ.m * gwas$N.m) / (gwas$N.f + gwas$N.m),
    info = pmin(gwas$INFO.f, gwas$INFO.m),
    beta = (gwas$BETA.f * w_f + gwas$BETA.m * w_m) / (w_f + w_m),
    beta_se = sqrt(1 / (w_f + w_m)),
    n_eff = neff_ldpred2,
    stringsAsFactors = FALSE
)
sumstats$p <- 2 * pnorm(-abs(sumstats$beta / sumstats$beta_se))
rm(gwas); gc()

sumstats <- sumstats[sumstats$info >= 0.8 & sumstats$freq >= 0.01 & sumstats$freq <= 0.99 &
                     !is.na(sumstats$beta) & !is.na(sumstats$beta_se) & sumstats$beta_se > 0, ]

sumstats$chr <- as.character(sumstats$chr)
sumstats$chr[sumstats$chr %in% c("X", "chrX")] <- "23"
sumstats$chr <- as.integer(sumstats$chr)

cat(sprintf("  GWAS tras QC: %s variantes\n", format(nrow(sumstats), big.mark = ",")))

# ================================================================
# [2] Filtrar GWAS a HapMap3 por posición + alelos
# ================================================================
cat("\n[2] Filtrando a HapMap3...\n")

pos_col <- NULL
for (cn in c("pos_hg38", "pos", "pos_hg19")) {
    if (cn %in% names(map_hm3)) {
        key_hm3 <- paste(map_hm3$chr, map_hm3[[cn]], sep = ":")
        key_gwas <- paste(sumstats$chr, sumstats$pos, sep = ":")
        n <- sum(key_gwas %in% key_hm3)
        if (n > 100000) {
            pos_col <- cn
            cat(sprintf("  Usando columna '%s': %s matches con GWAS\n", cn, format(n, big.mark = ",")))
            break
        }
    }
}

if (is.null(pos_col)) {
    cat("\n  ⚠️o se encontró match posicional HapMap3↔GWAS (hg19 vs hg38)\n")
    cat("  Fallback: LD-pruning del GWAS a ~1M variantes con PLINK\n\n")

    gwas_ids_file <- file.path(tmp_dir, "gwas_ids_for_pruning.txt")
    gwas_extract <- paste0("chr", sumstats$chr, ":", sumstats$pos, ":",
                           sumstats$a1, ":", sumstats$a0)
    gwas_extract2 <- paste0("chr", sumstats$chr, ":", sumstats$pos, ":",
                            sumstats$a0, ":", sumstats$a1)
    writeLines(c(gwas_extract, gwas_extract2), gwas_ids_file)

    cat("  LD-pruning por cromosoma (r2 < 0.2, ventana 1Mb)...\n")
    pruned_ids_all <- c()
    for (chr in c(1:22, 23)) {
        chr_label <- ifelse(chr == 23, "X", as.character(chr))
        prefix <- file.path(plink_dir, sprintf("mcps-freeze150k_qcd_chr%s_ivs", chr_label))
        if (!file.exists(paste0(prefix, ".bed"))) next

        prune_out <- file.path(tmp_dir, sprintf("prune_chr%d", chr))
        system(sprintf(
            "plink2 --bfile %s --extract %s --indep-pairwise 1000 100 0.2 --out %s --threads 4 --memory 4000 2>/dev/null",
            prefix, gwas_ids_file, prune_out
        ))

        prune_in <- paste0(prune_out, ".prune.in")
        if (file.exists(prune_in)) {
            ids <- readLines(prune_in)
            pruned_ids_all <- c(pruned_ids_all, ids)
            cat(sprintf("    chr%s: %s variantes tras pruning\n", chr_label, format(length(ids), big.mark = ",")))
        }
    }
    cat(sprintf("  Total tras LD-pruning: %s variantes\n", format(length(pruned_ids_all), big.mark = ",")))

    pruned_parts <- strsplit(pruned_ids_all, ":")
    pruned_chr <- as.integer(gsub("chr", "", sapply(pruned_parts, `[`, 1)))
    pruned_chr[is.na(pruned_chr)] <- 23
    pruned_pos <- as.integer(sapply(pruned_parts, `[`, 2))
    pruned_key <- paste(pruned_chr, pruned_pos, sep = ":")
    gwas_key <- paste(sumstats$chr, sumstats$pos, sep = ":")

    sumstats <- sumstats[gwas_key %in% pruned_key, ]
    cat(sprintf("  GWAS filtrado: %s variantes\n", format(nrow(sumstats), big.mark = ",")))

} else {
    key_hm3 <- paste(map_hm3$chr, map_hm3[[pos_col]], sep = ":")
    key_gwas <- paste(sumstats$chr, sumstats$pos, sep = ":")
    sumstats <- sumstats[key_gwas %in% key_hm3, ]
    cat(sprintf("  GWAS filtrado a HapMap3: %s variantes\n", format(nrow(sumstats), big.mark = ",")))
}

cat(sprintf("  Variantes por cromosoma (muestra): chr1=%d, chr22=%d\n", sum(sumstats$chr == 1), sum(sumstats$chr == 22)))

# ================================================================
# [3] LDpred2-auto por cromosoma
# ================================================================
cat("\n[3] LDpred2-auto por cromosoma...\n")

prs_total <- NULL
fam_global <- NULL
matched_total <- 0
results_log <- list()

for (chr in c(1:22, 23)) {
    chr_label <- ifelse(chr == 23, "X", as.character(chr))
    prefix <- file.path(plink_dir, sprintf("mcps-freeze150k_qcd_chr%s_ivs", chr_label))
    bed_file <- paste0(prefix, ".bed")

    if (!file.exists(bed_file)) next

    sumstats_chr <- sumstats[sumstats$chr == chr, ]
    if (nrow(sumstats_chr) < 100) {
        cat(sprintf("  chr%s: solo %d variantes, saltando\n", chr_label, nrow(sumstats_chr)))
        next
    }

    cat(sprintf("\n  >>> chr%s: %s variantes GWAS <<<\n", chr_label, format(nrow(sumstats_chr), big.mark = ",")))
    t0 <- Sys.time()

    bk_file <- file.path(tmp_dir, sprintf("chr%s", chr_label))
    snp_readBed(bed_file, backingfile = bk_file)
    obj <- snp_attach(paste0(bk_file, ".rds"))
    G <- obj$genotypes

    if (is.null(prs_total)) {
        prs_total <- rep(0, nrow(G))
        fam_global <- obj$fam
    }

    map_df <- data.frame(
        chr = obj$map$chromosome, pos = obj$map$physical.pos,
        a1 = obj$map$allele1, a0 = obj$map$allele2,
        stringsAsFactors = FALSE
    )
    info_snp <- snp_match(sumstats_chr, map_df, strand_flip = FALSE, join_by_pos = TRUE, match.min.prop = 0.01)

    cat(sprintf("      Matched: %s variantes\n", format(nrow(info_snp), big.mark = ",")))
    matched_total <- matched_total + nrow(info_snp)

    if (nrow(info_snp) < 50) {
        cat("      Muy pocas variantes, saltando\n")
        rm(obj, G); gc()
        next
    }

    # FIX CRÍTICO FINAL: Filtrar monomórficas + Guardado Compacto
    ind_col <- info_snp$`_NUM_ID_`
    sd_vals <- big_scale()(G, ind.col = ind_col)$scale
    keep_var <- which(sd_vals > 1e-6)
    
    if (length(keep_var) < nrow(info_snp)) {
        cat(sprintf("      Removidas %d variantes con varianza ~0\n", nrow(info_snp) - length(keep_var)))
        info_snp <- info_snp[keep_var, ]
        ind_col <- info_snp$`_NUM_ID_`
    }
    
    if (nrow(info_snp) < 50) {
        cat("      Muy pocas variantes tras filtro varianza, saltando\n")
        rm(obj, G); gc()
        next
    }
    
    cat(sprintf("      Calculando LD matrix (%s variantes)...\n", format(length(ind_col), big.mark=",")))
    
    corr0 <- snp_cor(
        Gna = G, ind.col = ind_col,
        infos.pos = info_snp$pos / 1000,
        size = 3000, ncores = ncores,
        thr_r2 = 0.005
    )
    
    cat(sprintf("      LD matrix: %s x %s, nnz=%s, clase=%s\n",
        format(nrow(corr0), big.mark=","), format(ncol(corr0), big.mark=","),
        format(length(corr0@x), big.mark=","), class(corr0)))

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

    h2_est <- sapply(multi_auto, function(x) x$h2_est)
    p_est  <- sapply(multi_auto, function(x) x$p_est)
    keep <- which(h2_est > 0.0005 & h2_est < 0.5 & p_est > 1e-6 & p_est < 0.9)
    if (length(keep) < 3) keep <- seq_along(multi_auto)
    beta_final <- rowMeans(sapply(multi_auto[keep], function(x) x$beta_est))

    cat(sprintf("      h2=%.4f, p=%.2e, chains=%d/%d, nonzero=%s\n",
        median(h2_est[keep]), median(p_est[keep]), length(keep),
        length(multi_auto), format(sum(beta_final != 0), big.mark = ",")))

    # Scoring
    prs_chr <- big_prodVec(X = G, y.col = beta_final, ind.col = info_snp$`_NUM_ID_`, ncores = ncores)
    prs_total <- prs_total + prs_chr

    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    cat(sprintf("      Tiempo: %.1f min\n", elapsed))

    results_log[[chr_label]] <- list(
        n_matched = nrow(info_snp), h2 = median(h2_est[keep]),
        p = median(p_est[keep]), elapsed_min = elapsed
    )

    rm(obj, G, corr, multi_auto, info_snp, beta_final, prs_chr); gc()
}

# ================================================================
# [4] Exportar
# ================================================================
cat(sprintf("\n[4] Exportando (%s variantes totales)...\n", format(matched_total, big.mark = ",")))

cat("\n  Resumen por cromosoma:\n")
for (cn in names(results_log)) {
    r <- results_log[[cn]]
    cat(sprintf("    chr%s: %s vars, h2=%.4f, p=%.2e, %.1f min\n",
        cn, format(r$n_matched, big.mark = ","), r$h2, r$p, r$elapsed_min))
}

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

df_out <- data.frame(PATID = patids, IID = fam_global$sample.ID,
                     SEX = sexes, PRS_TOTAL = prs_total, stringsAsFactors = FALSE)
df_out <- df_out[df_out$SEX >= 0, ]

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
out_file <- file.path(out_dir, "prs_diabetes_ldpred2_final.tsv")
write.table(df_out, out_file, sep = "\t", row.names = FALSE, quote = FALSE)
saveRDS(results_log, file.path(out_dir, "ldpred2_hm3_results_log.rds"))

cat(sprintf("\n Pipeline completado\n"))
cat(sprintf("   Output: %s\n", out_file))
cat(sprintf("   Muestras: %s\n", format(nrow(df_out), big.mark = ",")))
cat(sprintf("   Variantes: %s\n", format(matched_total, big.mark = ",")))
REOF

apptainer exec --bind /mnt/cephfs:/mnt/cephfs --bind "$TMPDIR:/tmp,$TMPDIR:/var/tmp" \
    "$SIF" /opt/miniforge3/envs/r-bio/bin/Rscript "$R_SCRIPT"

rm -rf "$TMPDIR"

echo ""
echo "============================================================"
echo "  PIPELINE B v7 HapMap3 COMPLETADO"
echo "============================================================"
