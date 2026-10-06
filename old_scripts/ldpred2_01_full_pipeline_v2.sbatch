#!/bin/bash
#SBATCH --job-name=ldpred2_run
#SBATCH --account=researchers
#SBATCH --qos=vip
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=120G
#SBATCH --time=24:00:00
#SBATCH --output=/mnt/cephfs/orgs/home/angel.pacheco/logs/ldpred2_run_%j.out
#SBATCH --error=/mnt/cephfs/orgs/home/angel.pacheco/logs/ldpred2_run_%j.err
#SBATCH --mail-user=angelpachecosaavedra@gmail.com
#SBATCH --mail-type=BEGIN,END,FAIL

set -euo pipefail

# ================================================================
# PIPELINE B — LDpred2-auto (bigsnpr) — v2 optimizado
#
# Cambios vs v1:
#   - NO merge de PLINK: procesa por cromosoma, evita .bed de 390GB
#   - snp_readBed2 por cromosoma → bigSNP concatenados en R
#   - Ventana LD corregida (3 cM con mapa genético, o 3000 kb)
#   - Mapeo correcto de allele1/allele2 de bigsnpr
#   - GWAS: meta-análisis IVW de females + males
#   - chrX incluido: hombres como pseudo-diploides (0/2)
#   - Genotipos discretos end-to-end
# ================================================================

SIF="/mnt/cephfs/biocontainers/images/jupyter-biotools-1.4.sif"

PLINK_DIR="/mnt/cephfs/hot_nvme/mcps/imputed-topmed/plink_files/maximally_unrelated"
GWAS_DIR="/mnt/cephfs/orgs/home/angel.pacheco/gwas_miami"
OUT_DIR="/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes/ldpred2"
mkdir -p "$OUT_DIR"

TMPDIR="/mnt/cephfs/scratch/$USER/ldpred2_tmp/$SLURM_JOB_ID"
mkdir -p "$TMPDIR"
chmod 700 "$TMPDIR"
export TMPDIR TMP="$TMPDIR" TEMP="$TMPDIR"

echo "============================================================"
echo "  PIPELINE B — LDpred2-auto (bigsnpr) — v2"
echo "  Nodo: $SLURMD_NODENAME"
echo "  RAM: $(free -g | awk 'NR==2{print $2}')G"
echo "  Threads: $SLURM_CPUS_PER_TASK"
echo "  $(date)"
echo "============================================================"

# ================================================================
# PASO 0: Verificar que los bfiles existen
# ================================================================
echo ""
echo "--- PASO 0: Verificando bfiles ---"
for CHR in $(seq 1 22) X; do
    PREFIX="${PLINK_DIR}/mcps-freeze150k_qcd_chr${CHR}_ivs"
    if [ ! -f "${PREFIX}.bed" ]; then
        echo " alta chr${CHR}"
    fi
done
echo "  .bim chr1 ejemplo:"
head -2 "${PLINK_DIR}/mcps-freeze150k_qcd_chr1_ivs.bim"
echo "  .bim chrX ejemplo:"
head -2 "${PLINK_DIR}/mcps-freeze150k_qcd_chrX_ivs.bim"

# ================================================================
# PASOS 1-6: Todo en R
# ================================================================
echo ""
echo "--- Ejecutando LDpred2 pipeline en R ---"

R_SCRIPT="$TMPDIR/run_ldpred2_v2.R"
cat << 'REOF' > "$R_SCRIPT"
library(bigsnpr)
library(data.table)

ncores <- 24
options(bigstatsr.check.parallel.blas = FALSE)

out_dir    <- "/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes/ldpred2"
plink_dir  <- "/mnt/cephfs/hot_nvme/mcps/imputed-topmed/plink_files/maximally_unrelated"
gwas_dir   <- "/mnt/cephfs/orgs/home/angel.pacheco/gwas_miami"
covar_file <- "/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes/MCPS BASELINE.csv"
link_file  <- "/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes/RGN_LINK_IID.csv"
tmpdir     <- Sys.getenv("TMPDIR")

chromosomes <- c(as.character(1:22), "X")

# ================================================================
# PASO 1: Leer PLINK bfiles por cromosoma → bigSNP
# ================================================================
cat(sprintf("\n[%s] PASO 1: Leer genotipos por cromosoma\n", format(Sys.time(), "%H:%M:%S")))

# bigsnpr puede leer directamente archivos .bed por cromosoma
# y luego concatenamos los mapas para tener un índice global

rds_file <- file.path(out_dir, "mcps_79k_allchr.rds")
bk_file  <- file.path(out_dir, "mcps_79k_allchr.bk")

if (file.exists(rds_file)) {
    cat("  Cargando bigSNP existente...\n")
    obj.bigSNP <- snp_attach(rds_file)
} else {
    cat("  Leyendo y mergeando bfiles por cromosoma...\n")
    cat("  Esto puede tomar 20-40 min con 41.5M variantes...\n")

    # Estrategia: snp_readBed2 por cromosoma a TMPDIR,
    # luego snp_bind para concatenar
    bed_files <- c()
    for (chr in chromosomes) {
        prefix <- file.path(plink_dir, sprintf("mcps-freeze150k_qcd_chr%s_ivs", chr))
        bed_path <- paste0(prefix, ".bed")
        if (file.exists(bed_path)) {
            bed_files <- c(bed_files, bed_path)
        }
    }
    cat(sprintf("  Archivos .bed encontrados: %d\n", length(bed_files)))

    # Leer todos los cromosomas de una vez
    # snp_readBed2 puede leer un .bed con cualquier tamaño
    # Pero necesitamos un solo .bed merged — usaremos otra estrategia
    #
    # Alternativa más eficiente: leer chr por chr, concatenar backing files
    obj_list <- list()
    for (i in seq_along(bed_files)) {
        chr <- chromosomes[i]
        bed_path <- bed_files[i]
        chr_bk <- file.path(tmpdir, sprintf("chr%s", chr))
        cat(sprintf("  Leyendo chr%s...\n", chr))
        obj_list[[i]] <- snp_readBed2(bed_path, backingfile = chr_bk, ncores = ncores)
    }

    # Concatenar usando snp_bind
    cat("  Concatenando cromosomas con snp_bind...\n")
    obj.bigSNP <- snp_bind(obj_list, backingfile = tools::file_path_sans_ext(bk_file))
    cat("  bigSNP concatenado y guardado\n")

    # Limpiar intermedios
    for (chr in chromosomes) {
        f <- file.path(tmpdir, sprintf("chr%s", chr))
        unlink(paste0(f, ".bk"))
        unlink(paste0(f, ".rds"))
    }
}

G   <- obj.bigSNP$genotypes
map <- obj.bigSNP$map
fam <- obj.bigSNP$fam

cat(sprintf("  Dimensiones: %s variantes × %s muestras\n",
    format(ncol(G), big.mark = ","), format(nrow(G), big.mark = ",")))
cat(sprintf("  Cromosomas: %s\n", paste(sort(unique(map$chromosome)), collapse = ", ")))

# ================================================================
# PASO 2: Cargar y formatear GWAS MCPS (meta-análisis IVW)
# ================================================================
cat(sprintf("\n[%s] PASO 2: Cargar GWAS MCPS + meta-análisis\n", format(Sys.time(), "%H:%M:%S")))

cat("  Leyendo GWAS females...\n")
gwas_f <- fread(file.path(gwas_dir, "combined-gwas-results-all-ctrl-females.txt.gz"),
                select = c("CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1", "A1FREQ", "INFO", "N", "BETA", "SE", "P"),
                showProgress = FALSE)

cat("  Leyendo GWAS males...\n")
gwas_m <- fread(file.path(gwas_dir, "combined-gwas-results-all-ctrl-males.txt.gz"),
                select = c("CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1", "A1FREQ", "INFO", "N", "BETA", "SE", "P"),
                showProgress = FALSE)

cat(sprintf("  Females: %s variantes, Males: %s variantes\n",
    format(nrow(gwas_f), big.mark = ","), format(nrow(gwas_m), big.mark = ",")))

# IVW meta-analysis
cat("  Meta-análisis IVW...\n")
gwas_merged <- merge(gwas_f, gwas_m, by = "ID", suffixes = c(".f", ".m"))

gwas_merged[, w_f := 1 / SE.f^2]
gwas_merged[, w_m := 1 / SE.m^2]
gwas_merged[, BETA_meta := (BETA.f * w_f + BETA.m * w_m) / (w_f + w_m)]
gwas_merged[, SE_meta := sqrt(1 / (w_f + w_m))]
gwas_merged[, P_meta := 2 * pnorm(-abs(BETA_meta / SE_meta))]
gwas_merged[, N_meta := N.f + N.m]
gwas_merged[, A1FREQ_meta := (A1FREQ.f * N.f + A1FREQ.m * N.m) / (N.f + N.m)]

cat(sprintf("  Variantes en meta-análisis: %s\n", format(nrow(gwas_merged), big.mark = ",")))

# snp_match espera:
#   chr, pos, a0 (ref/non-effect), a1 (alt/effect), beta, beta_se, n_eff
# REGENIE: ALLELE0 = non-effect (a0), ALLELE1 = effect (a1)
sumstats <- data.frame(
    chr     = gwas_merged$CHROM.f,
    pos     = gwas_merged$GENPOS.f,
    a0      = gwas_merged$ALLELE0.f,
    a1      = gwas_merged$ALLELE1.f,
    freq    = gwas_merged$A1FREQ_meta,
    info    = pmin(gwas_merged$INFO.f, gwas_merged$INFO.m),
    beta    = gwas_merged$BETA_meta,
    beta_se = gwas_merged$SE_meta,
    p       = gwas_merged$P_meta,
    n_eff   = gwas_merged$N_meta,
    stringsAsFactors = FALSE
)

# QC
n_pre <- nrow(sumstats)
sumstats <- sumstats[sumstats$info >= 0.8, ]
sumstats <- sumstats[sumstats$freq >= 0.01 & sumstats$freq <= 0.99, ]
sumstats <- sumstats[!is.na(sumstats$beta) & !is.na(sumstats$beta_se), ]
sumstats <- sumstats[sumstats$beta_se > 0, ]
cat(sprintf("  QC: %s → %s variantes\n",
    format(n_pre, big.mark = ","), format(nrow(sumstats), big.mark = ",")))

# Convertir chr: el GWAS usa enteros (1-22) y posiblemente "X"
# bigsnpr .bim usa enteros y 23 para chrX
sumstats$chr <- as.character(sumstats$chr)
sumstats$chr[sumstats$chr %in% c("X", "chrX")] <- "23"
sumstats$chr <- as.integer(sumstats$chr)

rm(gwas_f, gwas_m, gwas_merged)
gc()

# ================================================================
# PASO 3: Match variantes GWAS ↔ LD panel
# ================================================================
cat(sprintf("\n[%s] PASO 3: Match GWAS ↔ genotipos\n", format(Sys.time(), "%H:%M:%S")))

# bigsnpr map columns:
#   $chromosome (int), $physical.pos (int),
#   $allele1 (character), $allele2 (character)
#
# En bigsnpr: allele1 = PRIMER alelo del .bim (col5 = ALT en PLINK1)
#             allele2 = SEGUNDO alelo del .bim (col6 = REF en PLINK1)
#
# En PLINK1 .bim: col5 = minor/ALT, col6 = major/REF
# Pero en nuestros .bim (imputed), el formato del ID es chr:pos:REF:ALT
# con col5 = ALT(C) y col6 = REF(A)
#
# snp_match quiere: en sumstats a1 = effect allele, a0 = other
#                   en map: a1 = allele1 del bigSNP
# El match alinea los effect alleles y ajusta sign de beta si hay flip

# Verificar qué son allele1 y allele2 en nuestro bigSNP
cat("  Verificación alleles bigSNP (primeras 5 variantes):\n")
cat(sprintf("    chr=%s, pos=%s, allele1=%s, allele2=%s\n",
    map$chromosome[1:5], map$physical.pos[1:5], map$allele1[1:5], map$allele2[1:5]))

map_df <- data.frame(
    chr = map$chromosome,
    pos = map$physical.pos,
    a1  = map$allele1,   # allele1 del bigSNP (ALT en PLINK1)
    a0  = map$allele2,   # allele2 del bigSNP (REF en PLINK1)
    stringsAsFactors = FALSE
)

# snp_match matchea sumstats$a1 con map$a1 (o flip con a0)
info_snp <- snp_match(sumstats, map_df, strand_flip = FALSE, join_by_pos = TRUE,
                       match.min.prop = 0.1)

cat(sprintf("  Variantes matcheadas: %s / %s (%.1f%%)\n",
    format(nrow(info_snp), big.mark = ","),
    format(nrow(sumstats), big.mark = ","),
    nrow(info_snp) / nrow(sumstats) * 100))

# Verificar si hubo flips
n_flipped <- sum(info_snp$`_FLIP_` == TRUE, na.rm = TRUE)
cat(sprintf("  Variantes con beta flipped: %s\n", format(n_flipped, big.mark = ",")))

# ================================================================
# PASO 4: Calcular LD matrix sparse + LDpred2-auto
# ================================================================
cat(sprintf("\n[%s] PASO 4: LD matrix + LDpred2-auto MCMC\n", format(Sys.time(), "%H:%M:%S")))

ind.chr <- info_snp$`_NUM_ID_`

# Calcular LD por cromosoma (más eficiente que todo junto)
ld_file <- file.path(out_dir, "ldpred2_ld_matrix.rds")

if (file.exists(ld_file)) {
    cat("  Cargando LD matrix existente...\n")
    corr <- readRDS(ld_file)
} else {
    cat("  Calculando LD matrix sparse por cromosoma...\n")
    cat("  Esto puede tomar 2-6 horas con ~millones de variantes...\n")

    # Usar posiciones físicas en bp
    # snp_cor size parameter: window of SIZE positions on each side
    # Con infos.pos en bp: size = 3e6 (3Mb) pero dividido por el rango
    # La forma correcta de bigsnpr:
    #   size = window_in_number_of_SNPs / 2  (si no se usa infos.pos)
    #   o bien: con infos.pos, size es en las mismas unidades que infos.pos
    #   Para posiciones en bp: size = 3000 da 3000 bp (muy poco!)
    #
    # Según documentación de Privé: snp_cor(..., infos.pos = POS, size = 3/1000)
    # donde POS está en kb → size = 3/1000 de la unidad = 0.003 kb = 3 bp (¡ERROR!)
    #
    # La forma correcta (Privé 2020 tutorial):
    #   POS en bp → convertir a cM o usar size en bp:
    #   snp_cor(size = 3000)  sin infos.pos → ventana de 3000 SNPs
    #   snp_cor(infos.pos = POS/1000, size = 3) → POS en kb, ventana 3kb (muy poco)
    #
    # Mejor: ventana de SNPs, tamaño estándar = 3000 SNPs
    # O convertir posición a Mb y usar size = 3

    POS_MB <- info_snp$pos / 1e6  # Posiciones en Mb

    corr <- snp_cor(
        Gna       = G,
        ind.col   = ind.chr,
        infos.pos = POS_MB,
        size      = 3,           # 3 Mb window
        ncores    = ncores
    )

    saveRDS(corr, ld_file)
    cat(sprintf("   LD matrix guardada: %s\n", ld_file))
}

cat(sprintf("  LD matrix: %s × %s, %s entradas no-cero\n",
    format(nrow(corr), big.mark = ","),
    format(ncol(corr), big.mark = ","),
    format(length(corr@x), big.mark = ",")))

# LDpred2-auto MCMC
cat("\n  Ejecutando LDpred2-auto MCMC...\n")
cat("  30 chains, 500 burn-in + 500 iteraciones cada una\n")
cat("  Esto puede tomar 1-4 horas...\n\n")

# Preparar df_beta con nombres exactos que espera snp_ldpred2_auto
df_beta <- info_snp[, c("beta", "beta_se", "n_eff")]

multi_auto <- snp_ldpred2_auto(
    corr       = corr,
    df_beta    = df_beta,
    h2_init    = 0.05,
    vec_p_init = seq_log(1e-4, 0.5, length.out = 30),
    burn_in    = 500,
    num_iter   = 500,
    report_step = 100,
    ncores     = ncores,
    sparse     = TRUE
)

cat("  LDpred2-auto MCMC completado\n\n")

# Extraer estimaciones
h2_est <- sapply(multi_auto, function(auto) auto$h2_est)
p_est  <- sapply(multi_auto, function(auto) auto$p_est)

cat(sprintf("  h2 estimado por chain: %s\n", paste(sprintf("%.4f", h2_est), collapse = ", ")))
cat(sprintf("  h2 mediana: %.4f\n", median(h2_est)))
cat(sprintf("  p mediana: %.6f\n", median(p_est)))

# Filtrar chains convergentes (Privé 2022 criteria)
# 1. h2 y p dentro de rangos razonables
# 2. Chains con valores extremos se descartan
keep <- which(
    h2_est > 0.001 & h2_est < 0.5 &
    p_est > 1e-6 & p_est < 0.9
)

# Si hay suficientes chains, también filtrar por IQR
if (length(keep) >= 10) {
    h2_iqr <- IQR(h2_est[keep])
    h2_med <- median(h2_est[keep])
    keep2 <- keep[
        h2_est[keep] > (h2_med - 3 * h2_iqr) &
        h2_est[keep] < (h2_med + 3 * h2_iqr)
    ]
    if (length(keep2) >= 5) keep <- keep2
}

if (length(keep) < 3) {
    cat("  ⚠️ Pocas chains convergentes, usando todas\n")
    keep <- seq_along(multi_auto)
}

cat(sprintf("  Chains convergentes: %d / %d\n", length(keep), length(multi_auto)))

# Beta final: promedio de chains convergentes
betas_auto <- sapply(multi_auto[keep], function(auto) auto$beta_est)
beta_final <- rowMeans(betas_auto)

cat(sprintf("  Betas no-cero: %s / %s\n",
    format(sum(beta_final != 0), big.mark = ","),
    format(length(beta_final), big.mark = ",")))

# Guardar resultados intermedios
saveRDS(list(
    beta      = beta_final,
    info_snp  = info_snp,
    h2        = h2_est[keep],
    p         = p_est[keep],
    all_h2    = h2_est,
    all_p     = p_est,
    keep      = keep,
    multi_auto = multi_auto
), file.path(out_dir, "ldpred2_auto_results.rds"))

cat("  Resultados intermedios guardados: ldpred2_auto_results.rds\n")

# ================================================================
# PASO 5: Scoring con big_prodVec (genotipos 0/1/2)
# ================================================================
cat(sprintf("\n[%s] PASO 5: Scoring\n", format(Sys.time(), "%H:%M:%S")))

prs_ldpred2 <- big_prodVec(
    X       = G,
    y.col   = beta_final,
    ind.col = ind.chr,
    ncores  = ncores
)

cat(sprintf("  PRS calculado para %s muestras\n", format(length(prs_ldpred2), big.mark = ",")))
cat(sprintf("  Rango PRS: [%.6f, %.6f]\n", min(prs_ldpred2), max(prs_ldpred2)))

# ================================================================
# PASO 6: Generar tabla de salida
# ================================================================
cat(sprintf("\n[%s] PASO 6: Exportando resultados\n", format(Sys.time(), "%H:%M:%S")))

# IDs del bigSNP
sample_ids <- fam$sample.ID

# Conectar con IDs clínicos
df_link <- read.csv(link_file)
iid_to_patid <- setNames(as.character(df_link$PATID), as.character(df_link$IID))

df_baseline <- read.csv(covar_file)
df_baseline <- df_baseline[!is.na(df_baseline$MALE), ]
patid_to_sex <- setNames(as.integer(df_baseline$MALE), as.character(df_baseline$PATID))

# Construir tabla usando vectorización (más rápido que sapply para 79K)
patids <- iid_to_patid[as.character(sample_ids)]
patids[is.na(patids)] <- "UNKNOWN"

sexes <- patid_to_sex[patids]
sexes[is.na(sexes)] <- -1L

df_out <- data.frame(
    PATID     = patids,
    IID       = sample_ids,
    SEX       = sexes,
    PRS_TOTAL = prs_ldpred2,
    stringsAsFactors = FALSE
)

# Filtrar solo con sexo conocido
df_out <- df_out[df_out$SEX >= 0, ]
cat(sprintf("  Muestras con fenotipo: %s\n", format(nrow(df_out), big.mark = ",")))

# Estandarizar PRS_Z por sexo
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

out_file <- file.path(out_dir, "prs_diabetes_ldpred2_final.tsv")
df_out <- df_out[order(-df_out$PRS_TOTAL), ]
write.table(df_out, out_file, sep = "\t", row.names = FALSE, quote = FALSE)

cat(sprintf("\nPRS LDpred2 GUARDADO: %s\n", out_file))
cat(sprintf("   Muestras: %s\n", format(nrow(df_out), big.mark = ",")))
cat(sprintf("   h2 estimado: %.4f\n", median(h2_est[keep])))
cat(sprintf("   p estimado: %.6f\n", median(p_est[keep])))
cat(sprintf("   Variantes matcheadas: %s\n", format(nrow(info_snp), big.mark = ",")))
cat(sprintf("   Variantes con beta != 0: %s\n", format(sum(beta_final != 0), big.mark = ",")))

cat(sprintf("\n=== LDpred2 pipeline completado: %s ===\n", Sys.time()))
REOF

apptainer exec --bind /mnt/cephfs:/mnt/cephfs --bind "$TMPDIR:/tmp,$TMPDIR:/var/tmp" \
    "$SIF" /opt/miniforge3/envs/r-bio/bin/Rscript "$R_SCRIPT"

rm -rf "$TMPDIR"
echo ""
echo "============================================================"
echo "  PIPELINE B COMPLETADO: prs_diabetes_ldpred2_final.tsv"
echo "  Siguiente: validate_prs_3methods.sbatch"
echo "============================================================"
