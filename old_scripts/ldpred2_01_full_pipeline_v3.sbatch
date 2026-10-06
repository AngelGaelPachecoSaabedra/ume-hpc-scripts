#!/bin/bash
#SBATCH --job-name=ldpred2_v3
#SBATCH --account=researchers
#SBATCH --qos=vip
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=120G
#SBATCH --time=48:00:00
#SBATCH --output=/mnt/cephfs/orgs/home/angel.pacheco/logs/ldpred2_v3_%j.out
#SBATCH --error=/mnt/cephfs/orgs/home/angel.pacheco/logs/ldpred2_v3_%j.err
#SBATCH --mail-user=angelpachecosaavedra@gmail.com
#SBATCH --mail-type=BEGIN,END,FAIL

set -euo pipefail

# ================================================================
# PIPELINE B — LDpred2-auto (bigsnpr) — v3
#
# Fix v2→v3: snp_bind no existe en bigsnpr.
# Estrategia correcta:
#   1. Extraer IDs del GWAS que pasan QC
#   2. PLINK merge chr1-22+X filtrando SOLO esas variantes
#      (~15-20M → .bed de ~170GB en disco persistente)
#   3. snp_readBed2 del merged
#   4. snp_match, snp_cor, LDpred2-auto, scoring
#
# Archivos .bk/.rds se cachean en out_dir para re-runs
# ================================================================

SIF="/mnt/cephfs/biocontainers/images/jupyter-biotools-1.4.sif"

PLINK_DIR="/mnt/cephfs/hot_nvme/mcps/imputed-topmed/plink_files/maximally_unrelated"
GWAS_DIR="/mnt/cephfs/orgs/home/angel.pacheco/gwas_miami"
OUT_DIR="/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes/ldpred2"
mkdir -p "$OUT_DIR"

TMPDIR="/mnt/cephfs/scratch/$USER/ldpred2_v3_tmp/$SLURM_JOB_ID"
mkdir -p "$TMPDIR"
chmod 700 "$TMPDIR"
export TMPDIR TMP="$TMPDIR" TEMP="$TMPDIR"

echo "============================================================"
echo "  PIPELINE B — LDpred2-auto (bigsnpr) — v3"
echo "  Nodo: $SLURMD_NODENAME"
echo "  RAM: $(free -g | awk 'NR==2{print $2}')G"
echo "  Threads: $SLURM_CPUS_PER_TASK"
echo "  $(date)"
echo "============================================================"

# ================================================================
# PASO 0: Extraer IDs de variantes del GWAS para pre-filtro
# ================================================================
MERGED_BED="${OUT_DIR}/mcps_79k_gwas_variants.bed"

if [ -f "$MERGED_BED" ] && [ -f "${OUT_DIR}/mcps_79k_gwas_variants.rds" ]; then
    echo ""
    echo "--- Merged .bed y .rds ya existen, saltando PLINK merge ---"
    N_VARS=$(wc -l < "${OUT_DIR}/mcps_79k_gwas_variants.bim")
    N_SAMPS=$(wc -l < "${OUT_DIR}/mcps_79k_gwas_variants.fam")
    echo "  Cached: ${N_VARS} variantes × ${N_SAMPS} muestras"
else
    echo ""
    echo "--- PASO 0A: Extrayendo IDs de variantes del GWAS ---"

    # Extraer IDs únicos de ambos GWAS (females + males) que pasan INFO >= 0.8
    # ID format en GWAS: "1:598941:G:A" → necesita "chr" prefix para match con .bim
    GWAS_IDS="$TMPDIR/gwas_variant_ids.txt"

    echo "  Extrayendo IDs de females GWAS (INFO >= 0.8)..."
    zcat "${GWAS_DIR}/combined-gwas-results-all-ctrl-females.txt.gz" | \
        awk -F'\t' 'NR > 1 && $7 >= 0.8 {print "chr"$3}' > "$TMPDIR/ids_f.txt"

    echo "  Extrayendo IDs de males GWAS (INFO >= 0.8)..."
    zcat "${GWAS_DIR}/combined-gwas-results-all-ctrl-males.txt.gz" | \
        awk -F'\t' 'NR > 1 && $7 >= 0.8 {print "chr"$3}' > "$TMPDIR/ids_m.txt"

    # Intersección: solo variantes presentes en AMBOS sexos (para meta-análisis)
    echo "  Intersectando..."
    sort -S 8G --parallel=8 "$TMPDIR/ids_f.txt" > "$TMPDIR/ids_f_sorted.txt"
    sort -S 8G --parallel=8 "$TMPDIR/ids_m.txt" > "$TMPDIR/ids_m_sorted.txt"
    comm -12 "$TMPDIR/ids_f_sorted.txt" "$TMPDIR/ids_m_sorted.txt" > "$GWAS_IDS"
    rm -f "$TMPDIR/ids_f.txt" "$TMPDIR/ids_m.txt" "$TMPDIR/ids_f_sorted.txt" "$TMPDIR/ids_m_sorted.txt"

    N_GWAS=$(wc -l < "$GWAS_IDS")
    echo "  IDs GWAS con INFO >= 0.8 en ambos sexos: ${N_GWAS}"
    echo "  Ejemplos: $(head -3 "$GWAS_IDS" | tr '\n' ' ')"

    # ================================================================
    # PASO 0B: PLINK merge chr1-22+X filtrando por IDs GWAS
    # ================================================================
    echo ""
    echo "--- PASO 0B: PLINK merge con pre-filtro ---"

    # Extraer variantes por cromosoma y mergear
    FILTERED_DIR="$TMPDIR/filtered_chr"
    mkdir -p "$FILTERED_DIR"

    FIRST_CHR=""
    MERGE_LIST="$TMPDIR/merge_list.txt"
    > "$MERGE_LIST"

    for CHR in $(seq 1 22) X; do
        PREFIX="${PLINK_DIR}/mcps-freeze150k_qcd_chr${CHR}_ivs"
        if [ ! -f "${PREFIX}.bed" ]; then continue; fi

        FILT_PREFIX="${FILTERED_DIR}/chr${CHR}_filt"

        # Filtrar variantes que están en el GWAS
        apptainer exec --bind /mnt/cephfs:/mnt/cephfs --bind "$TMPDIR:/tmp" "$SIF" \
            plink2 --bfile "$PREFIX" \
            --extract "$GWAS_IDS" \
            --make-bed --out "$FILT_PREFIX" \
            --threads 4 --memory 8000 2>/dev/null

        if [ -f "${FILT_PREFIX}.bed" ]; then
            N_CHR=$(wc -l < "${FILT_PREFIX}.bim")
            echo "  chr${CHR}: ${N_CHR} variantes tras filtro"

            if [ -z "$FIRST_CHR" ]; then
                FIRST_CHR="$FILT_PREFIX"
            else
                echo "$FILT_PREFIX" >> "$MERGE_LIST"
            fi
        fi
    done

    echo ""
    echo "  Mergeando cromosomas filtrados..."

    MERGED_PREFIX="${OUT_DIR}/mcps_79k_gwas_variants"

    if [ -s "$MERGE_LIST" ]; then
        apptainer exec --bind /mnt/cephfs:/mnt/cephfs --bind "$TMPDIR:/tmp" "$SIF" \
            plink2 --bfile "$FIRST_CHR" \
            --pmerge-list "$MERGE_LIST" bfile \
            --make-bed --out "$MERGED_PREFIX" \
            --threads $SLURM_CPUS_PER_TASK --memory 100000 2>&1 | tail -5

        # Si pmerge falla (triallelic), intentar plink1.9
        if [ ! -f "${MERGED_PREFIX}.bed" ]; then
            echo "  plink2 --pmerge falló, intentando plink1.9..."
            MERGE_LIST19="$TMPDIR/merge_list_p19.txt"
            > "$MERGE_LIST19"
            while read -r line; do
                echo "${line}.bed ${line}.bim ${line}.fam" >> "$MERGE_LIST19"
            done < "$MERGE_LIST"

            apptainer exec --bind /mnt/cephfs:/mnt/cephfs --bind "$TMPDIR:/tmp" "$SIF" \
                plink --bfile "$FIRST_CHR" \
                --merge-list "$MERGE_LIST19" \
                --make-bed --out "$MERGED_PREFIX" \
                --threads $SLURM_CPUS_PER_TASK --memory 100000 2>&1 | tail -5
        fi
    fi

    if [ ! -f "${MERGED_PREFIX}.bed" ]; then
        echo " Merge falló"
        exit 1
    fi

    N_VARS=$(wc -l < "${MERGED_PREFIX}.bim")
    N_SAMPS=$(wc -l < "${MERGED_PREFIX}.fam")
    echo "   Merged: ${N_VARS} variantes × ${N_SAMPS} muestras"
    echo "  Archivo: ${MERGED_PREFIX}.bed ($(du -sh "${MERGED_PREFIX}.bed" | awk '{print $1}'))"

    # Limpiar filtrados temporales
    rm -rf "$FILTERED_DIR"
fi

# ================================================================
# PASOS 1-6: Todo en R
# ================================================================
echo ""
echo "--- Ejecutando LDpred2 en R ---"

R_SCRIPT="$TMPDIR/run_ldpred2_v3.R"
cat << 'REOF' > "$R_SCRIPT"
library(bigsnpr)
library(data.table)

ncores <- 24
options(bigstatsr.check.parallel.blas = FALSE)

out_dir    <- "/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes/ldpred2"
gwas_dir   <- "/mnt/cephfs/orgs/home/angel.pacheco/gwas_miami"
covar_file <- "/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes/MCPS BASELINE.csv"
link_file  <- "/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes/RGN_LINK_IID.csv"

merged_bed <- file.path(out_dir, "mcps_79k_gwas_variants.bed")

# ================================================================
# PASO 1: Convertir PLINK merged → bigSNP
# ================================================================
cat(sprintf("\n[%s] PASO 1: Convertir PLINK → bigSNP\n", format(Sys.time(), "%H:%M:%S")))

rds_file <- file.path(out_dir, "mcps_79k_gwas_variants.rds")
bk_file  <- file.path(out_dir, "mcps_79k_gwas_variants.bk")

if (file.exists(rds_file)) {
    cat("  Cargando bigSNP existente...\n")
    obj.bigSNP <- snp_attach(rds_file)
} else {
    cat(sprintf("  Leyendo %s...\n", merged_bed))
    cat("  Esto puede tomar 15-30 min...\n")
    obj.bigSNP <- snp_readBed2(
        merged_bed,
        backingfile = tools::file_path_sans_ext(bk_file),
        ncores = ncores
    )
    cat("   bigSNP creado\n")
}

G   <- obj.bigSNP$genotypes
map <- obj.bigSNP$map
fam <- obj.bigSNP$fam

cat(sprintf("  Dimensiones: %s variantes × %s muestras\n",
    format(ncol(G), big.mark = ","), format(nrow(G), big.mark = ",")))
cat(sprintf("  Cromosomas: %s\n", paste(sort(unique(map$chromosome)), collapse = ", ")))

# Verificar estructura del map
cat("\n  Estructura map (primeras 5):\n")
cat(sprintf("  Columnas: %s\n", paste(names(map), collapse = ", ")))
print(head(map, 5))

# ================================================================
# PASO 2: Cargar GWAS MCPS + meta-análisis IVW
# ================================================================
cat(sprintf("\n[%s] PASO 2: Cargar GWAS + meta-análisis\n", format(Sys.time(), "%H:%M:%S")))

cat("  Leyendo GWAS females...\n")
gwas_f <- fread(file.path(gwas_dir, "combined-gwas-results-all-ctrl-females.txt.gz"),
                select = c("CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1", "A1FREQ", "INFO", "N", "BETA", "SE", "P"),
                showProgress = FALSE)

cat("  Leyendo GWAS males...\n")
gwas_m <- fread(file.path(gwas_dir, "combined-gwas-results-all-ctrl-males.txt.gz"),
                select = c("CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1", "A1FREQ", "INFO", "N", "BETA", "SE", "P"),
                showProgress = FALSE)

cat(sprintf("  Females: %s, Males: %s\n",
    format(nrow(gwas_f), big.mark = ","), format(nrow(gwas_m), big.mark = ",")))

# IVW meta-analysis
cat("  Meta-análisis IVW...\n")
gwas_merged <- merge(gwas_f, gwas_m, by = "ID", suffixes = c(".f", ".m"))
rm(gwas_f, gwas_m); gc()

gwas_merged[, w_f := 1 / SE.f^2]
gwas_merged[, w_m := 1 / SE.m^2]
gwas_merged[, BETA_meta := (BETA.f * w_f + BETA.m * w_m) / (w_f + w_m)]
gwas_merged[, SE_meta := sqrt(1 / (w_f + w_m))]
gwas_merged[, P_meta := 2 * pnorm(-abs(BETA_meta / SE_meta))]
gwas_merged[, N_meta := N.f + N.m]
gwas_merged[, A1FREQ_meta := (A1FREQ.f * N.f + A1FREQ.m * N.m) / (N.f + N.m)]

# Formatear para snp_match
# REGENIE: ALLELE0 = non-effect, ALLELE1 = effect
# snp_match: a1 = effect, a0 = reference (non-effect)
sumstats <- data.frame(
    chr     = gwas_merged$CHROM.f,
    pos     = gwas_merged$GENPOS.f,
    a1      = gwas_merged$ALLELE1.f,   # effect allele
    a0      = gwas_merged$ALLELE0.f,   # non-effect allele
    freq    = gwas_merged$A1FREQ_meta,
    info    = pmin(gwas_merged$INFO.f, gwas_merged$INFO.m),
    beta    = gwas_merged$BETA_meta,
    beta_se = gwas_merged$SE_meta,
    p       = gwas_merged$P_meta,
    n_eff   = gwas_merged$N_meta,
    stringsAsFactors = FALSE
)
rm(gwas_merged); gc()

# QC
n_pre <- nrow(sumstats)
sumstats <- sumstats[sumstats$info >= 0.8, ]
sumstats <- sumstats[sumstats$freq >= 0.01 & sumstats$freq <= 0.99, ]
sumstats <- sumstats[!is.na(sumstats$beta) & !is.na(sumstats$beta_se), ]
sumstats <- sumstats[sumstats$beta_se > 0, ]

# Convertir chr: GWAS usa enteros (1-22), PLINK .bim usa enteros y 23 para chrX
sumstats$chr <- as.character(sumstats$chr)
sumstats$chr[sumstats$chr %in% c("X", "chrX")] <- "23"
sumstats$chr <- as.integer(sumstats$chr)

cat(sprintf("  QC: %s → %s variantes\n",
    format(n_pre, big.mark = ","), format(nrow(sumstats), big.mark = ",")))

# ================================================================
# PASO 3: Match variantes GWAS ↔ genotipos
# ================================================================
cat(sprintf("\n[%s] PASO 3: snp_match\n", format(Sys.time(), "%H:%M:%S")))

# bigsnpr map tras snp_readBed2:
#   chromosome, marker.ID, genetic.dist, physical.pos, allele1, allele2
# En PLINK1 .bim: col1=chr, col2=ID, col3=cM, col4=bp, col5=A1(ALT), col6=A2(REF)
# bigsnpr lee: allele1 = col5 (ALT), allele2 = col6 (REF)

cat("  Verificación alleles bigSNP:\n")
for (i in 1:min(5, nrow(map))) {
    cat(sprintf("    chr%s:%s allele1(ALT)=%s allele2(REF)=%s  ID=%s\n",
        map$chromosome[i], map$physical.pos[i], map$allele1[i], map$allele2[i], map$marker.ID[i]))
}

# snp_match: matchea sumstats$a1 con map_df$a1
# En nuestro caso:
#   sumstats$a1 = ALLELE1 de REGENIE = effect allele = ALT
#   map$allele1 = col5 del .bim = ALT
# → Direct match cuando GWAS y genotipos reportan el mismo ALT

map_df <- data.frame(
    chr = map$chromosome,
    pos = map$physical.pos,
    a1  = map$allele1,    # ALT del .bim
    a0  = map$allele2,    # REF del .bim
    stringsAsFactors = FALSE
)

info_snp <- snp_match(sumstats, map_df, strand_flip = FALSE, join_by_pos = TRUE,
                       match.min.prop = 0.05)

n_matched <- nrow(info_snp)
cat(sprintf("  Variantes matcheadas: %s / %s (%.1f%%)\n",
    format(n_matched, big.mark = ","),
    format(nrow(sumstats), big.mark = ","),
    n_matched / nrow(sumstats) * 100))

if ("_FLIP_" %in% names(info_snp)) {
    n_flip <- sum(info_snp$`_FLIP_`, na.rm = TRUE)
    cat(sprintf("  Flipped (beta * -1): %s\n", format(n_flip, big.mark = ",")))
}

if (n_matched < 100000) {
    cat("\n⚠️uy pocas variantes matcheadas. Posible problema de formato.\n")
    cat("  Verificando muestra de IDs...\n")
    cat("  Sumstats (chr, pos, a1, a0):\n")
    print(head(sumstats[, c("chr", "pos", "a1", "a0")], 5))
    cat("  Map (chr, pos, a1, a0):\n")
    print(head(map_df, 5))
}

rm(sumstats); gc()

# ================================================================
# PASO 4: LD matrix + LDpred2-auto MCMC
# ================================================================
cat(sprintf("\n[%s] PASO 4: LD matrix + LDpred2-auto\n", format(Sys.time(), "%H:%M:%S")))

ind.chr <- info_snp$`_NUM_ID_`

ld_file <- file.path(out_dir, "ldpred2_ld_matrix.rds")

if (file.exists(ld_file)) {
    cat("  Cargando LD matrix existente...\n")
    corr <- readRDS(ld_file)
} else {
    cat(sprintf("  Calculando LD matrix (%s variantes)...\n",
        format(length(ind.chr), big.mark = ",")))
    cat("  Ventana: 3000 kb (~3 Mb), con infos.chr para evitar LD trans-cromosómico\n")
    cat("  Esto puede tomar 2-8 horas...\n")

    # Según tutorial oficial de Privé (LDpred2):
    #   - infos.pos en kb (pos / 1000)
    #   - size = 3000 (3000 kb = 3 Mb)
    #   - infos.chr OBLIGATORIO para no calcular LD entre cromosomas
    corr <- snp_cor(
        Gna       = G,
        ind.col   = ind.chr,
        infos.chr = info_snp$chr,
        infos.pos = info_snp$pos / 1000,
        size      = 3000,
        ncores    = ncores
    )

    saveRDS(corr, ld_file)
    cat(sprintf("  LD matrix guardada: %s (%.1f GB)\n",
        ld_file, file.size(ld_file) / 1e9))
}

cat(sprintf("  LD matrix: %s × %s\n",
    format(nrow(corr), big.mark = ","), format(ncol(corr), big.mark = ",")))

# LDpred2-auto MCMC
cat("\n  Lanzando LDpred2-auto (30 chains, 500+500 iter)...\n")
cat("  Estimación: 1-6 horas...\n\n")

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

cat("  LDpred2-auto completado\n\n")

# Extraer estimaciones de todas las chains
h2_est <- sapply(multi_auto, function(auto) auto$h2_est)
p_est  <- sapply(multi_auto, function(auto) auto$p_est)

cat(sprintf("  h2 por chain: min=%.4f, median=%.4f, max=%.4f\n",
    min(h2_est), median(h2_est), max(h2_est)))
cat(sprintf("  p  por chain: min=%.2e, median=%.2e, max=%.2e\n",
    min(p_est), median(p_est), max(p_est)))

# Filtrar chains convergentes
keep <- which(h2_est > 0.001 & h2_est < 0.5 & p_est > 1e-6 & p_est < 0.9)

if (length(keep) >= 10) {
    h2_med <- median(h2_est[keep])
    h2_iqr <- IQR(h2_est[keep])
    keep2 <- keep[h2_est[keep] > (h2_med - 3*h2_iqr) & h2_est[keep] < (h2_med + 3*h2_iqr)]
    if (length(keep2) >= 5) keep <- keep2
}

if (length(keep) < 3) {
    cat("  ⚠️ocas chains convergentes, usando todas\n")
    keep <- seq_along(multi_auto)
}

cat(sprintf("  Chains usadas: %d / %d\n", length(keep), length(multi_auto)))

# Beta final: promedio de chains convergentes
betas_auto <- sapply(multi_auto[keep], function(auto) auto$beta_est)
beta_final <- rowMeans(betas_auto)

cat(sprintf("  Betas no-cero: %s / %s\n",
    format(sum(beta_final != 0), big.mark = ","),
    format(length(beta_final), big.mark = ",")))

# Guardar
saveRDS(list(
    beta = beta_final, info_snp = info_snp,
    h2 = h2_est, p = p_est, keep = keep,
    multi_auto = multi_auto
), file.path(out_dir, "ldpred2_auto_results.rds"))
cat("  Resultados guardados: ldpred2_auto_results.rds\n")

# ================================================================
# PASO 5: Scoring
# ================================================================
cat(sprintf("\n[%s] PASO 5: Scoring\n", format(Sys.time(), "%H:%M:%S")))

prs_ldpred2 <- big_prodVec(
    X       = G,
    y.col   = beta_final,
    ind.col = ind.chr,
    ncores  = ncores
)

cat(sprintf("  PRS para %s muestras\n", format(length(prs_ldpred2), big.mark = ",")))
cat(sprintf("  Rango: [%.6f, %.6f]\n", min(prs_ldpred2), max(prs_ldpred2)))

# ================================================================
# PASO 6: Exportar
# ================================================================
cat(sprintf("\n[%s] PASO 6: Exportando\n", format(Sys.time(), "%H:%M:%S")))

sample_ids <- fam$sample.ID

df_link <- read.csv(link_file)
iid_to_patid <- setNames(as.character(df_link$PATID), as.character(df_link$IID))

df_baseline <- read.csv(covar_file)
df_baseline <- df_baseline[!is.na(df_baseline$MALE), ]
patid_to_sex <- setNames(as.integer(df_baseline$MALE), as.character(df_baseline$PATID))

patids <- iid_to_patid[as.character(sample_ids)]
patids[is.na(patids)] <- "UNKNOWN"
sexes <- patid_to_sex[patids]
sexes[is.na(sexes)] <- -1L

df_out <- data.frame(
    PATID = patids, IID = sample_ids, SEX = sexes,
    PRS_TOTAL = prs_ldpred2, stringsAsFactors = FALSE
)
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

out_file <- file.path(out_dir, "prs_diabetes_ldpred2_final.tsv")
df_out <- df_out[order(-df_out$PRS_TOTAL), ]
write.table(df_out, out_file, sep = "\t", row.names = FALSE, quote = FALSE)

cat(sprintf("\n PRS LDpred2 GUARDADO: %s\n", out_file))
cat(sprintf("   Muestras: %s\n", format(nrow(df_out), big.mark = ",")))
cat(sprintf("   h2 estimado: %.4f\n", median(h2_est[keep])))
cat(sprintf("   p estimado: %.2e\n", median(p_est[keep])))
cat(sprintf("   Variantes matcheadas: %s\n", format(n_matched, big.mark = ",")))
cat(sprintf("   Variantes con beta != 0: %s\n", format(sum(beta_final != 0), big.mark = ",")))
cat(sprintf("\n=== LDpred2 pipeline completado: %s ===\n", Sys.time()))
REOF

apptainer exec --bind /mnt/cephfs:/mnt/cephfs --bind "$TMPDIR:/tmp,$TMPDIR:/var/tmp" \
    "$SIF" /opt/miniforge3/envs/r-bio/bin/Rscript "$R_SCRIPT"

rm -rf "$TMPDIR"
echo ""
echo "============================================================"
echo "  PIPELINE B v3 COMPLETADO"
echo "  Output: ${OUT_DIR}/prs_diabetes_ldpred2_final.tsv"
echo "============================================================"
