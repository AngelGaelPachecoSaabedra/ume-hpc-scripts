#!/bin/bash
#SBATCH --job-name=ldpred2_uni
#SBATCH --account=researchers
#SBATCH --qos=vip
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=120G
#SBATCH --time=48:00:00
#SBATCH --output=/mnt/cephfs/orgs/home/angel.pacheco/logs/ldpred2_uni_%j.out
#SBATCH --error=/mnt/cephfs/orgs/home/angel.pacheco/logs/ldpred2_uni_%j.err
#SBATCH --mail-user=angelpachecosaavedra@gmail.com
#SBATCH --mail-type=BEGIN,END,FAIL

set -euo pipefail

# ================================================================
# PIPELINE — LDpred2-auto (bigsnpr) —
# 
# Flujo Inteligente:
# 1. Comprueba si el merged .bed gigante existe.
# 2. Si no existe, lo crea extrayendo IDs y usando PLINK pmerge.
# 3. Divide el archivo gigante por cromosomas en el directorio Scratch.
# 4. R lee los cromosomas secuencialmente desde Scratch (evita límite 1.1TB).
# 5. MCMC, Scoring, y limpieza final automática.
# ================================================================

SIF="/mnt/cephfs/biocontainers/images/jupyter-biotools-1.4.sif"

PLINK_DIR="/mnt/cephfs/hot_nvme/mcps/imputed-topmed/plink_files/maximally_unrelated"
GWAS_DIR="/mnt/cephfs/orgs/home/angel.pacheco/gwas_miami"
OUT_DIR="/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes/ldpred2"
mkdir -p "$OUT_DIR"

# Entorno Scratch para archivos temporales masivos
TMPDIR="/mnt/cephfs/scratch/$USER/ldpred2_uni_tmp/$SLURM_JOB_ID"
mkdir -p "$TMPDIR"
chmod 700 "$TMPDIR"
export TMPDIR TMP="$TMPDIR" TEMP="$TMPDIR"

echo "============================================================"
echo "  PIPELINE B — LDpred2-auto — MAESTRO UNIFICADO"
echo "  Nodo: $SLURMD_NODENAME"
echo "  RAM: $(free -g | awk 'NR==2{print $2}')G"
echo "  Workspace: $TMPDIR"
echo "  $(date)"
echo "============================================================"

# ================================================================
# PASO 0: Gestión del Archivo Gigante y Scratch
# ================================================================
MERGED_PREFIX="${OUT_DIR}/mcps_79k_gwas_variants"

if [ -f "${MERGED_PREFIX}.bed" ]; then
    echo ""
    echo "--- Merged .bed de 258GB detectado, saltando creación ---"
    N_VARS=$(wc -l < "${MERGED_PREFIX}.bim")
    N_SAMPS=$(wc -l < "${MERGED_PREFIX}.fam")
    echo "  Cached: ${N_VARS} variantes × ${N_SAMPS} muestras"
    
    echo "  Dividiendo el archivo gigante en cromosomas hacia Scratch para R..."
    for CHR in $(seq 1 22); do
        OUT_PREFIX="${TMPDIR}/chr${CHR}_filt"
        apptainer exec --bind /mnt/cephfs:/mnt/cephfs --bind "$TMPDIR:/tmp" "$SIF" \
            plink2 --bfile "${MERGED_PREFIX}" --chr $CHR \
            --make-bed --out "$OUT_PREFIX" --threads 4 --memory 8000 2>/dev/null
    done
else
    echo ""
    echo "--- Archivo gigante no detectado. Creándolo desde cero ---"
    GWAS_IDS="$TMPDIR/gwas_variant_ids.txt"

    echo "  Extrayendo IDs de GWAS (INFO >= 0.8)..."
    zcat "${GWAS_DIR}/combined-gwas-results-all-ctrl-females.txt.gz" | awk -F'\t' 'NR>1 && $7>=0.8 {print "chr"$3}' > "$TMPDIR/ids_f.txt"
    zcat "${GWAS_DIR}/combined-gwas-results-all-ctrl-males.txt.gz" | awk -F'\t' 'NR>1 && $7>=0.8 {print "chr"$3}' > "$TMPDIR/ids_m.txt"
    sort -S 8G --parallel=8 "$TMPDIR/ids_f.txt" > "$TMPDIR/ids_f_s.txt"
    sort -S 8G --parallel=8 "$TMPDIR/ids_m.txt" > "$TMPDIR/ids_m_s.txt"
    comm -12 "$TMPDIR/ids_f_s.txt" "$TMPDIR/ids_m_s.txt" > "$GWAS_IDS"

    echo "  Filtrando PLINK por cromosoma hacia Scratch..."
    MERGE_LIST="$TMPDIR/merge_list.txt"
    > "$MERGE_LIST"
    FIRST_CHR=""
    
    for CHR in $(seq 1 22) X; do
        PREFIX="${PLINK_DIR}/mcps-freeze150k_qcd_chr${CHR}_ivs"
        OUT_PREFIX="${TMPDIR}/chr${CHR}_filt"
        
        if [ -f "${PREFIX}.bed" ]; then
            apptainer exec --bind /mnt/cephfs:/mnt/cephfs --bind "$TMPDIR:/tmp" "$SIF" \
                plink2 --bfile "$PREFIX" --extract "$GWAS_IDS" \
                --make-bed --out "$OUT_PREFIX" --threads 4 --memory 8000 2>/dev/null
            
            if [ -f "${OUT_PREFIX}.bed" ]; then
                if [ -z "$FIRST_CHR" ]; then
                    FIRST_CHR="$OUT_PREFIX"
                else
                    echo "$OUT_PREFIX" >> "$MERGE_LIST"
                fi
            fi
        fi
    done

    echo "  Mergeando cromosomas para crear archivo maestro consolidado..."
    if [ -s "$MERGE_LIST" ]; then
        apptainer exec --bind /mnt/cephfs:/mnt/cephfs --bind "$TMPDIR:/tmp" "$SIF" \
            plink2 --bfile "$FIRST_CHR" --pmerge-list "$MERGE_LIST" bfile \
            --make-bed --out "$MERGED_PREFIX" --threads $SLURM_CPUS_PER_TASK --memory 100000 2>/dev/null
        
        # Fallback a plink 1.9 si pmerge falla
        if [ ! -f "${MERGED_PREFIX}.bed" ]; then
            echo "  plink2 --pmerge falló, intentando plink1.9..."
            MERGE_LIST19="$TMPDIR/merge_list_p19.txt"
            > "$MERGE_LIST19"
            while read -r line; do echo "${line}.bed ${line}.bim ${line}.fam" >> "$MERGE_LIST19"; done < "$MERGE_LIST"
            apptainer exec --bind /mnt/cephfs:/mnt/cephfs --bind "$TMPDIR:/tmp" "$SIF" \
                plink --bfile "$FIRST_CHR" --merge-list "$MERGE_LIST19" \
                --make-bed --out "$MERGED_PREFIX" --threads $SLURM_CPUS_PER_TASK --memory 100000 2>/dev/null
        fi
    fi
fi

# ================================================================
# PASOS 1-6: LDpred2 Per-Chromosome en R
# ================================================================
echo ""
echo "--- Ejecutando LDpred2 Chunked en R ---"

R_SCRIPT="$TMPDIR/run_ldpred2_uni.R"
cat << 'REOF' > "$R_SCRIPT"
library(bigsnpr)
library(data.table)

ncores <- 24
options(bigstatsr.check.parallel.blas = FALSE)
options(bigstatsr.ncores.max = 24) # FIX para el warning de Slurm

out_dir    <- "/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes/ldpred2"
gwas_dir   <- "/mnt/cephfs/orgs/home/angel.pacheco/gwas_miami"
tmp_dir    <- Sys.getenv("TMPDIR")

cat("\n[1] Preparando GWAS Meta-análisis...\n")
gf <- fread(file.path(gwas_dir, "combined-gwas-results-all-ctrl-females.txt.gz"), select=c("CHROM","GENPOS","ID","ALLELE0","ALLELE1","A1FREQ","INFO","N","BETA","SE","P"), showProgress=F)
gm <- fread(file.path(gwas_dir, "combined-gwas-results-all-ctrl-males.txt.gz"), select=c("CHROM","GENPOS","ID","ALLELE0","ALLELE1","A1FREQ","INFO","N","BETA","SE","P"), showProgress=F)
gwas <- merge(gf, gm, by="ID", suffixes=c(".f", ".m"))
rm(gf, gm); gc()

w_f <- 1/gwas$SE.f^2; w_m <- 1/gwas$SE.m^2
sumstats <- data.frame(
    chr = as.integer(gwas$CHROM.f), pos = gwas$GENPOS.f,
    a1 = gwas$ALLELE1.f, a0 = gwas$ALLELE0.f, freq = (gwas$A1FREQ.f*gwas$N.f + gwas$A1FREQ.m*gwas$N.m)/(gwas$N.f+gwas$N.m),
    info = pmin(gwas$INFO.f, gwas$INFO.m),
    beta = (gwas$BETA.f*w_f + gwas$BETA.m*w_m)/(w_f+w_m),
    beta_se = sqrt(1/(w_f+w_m)), n_eff = gwas$N.f + gwas$N.m,
    stringsAsFactors = FALSE
)
sumstats$p <- 2 * pnorm(-abs(sumstats$beta / sumstats$beta_se))
sumstats <- sumstats[sumstats$info >= 0.8 & sumstats$freq >= 0.01 & sumstats$freq <= 0.99 & !is.na(sumstats$beta) & !is.na(sumstats$beta_se) & sumstats$beta_se > 0, ]
rm(gwas); gc()

prs_total <- NULL
fam_global <- NULL
matched_total <- 0

cat("\n[2] MCMC y Scoring por Cromosoma...\n")
for (chr in 1:22) {
    bk_file <- file.path(tmp_dir, sprintf("chr%d_filt", chr))
    bed_file <- paste0(bk_file, ".bed")
    
    if (!file.exists(bed_file)) next
    cat(sprintf("\n  >>> Procesando chr%d <<<\n", chr))
    
    # Lectura secuencial estable
    snp_readBed(bed_file, backingfile = bk_file)
    obj.bigSNP <- snp_attach(paste0(bk_file, ".rds"))
    G <- obj.bigSNP$genotypes
    
    if (is.null(prs_total)) {
        prs_total <- rep(0, nrow(G))
        fam_global <- obj.bigSNP$fam
    }
    
    # Match
    map_df <- data.frame(chr=obj.bigSNP$map$chromosome, pos=obj.bigSNP$map$physical.pos, a1=obj.bigSNP$map$allele1, a0=obj.bigSNP$map$allele2, stringsAsFactors=F)
    sumstats_chr <- sumstats[sumstats$chr == chr, ]
    info_snp <- snp_match(sumstats_chr, map_df, strand_flip=FALSE, join_by_pos=TRUE, match.min.prop=0.01)
    
    cat(sprintf("      Matched: %d variantes\n", nrow(info_snp)))
    matched_total <- matched_total + nrow(info_snp)
    if (nrow(info_snp) == 0) next
    
    # LD Matrix
    # LD Matrix
    corr <- snp_cor(Gna = G, ind.col = info_snp$`_NUM_ID_`, infos.pos = info_snp$pos / 1000, size = 3000, ncores = ncores)
    #corr <- snp_cor(Gna = G, ind.col = info_snp$`_NUM_ID_`, infos.chr = info_snp$chr, infos.pos = info_snp$pos / 1000, size = 3000, ncores = ncores)
    
    # MCMC LDpred2-auto
    df_beta <- info_snp[, c("beta", "beta_se", "n_eff")]
    multi_auto <- snp_ldpred2_auto(
        corr = corr, df_beta = df_beta, h2_init = 0.01,
        vec_p_init = seq_log(1e-4, 0.5, length.out = 30),
        burn_in = 500, num_iter = 500, report_step = 100, ncores = ncores, sparse = TRUE
    )
    
    keep <- which(sapply(multi_auto, function(x) x$h2_est) > 0.001)
    if (length(keep) < 3) keep <- seq_along(multi_auto)
    beta_final <- rowMeans(sapply(multi_auto[keep], function(x) x$beta_est))
    
    # Scoring Inmediato
    prs_chr <- big_prodVec(X = G, y.col = beta_final, ind.col = info_snp$`_NUM_ID_`, ncores = ncores)
    prs_total <- prs_total + prs_chr
    
    # Limpieza agresiva de RAM
    rm(obj.bigSNP, G, corr, multi_auto); gc()
}

cat(sprintf("\n[3] Exportando (Total variantes usadas: %d)...\n", matched_total))
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

df_out <- data.frame(PATID = patids, IID = fam_global$sample.ID, SEX = sexes, PRS_TOTAL = prs_total, stringsAsFactors = FALSE)
df_out <- df_out[df_out$SEX >= 0, ]

# Estandarizar Z-score y percentiles por sexo
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
cat(sprintf("Pipeline completado. Resultados en: %s\n", out_file))
REOF

apptainer exec --bind /mnt/cephfs:/mnt/cephfs --bind "$TMPDIR:/tmp,$TMPDIR:/var/tmp" \
    "$SIF" /opt/miniforge3/envs/r-bio/bin/Rscript "$R_SCRIPT"

# ================================================================
# Limpieza Final
# ================================================================
echo ""
echo "  Limpiando directorio Scratch temporal ($TMPDIR)..."
rm -rf "$TMPDIR"

echo "============================================================"
echo "  PIPELINE B COMPLETADO CON ÉXITO"
echo "============================================================"
