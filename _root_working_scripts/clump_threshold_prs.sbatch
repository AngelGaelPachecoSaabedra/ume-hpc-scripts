#!/bin/bash
#SBATCH --job-name=ct_prs
#SBATCH --account=researchers
#SBATCH --qos=vip
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=/mnt/cephfs/orgs/home/angel.pacheco/logs/ct_prs_%j.out
#SBATCH --error=/mnt/cephfs/orgs/home/angel.pacheco/logs/ct_prs_%j.err

set -euo pipefail

SIF="/mnt/cephfs/biocontainers/images/jupyter-biotools-1.2.sif"

# Temporales en Scratch
TMPDIR="/mnt/cephfs/scratch/$USER/ct_tmp/$SLURM_JOB_ID"
mkdir -p "$TMPDIR"
chmod 700 "$TMPDIR"
export TMPDIR TMP="$TMPDIR" TEMP="$TMPDIR"

# Rutas
GWAS_DIR="/mnt/cephfs/orgs/home/angel.pacheco/gwas_miami"
PLINK_DIR="/mnt/cephfs/hot_nvme/mcps/imputed-topmed/plink-prs"
OUT_DIR="/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes/ct_gwas"
mkdir -p "$OUT_DIR"

echo "=============================================="
echo "CLUMPING + THRESHOLDING - PRS MCPS GWAS"
echo "=============================================="

# ================================================================
# PASO 1: Preparar archivos GWAS en formato plink2 --clump
#
# plink2 --clump espera un archivo con columnas:
#   ID (variant ID), P (p-value)
#   opcionalmente: BETA, A1
#
# Nuestro GWAS de REGENIE tiene:
#   $3=ID, $4=ALLELE0(REF), $5=ALLELE1(ALT/effect), $10=BETA, $15=P
# ================================================================

echo ">> Paso 1: Preparando summary statistics..."

for SEX in females males; do
    echo "  Procesando ${SEX}..."
    
    GWAS_FILE="${GWAS_DIR}/combined-gwas-results-all-ctrl-${SEX}.txt.gz"
    CLEAN_FILE="${OUT_DIR}/gwas_clean_${SEX}.tsv"
    
    # Extraer, filtrar p < 5e-8, solo SNVs (1 carácter por alelo), solo autosomas
    # Agregar prefijo "chr" al ID para que coincida con los PLINK
    # GWAS: "1:598941:G:A" → "chr1:598941:G:A"
    # Excluir chrX (CHROM == 23 o X)
    zcat "$GWAS_FILE" | awk -F'\t' '
        NR == 1 {
            print "SNP\tCHR\tBP\tA1\tA2\tBETA\tP"
            next
        }
        $15 < 5e-8 && length($4) == 1 && length($5) == 1 && $1 != "23" && $1 != "X" {
            snp = "chr"$3
            print snp"\t"$1"\t"$2"\t"$5"\t"$4"\t"$10"\t"$15
        }
    ' > "$CLEAN_FILE"
    
    N_VARS=$(tail -n +2 "$CLEAN_FILE" | wc -l)
    echo "    Variantes significativas (SNVs): ${N_VARS}"
done

# ================================================================
# PASO 2: Clumping por cromosoma con plink2
#
# Parámetros C+T estándar:
#   --clump-p1 5e-8     : p-value threshold para variantes índice
#   --clump-r2 0.1      : r² threshold para declarar LD
#   --clump-kb 250      : ventana de 250kb alrededor de cada índice
#
# El resultado conserva solo la variante más significativa
# por cada bloque de LD independiente.
# ================================================================

echo ""
echo ">> Paso 2: Clumping por cromosoma..."

for SEX in females males; do
    echo "  === ${SEX} ==="
    
    CLEAN_FILE="${OUT_DIR}/gwas_clean_${SEX}.tsv"
    CLUMP_MERGED="${OUT_DIR}/clumped_${SEX}_all.tsv"
    
    # Header del archivo final
    echo -e "CHR\tBP\tSNP\tA1\tA2\tBETA\tP" > "$CLUMP_MERGED"
    
    for CHR in $(seq 1 22); do
        PLINK_PREFIX="${PLINK_DIR}/MCPS.HIS.TOPMED.b38.FR.2602.TARGET.imputed.chr${CHR}.plink"
        CHR_CODE="$CHR"
        CLUMP_OUT="${TMPDIR}/clump_${SEX}_chr${CHR}"
        
        # Verificar que el archivo PLINK existe
        if [ ! -f "${PLINK_PREFIX}.bed" ]; then
            echo "    chr${CHR}: PLINK no encontrado, saltando..."
            continue
        fi
        
        # Contar variantes de este cromosoma en el GWAS filtrado
        N_CHR=$(awk -F'\t' -v c="$CHR_CODE" 'NR>1 && $2==c' "$CLEAN_FILE" | wc -l)
        if [ "$N_CHR" -eq 0 ]; then
            echo "    chr${CHR}: 0 variantes significativas, saltando..."
            continue
        fi
        
        # Ejecutar clumping
        apptainer exec \
            --bind /mnt/cephfs:/mnt/cephfs \
            --bind "$TMPDIR:/tmp,$TMPDIR:/var/tmp" \
            "$SIF" \
            plink2 \
                --bfile "$PLINK_PREFIX" \
                --clump "$CLEAN_FILE" \
                --clump-p1 5e-8 \
                --clump-r2 0.1 \
                --clump-kb 250 \
                --clump-snp-field SNP \
                --clump-p-field P \
                --threads $SLURM_CPUS_PER_TASK \
                --out "$CLUMP_OUT" \
                2>/dev/null || true
        
        # Extraer variantes que sobrevivieron el clumping
        if [ -f "${CLUMP_OUT}.clumps" ]; then
            # plink2 .clumps tiene las variantes índice
            # Cruzar con el GWAS original para recuperar BETA y alelos
            awk 'NR>1 {print $3}' "${CLUMP_OUT}.clumps" > "${TMPDIR}/clump_ids_${SEX}_chr${CHR}.txt"
            
            N_CLUMPED=$(wc -l < "${TMPDIR}/clump_ids_${SEX}_chr${CHR}.txt")
            echo "    chr${CHR}: ${N_CHR} → ${N_CLUMPED} variantes independientes"
            
            # Recuperar info completa de las variantes clumped
            grep -F -f "${TMPDIR}/clump_ids_${SEX}_chr${CHR}.txt" "$CLEAN_FILE" | \
                awk -F'\t' '{print $2"\t"$3"\t"$1"\t"$4"\t"$5"\t"$6"\t"$7}' >> "$CLUMP_MERGED"
        else
            echo "    chr${CHR}: clumping no produjo resultados"
        fi
    done
    
    N_FINAL=$(tail -n +2 "$CLUMP_MERGED" | wc -l)
    echo "  Total ${SEX} después de C+T: ${N_FINAL} variantes independientes"
    echo ""
done

# ================================================================
# PASO 3: Resumen final
# ================================================================
echo ">> Paso 3: Resumen"
echo ""

for SEX in females males; do
    FILE="${OUT_DIR}/clumped_${SEX}_all.tsv"
    echo "=== ${SEX} ==="
    echo "Archivo: ${FILE}"
    echo "Variantes: $(tail -n +2 "$FILE" | wc -l)"
    echo "Cromosomas representados: $(tail -n +2 "$FILE" | cut -f1 | sort -u | tr '\n' ' ')"
    echo "Top 5 variantes más significativas:"
    tail -n +2 "$FILE" | sort -t$'\t' -k7 -g | head -5
    echo ""
done

rm -rf "$TMPDIR"
echo "=============================================="
echo "C+T COMPLETADO"
echo "Archivos finales en: ${OUT_DIR}/"
echo "  clumped_females_all.tsv"
echo "  clumped_males_all.tsv"
echo "=============================================="
