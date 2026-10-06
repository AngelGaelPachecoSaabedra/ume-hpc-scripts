#!/bin/bash
#SBATCH --job-name=sbrc_liftover
#SBATCH --account=researchers
#SBATCH --qos=vip
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/mnt/cephfs/orgs/home/angel.pacheco/logs/sbrc_liftover_%j.out
#SBATCH --error=/mnt/cephfs/orgs/home/angel.pacheco/logs/sbrc_liftover_%j.err
#SBATCH --mail-user=angelpachecosaavedra@gmail.com
#SBATCH --mail-type=END,FAIL

set -euo pipefail

# ================================================================
# FIX v3: LiftOver T2DGGI hg19 → hg38
#
# Problema: T2DGGI usa GRCh37/hg19, MCPS-TOPMed usa GRCh38/hg38
# Solo 1.2% matcheaba porque las coordenadas son de ensamblajes
# diferentes.
#
# Pipeline:
#   1. Descargar liftOver binary + chain file hg19→hg38
#   2. Convertir .ma → BED (hg19)
#   3. liftOver BED hg19 → BED hg38
#   4. Reconstruir .ma con coordenadas hg38
#   5. Rebuild SNP IDs con ambas orientaciones de alelos vs .bim
#   6. Validar match rate
# ================================================================

SIF="/mnt/cephfs/biocontainers/images/jupyter-biotools-1.4.sif"

PLINK_DIR="/mnt/cephfs/hot_nvme/mcps/imputed-topmed/plink_files/maximally_unrelated"
OUT_BASE="/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes/sbayesrc"
MA_ORIG="${OUT_BASE}/T2DGGI_ALL_autosomas.ma"
#MA_LIFTED=${OUT_BASE}/T2DGGI_ALL_hg38.ma"
MA_FINAL="${OUT_BASE}/T2DGGI_ALL_hg38_matched.ma"

TOOLS_DIR="/mnt/cephfs/orgs/home/angel.pacheco/tools"
mkdir -p "$TOOLS_DIR"

TMPDIR="/mnt/cephfs/scratch/$USER/sbrc_liftover_tmp/$SLURM_JOB_ID"
mkdir -p "$TMPDIR"
chmod 700 "$TMPDIR"
export TMPDIR TMP="$TMPDIR" TEMP="$TMPDIR"

echo "============================================================"
echo "  FIX v3: LiftOver T2DGGI de hg19 → hg38"
echo "  $(date)"
echo "============================================================"

# ================================================================
# 1. Descargar liftOver + chain file
# ================================================================
echo ""
echo "--- 1. Obteniendo liftOver y chain file ---"

LIFTOVER="$TOOLS_DIR/liftOver"
CHAIN="$TOOLS_DIR/hg19ToHg38.over.chain.gz"

if [ ! -x "$LIFTOVER" ]; then
    echo "  Descargando liftOver binary..."
    wget -q "https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver" -O "$LIFTOVER" 2>&1 || {
        echo "  ⚠️get directo falló, intentando con curl..."
        curl -sL "https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver" -o "$LIFTOVER"
    }
    chmod +x "$LIFTOVER"
fi

if [ ! -f "$CHAIN" ]; then
    echo "  Descargando chain file hg19→hg38..."
    wget -q "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz" -O "$CHAIN" 2>&1 || {
        echo "  ⚠️wget directo falló, intentando con curl..."
        curl -sL "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz" -o "$CHAIN"
    }
fi

echo "  liftOver: $(ls -lh "$LIFTOVER" | awk '{print $5}')"
echo "  Chain:    $(ls -lh "$CHAIN" | awk '{print $5}')"

# ================================================================
# 2. Convertir .ma → BED (hg19)
# ================================================================
echo ""
echo "--- 2. Convirtiendo .ma a BED format (hg19) ---"

BED_HG19="$TMPDIR/t2dggi_hg19.bed"
METADATA="$TMPDIR/t2dggi_metadata.tsv"

# BED format: chr  start(0-based)  end(1-based)  name
# Guardar metadata aparte para reconstruir después del liftOver
# El SNP ID actual en .ma es: CHR:POS:A1:A2 (ya sin chr prefix en original)
# Necesitamos extraer CHR y POS del ID

awk -F'\t' '
NR == 1 { next }
{
    snp = $1
    # Parsear SNP ID: puede ser "3:5484355:T:C" o ya "chr3:..."
    split(snp, parts, ":")
    chr = parts[1]
    pos = parts[2]

    # Asegurar prefijo chr para BED
    if (chr !~ /^chr/) chr = "chr" chr

    # BED es 0-based start
    start = pos - 1
    end = pos

    # Name = índice de línea para reconnectar después
    print chr "\t" start "\t" end "\t" NR-1
}
' "$MA_ORIG" > "$BED_HG19"

# Guardar todas las columnas originales indexadas
awk -F'\t' 'NR > 1 { print NR-1 "\t" $0 }' "$MA_ORIG" > "$METADATA"

N_BED=$(wc -l < "$BED_HG19")
echo "  Variantes en BED: ${N_BED}"
echo "  Primeras líneas:"
head -3 "$BED_HG19"

# ================================================================
# 3. Ejecutar liftOver hg19 → hg38
# ================================================================
echo ""
echo "--- 3. Ejecutando liftOver ---"

BED_HG38="$TMPDIR/t2dggi_hg38.bed"
BED_UNMAPPED="$TMPDIR/t2dggi_unmapped.bed"

"$LIFTOVER" "$BED_HG19" "$CHAIN" "$BED_HG38" "$BED_UNMAPPED"

N_LIFTED=$(wc -l < "$BED_HG38")
N_UNMAPPED=$(grep -c "^[^#]" "$BED_UNMAPPED" 2>/dev/null || echo 0)

echo "  Resultados:"
echo "    Lifted exitosamente: ${N_LIFTED}"
echo "    No mapeables:        ${N_UNMAPPED}"
echo "    Tasa de éxito:       $(echo "scale=1; $N_LIFTED * 100 / $N_BED" | bc)%"

echo "  Primeras líneas hg38:"
head -3 "$BED_HG38"

# ================================================================
# 4. Reconstruir .ma con coordenadas hg38
# ================================================================
echo ""
echo "--- 4. Reconstruyendo .ma con coordenadas hg38 ---"

PYTHON_SCRIPT="$TMPDIR/rebuild_ma.py"
cat << 'PYEOF' > "$PYTHON_SCRIPT"
import time
import sys

bed_hg38_file = sys.argv[1]
metadata_file = sys.argv[2]
bim_ids_file  = sys.argv[3]
ma_output     = sys.argv[4]

# Paso 1: Cargar posiciones hg38 del liftOver
print(f"[{time.strftime('%H:%M:%S')}] Cargando coordenadas hg38...")
hg38_coords = {}  # idx → (chr, pos_hg38)
with open(bed_hg38_file) as f:
    for line in f:
        cols = line.strip().split('\t')
        chrom = cols[0]        # chr20
        pos_hg38 = int(cols[2])  # end = 1-based position
        idx = int(cols[3])       # índice original
        hg38_coords[idx] = (chrom, pos_hg38)
print(f"  Coordenadas hg38 cargadas: {len(hg38_coords):,}")

# Paso 2: Cargar IDs del .bim
print(f"[{time.strftime('%H:%M:%S')}] Cargando IDs del .bim...")
bim_ids = set()
with open(bim_ids_file) as f:
    for line in f:
        bim_ids.add(line.strip())
print(f"  IDs .bim cargados: {len(bim_ids):,}")

# Paso 3: Reconstruir .ma con IDs que matcheen el .bim
print(f"[{time.strftime('%H:%M:%S')}] Reconstruyendo .ma...")

matched_fwd = 0
matched_rev = 0
unmatched_no_lift = 0
unmatched_no_bim = 0
total = 0

with open(metadata_file) as fin, open(ma_output, 'w') as fout:
    fout.write("SNP\tA1\tA2\tfreq\tb\tse\tp\tN\n")

    for line in fin:
        cols = line.strip().split('\t')
        idx = int(cols[0])
        # cols[1:] = SNP_orig, A1, A2, freq, b, se, p, N
        snp_orig = cols[1]
        a1 = cols[2]   # Effect allele
        a2 = cols[3]   # Non-effect allele
        freq = cols[4]
        beta = cols[5]
        se = cols[6]
        p = cols[7]
        n = cols[8]
        total += 1

        # ¿Tiene coordenada hg38?
        if idx not in hg38_coords:
            unmatched_no_lift += 1
            continue

        chrom, pos38 = hg38_coords[idx]

        # Construir IDs candidatos (formato .bim: chr20:60705:A:C)
        id_fwd = f"{chrom}:{pos38}:{a1}:{a2}"
        id_rev = f"{chrom}:{pos38}:{a2}:{a1}"

        if id_fwd in bim_ids:
            new_id = id_fwd
            matched_fwd += 1
        elif id_rev in bim_ids:
            new_id = id_rev
            matched_rev += 1
        else:
            unmatched_no_bim += 1
            continue

        fout.write(f"{new_id}\t{a1}\t{a2}\t{freq}\t{beta}\t{se}\t{p}\t{n}\n")

        if total % 2000000 == 0:
            print(f"  ... {total:,} procesadas, {matched_fwd + matched_rev:,} matched")

matched_total = matched_fwd + matched_rev
print(f"\n[{time.strftime('%H:%M:%S')}] RESULTADO FINAL:")
print(f"  Total variantes en .ma:     {total:,}")
print(f"  Lifted a hg38:              {len(hg38_coords):,} ({len(hg38_coords)/total*100:.1f}%)")
print(f"  Match forward (A1:A2):      {matched_fwd:,}")
print(f"  Match reverse (A2:A1):      {matched_rev:,}")
print(f"  Total matched con .bim:     {matched_total:,} ({matched_total/total*100:.1f}%)")
print(f"  Sin liftOver:               {unmatched_no_lift:,}")
print(f"  Sin match en .bim (hg38):   {unmatched_no_bim:,}")
PYEOF

# Primero regenerar la lista de IDs del .bim
BIM_IDS="$TMPDIR/all_bim_ids.txt"
echo "  Extrayendo IDs del .bim (chr1-22)..."
for CHR in $(seq 1 22); do
    BIM="${PLINK_DIR}/mcps-freeze150k_qcd_chr${CHR}_ivs.bim"
    if [ -f "$BIM" ]; then
        awk '{print $2}' "$BIM"
    fi
done | sort -S 8G --parallel=8 > "$BIM_IDS"
echo "  IDs .bim: $(wc -l < "$BIM_IDS")"

apptainer exec --bind /mnt/cephfs:/mnt/cephfs --bind "$TMPDIR:/tmp,$TMPDIR:/var/tmp" \
    "$SIF" /opt/miniforge3/envs/vcf2zarr/bin/python "$PYTHON_SCRIPT" \
    "$BED_HG38" "$METADATA" "$BIM_IDS" "$MA_FINAL"

# ================================================================
# 5. Validación
# ================================================================
echo ""
echo "--- 5. Validación ---"

echo "  Archivo final: $MA_FINAL"
echo "  Primeras líneas:"
head -5 "$MA_FINAL" | column -t

N_FINAL=$(tail -n +2 "$MA_FINAL" | wc -l)
echo ""
echo "  Variantes totales: ${N_FINAL}"

# Quick check chr20
echo ""
echo "  Validación chr20:"
awk '{print $2}' "${PLINK_DIR}/mcps-freeze150k_qcd_chr20_ivs.bim" | sort > "$TMPDIR/bim20.txt"
awk 'NR>1 && $1 ~ /^chr20:/' "$MA_FINAL" | awk '{print $1}' | sort > "$TMPDIR/ma20.txt"
N20_MA=$(wc -l < "$TMPDIR/ma20.txt")
N20_MATCH=$(comm -12 "$TMPDIR/bim20.txt" "$TMPDIR/ma20.txt" | wc -l)
echo "  .ma chr20: ${N20_MA} variantes"
echo "  Match chr20: ${N20_MATCH} / ${N20_MA} ($(echo "scale=1; $N20_MATCH * 100 / $N20_MA" | bc)%)"

# ================================================================
# 6. Limpiar LD directory viejo
# ================================================================
echo ""
echo "--- 6. Limpiando LD directory ---"
LD_DIR="${OUT_BASE}/ldm_79k"
if [ -d "$LD_DIR" ]; then
    rm -rf "$LD_DIR"
    echo "  ldm_79k limpiado"
fi

rm -rf "$TMPDIR"

echo ""
echo "============================================================"
echo "  LIFTOVER COMPLETADO"
echo ""
echo "  Archivo .ma final (hg38): $MA_FINAL"
echo "  Variantes: ${N_FINAL}"
echo ""
echo "  Actualizar sbrc_02_ldm_sbayesrc.sbatch:"
echo "    MA_FILE=\"${MA_FINAL}\""
echo ""
echo "  SBayesRC correrá SIN anotaciones (SBayesR mode)."
echo "  Cambiar en el R script: annot = NULL"
echo ""
echo "  Luego: sbatch sbrc_02_ldm_sbayesrc.sbatch"
echo "=======================================
