#!/bin/bash
#SBATCH --job-name=sbrc_fix2
#SBATCH --account=researchers
#SBATCH --qos=vip
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=/mnt/cephfs/orgs/home/angel.pacheco/logs/sbrc_fix2_%j.out
#SBATCH --error=/mnt/cephfs/orgs/home/angel.pacheco/logs/sbrc_fix2_%j.err

set -euo pipefail

# ================================================================
# FIX v2: SNP ID mismatch — problema de ORDEN DE ALELOS
#
# El prefijo chr ya se arregló, pero solo 695/448K matchean en chr20.
# Causa: el .bim usa CHR:POS:REF:ALT y el .ma usa CHR:POS:EA:NEA.
# Cuando EA != REF, el ID tiene alelos invertidos y no matchea.
#
# Solución: Reconstruir el SNP ID del .ma probando ambas orientaciones
# de alelos contra un set de IDs del .bim.
#
# Estrategia: Para cada variante del .ma, generar 2 IDs candidatos:
#   1. chr{CHR}:{POS}:{A1}:{A2}  (original)
#   2. chr{CHR}:{POS}:{A2}:{A1}  (flipped)
# y usar el que exista en el .bim.
# ================================================================

SIF="/mnt/cephfs/biocontainers/images/jupyter-biotools-1.4.sif"
SBRC_SIF="/mnt/cephfs/biocontainers/images/sbayesrc-0.2.6.sif"

PLINK_DIR="/mnt/cephfs/hot_nvme/mcps/imputed-topmed/plink_files/maximally_unrelated"
OUT_BASE="/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes/sbayesrc"
MA_CHRFIX="${OUT_BASE}/T2DGGI_ALL_autosomas_chrfix.ma"
MA_FINAL="${OUT_BASE}/T2DGGI_ALL_autosomas_matched.ma"

TMPDIR="/mnt/cephfs/scratch/$USER/sbrc_fix2_tmp/$SLURM_JOB_ID"
mkdir -p "$TMPDIR"
chmod 700 "$TMPDIR"
export TMPDIR TMP="$TMPDIR" TEMP="$TMPDIR"

echo "============================================================"
echo "  FIX v2: Reconstruir SNP IDs para match con .bim"
echo "  $(date)"
echo "============================================================"

# ================================================================
# 1. Extraer TODOS los SNP IDs del .bim (chr1-22)
# ================================================================
echo ""
echo "--- 1. Construyendo set de IDs del .bim ---"

BIM_IDS="$TMPDIR/all_bim_ids.txt"

for CHR in $(seq 1 22); do
    BIM="${PLINK_DIR}/mcps-freeze150k_qcd_chr${CHR}_ivs.bim"
    if [ -f "$BIM" ]; then
        awk '{print $2}' "$BIM"
    fi
done | sort -S 8G --parallel=8 > "$BIM_IDS"

N_BIM=$(wc -l < "$BIM_IDS")
echo "  Total IDs en .bim: ${N_BIM}"
echo "  Ejemplos:"
head -3 "$BIM_IDS"

# ================================================================
# 2. Reconstruir IDs del .ma con match bidireccional
# ================================================================
echo ""
echo "--- 2. Reconstruyendo IDs del .ma ---"

PYTHON_SCRIPT="$TMPDIR/fix_ma_ids.py"
cat << 'PYEOF' > "$PYTHON_SCRIPT"
import time
import sys

bim_ids_file = sys.argv[1]
ma_input     = sys.argv[2]
ma_output    = sys.argv[3]

print(f"[{time.strftime('%H:%M:%S')}] Cargando IDs del .bim...")
bim_ids = set()
with open(bim_ids_file) as f:
    for line in f:
        bim_ids.add(line.strip())
print(f"  IDs cargados: {len(bim_ids):,}")

print(f"[{time.strftime('%H:%M:%S')}] Procesando .ma y reconstruyendo IDs...")

matched_fwd = 0    # chr:pos:A1:A2 matchea
matched_rev = 0    # chr:pos:A2:A1 matchea
matched_pos = 0    # match por chr:pos solo (con alelos del .bim)
unmatched = 0
total = 0

with open(ma_input) as fin, open(ma_output, 'w') as fout:
    header = fin.readline()
    fout.write(header)
    
    for line in fin:
        total += 1
        cols = line.strip().split('\t')
        snp_id = cols[0]   # chr3:5484355:T:C
        a1 = cols[1]       # Effect allele
        a2 = cols[2]       # Non-effect allele
        
        # Parsear el SNP ID actual
        parts = snp_id.split(':')
        if len(parts) >= 4:
            chrom = parts[0]   # chr3
            pos   = parts[1]   # 5484355
        else:
            # ID no tiene formato chr:pos:a1:a2 — skip
            unmatched += 1
            continue
        
        # Candidato 1: ID original (ya con chr prefix)
        id_fwd = f"{chrom}:{pos}:{a1}:{a2}"
        
        # Candidato 2: alelos invertidos
        id_rev = f"{chrom}:{pos}:{a2}:{a1}"
        
        if id_fwd in bim_ids:
            cols[0] = id_fwd
            matched_fwd += 1
        elif id_rev in bim_ids:
            cols[0] = id_rev
            matched_rev += 1
        else:
            # Ni forward ni reverse matchea — descartar
            unmatched += 1
            continue
        
        fout.write('\t'.join(cols) + '\n')
        
        if total % 2000000 == 0:
            print(f"  ... {total:,} procesadas, {matched_fwd + matched_rev:,} matched")

print(f"\n[{time.strftime('%H:%M:%S')}] RESULTADO:")
print(f"  Total variantes: {total:,}")
print(f"  Match forward (A1:A2): {matched_fwd:,}")
print(f"  Match reverse (A2:A1): {matched_rev:,}")
print(f"  Total matched: {matched_fwd + matched_rev:,} ({(matched_fwd + matched_rev)/total*100:.1f}%)")
print(f"  Sin match: {unmatched:,}")
PYEOF

apptainer exec --bind /mnt/cephfs:/mnt/cephfs --bind "$TMPDIR:/tmp,$TMPDIR:/var/tmp" \
    "$SIF" /opt/miniforge3/envs/vcf2zarr/bin/python "$PYTHON_SCRIPT" \
    "$BIM_IDS" "$MA_CHRFIX" "$MA_FINAL"

echo ""
echo "  Verificación del archivo final:"
head -3 "$MA_FINAL" | column -t
echo ""
N_FINAL=$(tail -n +2 "$MA_FINAL" | wc -l)
echo "  Variantes en .ma final: ${N_FINAL}"

# Quick check chr20
echo ""
echo "  Quick match validation (chr20):"
awk '{print $2}' "${PLINK_DIR}/mcps-freeze150k_qcd_chr20_ivs.bim" | sort > "$TMPDIR/bim20.txt"
awk 'NR>1 && $1 ~ /^chr20:/' "$MA_FINAL" | awk '{print $1}' | sort > "$TMPDIR/ma20.txt"
N20=$(comm -12 "$TMPDIR/bim20.txt" "$TMPDIR/ma20.txt" | wc -l)
echo "  chr20: ${N20} variantes en común (de $(wc -l < "$TMPDIR/ma20.txt") en .ma)"

# ================================================================
# 3. Buscar anotaciones en el contenedor SBayesRC (timeout 60s)
# ================================================================
echo ""
echo "--- 3. Buscando anotaciones BaselineLD ---"

timeout 60 apptainer exec "$SBRC_SIF" bash -c '
# Buscar solo en ubicaciones probables (no find / completo)
echo "  Checando R package extdata..."
Rscript -e "
pkg_path <- system.file(package=\"SBayesRC\")
cat(sprintf(\"  Package path: %s\n\", pkg_path))
extdata <- system.file(\"extdata\", package=\"SBayesRC\")
cat(sprintf(\"  extdata path: %s\n\", extdata))
if (nchar(extdata) > 0) {
    files <- list.files(extdata, full.names=TRUE)
    cat(sprintf(\"  Files: %s\n\", paste(files, collapse=\", \")))
} else {
    cat(\"  No extdata found in package\n\")
}
" 2>/dev/null

echo ""
echo "  Checando paths comunes..."
for d in /opt /usr/local/lib/R /usr/lib/R /home; do
    f=$(find "$d" -maxdepth 5 -name "annot_baseline*" -type f 2>/dev/null | head -3)
    if [ -n "$f" ]; then echo "  FOUND: $f"; fi
done

echo ""
echo "  Checando GCTB..."
which gctb 2>/dev/null || find /usr /opt -maxdepth 3 -name "gctb" -type f 2>/dev/null | head -3
' 2>/dev/null || echo "  (timeout alcanzado buscando anotaciones)"

# ================================================================
# 4. Limpiar LD directory viejo para re-run
# ================================================================
echo ""
echo "--- 4. Limpiando LD directory anterior ---"
LD_DIR="${OUT_BASE}/ldm_79k"
if [ -d "$LD_DIR" ]; then
    echo "  Removiendo $LD_DIR..."
    rm -rf "$LD_DIR"
    echo "  ✅ Limpiado"
fi

rm -rf "$TMPDIR"

echo ""
echo "============================================================"
echo "  FIX v2 COMPLETADO"
echo ""
echo "  Archivo .ma corregido: $MA_FINAL"
echo ""
echo "  Para re-ejecutar SBayesRC, cambiar en sbrc_02:"
echo "    MA_FILE=\"${MA_FINAL}\""
echo ""
echo "  Luego: sbatch sbrc_02_ldm_sbayesrc.sbatch"
echo "============================================================"
