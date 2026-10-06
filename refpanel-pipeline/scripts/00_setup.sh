#!/bin/bash
set -euo pipefail
source /mnt/cephfs/orgs/home/angel.pacheco/refpanel-pipeline/scripts/config.sh

echo "=== Setup: MCPS10k Reference Panel (chr1-22 + chrX) ==="

mkdir -p "${SCRATCH}"/{01_concat,02_filter,03_phase}
mkdir -p "${PANEL_DIR}"/{msav,bref3,vcf,imp5_chunks,qc}
mkdir -p "${LOG_DIR}"

echo "  Scratch: ${SCRATCH}/"
echo "  Output:  ${PANEL_DIR}/"
echo ""

# Containers
for SIF in "$BIO_SIF" "$IMP_SIF"; do
    [ -f "$SIF" ] && echo "  ✓ $(basename $SIF)" || echo "  ✗ FALTA: $SIF"
done

# pVCFs (chr1-22 + chrX = 46 chunks)
echo ""
MISSING=0
for CHR in $(seq 1 22) X; do
    IFS='|' read -ra NAMES <<< "${CHR_CHUNKS[$CHR]}"
    for NAME in "${NAMES[@]}"; do
        [ -f "${PVCF_DIR}/${NAME}" ] || { echo "  ✗ ${NAME}"; ((MISSING++)); }
    done
done
echo "  pVCF chunks: $((46 - MISSING))/46 presentes"

# Reference
echo ""
[ -f "$REF_FASTA" ] && echo "  ✓ hg38.fa" || echo "  ✗ FALTA: $REF_FASTA"
[ -f "${REF_FASTA}.fai" ] && echo "  ✓ hg38.fa.fai" || echo "  ⚠ .fai faltante"

# Sex file for chrX
echo ""
if [ -f "$SEX_FILE" ]; then
    N=$(wc -l < "$SEX_FILE")
    echo "  ✓ wgs_sex.txt (${N} lines)"
else
    echo "  ✗ FALTA: ${SEX_FILE} (necesario para chrX)"
fi

echo ""
df -h /mnt/cephfs/ | tail -1 | awk '{printf "  Espacio: %s disponible de %s\n", $4, $2}'

echo ""
echo "=== Listo. Lanza con: bash submit_all.sh ==="
