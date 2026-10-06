#!/bin/bash
set -euo pipefail
source /mnt/cephfs/orgs/home/angel.pacheco/refpanel-pipeline/scripts/config.sh

echo "=== Validación: MCPS10k Reference Panel (chr1-22 + chrX) ==="
echo ""

PASS=0; FAIL=0; TOTAL_V=0

# 1. Archivos
echo "=== 1. Archivos ==="
for CHR in $(seq 1 22) X; do
    OK=true
    for F in "vcf/chr${CHR}.phased.vcf.gz" "msav/chr${CHR}.msav" "bref3/chr${CHR}.phased.bref3"; do
        [ -f "${PANEL_DIR}/${F}" ] || { echo "  ✗ FALTA: ${F}"; OK=false; ((FAIL++)); }
    done
    $OK && ((PASS++))
done
echo "  ${PASS}/23 cromosomas completos"
echo ""

# 2. Faseo (spot check)
echo "=== 2. Faseo ==="
for CHR in 1 11 22 X; do
    VCF="${PANEL_DIR}/vcf/chr${CHR}.phased.vcf.gz"
    [ -f "$VCF" ] || continue
    UNP=$(apptainer exec ${BIND} "${BIO_SIF}" bcftools query -f '[%GT\t]\n' "$VCF" | head -3000 | grep -c '/' || true)
    [ "$UNP" -eq 0 ] && echo "  chr${CHR}: ✓" || echo "  chr${CHR}: ✗ ${UNP} unphased"
done
echo ""

# 3. Conteos
echo "=== 3. Variantes ==="
printf "  %-6s %12s\n" "Chr" "Variants"
for CHR in $(seq 1 22) X; do
    VCF="${PANEL_DIR}/vcf/chr${CHR}.phased.vcf.gz"
    [ -f "$VCF" ] || continue
    N=$(apptainer exec ${BIND} "${BIO_SIF}" bcftools view -H "$VCF" | wc -l)
    printf "  chr%-3s %'12d\n" "$CHR" "$N"
    TOTAL_V=$((TOTAL_V + N))
done
printf "\n  TOTAL  %'12d\n" "$TOTAL_V"

echo ""
echo "=== Panel: ${PANEL_DIR}/ ==="
echo "  Samples: 10,008 | Build: GRCh38 | Passed: ${PASS} | Failed: ${FAIL}"
