#!/usr/bin/env bash
# download_clinical_genes.sh
# ==========================
# Download external clinical gene sources for the clinical genes pipeline.
#
# Downloads:
#   1. GenCC submissions CSV (primary: curated gene-disease associations)
#   2. ClinVar gene-condition-source mapping (optional, public NCBI FTP)
#
# Usage:
#   bash scripts/download_clinical_genes.sh
#
# Output:
#   resources/clinical_genes/raw/gencc_submissions.csv
#   resources/clinical_genes/raw/clinvar_gene_condition.txt  (optional)

set -euo pipefail

# Resolve paths relative to the script location (app root)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
APP_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
RAW_DIR="$APP_ROOT/resources/clinical_genes/raw"

mkdir -p "$RAW_DIR"

echo "=================================================================="
echo "  Clinical Genes — Source Download"
echo "  App root : $APP_ROOT"
echo "  Output   : $RAW_DIR"
echo "=================================================================="
echo ""

# ── 1. GenCC submissions CSV ──────────────────────────────────────────────────
GENCC_URL="https://search.thegencc.org/download/action/submissions-export-csv"
GENCC_OUT="$RAW_DIR/gencc_submissions.csv"

echo "[1/2] GenCC submissions CSV"
echo "  URL : $GENCC_URL"
echo "  Out : $GENCC_OUT"

if curl -fL --progress-bar \
        --connect-timeout 30 \
        --max-time 120 \
        "$GENCC_URL" \
        -o "$GENCC_OUT"; then
    NLINES=$(wc -l < "$GENCC_OUT" 2>/dev/null || echo "?")
    SIZE=$(du -sh "$GENCC_OUT" 2>/dev/null | cut -f1 || echo "?")
    echo "  OK — $SIZE, ~$NLINES lines"
else
    echo "  ERROR: GenCC download failed. Check connectivity."
    exit 1
fi

echo ""

# ── 2. ClinVar gene-condition source mapping (optional) ───────────────────────
# Public NCBI FTP — provides gene→condition mappings that complement GenCC.
CLINVAR_URL="https://ftp.ncbi.nlm.nih.gov/pub/clinvar/gene_condition_source_id"
CLINVAR_OUT="$RAW_DIR/clinvar_gene_condition.txt"

echo "[2/2] ClinVar gene-condition mapping (optional)"
echo "  URL : $CLINVAR_URL"
echo "  Out : $CLINVAR_OUT"

if curl -fL --progress-bar \
        --connect-timeout 30 \
        --max-time 180 \
        "$CLINVAR_URL" \
        -o "$CLINVAR_OUT" 2>/dev/null; then
    NLINES=$(wc -l < "$CLINVAR_OUT" 2>/dev/null || echo "?")
    SIZE=$(du -sh "$CLINVAR_OUT" 2>/dev/null | cut -f1 || echo "?")
    echo "  OK — $SIZE, ~$NLINES lines"
else
    echo "  SKIP — ClinVar download unavailable (network or file moved)"
    echo "         Builder will proceed without ClinVar source."
    rm -f "$CLINVAR_OUT"
fi

echo ""
echo "=================================================================="
echo "  Download complete."
echo ""
echo "  Files in $RAW_DIR:"
ls -lh "$RAW_DIR/" 2>/dev/null || echo "  (empty)"
echo ""
echo "  Next step:"
echo "    python3 src/services/clinical_genes_builder.py"
echo "=================================================================="
