#!/bin/bash
# ================================================================
# Verifica qué panel LD usaste en el job SBayesRC 2932
# ¿UKBB descargado? ¿MCPS generado in-house? ¿Cuántos individuos?
# ================================================================

OUT_BASE="/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes/sbayesrc"
LD_DIR="${OUT_BASE}/ldm_79k"

echo "============================================================"
echo "  Verificando panel LD usado en SBayesRC"
echo "============================================================"

# 1. Existencia y tamaño
echo ""
echo "--- 1. Panel LD: estructura general ---"
if [[ -d "$LD_DIR" ]]; then
  du -sh "$LD_DIR"
  echo "  Total archivos: $(ls "$LD_DIR" | wc -l)"
  echo ""
  echo "  Primeros 10 archivos:"
  ls "$LD_DIR" | head -10
else
  echo "  ⚠ NO existe $LD_DIR"
fi

# 2. ldm.info: nos dice cuántos bloques y de dónde son
echo ""
echo "--- 2. ldm.info (metadata de bloques) ---"
if [[ -f "$LD_DIR/ldm.info" ]]; then
  echo "  Header:"
  head -1 "$LD_DIR/ldm.info"
  echo ""
  echo "  Primeros 3 bloques:"
  head -4 "$LD_DIR/ldm.info" | tail -3 | column -t
  echo ""
  echo "  Total de bloques: $(tail -n +2 "$LD_DIR/ldm.info" | wc -l)"
else
  echo "  No existe ldm.info — SBayesRC v0.2.6 puede no generarlo"
fi

# 3. snp.info: nos dice cuántas variantes
echo ""
echo "--- 3. snp.info (variantes en el panel) ---"
if [[ -f "$LD_DIR/snp.info" ]]; then
  echo "  Header:"
  head -1 "$LD_DIR/snp.info"
  echo ""
  echo "  Primeras 3 variantes:"
  head -4 "$LD_DIR/snp.info" | tail -3 | column -t
  echo ""
  N_SNP=$(tail -n +2 "$LD_DIR/snp.info" | wc -l)
  echo "  Total variantes: $(printf "%'d" $N_SNP)"
fi

# 4. Tamaño de un bloque de ejemplo
echo ""
echo "--- 4. Tamaño bloque de ejemplo (block1.eigen.bin) ---"
if [[ -f "$LD_DIR/block1.eigen.bin" ]]; then
  ls -lh "$LD_DIR/block1.eigen.bin"
fi

# 5. ld.sh — el script que generó el panel
echo ""
echo "--- 5. ld.sh (script de generación) ---"
LD_SH="${OUT_BASE}/ld.sh"
if [[ -f "$LD_SH" ]]; then
  echo "  Existe: $LD_SH"
  echo "  Líneas totales: $(wc -l < "$LD_SH")"
  echo ""
  echo "  Primeras 30 líneas (mostrará si menciona MCPS, UKBB o panel descargado):"
  echo "  ───────────────────────────────────────"
  head -30 "$LD_SH"
  echo "  ───────────────────────────────────────"
  echo ""
  echo "  Búsqueda de keywords:"
  grep -i "ukb\|biobank\|download\|wget\|curl" "$LD_SH" | head -5 || echo "    (sin menciones de UKBB/download)"
  echo ""
  grep -i "mcps\|maximally_unrelated\|freeze" "$LD_SH" | head -5 || echo "    (sin menciones de MCPS)"
else
  echo "  ⚠ NO existe $LD_SH — esto es raro, debería haberse generado"
fi

# 6. Buscar genotipos referenciados (probable .bfile)
echo ""
echo "--- 6. ¿Qué .bfile usó el panel LD? ---"
if [[ -f "$LD_SH" ]]; then
  grep -iE "bfile|--bfile" "$LD_SH" | head -3
fi

# 7. Number of individuals — esto es lo que define el "79k"
echo ""
echo "--- 7. Cuántos individuos tiene el panel LD ---"
# Si el panel es de MCPS, debería coincidir con el .fam de los maximally_unrelated
echo "  Individuos en MCPS (maximally_unrelated chr1):"
MCPS_FAM="/mnt/cephfs/hot_nvme/mcps/imputed-topmed/plink_files/maximally_unrelated/mcps-freeze150k_qcd_chr1_ivs.fam"
if [[ -f "$MCPS_FAM" ]]; then
  N_MCPS=$(wc -l < "$MCPS_FAM")
  echo "    $(printf "%'d" $N_MCPS) individuos en .fam de MCPS"
fi

# Si el panel UKBB descargado fue, los archivos típicos son ~22 MB para block1
# Si es MCPS con 79k indivs, los bloques serían más grandes

# 8. Buscar logs antiguos de generación del panel
echo ""
echo "--- 8. Logs antiguos de generación del panel ---"
find /mnt/cephfs/orgs/home/angel.pacheco/logs/ -name "sbrc_*" -mtime +7 2>/dev/null | head -5
find "$OUT_BASE" -name "*.log" 2>/dev/null | head -5

echo ""
echo "============================================================"
echo "  VERIFICACIÓN COMPLETA"
echo "============================================================"
