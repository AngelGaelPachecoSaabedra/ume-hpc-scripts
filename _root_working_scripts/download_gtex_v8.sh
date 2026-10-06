#!/bin/bash
# =============================================================================
# Descarga completa de GTEx Analysis v8 QTLs
# Bucket: gs://gtex-resources/GTEx_Analysis_v8_QTLs/
# Tamaño total: ~2.23 TiB
# Costo estimado de egress: ~$270 USD
# =============================================================================

PROJECT_ID="project-e5e73628-a322-418b-aee"
BASE_GCS="gs://gtex-resources/GTEx_Analysis_v8_QTLs"
BASE_LOCAL="/mnt/cephfs/hot_nvme/gtex/GTEx_v8_QTLs"

# Crear estructura de carpetas
CARPETAS=(
  "eQTL_all_associations"
  "EUR_eQTL_all_associations"
  "EUR_eQTL_covariates"
  "EUR_eQTL_expression_matrices"
  "EUR_sQTL"
  "EUR_sQTL_all_associations"
  "EUR_sQTL_covariates"
  "EUR_sQTL_splicing_phenotype_groups"
  "EUR_sQTL_splicing_phenotypes"
  "sQTL"
  "sQTL_all_associations"
)

echo "============================================="
echo " GTEx v8 QTLs - Descarga completa"
echo " Destino: $BASE_LOCAL"
echo " Costo estimado: ~\$270 USD (egress)"
echo "============================================="
echo ""

# Crear directorio base
mkdir -p "$BASE_LOCAL"

# 1. Descargar el archivo tar suelto
echo "[1/12] Descargando archivo suelto: GTEx_Analysis_v8_sQTL_leafcutter_counts.tar"
gsutil -m -u "$PROJECT_ID" cp \
  "${BASE_GCS}/GTEx_Analysis_v8_sQTL_leafcutter_counts.tar" \
  "$BASE_LOCAL/"
echo ""

# 2. Descargar cada carpeta
COUNTER=2
for carpeta in "${CARPETAS[@]}"; do
  GCS_PATH="${BASE_GCS}/GTEx_Analysis_v8_${carpeta}/"
  LOCAL_PATH="${BASE_LOCAL}/${carpeta}"

  echo "[${COUNTER}/12] Descargando: ${carpeta}"
  echo "  Origen:  ${GCS_PATH}"
  echo "  Destino: ${LOCAL_PATH}"

  mkdir -p "$LOCAL_PATH"
  gsutil -m -u "$PROJECT_ID" cp -r "${GCS_PATH}*" "$LOCAL_PATH/"

  if [ $? -eq 0 ]; then
    echo "  ✓ Completado: ${carpeta}"
  else
    echo "  ✗ ERROR en: ${carpeta}"
  fi
  echo ""
  COUNTER=$((COUNTER + 1))
done

# 3. Resumen final
echo "============================================="
echo " Descarga finalizada"
echo " Ubicación: $BASE_LOCAL"
echo "============================================="
echo ""
echo "Estructura:"
du -sh "$BASE_LOCAL"/*/ 2>/dev/null
du -sh "$BASE_LOCAL"/*.tar 2>/dev/null
echo ""
echo "Total:"
du -sh "$BASE_LOCAL"
