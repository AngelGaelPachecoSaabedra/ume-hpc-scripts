#!/bin/bash
#SBATCH --job-name=index_bgen
#SBATCH --account=researchers
#SBATCH --qos=vip
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=/mnt/cephfs/orgs/home/angel.pacheco/logs/index_bgen_%j.out
#SBATCH --error=/mnt/cephfs/orgs/home/angel.pacheco/logs/index_bgen_%j.err

set -euo pipefail

# ============================================================
# RUTAS
# ============================================================
CONTAINER_DIR="/mnt/cephfs/orgs/home/angel.pacheco/containers"
SIF="${CONTAINER_DIR}/bgenix.sif"
BGEN_FILE="/mnt/cephfs/hot_nvme/mcps/imputed-topmed/bgen_files/MCPS_Freeze_150.GT_hg38.pVCF.rgcpid.QC2.TOPMED_dosages.bgen"

echo "============================================================"
echo "INICIANDO INDEXACIÓN CON BGENIX 1.1.7"
echo "Archivo: $BGEN_FILE"
echo "============================================================"

# Comprobación de seguridad (Ruta de micromamba)
echo "Comprobando dependencias en el contenedor..."
if ! apptainer exec "$SIF" test -x /opt/condaenv/bin/bgenix; then
    echo "ERROR: bgenix no encontrado en /opt/condaenv/bin/"
    exit 1
fi

# Ejecutar la indexación
# Usamos 'run' porque ya definiste el comando por defecto en el contenedor
apptainer run --bind /mnt/cephfs:/mnt/cephfs "$SIF" \
    -g "$BGEN_FILE" -index

echo "============================================================"
echo "INDEXACIÓN COMPLETADA"
echo "Verifica el archivo: ${BGEN_FILE}.bgi"
echo "============================================================"
