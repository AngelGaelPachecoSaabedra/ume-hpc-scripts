#!/bin/bash
#SBATCH --job-name=reindex_vcf
#SBATCH --account=researchers
#SBATCH --qos=vip
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4          # 4 CPUs son perfectos para bcftools index
#SBATCH --mem=16G                  # No necesitamos 128G para indexar
#SBATCH --time=04:00:00            # Tiempo de sobra para 150k muestras
#SBATCH --array=1-23               # 1-22 para autosomas, 23 será para el Chr X
#SBATCH --output=/mnt/cephfs/hot_nvme/mcps/imputed-topmed/vcf_files/logs/reindex_%a_%j.out
#SBATCH --error=/mnt/cephfs/hot_nvme/mcps/imputed-topmed/vcf_files/logs/reindex_%a_%j.err

set -euo pipefail

# --- CONFIGURACIÓN ---
SIF="/mnt/cephfs/biocontainers/images/jupyter-biotools-1.2.sif"
VCF_DIR="/mnt/cephfs/hot_nvme/mcps/imputed-topmed/vcf_files"

# Crear carpeta de logs si no existe
mkdir -p "${VCF_DIR}/logs"

# --- LÓGICA DEL CROMOSOMA ---
# Si el array es 23, lo cambiamos a "X"
CHR=$SLURM_ARRAY_TASK_ID
if [ "$CHR" -eq 23 ]; then
    CHR="X"
fi

VCF_FILE="${VCF_DIR}/mcps-freeze150k_qcd_chr${CHR}.vcf.gz"

echo "=========================================================="
echo " Re-indexando Cromosoma: $CHR"
echo " Archivo: $VCF_FILE"
echo " CPUs asignados: $SLURM_CPUS_PER_TASK"
echo "=========================================================="

# Validar que el archivo existe
if [ ! -f "$VCF_FILE" ]; then
    echo "Error: No se encontró el archivo $VCF_FILE"
    exit 1
fi

# --- EJECUCIÓN ---
# Usamos --threads para que bcftools aproveche los 4 CPUs que pedimos
# -t crea el índice .tbi
# -f fuerza la sobrescritura del índice viejo o corrupto

apptainer exec --bind /mnt/cephfs:/mnt/cephfs "$SIF" \
    bcftools index --threads "$SLURM_CPUS_PER_TASK" -t -f "$VCF_FILE"

echo "=========================================================="
echo " Éxito: Índice reconstruido para Chr $CHR"
echo "=========================================================="
