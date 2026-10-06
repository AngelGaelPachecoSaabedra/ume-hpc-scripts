#!/bin/bash
#SBATCH --job-name=reindex_fix
#SBATCH --account=researchers
#SBATCH --qos=vip
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8        # Subimos a 8 CPUs para optimizar velocidad sin cuello de botella
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --array=1,2,21,22      # Freno de disco: Solo los 4 que fallaron, máximo 2 a la vez
#SBATCH --output=/mnt/cephfs/hot_nvme/mcps/imputed-topmed/vcf_files/logs/fix_%a_%j.out
#SBATCH --error=/mnt/cephfs/hot_nvme/mcps/imputed-topmed/vcf_files/logs/fix_%a_%j.err

# Modo estricto para evitar fallos silenciosos
set -euo pipefail

# --- SISTEMA DE LOGS BLINDADO ---
LOGDIR="/mnt/cephfs/orgs/home/$USER/logs"
mkdir -p "$LOGDIR"

# Redirigir stdout y stderr simultáneamente al archivo de log personal
# (Le agregué unos guiones bajos para que el nombre del archivo sea más legible)
exec > >(tee -a "$LOGDIR/${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID:-na}_${SLURMD_NODENAME}.run.log") 2>&1

# --- CONFIGURACIÓN ---
SIF="/mnt/cephfs/biocontainers/images/jupyter-biotools-1.2.sif"
VCF_DIR="/mnt/cephfs/hot_nvme/mcps/imputed-topmed/vcf_files"
CHR=$SLURM_ARRAY_TASK_ID
VCF_FILE="${VCF_DIR}/mcps-freeze150k_qcd_chr${CHR}.vcf.gz"

echo "=========================================================="
echo " Recuperando Índice - Cromosoma: $CHR"
echo " Archivo: $VCF_FILE"
echo " CPUs: $SLURM_CPUS_PER_TASK"
echo " Nodo: $SLURMD_NODENAME"
echo " Log local: $LOGDIR/${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID:-na}_${SLURMD_NODENAME}.run.log"
echo "=========================================================="

# Validar que el archivo existe
if [ ! -f "$VCF_FILE" ]; then
    echo "Error: No se encontró el archivo $VCF_FILE"
    exit 1
fi

# --- EJECUCIÓN ---
echo "Iniciando bcftools index..."
apptainer exec --bind /mnt/cephfs:/mnt/cephfs "$SIF" \
    bcftools index --threads "$SLURM_CPUS_PER_TASK" -t -f "$VCF_FILE"

echo "=========================================================="
echo " Éxito: Índice recuperado para Chr $CHR"
echo "=========================================================="
