#!/bin/bash
#SBATCH --job-name=reindex_bp
#SBATCH --account=researchers
#SBATCH --qos=vip
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --array=1,2,21,22%2  # Freno de disco Ceph: Máximo 2 a la vez

############################
# NIVEL 1 – OBLIGATORIO
# stdout / stderr completos por array
############################
#SBATCH --output=/mnt/cephfs/hot_nvme/mcps/imputed-topmed/vcf_files/logs/%x_%A_%a.out
#SBATCH --error=/mnt/cephfs/hot_nvme/mcps/imputed-topmed/vcf_files/logs/%x_%A_%a.err

############################
# NIVEL 2 – CAPTURAR TODO (Modo Paranoico)
############################
set -euo pipefail

LOGDIR="/mnt/cephfs/orgs/home/$USER/logs"
mkdir -p "$LOGDIR"

# Captura ABSOLUTAMENTE TODO lo que ocurra en tiempo real al log personal
exec > >(tee -a "$LOGDIR/${SLURM_JOB_NAME}_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID:-na}_${SLURMD_NODENAME}.run.log") 2>&1

############################
# NIVEL 3 – AISLAMIENTO DE TMP (Evita colapsos de /tmp local)
############################
# Crear un TMP único para este trabajo dentro del nodo
BASE_TMP="${SLURM_TMPDIR:-/tmp}"
TMPDIR="$BASE_TMP/$USER/$SLURM_JOB_ID"
mkdir -p "$TMPDIR"
chmod 700 "$TMPDIR"

export TMPDIR
export TMP="$TMPDIR"
export TEMP="$TMPDIR"

# Aislamiento agresivo de Apptainer para que no toque tu home ni el /tmp general
export APPTAINER_TMPDIR="$TMPDIR/apptainer-tmp"
export APPTAINER_CACHEDIR="$TMPDIR/apptainer-cache"
mkdir -p "$APPTAINER_TMPDIR" "$APPTAINER_CACHEDIR"

############################
# NIVEL 4 – CONTEXTO FORENSE
############################
echo "==== SLURM CONTEXT ===="
echo "JOBID=$SLURM_JOB_ID"
echo "ARRAY_TASK=${SLURM_ARRAY_TASK_ID:-none}"
echo "USER=$USER"
echo "HOST=$(hostname)"
echo "START=$(date -Is)"
echo "PWD=$PWD"
echo "SCRIPT=$0"
echo ""

echo "==== ENV (relevante) ===="
env | grep -E '^(PATH=|TMP|TEMP|APPTAINER|SINGULARITY|CONDA|PYTHON|LD_LIBRARY_PATH)' || true
echo ""

echo "==== LIMITS ===="
ulimit -a
echo ""

echo "==== FILESYSTEM ANTES DE EMPEZAR ===="
df -h /tmp /mnt/cephfs "$TMPDIR" | cat
echo ""

############################
# EJECUCIÓN DEL TRABAJO
############################
echo "==== JOB START ===="

SIF="/mnt/cephfs/biocontainers/images/jupyter-biotools-1.2.sif"
VCF_DIR="/mnt/cephfs/hot_nvme/mcps/imputed-topmed/vcf_files"
CHR=$SLURM_ARRAY_TASK_ID
VCF_FILE="${VCF_DIR}/mcps-freeze150k_qcd_chr${CHR}.vcf.gz"

if [ ! -f "$VCF_FILE" ]; then
    echo "Error forense: No se encontró el archivo $VCF_FILE"
    exit 1
fi

echo "Procesando Cromosoma $CHR..."

# Magia de Apptainer: Inyectamos nuestro TMPDIR seguro como si fuera el /tmp y /var/tmp del contenedor
apptainer exec \
  --bind /mnt/cephfs:/mnt/cephfs \
  --bind "$TMPDIR:/tmp,$TMPDIR:/var/tmp" \
  "$SIF" \
  bcftools index --threads "$SLURM_CPUS_PER_TASK" -t -f "$VCF_FILE"

echo "==== JOB END ===="
echo "END=$(date -Is)"

# Limpieza del directorio temporal para no dejar basura en el nodo
rm -rf "$TMPDIR"
echo "Limpieza de TMP completada."
