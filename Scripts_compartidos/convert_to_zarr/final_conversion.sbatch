#!/bin/bash
#SBATCH --job-name=vcf2zarr_prod
#SBATCH --account=researchers
#SBATCH --qos=vip
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --array=1,2,21,22  # El 23 será el ChrX. Límite de 3 simultáneos para cuidar CephFS
#SBATCH --output=/mnt/cephfs/hot_nvme/mcps/imputed-topmed/zar_files/logs/chr_%a.out
#SBATCH --error=/mnt/cephfs/hot_nvme/mcps/imputed-topmed/zar_files/logs/chr_%a.err

######################################################################
# NIVEL 1: CAPTURA FORENSE (Modo Paranoico)
######################################################################
set -euo pipefail

# Captura ABSOLUTAMENTE TODO lo que ocurra en tiempo real al log personal
LOGDIR_FORENSIC="/mnt/cephfs/hot_nvme/mcps/imputed-topmed/zar_files/logs/forensics"
mkdir -p "$LOGDIR_FORENSIC"
exec > >(tee -a "$LOGDIR_FORENSIC/${SLURM_JOB_NAME}_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID:-na}_${SLURMD_NODENAME}.run.log") 2>&1

######################################################################
# NIVEL 2: AISLAMIENTO DE TMP (Evita colapsos en la RAM/Disco local)
######################################################################
BASE_TMP="${SLURM_TMPDIR:-/tmp}"
TMPDIR="$BASE_TMP/$USER/$SLURM_JOB_ID"
mkdir -p "$TMPDIR"
chmod 700 "$TMPDIR"

export TMPDIR
export TMP="$TMPDIR"
export TEMP="$TMPDIR"

# Aislamiento de Apptainer
export APPTAINER_TMPDIR="$TMPDIR/apptainer-tmp"
export APPTAINER_CACHEDIR="$TMPDIR/apptainer-cache"
mkdir -p "$APPTAINER_TMPDIR" "$APPTAINER_CACHEDIR"

######################################################################
# NIVEL 3: CONFIGURACIÓN Y MAPEO DE VARIABLES
######################################################################
SIF="/mnt/cephfs/biocontainers/images/jupyter-biotools-1.2.sif"
BASE_DIR="/mnt/cephfs/hot_nvme/mcps/imputed-topmed"
VCF_DIR="${BASE_DIR}/vcf_files"
OUT_DIR="${BASE_DIR}/zar_files"

# Lógica para convertir el Job 23 en Cromosoma X
CHR_NUM=$SLURM_ARRAY_TASK_ID
if [ "$CHR_NUM" -eq 23 ]; then
    CHR_NAME="X"
else
    CHR_NAME="$CHR_NUM"
fi

VCF="${VCF_DIR}/mcps-freeze150k_qcd_chr${CHR_NAME}.vcf.gz"
ICF="${OUT_DIR}/intermediate/chr${CHR_NAME}_icf"  # Carpeta temporal Zarr
ZARR="${OUT_DIR}/chr${CHR_NAME}.zarr"             # Resultado final

mkdir -p "${OUT_DIR}/intermediate"
mkdir -p "${OUT_DIR}/logs"

echo "==== SLURM CONTEXT ===="
echo "Host: $(hostname)"
echo "Memoria asignada: ${SLURM_MEM_PER_NODE:-unknown}M"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Directorio temporal seguro: $TMPDIR"
df -h /tmp /mnt/cephfs "$TMPDIR" | cat
echo "======================="

echo "=== Procesando Cromosoma $CHR_NAME ==="
echo "Input: $VCF"
echo "Temp:  $ICF"
echo "Final: $ZARR"

# Validar que existe el VCF
if [ ! -f "$VCF" ]; then
    echo "Error forense: No encuentro el archivo $VCF"
    exit 1
fi

# Limpieza previa por si hubo intentos fallidos anteriores
rm -rf "$ICF" "$ZARR"

######################################################################
# NIVEL 4: EJECUCIÓN DEL PIPELINE
######################################################################

# --- PASO 1: EXPLODE ---
echo "[$(date -Is)] --> Iniciando EXPLODE..."
apptainer exec \
    --bind /mnt/cephfs:/mnt/cephfs \
    --bind "$TMPDIR:/tmp,$TMPDIR:/var/tmp" \
    "$SIF" \
    vcf2zarr explode "$VCF" "$ICF" \
    --worker-processes "$SLURM_CPUS_PER_TASK" \
    --progress

# --- PASO 2: ENCODE ---
echo "[$(date -Is)] --> Iniciando ENCODE..."
apptainer exec \
    --bind /mnt/cephfs:/mnt/cephfs \
    --bind "$TMPDIR:/tmp,$TMPDIR:/var/tmp" \
    "$SIF" \
    vcf2zarr encode "$ICF" "$ZARR" \
    --worker-processes "$SLURM_CPUS_PER_TASK" \
    --progress

# --- LIMPIEZA FINAL ---
echo "[$(date -Is)] --> Borrando carpetas intermedias y TMP local..."
rm -rf "$ICF"
rm -rf "$TMPDIR"

echo "=== COMPLETADO EXITOSAMENTE CHR $CHR_NAME ==="
