#!/bin/bash
#SBATCH --job-name=vcf2zarr_wgs_wes
#SBATCH --account=researchers
#SBATCH --qos=vip
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=150G
#SBATCH --array=1      # Array de 1 tarea (WES)
#SBATCH --output=/mnt/cephfs/hot_nvme/mcps/logs/vcf2zarr_complete_%a.out
#SBATCH --error=/mnt/cephfs/hot_nvme/mcps/logs/vcf2zarr_complete_%a.err

######################################################################
# NIVEL 1: CAPTURA FORENSE (Modo Paranoico)
######################################################################
set -euo pipefail

# Aseguramos que exista un directorio global para los logs forenses
LOGDIR_FORENSIC="/mnt/cephfs/hot_nvme/mcps/logs/forensics"
mkdir -p "$LOGDIR_FORENSIC"
exec > >(tee -a "$LOGDIR_FORENSIC/${SLURM_JOB_NAME}_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID:-na}_${SLURMD_NODENAME}.run.log") 2>&1

######################################################################
# NIVEL 2: AISLAMIENTO DE TMP EN SCRATCH NVME
######################################################################
# Definimos tu nueva ruta de scratch ultrarrápida
FAST_SCRATCH="/mnt/cephfs/m.2/scratch/angel.pacheco"
mkdir -p "$FAST_SCRATCH"

# Redirigimos todo el TMP a los discos rápidos en lugar del /tmp local del nodo
TMPDIR="$FAST_SCRATCH/tmp_jobs/$SLURM_JOB_ID"
mkdir -p "$TMPDIR"
chmod 700 "$TMPDIR"

export TMPDIR
export TMP="$TMPDIR"
export TEMP="$TMPDIR"

# Aislamiento de Apptainer dentro del NVMe
export APPTAINER_TMPDIR="$TMPDIR/apptainer-tmp"
export APPTAINER_CACHEDIR="$TMPDIR/apptainer-cache"
mkdir -p "$APPTAINER_TMPDIR" "$APPTAINER_CACHEDIR"

######################################################################
# NIVEL 3: CONFIGURACIÓN Y MAPEO DE VARIABLES (WES vs WGS)
######################################################################
SIF="/mnt/cephfs/biocontainers/images/jupyter-biotools-1.2.sif"

# Lógica para elegir entre WES y WGS basándose en el SLURM_ARRAY_TASK_ID
if [ "$SLURM_ARRAY_TASK_ID" -eq 1 ]; then
    DATA_TYPE="WES"
    BASE_DIR="/mnt/cephfs/hot_nvme/mcps/whole_exome_sequencing/oqfe-svm"
    VCF="${BASE_DIR}/wes-svm-filtered.vcf.gz"
elif [ "$SLURM_ARRAY_TASK_ID" -eq 2 ]; then
    DATA_TYPE="WGS"
    BASE_DIR="/mnt/cephfs/hot_nvme/mcps/whole_genome_sequencing/oqfe-svm"
    VCF="${BASE_DIR}/wgs-svm-filtered.vcf.gz"
else
    echo "Error forense: Índice de array $SLURM_ARRAY_TASK_ID no es válido (solo 1 o 2)."
    exit 1
fi

OUT_DIR="${BASE_DIR}/zar_files"
# MODIFICACIÓN CLAVE: La carpeta intermedia ahora vive en el NVMe
ICF="${FAST_SCRATCH}/intermediate/${DATA_TYPE}_icf"      
ZARR="${OUT_DIR}/${DATA_TYPE}-svm-filtered.zarr"    # Resultado final se queda en su ruta original

# Crear directorios de salida
mkdir -p "${FAST_SCRATCH}/intermediate"
mkdir -p "${OUT_DIR}/logs"
mkdir -p "${OUT_DIR}/intermediate" # Por si acaso otros scripts lo buscan

echo "==== SLURM CONTEXT ===="
echo "Host: $(hostname)"
echo "Memoria asignada: ${SLURM_MEM_PER_NODE:-unknown}M"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Directorio Scratch NVMe: $FAST_SCRATCH"
echo "Directorio temporal Apptainer: $TMPDIR"
echo "======================="

echo "=== Procesando dataset: $DATA_TYPE ==="
echo "Input: $VCF"
echo "Temp (NVMe): $ICF"
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
echo "[$(date -Is)] --> Iniciando EXPLODE para $DATA_TYPE..."
apptainer exec \
    --bind /mnt/cephfs:/mnt/cephfs \
    --bind "$TMPDIR:/tmp,$TMPDIR:/var/tmp" \
    "$SIF" \
    vcf2zarr explode "$VCF" "$ICF" \
    --worker-processes "$SLURM_CPUS_PER_TASK" \
    --progress

# --- PASO 2: ENCODE ---
echo "[$(date -Is)] --> Iniciando ENCODE para $DATA_TYPE..."
apptainer exec \
    --bind /mnt/cephfs:/mnt/cephfs \
    --bind "$TMPDIR:/tmp,$TMPDIR:/var/tmp" \
    "$SIF" \
    vcf2zarr encode "$ICF" "$ZARR" \
    --worker-processes "$SLURM_CPUS_PER_TASK" \
    --progress

# --- LIMPIEZA FINAL ---
echo "[$(date -Is)] --> Borrando carpetas intermedias en NVMe..."
rm -rf "$ICF"
rm -rf "$TMPDIR"

echo "=== COMPLETADO EXITOSAMENTE DATASET $DATA_TYPE ==="
