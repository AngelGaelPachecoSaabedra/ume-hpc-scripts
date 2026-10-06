#!/bin/bash
#SBATCH --job-name=plink2zarr_up2
#SBATCH --account=researchers
#SBATCH --qos=vip
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=110G
#SBATCH --time=06:00:00
#SBATCH --array=2
#SBATCH --output=/mnt/cephfs/hot_nvme/mcps/imputed-topmed/plink_to_zar/logs/plink_update_%a.out
#SBATCH --error=/mnt/cephfs/hot_nvme/mcps/imputed-topmed/plink_to_zar/logs/plink_update_%a.err

######################################################################
# NIVEL 1 & 2: FORENSE Y AISLAMIENTO DE TMP (VERSIÓN SCRATCH HDD)
######################################################################
set -euo pipefail

LOGDIR_FORENSIC="/mnt/cephfs/hot_nvme/mcps/imputed-topmed/plink_to_zar/logs/forensics"
mkdir -p "$LOGDIR_FORENSIC"
exec > >(tee -a "$LOGDIR_FORENSIC/${SLURM_JOB_NAME}_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}_${SLURMD_NODENAME}.run.log") 2>&1

# Movemos los Terabytes temporales a los discos mecánicos
BASE_TMP="/mnt/cephfs/scratch/$USER/plink_tmp"
TMPDIR="$BASE_TMP/$SLURM_JOB_ID"
mkdir -p "$TMPDIR"
chmod 700 "$TMPDIR"
export TMPDIR TMP="$TMPDIR" TEMP="$TMPDIR"

export APPTAINER_TMPDIR="$TMPDIR/apptainer-tmp"
export APPTAINER_CACHEDIR="$TMPDIR/apptainer-cache"
mkdir -p "$APPTAINER_TMPDIR" "$APPTAINER_CACHEDIR"

######################################################################
# NIVEL 3: CONFIGURACIÓN DE RUTAS Y ENTORNO
######################################################################
SIF="/mnt/cephfs/biocontainers/images/jupyter-biotools-1.2.sif"

BASE_DIR="/mnt/cephfs/hot_nvme/mcps/imputed-topmed"
PLINK_DIR="$BASE_DIR/plink-prs"
ZARR_DIR="${BASE_DIR}/plink_to_zar"

mkdir -p "$ZARR_DIR"

UPDATE_NUM=$SLURM_ARRAY_TASK_ID
PREFIX="MCPS.HIS.TOPMED.b38.FR.2602.TARGET.imputed.plink.update${UPDATE_NUM}"

PLINK_PREFIX="${PLINK_DIR}/${PREFIX}"
ZARR_OUT="${ZARR_DIR}/${PREFIX}.zarr"

echo "==== INICIANDO CONVERSIÓN PLINK -> ZARR ===="
echo "Fecha inicio:  $(date -Is)"
echo "Update Num:    $UPDATE_NUM"
echo "Input PREFIX:  $PLINK_PREFIX"
echo "Output Zarr:   $ZARR_OUT"
echo "Directorio TMP: $TMPDIR"
echo "Nodo:          ${SLURMD_NODENAME}"
echo "CPUs asignados: ${SLURM_CPUS_PER_TASK}"

# Validar que los tres archivos plink existan
for ext in bed bim fam; do
    if [ ! -f "${PLINK_PREFIX}.${ext}" ]; then
        echo "ERROR: No se encontró ${PLINK_PREFIX}.${ext}"
        exit 1
    fi
done
echo "Archivos .bed/.bim/.fam verificados"

# Borrar intentos viejos si existen
rm -rf "$ZARR_OUT"

######################################################################
# NIVEL 4: EJECUCIÓN CON HERRAMIENTA NATIVA DE BIO2ZARR
######################################################################

# --- PASO ÚNICO: CONVERT ---
echo "[$(date -Is)] --> Iniciando CONVERT directo de PLINK a Zarr..."
apptainer exec \
    --bind /mnt/cephfs:/mnt/cephfs \
    --bind "$TMPDIR:/tmp,$TMPDIR:/var/tmp" \
    "$SIF" \
    /opt/miniforge3/envs/vcf2zarr/bin/plink2zarr convert "$PLINK_PREFIX" "$ZARR_OUT" \
    --worker-processes "$SLURM_CPUS_PER_TASK"

EXIT_CODE=$?

# Limpieza
echo "[$(date -Is)] --> Limpiando 5.5TB de basura en scratch..."
rm -rf "$TMPDIR"

if [ $EXIT_CODE -eq 0 ]; then
    echo "==== TRABAJO COMPLETADO EXITOSAMENTE ===="
    echo "Fecha fin: $(date -Is)"
    echo "Zarr disponible en: $ZARR_OUT"
else
    echo "==== TRABAJO FALLIDO (exit $EXIT_CODE) ===="
    exit $EXIT_CODE
fi
