#!/bin/bash
#SBATCH --job-name=vcf2zarr_chr10_multi
#SBATCH --account=researchers
#SBATCH --qos=vip
#SBATCH --partition=normal
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=50G
#SBATCH --time=24:00:00
#SBATCH --output=/mnt/cephfs/orgs/home/angel.pacheco/vcf2zarr/logs/vcf2zarr_multi_chr10_%j.out
#SBATCH --error=/mnt/cephfs/orgs/home/angel.pacheco/vcf2zarr/logs/vcf2zarr_multi_chr10_%j.err
#SBATCH -D /mnt/cephfs/orgs/home/angel.pacheco/vcf2zarr

set -euo pipefail

echo "======================================"
echo " VCF → Zarr con bio2zarr (multinodo)"
echo " Job ID:     $SLURM_JOB_ID"
echo " Nodos:      $SLURM_NNODES"
echo " Tareas:     $SLURM_NTASKS"
echo " CPUs/tarea: $SLURM_CPUS_PER_TASK"
echo " Inicio:     $(date)"
echo "======================================"
echo ""

# --- 1. Activar entorno ---
VENV_PATH=/mnt/cephfs/orgs/home/angel.pacheco/uv-envs/vcf2zarr
source "${VENV_PATH}/bin/activate"

echo "Verificando instalación..."
python - << 'EOF'
import bio2zarr
print(f"  ✓ bio2zarr importado, versión {bio2zarr.__version__}")
try:
    import cyvcf2
    print(f"  ✓ cyvcf2 importado, versión {cyvcf2.__version__}")
except Exception as e:
    print("  ✗ FALTA cyvcf2. Instala bio2zarr[vcf] en este entorno.")
    raise
EOF
echo ""

# --- 2. Rutas y parámetros ---
VCF_FILE="/mnt/cephfs/hot_nvme/mcps/imputed-topmed/vcf_files/mcps-freeze150k_qcd_chr10.vcf.gz"
ICF_PATH="/mnt/cephfs/orgs/home/angel.pacheco/vcf2zarr/icf_chr10_multi"
ZARR_OUTPUT="/mnt/cephfs/orgs/home/angel.pacheco/vcf2zarr/zarr_output_multi/mcps-freeze150k_qcd_chr10.zarr"

mkdir -p "$(dirname "$ICF_PATH")"
mkdir -p "$(dirname "$ZARR_OUTPUT")"

echo "Archivos:"
echo "  Input VCF: $VCF_FILE"
echo "  ICF:       $ICF_PATH"
echo "  Zarr:      $ZARR_OUTPUT"
echo ""

if [ ! -f "$VCF_FILE" ]; then
    echo "ERROR: VCF no encontrado: $VCF_FILE"
    exit 1
fi

echo "Tamaño VCF: $(du -h "$VCF_FILE" | cut -f1)"
echo ""

# Hilos por proceso (lo usa cyvcf2 / numexpr / OpenMP si aplican)
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export NUMEXPR_MAX_THREADS=$SLURM_CPUS_PER_TASK

TASKS=$SLURM_NTASKS

#######################################
# 3. Fase EXPLODE distribuida
#######################################

# Para explode usamos 8 particiones por tarea → 4*8 = 32 particiones
PARTITIONS_PER_TASK_EXP=8
NUM_PARTITIONS_EXP=$((TASKS * PARTITIONS_PER_TASK_EXP))

echo "Configuración distribuida (EXPLODE):"
echo "  Tareas              : $TASKS"
echo "  Particiones totales : $NUM_PARTITIONS_EXP"
echo "  Particiones/tarea   : $PARTITIONS_PER_TASK_EXP"
echo ""

echo "======================================"
echo " [EXPLODE] Fase distribuida"
echo "======================================"
echo ""

echo "[EXPLODE] Paso 1/3: dexplode-init"
vcf2zarr dexplode-init \
    "$VCF_FILE" \
    "$ICF_PATH" \
    --num-partitions "$NUM_PARTITIONS_EXP" \
    --force \
    --worker-processes "$SLURM_CPUS_PER_TASK" \
    --progress

echo ""
echo "[EXPLODE] Paso 2/3: dexplode-partition en paralelo (srun)..."
srun --ntasks=$TASKS --ntasks-per-node=1 bash -c '
set -euo pipefail
ICF_PATH="'"$ICF_PATH"'"
TASKS='"$TASKS"'
PARTITIONS_PER_TASK_EXP='"$PARTITIONS_PER_TASK_EXP"'

RANK=$SLURM_PROCID

START=$((RANK * PARTITIONS_PER_TASK_EXP))
END=$((START + PARTITIONS_PER_TASK_EXP - 1))

echo "[$(hostname)] Tarea $RANK procesando particiones $START..$END"

for P in $(seq $START $END); do
    echo "Tarea $RANK → dexplode-partition $P"
    vcf2zarr dexplode-partition "$ICF_PATH" "$P" --verbose
done
'

echo ""
echo "[EXPLODE] Paso 3/3: dexplode-finalise"
vcf2zarr dexplode-finalise "$ICF_PATH" --verbose

echo ""
echo "[EXPLODE] Completado."
echo ""

#######################################
# 4. Fase ENCODE distribuida
#######################################

echo "======================================"
echo " [ENCODE] Fase distribuida"
echo "======================================"
echo ""

# Parámetros de chunks
VARIANTS_CHUNK=10000
SAMPLES_CHUNK=1000

# Ejecutamos dencode-init y capturamos su salida para detectar num_partitions real
ENCODE_INIT_LOG=$(mktemp)

echo "[ENCODE] Paso 1/3: dencode-init"
vcf2zarr dencode-init \
    "$ICF_PATH" \
    "$ZARR_OUTPUT" \
    --num-partitions "$NUM_PARTITIONS_EXP" \
    --force \
    --variants-chunk-size "$VARIANTS_CHUNK" \
    --samples-chunk-size "$SAMPLES_CHUNK" \
    --progress | tee "$ENCODE_INIT_LOG"

# Extraer num_partitions de la salida
ENCODE_PARTITIONS=$(awk '/num_partitions/ {print $2}' "$ENCODE_INIT_LOG" | tail -n1)

if [ -z "${ENCODE_PARTITIONS:-}" ]; then
    echo "ERROR: no se pudo detectar num_partitions de dencode-init"
    cat "$ENCODE_INIT_LOG"
    rm -f "$ENCODE_INIT_LOG"
    exit 1
fi

rm -f "$ENCODE_INIT_LOG"

# Reparto de particiones de encode entre tareas
ENCODE_PARTITIONS_PER_TASK=$(( (ENCODE_PARTITIONS + TASKS - 1) / TASKS ))

echo ""
echo "Configuración distribuida (ENCODE):"
echo "  Particiones encode totales : $ENCODE_PARTITIONS"
echo "  Tareas                     : $TASKS"
echo "  Particiones/tarea (aprox.) : $ENCODE_PARTITIONS_PER_TASK"
echo ""

echo "[ENCODE] Paso 2/3: dencode-partition en paralelo (srun)..."
srun --ntasks=$TASKS --ntasks-per-node=1 bash -c '
set -euo pipefail
ZARR_OUTPUT="'"$ZARR_OUTPUT"'"
TASKS='"$TASKS"'
ENCODE_PARTITIONS='"$ENCODE_PARTITIONS"'
ENCODE_PARTITIONS_PER_TASK='"$ENCODE_PARTITIONS_PER_TASK"'

RANK=$SLURM_PROCID

START=$((RANK * ENCODE_PARTITIONS_PER_TASK))
if [ $START -ge $ENCODE_PARTITIONS ]; then
    echo "[$(hostname)] Tarea $RANK sin particiones de encode asignadas (start=$START >= $ENCODE_PARTITIONS)"
    exit 0
fi

END=$((START + ENCODE_PARTITIONS_PER_TASK - 1))
if [ $END -ge $ENCODE_PARTITIONS ]; then
    END=$((ENCODE_PARTITIONS - 1))
fi

echo "[$(hostname)] Tarea $RANK procesando particiones (encode) $START..$END"

for P in $(seq $START $END); do
    echo "Tarea $RANK → dencode-partition $P"
    vcf2zarr dencode-partition "$ZARR_OUTPUT" "$P" --verbose
done
'

echo ""
echo "[ENCODE] Paso 3/3: dencode-finalise"
vcf2zarr dencode-finalise "$ZARR_OUTPUT" --verbose --progress

echo ""
echo "======================================"
echo " Conversión multinodo COMPLETADA"
echo " Finalizado: $(date)"
echo "======================================"

deactivate
