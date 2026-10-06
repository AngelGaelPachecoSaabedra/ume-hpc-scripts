#!/bin/bash
#SBATCH --job-name=refpanel_01_concat
#SBATCH --account=researchers
#SBATCH --qos=vip
#SBATCH --partition=normal
#SBATCH --array=1-23
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=UNLIMITED
#SBATCH --output=/mnt/cephfs/orgs/home/angel.pacheco/logs/refpanel_01_concat_%A_%a.out
#SBATCH --error=/mnt/cephfs/orgs/home/angel.pacheco/logs/refpanel_01_concat_%A_%a.err
#SBATCH --mail-user=angelpachecosaavedra@gmail.com
#SBATCH --mail-type=BEGIN,END,FAIL

set -euo pipefail
source /mnt/cephfs/orgs/home/angel.pacheco/refpanel-pipeline/scripts/config.sh

CHR_NUM=${SLURM_ARRAY_TASK_ID}
CHR=$(get_chr_name ${CHR_NUM})
OUTDIR="${SCRATCH}/01_concat"
OUT_BCF="${OUTDIR}/chr${CHR}.norm.bcf"

echo "================================================================"
echo "  Paso 1: Concat + Normalize — chr${CHR}"
echo "  Job: ${SLURM_JOB_ID}, Node: ${SLURMD_NODENAME}"
echo "  Start: $(date)"
echo "================================================================"

mkdir -p "${OUTDIR}"

if [ -f "${OUT_BCF}" ] && [ -f "${OUT_BCF}.csi" ]; then
    echo "SKIP: ${OUT_BCF} ya existe"
    exit 0
fi

# Obtener chunks (usar X o número según corresponda)
IFS='|' read -ra CHUNK_NAMES <<< "${CHR_CHUNKS[$CHR]}"
CHUNKS=()
for NAME in "${CHUNK_NAMES[@]}"; do
    F="${PVCF_DIR}/${NAME}"
    [ -f "$F" ] || { echo "ERROR: No existe $F"; exit 1; }
    CHUNKS+=("$F")
    echo "  Chunk: $(basename $F) ($(ls -lh $F | awk '{print $5}'))"
done

T=${SLURM_CPUS_PER_TASK}

echo ""
echo "  concat → norm -m-any → norm --check-ref → BCF"
echo ""

apptainer exec ${BIND} "${BIO_SIF}" bash -c "
    set -euo pipefail
    bcftools concat --allow-overlaps --remove-duplicates --threads ${T} ${CHUNKS[*]} -Ou | \
    bcftools norm -m -any --threads ${T} -Ou | \
    bcftools norm --check-ref w -f ${REF_FASTA} --threads ${T} -Ob -o ${OUT_BCF}
    bcftools index --threads ${T} ${OUT_BCF}
"

echo ""
apptainer exec ${BIND} "${BIO_SIF}" bcftools stats ${OUT_BCF} | grep '^SN' | head -4
echo ""
echo "  Output: ${OUT_BCF} ($(ls -lh ${OUT_BCF} | awk '{print $5}'))"
echo "  Done: $(date)"
