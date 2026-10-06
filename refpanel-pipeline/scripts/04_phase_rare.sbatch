#!/bin/bash
#SBATCH --job-name=refpanel_04_phase_rare
#SBATCH --account=researchers
#SBATCH --qos=vip
#SBATCH --partition=normal
#SBATCH --array=1-23
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=UNLIMITED
#SBATCH --output=/mnt/cephfs/orgs/home/angel.pacheco/logs/refpanel_04_phase_rare_%A_%a.out
#SBATCH --error=/mnt/cephfs/orgs/home/angel.pacheco/logs/refpanel_04_phase_rare_%A_%a.err
#SBATCH --mail-user=angelpachecosaavedra@gmail.com
#SBATCH --mail-type=BEGIN,END,FAIL

set -euo pipefail
source /mnt/cephfs/orgs/home/angel.pacheco/refpanel-pipeline/scripts/config.sh

CHR=$(get_chr_name ${SLURM_ARRAY_TASK_ID})
IN_ALL="${SCRATCH}/02_filter/chr${CHR}.filtered.bcf"
IN_SCAFFOLD="${SCRATCH}/03_phase/chr${CHR}.common.phased.bcf"
OUTDIR="${SCRATCH}/03_phase"
OUT_BCF="${OUTDIR}/chr${CHR}.phased.bcf"
GMAP="${GMAP_DIR}/chr${CHR}.b38.gmap.gz"

echo "================================================================"
echo "  Paso 4: Phase Rare — chr${CHR}"
echo "  Job: ${SLURM_JOB_ID}, Node: ${SLURMD_NODENAME}"
echo "  Start: $(date)"
echo "================================================================"

[ -f "${IN_ALL}" ]      || { echo "ERROR: No existe ${IN_ALL}"; exit 1; }
[ -f "${IN_SCAFFOLD}" ] || { echo "ERROR: No existe ${IN_SCAFFOLD}"; exit 1; }

if [ -f "${OUT_BCF}" ] && [ -f "${OUT_BCF}.csi" ]; then
    echo "SKIP: ya existe"
    exit 0
fi

T=${SLURM_CPUS_PER_TASK}

echo -n "  Total variantes:  "; apptainer exec ${BIND} "${BIO_SIF}" bcftools view -H "${IN_ALL}" | wc -l
echo -n "  Scaffold comunes: "; apptainer exec ${BIND} "${BIO_SIF}" bcftools view -H "${IN_SCAFFOLD}" | wc -l

# --- chrX: haploides ---
EXTRA_ARGS=""
if [ "$CHR" == "X" ]; then
    HAPLOID_FILE="${OUTDIR}/chrX_haploid_males.txt"
    [ -f "${HAPLOID_FILE}" ] || awk -F'\t' 'NR>1 && $2==1 {print $1}' "${SEX_FILE}" > "${HAPLOID_FILE}"
    echo "  chrX: $(wc -l < ${HAPLOID_FILE}) haploid males"
    EXTRA_ARGS="--haploids ${HAPLOID_FILE}"
fi

echo ""

apptainer exec ${BIND} "${IMP_SIF}" \
    shapeit5_phase_rare \
        --input "${IN_ALL}" \
        --scaffold "${IN_SCAFFOLD}" \
        --map "${GMAP}" \
        --output "${OUT_BCF}" \
        --region "chr${CHR}" \
        --thread "${T}" \
        --log "${OUTDIR}/chr${CHR}.phased.log" \
        ${EXTRA_ARGS}


echo ""
echo -n "  Variantes faseadas: "
apptainer exec ${BIND} "${BIO_SIF}" bcftools view -H "${OUT_BCF}" | wc -l
echo "  Output: ${OUT_BCF} ($(ls -lh ${OUT_BCF} | awk '{print $5}'))"
echo "  Done: $(date)"
