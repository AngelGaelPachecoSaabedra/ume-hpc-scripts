#!/bin/bash
#SBATCH --job-name=refpanel_05_convert
#SBATCH --account=researchers
#SBATCH --qos=vip
#SBATCH --partition=normal
#SBATCH --array=1-23
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=UNLIMITED
#SBATCH --output=/mnt/cephfs/orgs/home/angel.pacheco/logs/refpanel_05_convert_%A_%a.out
#SBATCH --error=/mnt/cephfs/orgs/home/angel.pacheco/logs/refpanel_05_convert_%A_%a.err
#SBATCH --mail-user=angelpachecosaavedra@gmail.com
#SBATCH --mail-type=BEGIN,END,FAIL

set -euo pipefail
source /mnt/cephfs/orgs/home/angel.pacheco/refpanel-pipeline/scripts/config.sh

CHR=$(get_chr_name ${SLURM_ARRAY_TASK_ID})
IN_BCF="${SCRATCH}/03_phase/chr${CHR}.phased.bcf"

echo "================================================================"
echo "  Paso 5: Convert Formats — chr${CHR}"
echo "  Job: ${SLURM_JOB_ID}, Node: ${SLURMD_NODENAME}"
echo "  Start: $(date)"
echo "================================================================"

[ -f "${IN_BCF}" ] || { echo "ERROR: No existe ${IN_BCF}"; exit 1; }

mkdir -p "${PANEL_DIR}"/{vcf,msav,bref3,imp5_chunks}

T=${SLURM_CPUS_PER_TASK}

# 1. VCF.gz
VCF_OUT="${PANEL_DIR}/vcf/chr${CHR}.phased.vcf.gz"
if [ ! -f "${VCF_OUT}" ]; then
    echo "[1/4] BCF → VCF.gz..."
    apptainer exec ${BIND} "${BIO_SIF}" bash -c "
        bcftools view ${IN_BCF} --threads ${T} -Oz -o ${VCF_OUT}
        bcftools index --threads ${T} -t ${VCF_OUT}
    "
    echo "  ✓ $(ls -lh ${VCF_OUT} | awk '{print $5}')"
else
    echo "[1/4] SKIP VCF.gz"
fi

# 2. MSAV (Minimac4)
MSAV_OUT="${PANEL_DIR}/msav/chr${CHR}.msav"
if [ ! -f "${MSAV_OUT}" ]; then
    echo "[2/4] VCF → MSAV..."
    apptainer exec ${BIND} "${IMP_SIF}" \
        minimac4 --compress-reference "${VCF_OUT}" > "${MSAV_OUT}"
    echo "  ✓ $(ls -lh ${MSAV_OUT} | awk '{print $5}')"
else
    echo "[2/4] SKIP MSAV"
fi

# 3. bref3 (BEAGLE)
BREF3_OUT="${PANEL_DIR}/bref3/chr${CHR}.phased.bref3"
if [ ! -f "${BREF3_OUT}" ]; then
    echo "[3/4] VCF → bref3..."
    apptainer exec ${BIND} --env BREF3_MEM=16g "${IMP_SIF}" \
        bref3 "${VCF_OUT}" > "${BREF3_OUT}"
    echo "  ✓ $(ls -lh ${BREF3_OUT} | awk '{print $5}')"
else
    echo "[3/4] SKIP bref3"
fi

# 4. IMPUTE5 chunks
CHUNKS_OUT="${PANEL_DIR}/imp5_chunks/chr${CHR}.chunks.txt"
if [ ! -f "${CHUNKS_OUT}" ]; then
    echo "[4/4] IMPUTE5 chunks..."
    apptainer exec ${BIND} "${IMP_SIF}" \
        imp5Chunker \
            --h "${VCF_OUT}" \
            --r "chr${CHR}" \
            --o "${CHUNKS_OUT}" \
            --window-size 5000000 \
            --buffer-size 500000
    echo "  ✓ $(wc -l < ${CHUNKS_OUT}) chunks"
else
    echo "[4/4] SKIP chunks"
fi

echo ""
echo "=== chr${CHR} completo ==="
echo "  VCF:    $(ls -lh ${VCF_OUT} 2>/dev/null | awk '{print $5}')"
echo "  MSAV:   $(ls -lh ${MSAV_OUT} 2>/dev/null | awk '{print $5}')"
echo "  bref3:  $(ls -lh ${BREF3_OUT} 2>/dev/null | awk '{print $5}')"
echo "  Chunks: $(wc -l < ${CHUNKS_OUT} 2>/dev/null) regions"
echo "  Done: $(date)"
