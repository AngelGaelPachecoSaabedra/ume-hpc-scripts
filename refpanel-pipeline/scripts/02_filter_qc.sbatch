#!/bin/bash
#SBATCH --job-name=refpanel_02_filter
#SBATCH --account=researchers
#SBATCH --qos=vip
#SBATCH --partition=normal
#SBATCH --array=1-23
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=UNLIMITED
#SBATCH --output=/mnt/cephfs/orgs/home/angel.pacheco/logs/refpanel_02_filter_%A_%a.out
#SBATCH --error=/mnt/cephfs/orgs/home/angel.pacheco/logs/refpanel_02_filter_%A_%a.err
#SBATCH --mail-user=angelpachecosaavedra@gmail.com
#SBATCH --mail-type=BEGIN,END,FAIL

set -euo pipefail
source /mnt/cephfs/orgs/home/angel.pacheco/refpanel-pipeline/scripts/config.sh

CHR=$(get_chr_name ${SLURM_ARRAY_TASK_ID})
IN_BCF="${SCRATCH}/01_concat/chr${CHR}.norm.bcf"
OUTDIR="${SCRATCH}/02_filter"
OUT_BCF="${OUTDIR}/chr${CHR}.filtered.bcf"
QC_DIR="${PANEL_DIR}/qc"

echo "================================================================"
echo "  Paso 2: Filter + QC — chr${CHR}"
echo "  Job: ${SLURM_JOB_ID}, Node: ${SLURMD_NODENAME}"
echo "  Start: $(date)"
echo "================================================================"

mkdir -p "${OUTDIR}" "${QC_DIR}"
[ -f "${IN_BCF}" ] || { echo "ERROR: No existe ${IN_BCF}"; exit 1; }

if [ -f "${OUT_BCF}" ] && [ -f "${OUT_BCF}.csi" ]; then
    echo "SKIP: ya existe"
    exit 0
fi

T=${SLURM_CPUS_PER_TASK}

# Para chrX: no aplicar HWE
if [ "$CHR" == "X" ]; then
    HWE_EXPR="F_MISSING > ${MAX_MISS}"
    echo "  chrX: HWE filter DISABLED"
else
    HWE_EXPR="F_MISSING > ${MAX_MISS} || INFO/HWE < ${HWE_THRESH}"
fi

apptainer exec ${BIND} "${BIO_SIF}" bash -c "
    set -euo pipefail

    # 1) Set low-quality GTs to missing
    bcftools filter -S . -e 'FMT/GQ<${MIN_GQ} | FMT/DP<${MIN_DP}' --threads ${T} ${IN_BCF} -Ob -o ${OUTDIR}/chr${CHR}.tmp1.bcf

    # 2) Fill tags
    bcftools +fill-tags ${OUTDIR}/chr${CHR}.tmp1.bcf -Ob -o ${OUTDIR}/chr${CHR}.tmp2.bcf -- -t AN,AC,AF,F_MISSING,HWE

    # 3) Site filter
    bcftools view -e '${HWE_EXPR}' -m 2 -M 2 --threads ${T} ${OUTDIR}/chr${CHR}.tmp2.bcf -Ou | \
    bcftools view -e 'AC==0 || AC==AN' --threads ${T} -Ob -o ${OUT_BCF}

    # 4) Index
    bcftools index --threads ${T} ${OUT_BCF}

    # 5) QC stats
    bcftools query -f '%CHROM\t%POS\t%AF\t%F_MISSING\n' ${OUT_BCF} | gzip > ${QC_DIR}/chr${CHR}_stats.tsv.gz

    # Cleanup temps
    rm -f ${OUTDIR}/chr${CHR}.tmp1.bcf ${OUTDIR}/chr${CHR}.tmp2.bcf
"

N_IN=$(apptainer exec ${BIND} "${BIO_SIF}" bcftools view -H ${IN_BCF} | wc -l)
N_OUT=$(apptainer exec ${BIND} "${BIO_SIF}" bcftools view -H ${OUT_BCF} | wc -l)
echo ""
echo "  Input:  ${N_IN}"
echo "  Output: ${N_OUT} (removed $((N_IN - N_OUT)))"
echo "  Output: ${OUT_BCF} ($(ls -lh ${OUT_BCF} | awk '{print $5}'))"
echo "  Done: $(date)"
