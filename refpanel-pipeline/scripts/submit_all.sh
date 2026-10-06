#!/bin/bash
set -euo pipefail
DIR=/mnt/cephfs/orgs/home/angel.pacheco/refpanel-pipeline/scripts

echo "=== MCPS10k Reference Panel (chr1-22 + chrX) ==="
echo ""

JOB1=$(sbatch ${DIR}/01_concat_norm.sbatch | awk '{print $4}')
echo "  01 concat+norm:   Job ${JOB1}"

JOB2=$(sbatch --dependency=afterok:${JOB1} ${DIR}/02_filter_qc.sbatch | awk '{print $4}')
echo "  02 filter+qc:     Job ${JOB2} (after ${JOB1})"

JOB3=$(sbatch --dependency=afterok:${JOB2} ${DIR}/03_phase_common.sbatch | awk '{print $4}')
echo "  03 phase common:  Job ${JOB3} (after ${JOB2})"

JOB4=$(sbatch --dependency=afterok:${JOB3} ${DIR}/04_phase_rare.sbatch | awk '{print $4}')
echo "  04 phase rare:    Job ${JOB4} (after ${JOB3})"

JOB5=$(sbatch --dependency=afterok:${JOB4} ${DIR}/05_convert_formats.sbatch | awk '{print $4}')
echo "  05 convert:       Job ${JOB5} (after ${JOB4})"

echo ""
echo "  Monitor:  squeue -u \$USER"
echo "  Logs:     tail -f ~/logs/refpanel_01_concat_${JOB1}_*.out"
echo "  Validate: bash ${DIR}/06_validate.sh"
