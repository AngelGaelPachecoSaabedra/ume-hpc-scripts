#!/bin/bash
# ================================================================
# LDpred2 Recovery — LAUNCHER
# Lanza array 1-22 con 3 concurrentes (24 cores cada uno = 72 cores)
# ================================================================

#SBATCH --job-name=ldpred2_launcher
#SBATCH --account=researchers
#SBATCH --qos=vip
#SBATCH --partition=normal
#SBATCH --mem=1G
#SBATCH --cpus-per-task=1
#SBATCH --output=/mnt/cephfs/orgs/home/angel.pacheco/logs/ldpred2_launcher_%j.out

set -euo pipefail

# Verificar que los .sbk del 2931 existan antes de lanzar
PREV="/mnt/cephfs/scratch/angel.pacheco/ldpred2_hm3_tmp/2931"
missing=0
for chr in {1..22}; do
  if [[ ! -f "$PREV/corr_chr${chr}.sbk" ]]; then
    echo "FALTA: corr_chr${chr}.sbk"
    ((missing++)) || true
  fi
  if [[ ! -f "$PREV/chr${chr}.rds" ]]; then
    echo "FALTA: chr${chr}.rds"
    ((missing++)) || true
  fi
done

if [[ $missing -gt 0 ]]; then
  echo "ERROR: $missing archivos faltantes en $PREV"
  exit 1
fi

echo "✓ Los 22 .sbk + .rds existen en $PREV"

# Lanzar array 1-22, máximo 3 concurrentes (72 cores usados)
JOB_ARRAY=$(sbatch --parsable --array=1-22%3 \
    /mnt/cephfs/orgs/home/angel.pacheco/ldpred2_worker.sbatch)

echo "============================================================"
echo "  Array lanzado: job ID $JOB_ARRAY"
echo "  22 tareas, 3 concurrentes (24 cores × 3 = 72 cores)"
echo ""
echo "  Monitorear:"
echo "    squeue -u $USER"
echo "    ls /mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes/ldpred2/ldpred2_betas_chr*.rds"
echo ""
echo "  Logs:"
echo "    /mnt/cephfs/orgs/home/angel.pacheco/logs/ldpred2_worker_${JOB_ARRAY}_<chr>.out"
echo ""
echo "  Cuando todos terminen, mergear con:"
echo "    sbatch ldpred2_merge_final.sbatch"
echo "============================================================"
