#!/bin/bash
# ================================================================
# PIPELINE C v6 — LAUNCHER
# Lanza array 1-23 (22 autosomas + chrX) con 3 concurrentes
# REQUIERE: ldp2_custom_v6_prune.sbatch completado primero
# ================================================================

#SBATCH --job-name=ldp2_custom_launcher
#SBATCH --account=researchers
#SBATCH --qos=vip
#SBATCH --partition=normal
#SBATCH --mem=1G
#SBATCH --cpus-per-task=1
#SBATCH --output=/mnt/cephfs/orgs/home/angel.pacheco/logs/ldp2_custom_launcher_%j.out

set -euo pipefail

OUT_DIR="/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes/ldpred2_custom"
PRUNED_FILE="$OUT_DIR/sumstats_pruned.rds"

# Validar que el pruning se haya completado
if [[ ! -f "$PRUNED_FILE" ]]; then
  echo "ERROR: No existe $PRUNED_FILE"
  echo "Primero corre: sbatch ldp2_custom_v6_prune.sbatch"
  exit 1
fi
echo "✓ sumstats_pruned.rds existe"

# Lanzar array 1-23 (22 autosomas + chrX=23), máximo 3 concurrentes
JOB_ARRAY=$(sbatch --parsable --array=1-23%3 \
    /mnt/cephfs/orgs/home/angel.pacheco/LDPRED2_CUSTOM/ldp2_custom_v6_worker.sbatch)

echo "============================================================"
echo "  Array lanzado: job ID $JOB_ARRAY"
echo "  23 tareas (chr1-22 + chrX como chr23), 3 concurrentes"
echo "  24 cores × 3 = 72 cores usados"
echo ""
echo "  Monitorear:"
echo "    squeue -u $USER"
echo "    ls $OUT_DIR/ldpred2_custom_betas_chr*.rds | wc -l"
echo ""
echo "  Logs:"
echo "    /mnt/cephfs/orgs/home/angel.pacheco/logs/ldp2_custom_${JOB_ARRAY}_<chr>.out"
echo ""
echo "  Cuando todos terminen, mergear con:"
echo "    sbatch ldp2_custom_v6_merge.sbatch"
echo "============================================================"
