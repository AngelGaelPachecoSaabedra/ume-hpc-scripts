#!/bin/bash
#SBATCH --job-name=miami_R
#SBATCH --account=researchers
#SBATCH --partition=normal
#SBATCH --qos=vip
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=UNLIMITED
#SBATCH --output=/mnt/cephfs/orgs/home/%u/logs/miami_R_%j.out
#SBATCH --error=/mnt/cephfs/orgs/home/%u/logs/miami_R_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=angelpachecosaavedra@gmail.com

set -euo pipefail

SIF=/mnt/cephfs/biocontainers/images/jupyter-biotools-1.5.sif
WORKDIR=/mnt/cephfs/orgs/home/$USER/gwas_miami

mkdir -p /mnt/cephfs/orgs/home/$USER/logs

echo "============================================================"
echo "  Miami Plot R / Nature — Job $SLURM_JOB_ID"
echo "============================================================"
echo "Node:     $SLURMD_NODENAME"
echo "Start:    $(date)"
echo ""

cd "$WORKDIR"

apptainer exec --bind /mnt/cephfs:/mnt/cephfs "$SIF" \
  /opt/miniforge3/envs/r-bio/bin/Rscript miami_plot_nature.R

echo ""
ls -lh miami_t2d_nature.pdf top_snps_annotated.csv 2>/dev/null || true
echo "Completado: $(date)"
