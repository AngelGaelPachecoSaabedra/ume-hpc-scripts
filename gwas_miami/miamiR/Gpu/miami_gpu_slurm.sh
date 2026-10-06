#!/bin/bash
#SBATCH --job-name=miami_gpu
#SBATCH --account=researchers
#SBATCH --partition=gpu
#SBATCH --qos=vip
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --gres=gpu:1
#SBATCH --time=04:00:00
#SBATCH --output=/mnt/cephfs/orgs/home/%u/logs/miami_gpu_%j.out
#SBATCH --error=/mnt/cephfs/orgs/home/%u/logs/miami_gpu_%j.err

set -euo pipefail

SIF=/mnt/cephfs/biocontainers/images/jupyter-biotools-1.2.sif
WORKDIR=/mnt/cephfs/orgs/home/$USER/gwas_miami
OUTDIR=$WORKDIR/results_$SLURM_JOB_ID

mkdir -p /mnt/cephfs/orgs/home/$USER/logs
mkdir -p "$OUTDIR"

echo "============================================================"
echo "  Miami Plot GPU — Job $SLURM_JOB_ID"
echo "============================================================"
echo "Node:     $SLURMD_NODENAME"
echo "WORKDIR:  $WORKDIR"
echo "Start:    $(date)"
echo ""

cd "$WORKDIR"

# FASE 1: Cómputo GPU
echo "========== FASE 1: Cómputo GPU =========="
apptainer exec --nv --bind /mnt/cephfs:/mnt/cephfs "$SIF" \
  bash -c '
    JAX_LIBS=$(find /opt/miniforge3/envs/jax-gpu/lib/python3.10/site-packages/nvidia/ -type d -name lib 2>/dev/null | tr "\n" ":")
    JUP_LIBS=$(find /opt/miniforge3/envs/jupyter/lib/python3.11/site-packages/nvidia/ -type d -name lib 2>/dev/null | tr "\n" ":")
    export LD_LIBRARY_PATH="/.singularity.d/libs:${JAX_LIBS}${JUP_LIBS}"
    # Correr script fase 1
    /opt/miniforge3/envs/jax-gpu/bin/python miami_phase1_gpu.py
  '

echo ""
ls -lh miami_intermediate_*.npz 2>/dev/null

# FASE 2: Plotting (Con instalación automática de adjustText)
echo ""
echo "========== FASE 2: Generación del Plot =========="
apptainer exec --bind /mnt/cephfs:/mnt/cephfs "$SIF" \
  bash -c '
    echo "📦 Instalando adjustText en espacio de usuario..."
    /opt/miniforge3/envs/jupyter/bin/pip install --user adjustText --quiet
    
    echo "📊 Generando gráfico..."
    /opt/miniforge3/envs/jupyter/bin/python miami_phase2_plot.py
  '

# Copiar resultados
echo ""
cp -v miami_plot_final.png "$OUTDIR/" 2>/dev/null || true
cp -v top_snps_gpu_annotated.csv "$OUTDIR/" 2>/dev/null || true

echo "Completado: $(date)"
