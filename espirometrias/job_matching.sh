#!/bin/bash
#SBATCH --job-name=matching_espiro
#SBATCH --account=researchers
#SBATCH --qos=vip
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=/mnt/cephfs/orgs/home/angel.pacheco/espirometrias/matching_%j.out
#SBATCH --error=/mnt/cephfs/orgs/home/angel.pacheco/espirometrias/matching_%j.err
#SBATCH -D /mnt/cephfs/orgs/home/angel.pacheco/espirometrias

echo "=============================================="
echo "MATCHING AUTOMÁTICO DE ESPIROMETRÍAS"
echo "=============================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Nodo: $SLURM_NODELIST"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Memoria: $SLURM_MEM_PER_NODE MB"
echo "Inicio: $(date)"
echo "=============================================="

# Activar entorno virtual
source /home/angel.pacheco/.local/bin/env
cd /mnt/cephfs/orgs/home/angel.pacheco/espirometrias

# Activar el venv del proyecto
source .venv/bin/activate

# Ejecutar el script paralelizado usando todos los CPUs asignados
python matching_paralelo.py --cpus $SLURM_CPUS_PER_TASK

echo ""
echo "=============================================="
echo "Fin: $(date)"
echo "=============================================="
