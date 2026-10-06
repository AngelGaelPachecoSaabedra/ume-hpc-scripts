#!/bin/bash
#SBATCH --job-name=match_4fuentes
#SBATCH --account=researchers
#SBATCH --qos=vip
#SBATCH --nodes=4
#SBATCH --ntasks=4
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=48G
#SBATCH --time=02:00:00
#SBATCH --output=/mnt/cephfs/orgs/home/angel.pacheco/espirometrias/4fuentes_%j.out
#SBATCH --error=/mnt/cephfs/orgs/home/angel.pacheco/espirometrias/4fuentes_%j.err
#SBATCH -D /mnt/cephfs/orgs/home/angel.pacheco/espirometrias

echo "=============================================="
echo "MATCHING 4 FUENTES - 4 NODOS"
echo "=============================================="
echo "Job ID:        $SLURM_JOB_ID"
echo "Nodos:         $SLURM_JOB_NODELIST"
echo "CPUs/nodo:     $SLURM_CPUS_PER_TASK"
echo "Total CPUs:    $((SLURM_NTASKS * SLURM_CPUS_PER_TASK))"
echo "Inicio:        $(date)"
echo "=============================================="

cd /mnt/cephfs/orgs/home/angel.pacheco/espirometrias
source .venv/bin/activate

# Crear hostfile
scontrol show hostnames $SLURM_JOB_NODELIST > hostfile_$SLURM_JOB_ID
echo "Nodos:"
cat hostfile_$SLURM_JOB_ID

# Ejecutar con mpirun
mpirun --hostfile hostfile_$SLURM_JOB_ID \
       --map-by node \
       -np $SLURM_NTASKS \
       -x PATH \
       -x VIRTUAL_ENV \
       python matching_4fuentes_mpi.py

rm -f hostfile_$SLURM_JOB_ID

echo "=============================================="
echo "Fin: $(date)"
echo "=============================================="
