#!/bin/bash
#SBATCH --job-name=prs_orchestrator
#SBATCH --account=researchers
#SBATCH --qos=vip
#SBATCH --partition=normal
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=100G
#SBATCH --time=04:00:00
#SBATCH --output=/mnt/cephfs/orgs/home/gen/repos/pgscat/logs/slurm/prs_dask_%j.out

set -euo pipefail

: "${PGS_ID:?[ERROR] PGS_ID no está definido}"

# Activar tu entorno
source /home/gen/micromamba/etc/profile.d/mamba.sh
micromamba activate vcfarr-hpc

SCRIPT_DIR="/mnt/cephfs/orgs/home/gen/repos/pgscat/scripts"
PYTHON_EXEC="python"

# Configurar Dask
SCHEDULER_NODE=$(scontrol show hostnames "$SLURM_JOB_NODELIST" | head -n1)
SCHEDULER_PORT=$(( 10000 + $RANDOM % 50000 ))
SCHEDULER_ADDR="tcp://${SCHEDULER_NODE}:${SCHEDULER_PORT}"

echo "======================================================="
echo "  PRS Pipeline Dask Orchestrator — PGS_ID: ${PGS_ID}"
echo "======================================================="

# 1) Iniciar Scheduler
srun --overlap --nodes=1 --ntasks=1 --nodelist="$SCHEDULER_NODE" \
  $PYTHON_EXEC -m dask.distributed.cli.dask_scheduler --port "$SCHEDULER_PORT" &
sleep 5

# 2) Iniciar Workers
# Iniciar Workers
srun --overlap --nodes=$SLURM_NNODES --ntasks=$SLURM_NNODES --ntasks-per-node=1 \
  $PYTHON_EXEC -m dask.distributed.cli.dask_worker "$SCHEDULER_ADDR" \
    --nthreads 4 --memory-limit 90GB &

sleep 10

# 3) Script Orquestador en Python
# 3) Script Orquestador en Python
cat << 'EOF' > dask_runner.py
import sys
import subprocess
from dask.distributed import Client

def procesar_cromosoma(chrom, pgs_id, script_dir):
    # Usar sys.executable garantiza que estemos en el micromamba vcfarr-hpc
    cmd = [
        sys.executable, f"{script_dir}/compute_prs.py",
        "--pgs-id", pgs_id,
        "--chrom", chrom,
        "--window-size", "10000"
    ]
    res = subprocess.run(cmd, capture_output=True, text=True)
    
    if res.returncode == 0:
        return f" chr{chrom} OK"
    else:
        # compute_prs.py escupe sus errores en stdout, capturamos ambos:
        error_msg = f" chr{chrom} ERROR (Code {res.returncode}):\n"
        if res.stdout:
            error_msg += f"--- STDOUT ---\n{res.stdout}\n"
        if res.stderr:
            error_msg += f"--- STDERR ---\n{res.stderr}\n"
        return error_msg

if __name__ == "__main__":
    client = Client(sys.argv[1])
    chroms = [str(i) for i in range(1, 23)] + ["X"]
    
    print(f"Repartiendo 23 cromosomas a {len(client.scheduler_info()['workers'])} workers...")
    futuros = client.map(procesar_cromosoma, chroms, pgs_id=sys.argv[2], script_dir=sys.argv[3])
    
    for r in client.gather(futuros):
        print(r)
    client.close()
EOF

# 4) Ejecutar
$PYTHON_EXEC dask_runner.py "$SCHEDULER_ADDR" "$PGS_ID" "$SCRIPT_DIR"

jobs -p | xargs -r kill
echo "Fin del procesamiento."
