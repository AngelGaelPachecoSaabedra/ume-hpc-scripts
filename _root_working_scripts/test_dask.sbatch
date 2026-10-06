#!/bin/bash
#SBATCH --job-name=test_dask
#SBATCH --account=researchers
#SBATCH --qos=normal
#SBATCH --partition=normal
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --time=00:15:00
#SBATCH --output=/mnt/cephfs/orgs/home/%u/logs/test_dask_%j.out
#SBATCH --error=/mnt/cephfs/orgs/home/%u/logs/test_dask_%j.err

set -euo pipefail

SIF=/mnt/cephfs/biocontainers/images/jupyter-biotools-1.5.sif
DATA=/mnt/cephfs/hot_nvme
OUT=/mnt/cephfs/orgs/home/$USER/dask_output
LOG=/mnt/cephfs/orgs/home/$USER/logs/dask_$SLURM_JOB_ID

mkdir -p "$OUT" "$LOG"

# Apuntamos al entorno:
PYTHON_EXEC="/opt/miniforge3/envs/jupyter/bin/python"

# Generar puertos aleatorios para evitar colisiones con zombies
SCHEDULER_NODE=$(scontrol show hostnames "$SLURM_JOB_NODELIST" | head -n1)
SCHEDULER_PORT=$(( 10000 + $RANDOM % 50000 ))
DASHBOARD_PORT=$(( SCHEDULER_PORT + 1 ))
SCHEDULER_ADDR="tcp://${SCHEDULER_NODE}:${SCHEDULER_PORT}"

echo "=========================================="
echo "Iniciando despliegue de Dask en Slurm"
echo "Contenedor: jupyter-biotools-1.4.sif (Entorno jax-gpu)"
echo "Nodos asignados: $SLURM_JOB_NODELIST"
echo "Scheduler Address: $SCHEDULER_ADDR"
echo "Dashboard Port: $DASHBOARD_PORT"
echo "=========================================="

# 1) Scheduler (Escuchando en 0.0.0.0)
srun --nodes=1 --ntasks=1 --nodelist="$SCHEDULER_NODE" \
  apptainer exec --bind /mnt/cephfs:/mnt/cephfs --bind "$DATA:/data:ro" --bind "$OUT:/out:rw" "$SIF" \
  $PYTHON_EXEC -m dask.distributed.cli.dask_scheduler \
    --host 0.0.0.0 \
    --port "$SCHEDULER_PORT" \
    --dashboard-address ":$DASHBOARD_PORT" &

sleep 5

# 2) Workers
srun --nodes=$SLURM_NNODES --ntasks=$SLURM_NNODES --ntasks-per-node=1 \
  apptainer exec --bind /mnt/cephfs:/mnt/cephfs --bind "$DATA:/data:ro" --bind "$OUT:/out:rw" "$SIF" \
  $PYTHON_EXEC -m dask.distributed.cli.dask_worker \
    "$SCHEDULER_ADDR" \
    --nthreads "$SLURM_CPUS_PER_TASK" \
    --memory-limit 8GB \
    --name "worker-\$(hostname)" &

sleep 8

# 3) TRABAJO DE PRUEBA
apptainer exec --bind /mnt/cephfs:/mnt/cephfs --bind "$DATA:/data:ro" --bind "$OUT:/out:rw" "$SIF" \
  $PYTHON_EXEC - << EOF
import time
from dask.distributed import Client
import dask.array as da

# Inyectando la variable calculada por Bash:
sched = "${SCHEDULER_ADDR}"
print(f"\n[1] Intentando conectar al Scheduler en: {sched}")
client = Client(sched)

print("[2] Conexión exitosa. Esperando a que los workers se reporten...")
client.wait_for_workers(n_workers=1)
info = client.scheduler_info()
workers = info['workers']

print(f"[3] ¡Excelente! Hay {len(workers)} workers activos:")
for worker_addr, worker_data in workers.items():
    print(f"    -> Nombre: {worker_data['name']} | Cores: {worker_data['nthreads']} | RAM: {worker_data['memory_limit'] / (1024**3):.2f} GB")

print("\n[4] Iniciando prueba de cálculo distribuido...")
start_time = time.time()

# Matriz gigante para forzar a los nodos a trabajar por la red
x = da.random.random((20000, 20000), chunks=(2000, 2000))
resultado = x.mean().compute()

end_time = time.time()
print(f"[5] Cómputo finalizado con éxito.")
print(f"    -> Resultado (Media esperada ~0.5): {resultado:.5f}")
print(f"    -> Tiempo de ejecución: {end_time - start_time:.2f} segundos")

client.close()
EOF

echo -e "\n=========================================="
echo "Dask job terminado. Limpiando procesos..."
jobs -p | xargs -r kill
echo "Listo."
