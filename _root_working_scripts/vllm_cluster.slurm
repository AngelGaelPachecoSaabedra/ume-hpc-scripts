#!/bin/bash
#SBATCH --job-name=vllm_ray_cluster
#SBATCH --nodes=2
#SBATCH --nodelist=MCPS-HPC[04,05]   # <--- LA ORDEN ESTRICTA PARA SLURM
#SBATCH --account=researchers
#SBATCH --qos=vip
#SBATCH --partition=gpu     
#SBATCH --ntasks-per-node=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00

# 1. Variables Críticas
export NCCL_SOCKET_IFNAME=eth10g
export GLOO_SOCKET_IFNAME=$NCCL_SOCKET_IFNAME
export HF_HOME="/mnt/cephfs/orgs/home/angel.pacheco/.cache/huggingface"

# Librerías y Ejecutable de Ray
export APPTAINERENV_PYTHONPATH="/mnt/cephfs/orgs/home/angel.pacheco/ray_env"
RAY_BIN="/mnt/cephfs/orgs/home/angel.pacheco/ray_env/bin/ray"
IMAGE="/mnt/cephfs/orgs/home/angel.pacheco/vllm-openai.sif"

# 2. Identificar el nodo maestro
head_node=$(scontrol show hostnames $SLURM_JOB_NODELIST | head -n 1)

# Leemos la IP directamente de la interfaz de fibra
head_node_ip=$(srun --nodes=1 --ntasks=1 -w "$head_node" ip -4 addr show $NCCL_SOCKET_IFNAME | awk '/inet / {print $2}' | cut -d/ -f1)
port=6379

echo "Iniciando Ray Head Node en $head_node ($head_node_ip)"

# 3. Levantar el Head Node de Ray (FORZANDO 1 GPU)
apptainer exec --nv --bind /mnt/cephfs:/mnt/cephfs $IMAGE \
    python3 $RAY_BIN start --head --node-ip-address=$head_node_ip --port=$port --num-gpus=1 --block &

sleep 15 

# 4. Levantar los Worker Nodes (FORZANDO 1 GPU)
worker_nodes=$(scontrol show hostnames $SLURM_JOB_NODELIST | tail -n +2)
for node in $worker_nodes; do
    echo "Conectando Ray Worker Node en $node al Head Node..."
    srun --nodes=1 --ntasks=1 -w "$node" \
        apptainer exec --nv --bind /mnt/cephfs:/mnt/cephfs $IMAGE \
        python3 $RAY_BIN start --address="$head_node_ip:$port" --num-gpus=1 --block &
done

sleep 15

# 5. Lanzar el servidor de inferencia de vLLM
echo "Desplegando modelo en vLLM (Pipeline Parallel Size = $SLURM_NNODES)..."

apptainer exec --nv --bind /mnt/cephfs:/mnt/cephfs $IMAGE \
    python3 -m vllm.entrypoints.openai.api_server \
    --model "NousResearch/Meta-Llama-3-8B-Instruct" \
    --tensor-parallel-size 1 \
    --pipeline-parallel-size $SLURM_NNODES \
    --distributed-executor-backend ray \
    --host 0.0.0.0 \
    --port 8000 \
    --download-dir $HF_HOME

wait
