#!/bin/bash
#SBATCH --job-name=ddp_uv
#SBATCH --nodes=4
#SBATCH --ntasks=4
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --gres=gpu:1
#SBATCH --partition=gpu
#SBATCH --qos=vip
#SBATCH --time=24:00:00
#SBATCH --output=/mnt/cephfs/orgs/home/angel.pacheco/training_ddp/logs/ddp_%j.out
#SBATCH --error=/mnt/cephfs/orgs/home/angel.pacheco/training_ddp/logs/ddp_%j.err

TRAIN_DIR="/mnt/cephfs/orgs/home/angel.pacheco/training_ddp"
VENV="$TRAIN_DIR/.venv"
CSV_PATH="/mnt/cephfs/orgs/home/angel.pacheco/dataset_both_views2_clean.csv"
OUT_DIR="$TRAIN_DIR/resultados"
TAB_COLS="p_080607,p_080608,p_080609,p_080613,p_080614,p_080615,p_080616,p_080618,p_080619,p_080620,p_080624,p_080628,p_080629,cadera"

mkdir -p "$OUT_DIR" "$TRAIN_DIR/logs"

MASTER_ADDR=$(scontrol show hostnames $SLURM_JOB_NODELIST | head -n 1)

echo "Job: $SLURM_JOB_ID | Nodos: $SLURM_JOB_NODELIST | Master: $MASTER_ADDR"

srun --mpi=pmix bash -c '
export RANK=$SLURM_PROCID
export LOCAL_RANK=0
export WORLD_SIZE=4
export MASTER_ADDR='"$MASTER_ADDR"'
export MASTER_PORT=29500

# NCCL configuración conservadora para Ethernet
export NCCL_SOCKET_IFNAME=eno2
export NCCL_IB_DISABLE=1
export NCCL_P2P_DISABLE=1
export NCCL_NET_GDR_LEVEL=0
export NCCL_SHM_DISABLE=1
export NCCL_BUFFSIZE=2097152
export NCCL_PROTO=Simple
export NCCL_DEBUG=INFO

source '"$VENV"'/bin/activate

python3 '"$TRAIN_DIR"'/train_ddp_4gpu_real.py \
    --csv_path "'"$CSV_PATH"'" \
    --out_dir "'"$OUT_DIR"'" \
    --backbone efficientnet_b2 \
    --img_size 260 \
    --batch_size 8 \
    --epochs 35 \
    --lr 2e-4 \
    --dropout 0.35 \
    --workers 2 \
    --tab_cols "'"$TAB_COLS"'" \
    --w_peso 1 --w_altura 1 --w_ict 2.5
'
