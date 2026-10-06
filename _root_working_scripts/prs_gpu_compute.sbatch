#!/bin/bash
#SBATCH --job-name=prs_gpu
#SBATCH --account=researchers
#SBATCH --partition=gpu
#SBATCH --qos=vip
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=100G
#SBATCH --gres=gpu:1           
#SBATCH --time=04:00:00
#SBATCH --array=1-22           
#SBATCH --exclude=mcps-hpc-01  # <-- ¡CRÍTICO! Protege el nodo maestro
#SBATCH --output=/mnt/cephfs/orgs/home/angel.pacheco/logs/prs_gpu_%A_%a.out
#SBATCH --error=/mnt/cephfs/orgs/home/angel.pacheco/logs/prs_gpu_%A_%a.err

set -euo pipefail

SIF="/mnt/cephfs/biocontainers/images/jupyter-biotools-1.2.sif"

# Aislamiento TMP en Scratch
BASE_TMP="/mnt/cephfs/scratch/$USER/gpu_tmp"
TMPDIR="$BASE_TMP/${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$TMPDIR"
chmod 700 "$TMPDIR"
export TMPDIR TMP="$TMPDIR" TEMP="$TMPDIR"
export APPTAINER_TMPDIR="$TMPDIR/apptainer-tmp"

CHR=$SLURM_ARRAY_TASK_ID

echo "============================================================"
echo "  Motor PRS GPU - Cromosoma $CHR"
echo "  Job Array ID:   ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
echo "  Nodo:           $SLURMD_NODENAME"
echo "============================================================"

PYTHON_SCRIPT="$TMPDIR/compute_prs.py"

cat << 'EOF' > "$PYTHON_SCRIPT"
import sys
import os
import zarr
import time
import pandas as pd
import numpy as np
import cupy as cp

chrom = str(sys.argv[1])
print(f"[{time.strftime('%H:%M:%S')}] Iniciando Cromosoma {chrom} en GPU...")

# 1. Cargar Pesos de CADD
parquet_path = "/mnt/cephfs/orgs/home/angel.pacheco/matched_weights.parquet"
df_weights = pd.read_parquet(parquet_path)
df_w_chr = df_weights[df_weights['CHR'] == chrom]

# 2. Cargar Metadata Zarr
zarr_path = "/mnt/cephfs/hot_nvme/mcps/imputed-topmed/plink_to_zar/MCPS.HIS.TOPMED.b38.FR.2602.TARGET.imputed.plink.update2.zarr"
root = zarr.open_group(zarr_path, mode='r')

contig_id_list = list(root['contig_id'][:])
chrom_idx = contig_id_list.index(chrom)

z_chroms = root['variant_contig'][:]
idx_start = np.searchsorted(z_chroms, chrom_idx, side='left')
idx_end = np.searchsorted(z_chroms, chrom_idx, side='right')

z_pos = root['variant_position'][idx_start:idx_end]
z_alleles = root['variant_allele'][idx_start:idx_end]

df_zarr = pd.DataFrame({'POS': z_pos, 'REF': z_alleles[:, 0], 'ALT': z_alleles[:, 1]})

# 3. ALINEACIÓN MAESTRA
print(f"[{time.strftime('%H:%M:%S')}] Alineando pesos con la matriz genómica...")
df_merged = pd.merge(df_zarr, df_w_chr, on=['POS', 'REF', 'ALT'], how='left')
df_merged['PHRED'] = df_merged['PHRED'].fillna(0.0)
weights_array = df_merged['PHRED'].values.astype(np.float32)

# 4. PROCESAMIENTO POR LOTES (ANTI-OOM)
n_variants = idx_end - idx_start
n_samples = root['call_genotype'].shape[1]
print(f"[{time.strftime('%H:%M:%S')}] Extrayendo {n_variants:,} variantes y {n_samples:,} pacientes en lotes...")

# Vector en GPU lleno de ceros para acumular el PRS final
prs_gpu_total = cp.zeros(n_samples, dtype=cp.float32)

# Cucharadas de 10,000 variantes (Pesa ~5.6 GB, dejando la mitad de la GPU libre para cálculos)
CHUNK_SIZE = 10000 

for start_v in range(idx_start, idx_end, CHUNK_SIZE):
    end_v = min(start_v + CHUNK_SIZE, idx_end)
    
    # Extraer pesos de este bloque exacto
    local_start = start_v - idx_start
    local_end = end_v - idx_start
    W_chunk = cp.asarray(weights_array[local_start:local_end])
    
    # Extraer genotipos a RAM
    G_raw = root['call_genotype'][start_v:end_v, :, :]
    G_sum = np.sum(np.maximum(G_raw, 0), axis=2, dtype=np.float32)
    
    # Inyectar a VRAM y multiplicar
    G_gpu = cp.asarray(G_sum)
    prs_gpu_total += cp.dot(G_gpu.T, W_chunk)
    
    # ¡VACIAR LA MEMORIA! (Vital para no crashear)
    del G_raw, G_sum, G_gpu, W_chunk
    cp.get_default_memory_pool().free_all_blocks()
    
    print(f"[{time.strftime('%H:%M:%S')}] Progreso: {local_end:,} / {n_variants:,} variantes procesadas.")

# 6. Guardar Resultados
print(f"[{time.strftime('%H:%M:%S')}] Descargando resultados al CPU...")
prs_cpu = cp.asnumpy(prs_gpu_total)

out_dir = "/mnt/cephfs/hot_nvme/mcps/imputed-topmed/prs_results"
os.makedirs(out_dir, exist_ok=True)
out_file = f"{out_dir}/prs_chr_{chrom}.npy"
np.save(out_file, prs_cpu)

print(f"[{time.strftime('%H:%M:%S')}] ¡ÉXITO! Scores guardados en {out_file}")
EOF

# Ejecución aislando el entorno jax-gpu e inyectando las librerías CUDA/cuBLAS
apptainer exec --nv \
    --bind /mnt/cephfs:/mnt/cephfs \
    --bind "$TMPDIR:/tmp,$TMPDIR:/var/tmp" \
    "$SIF" \
    bash -lc "
        JAX_LIBS=\$(find /opt/miniforge3/envs/jax-gpu/lib/python3.10/site-packages/nvidia/ -type d -name lib 2>/dev/null | tr '\n' ':')
        export LD_LIBRARY_PATH=\"/.singularity.d/libs:/opt/miniforge3/envs/jax-gpu/lib:\${JAX_LIBS}\${LD_LIBRARY_PATH:-}\"
        
        /opt/miniforge3/envs/jax-gpu/bin/python \"$PYTHON_SCRIPT\" \"$CHR\"
    "

rm -rf "$TMPDIR"
echo "Trabajo del Cromosoma $CHR Terminado."
