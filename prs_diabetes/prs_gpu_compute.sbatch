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
#SBATCH --exclude=mcps-hpc-01
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
echo "  PRS Diabetes T2 (PGS000014) - GPU - Cromosoma $CHR"
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

chrom = sys.argv[1]
chrom_int = int(chrom)
chrom_str = str(chrom)
print(f"[{time.strftime('%H:%M:%S')}] Iniciando PRS Cromosoma {chrom} en GPU...")

# ================================================================
# 1. Cargar pesos PRS del Parquet
#    CHR en el Parquet es int32 (viene de Spark).
# ================================================================
parquet_path = "/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes/matched_prs_weights.parquet"
df_weights = pd.read_parquet(parquet_path)
df_w_chr = df_weights[df_weights['CHR'] == chrom_int].copy()
print(f"[{time.strftime('%H:%M:%S')}] Pesos PRS para chr{chrom}: {len(df_w_chr):,} variantes")

# ================================================================
# 2. Cargar metadata del Zarr para este cromosoma
# ================================================================
zarr_path = "/mnt/cephfs/hot_nvme/mcps/imputed-topmed/plink_to_zar/MCPS.HIS.TOPMED.b38.FR.2602.TARGET.imputed.plink.update2.zarr"
root = zarr.open_group(zarr_path, mode='r')

contig_id_list = list(root['contig_id'][:])
chrom_idx = contig_id_list.index(chrom_str)

z_chroms = root['variant_contig'][:]
idx_start = np.searchsorted(z_chroms, chrom_idx, side='left')
idx_end   = np.searchsorted(z_chroms, chrom_idx, side='right')

n_variants = idx_end - idx_start
n_samples  = root['call_genotype'].shape[1]
print(f"[{time.strftime('%H:%M:%S')}] Zarr chr{chrom}: {n_variants:,} variantes, {n_samples:,} muestras")

z_pos     = root['variant_position'][idx_start:idx_end]
z_alleles = root['variant_allele'][idx_start:idx_end]

df_zarr = pd.DataFrame({
    'POS': z_pos,
    'REF': z_alleles[:, 0].astype(str),
    'ALT': z_alleles[:, 1].astype(str)
})

# ================================================================
# 3. LEFT JOIN: Alinear pesos con el orden del Zarr
#    Preserva TODAS las filas del Zarr (tabla izquierda).
#    Variantes sin peso PRS reciben BETA=0 y FLIP=False.
# ================================================================
print(f"[{time.strftime('%H:%M:%S')}] Alineando pesos PRS con la matriz genómica (LEFT JOIN)...")

df_merged = pd.merge(
    df_zarr,
    df_w_chr[['POS', 'REF', 'ALT', 'BETA', 'FLIP']],
    on=['POS', 'REF', 'ALT'],
    how='left'
)

# Deduplicar: si una variante matcheó como directa y como invertida
# (variantes palindrómicas como A/T o C/G), conservar el match directo (FLIP=False)
if len(df_merged) > len(df_zarr):
    n_extra = len(df_merged) - len(df_zarr)
    print(f"[{time.strftime('%H:%M:%S')}] Detectados {n_extra} duplicados palindrómicos, conservando match directo...")
    df_merged['_orig_order'] = range(len(df_merged))
    df_merged = df_merged.sort_values(['_orig_order', 'FLIP'], ascending=[True, True]) \
                         .drop_duplicates(subset=['POS', 'REF', 'ALT'], keep='first') \
                         .sort_values('_orig_order') \
                         .drop(columns='_orig_order') \
                         .reset_index(drop=True)

# Validación: el merge debe tener el mismo tamaño que el Zarr
assert len(df_merged) == len(df_zarr), \
    f"ERROR: El merge cambió el tamaño! Zarr: {len(df_zarr)}, Merged: {len(df_merged)}. Posibles duplicados en pesos."

df_merged['BETA'] = df_merged['BETA'].fillna(0.0)
df_merged['FLIP'] = df_merged['FLIP'].fillna(False)

n_matched = (df_merged['BETA'] != 0.0).sum()
n_flipped = df_merged['FLIP'].sum()
print(f"[{time.strftime('%H:%M:%S')}] Variantes con peso PRS: {n_matched:,} / {n_variants:,}")
print(f"[{time.strftime('%H:%M:%S')}] Variantes con alelo invertido (FLIP): {n_flipped:,}")

weights_array = df_merged['BETA'].values.astype(np.float32)
flip_array    = df_merged['FLIP'].values

# ================================================================
# 4. CÁLCULO PRS EN GPU POR LOTES
#
# Para variantes normales:  PRS += G * BETA
# Para variantes FLIP:      PRS += (2 - G) * BETA
#
# FLIP significa que el effect_allele es el REF, no el ALT.
# El genotipo cuenta alelos ALT, así que para contar alelos
# del effect_allele (que aquí es REF) necesitamos invertir: 2 - G.
# ================================================================
print(f"[{time.strftime('%H:%M:%S')}] Calculando PRS en GPU por lotes...")

prs_gpu_total = cp.zeros(n_samples, dtype=cp.float32)
CHUNK_SIZE = 10000

for start_v in range(idx_start, idx_end, CHUNK_SIZE):
    end_v = min(start_v + CHUNK_SIZE, idx_end)
    local_start = start_v - idx_start
    local_end   = end_v - idx_start

    # Pesos de este chunk
    W_chunk    = cp.asarray(weights_array[local_start:local_end])
    flip_chunk = flip_array[local_start:local_end]

    # Genotipos: (variantes, muestras, 2 alelos) -> sumar alelos ALT
    G_raw = root['call_genotype'][start_v:end_v, :, :]
    G_sum = np.sum(np.maximum(G_raw, 0), axis=2, dtype=np.float32)

    # Invertir genotipos donde FLIP=True
    # G_sum cuenta alelos ALT. Si el effect_allele es REF,
    # necesitamos contar alelos REF = 2 - alelos ALT
    if flip_chunk.any():
        G_sum[flip_chunk] = 2.0 - G_sum[flip_chunk]

    # Subir a GPU y multiplicar
    G_gpu = cp.asarray(G_sum)
    prs_gpu_total += cp.dot(G_gpu.T, W_chunk)

    # Liberar memoria
    del G_raw, G_sum, G_gpu, W_chunk
    cp.get_default_memory_pool().free_all_blocks()

    print(f"[{time.strftime('%H:%M:%S')}] Progreso: {local_end:,} / {n_variants:,} variantes")

# ================================================================
# 5. Guardar resultado
# ================================================================
print(f"[{time.strftime('%H:%M:%S')}] Descargando resultados al CPU...")
prs_cpu = cp.asnumpy(prs_gpu_total)

out_dir = "/mnt/cephfs/hot_nvme/mcps/imputed-topmed/prs_diabetes"
os.makedirs(out_dir, exist_ok=True)
out_file = f"{out_dir}/prs_chr_{chrom}.npy"
np.save(out_file, prs_cpu)

print(f"[{time.strftime('%H:%M:%S')}] ÉXITO! PRS chr{chrom} guardado en {out_file}")
print(f"[{time.strftime('%H:%M:%S')}] Estadísticas: min={prs_cpu.min():.6f}, max={prs_cpu.max():.6f}, mean={prs_cpu.mean():.6f}")
EOF

# Ejecución con GPU
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
echo "PRS Cromosoma $CHR Terminado."
