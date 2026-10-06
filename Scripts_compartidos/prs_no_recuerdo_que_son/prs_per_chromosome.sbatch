#!/bin/bash
#SBATCH --job-name=prs_perchr
#SBATCH --account=researchers
#SBATCH --qos=vip
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:15:00
#SBATCH --output=/mnt/cephfs/orgs/home/angel.pacheco/logs/prs_perchr_%j.out
#SBATCH --error=/mnt/cephfs/orgs/home/angel.pacheco/logs/prs_perchr_%j.err

set -euo pipefail

SIF="/mnt/cephfs/biocontainers/images/jupyter-biotools-1.2.sif"

TMPDIR="/mnt/cephfs/scratch/$USER/merge_tmp/$SLURM_JOB_ID"
mkdir -p "$TMPDIR"
chmod 700 "$TMPDIR"
export TMPDIR TMP="$TMPDIR" TEMP="$TMPDIR"

echo "==== PRS POR CROMOSOMA (SIN SUMAR) ===="

PYTHON_SCRIPT="$TMPDIR/prs_perchr.py"
cat << 'EOF' > "$PYTHON_SCRIPT"
import numpy as np
import pandas as pd
import zarr
import time

print(f"[{time.strftime('%H:%M:%S')}] Generando tabla PRS por cromosoma...")

prs_dir = "/mnt/cephfs/hot_nvme/mcps/imputed-topmed/prs_diabetes"
out_dir = "/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes"

# Extraer IDs de muestras del Zarr
zarr_path = "/mnt/cephfs/hot_nvme/mcps/imputed-topmed/plink_to_zar/MCPS.HIS.TOPMED.b38.FR.2602.TARGET.imputed.plink.update2.zarr"
root = zarr.open_group(zarr_path, mode='r')

try:
    sample_ids = root['sample_id'][:]
    sample_ids = [s if isinstance(s, str) else s.decode() for s in sample_ids]
except KeyError:
    sample_ids = [f"sample_{i}" for i in range(140831)]

# Construir tabla con una columna por cromosoma
df = pd.DataFrame({'SAMPLE_ID': sample_ids})

for c in range(1, 23):
    prs_chr = np.load(f"{prs_dir}/prs_chr_{c}.npy")
    df[f'PRS_chr{c}'] = prs_chr
    print(f"  chr{c:2d}: min={prs_chr.min():.6f}, max={prs_chr.max():.6f}, mean={prs_chr.mean():.6f}")

# Agregar columna de PRS total y percentil
chr_cols = [f'PRS_chr{c}' for c in range(1, 23)]
df['PRS_TOTAL'] = df[chr_cols].sum(axis=1)
df['PERCENTILE'] = df['PRS_TOTAL'].rank(pct=True) * 100

# Ordenar por PRS total descendente
df = df.sort_values('PRS_TOTAL', ascending=False)

out_tsv = f"{out_dir}/prs_diabetes_per_chromosome.tsv"
df.to_csv(out_tsv, sep='\t', index=False, float_format='%.6f')

print(f"\n[{time.strftime('%H:%M:%S')}] Tabla guardada: {out_tsv}")
print(f"  Muestras: {len(df):,}")
print(f"  Columnas: SAMPLE_ID + 22 cromosomas + PRS_TOTAL")
EOF

apptainer exec \
    --bind /mnt/cephfs:/mnt/cephfs \
    --bind "$TMPDIR:/tmp,$TMPDIR:/var/tmp" \
    "$SIF" \
    /opt/miniforge3/envs/vcf2zarr/bin/python "$PYTHON_SCRIPT"

rm -rf "$TMPDIR"
echo "==== COMPLETADO ===="
