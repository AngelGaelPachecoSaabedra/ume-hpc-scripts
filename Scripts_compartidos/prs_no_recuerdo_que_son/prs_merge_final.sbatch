#!/bin/bash
#SBATCH --job-name=prs_merge
#SBATCH --account=researchers
#SBATCH --qos=vip
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:15:00
#SBATCH --output=/mnt/cephfs/orgs/home/angel.pacheco/logs/prs_merge_%j.out
#SBATCH --error=/mnt/cephfs/orgs/home/angel.pacheco/logs/prs_merge_%j.err

set -euo pipefail

SIF="/mnt/cephfs/biocontainers/images/jupyter-biotools-1.2.sif"

TMPDIR="/mnt/cephfs/scratch/$USER/merge_tmp/$SLURM_JOB_ID"
mkdir -p "$TMPDIR"
chmod 700 "$TMPDIR"
export TMPDIR TMP="$TMPDIR" TEMP="$TMPDIR"

echo "==== ETAPA 4: SUMA DE PRS POR CROMOSOMA ===="

PYTHON_SCRIPT="$TMPDIR/merge_prs.py"
cat << 'EOF' > "$PYTHON_SCRIPT"
import numpy as np
import pandas as pd
import zarr
import time
import os

print(f"[{time.strftime('%H:%M:%S')}] Sumando PRS de 22 cromosomas...")

prs_dir = "/mnt/cephfs/hot_nvme/mcps/imputed-topmed/prs_diabetes"
out_dir = "/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes"

# Verificar que existen los 22 archivos
missing = []
for c in range(1, 23):
    f = f"{prs_dir}/prs_chr_{c}.npy"
    if not os.path.exists(f):
        missing.append(c)

if missing:
    print(f"ERROR: Faltan cromosomas: {missing}")
    print("Asegúrate de que todos los jobs del array terminaron exitosamente.")
    exit(1)

# Sumar los 22 vectores
prs_total = None
for c in range(1, 23):
    prs_chr = np.load(f"{prs_dir}/prs_chr_{c}.npy")
    print(f"  chr{c:2d}: min={prs_chr.min():.6f}, max={prs_chr.max():.6f}, mean={prs_chr.mean():.6f}")
    if prs_total is None:
        prs_total = prs_chr.copy()
    else:
        prs_total += prs_chr

print(f"\n[{time.strftime('%H:%M:%S')}] PRS Total: min={prs_total.min():.6f}, max={prs_total.max():.6f}, mean={prs_total.mean():.6f}")

# Extraer IDs de muestras del Zarr
zarr_path = "/mnt/cephfs/hot_nvme/mcps/imputed-topmed/plink_to_zar/MCPS.HIS.TOPMED.b38.FR.2602.TARGET.imputed.plink.update2.zarr"
root = zarr.open_group(zarr_path, mode='r')

try:
    sample_ids = root['sample_id'][:]
    sample_ids = [s if isinstance(s, str) else s.decode() for s in sample_ids]
    print(f"[{time.strftime('%H:%M:%S')}] IDs de muestras extraídos del Zarr: {len(sample_ids)}")
except KeyError:
    print("Advertencia: No se encontró 'sample_id' en el Zarr, usando índices numéricos.")
    sample_ids = [f"sample_{i}" for i in range(len(prs_total))]

# Crear tabla final
df_result = pd.DataFrame({
    'SAMPLE_ID': sample_ids,
    'PRS_T2D': prs_total
})

# Ordenar por PRS descendente
df_result = df_result.sort_values('PRS_T2D', ascending=False)

# Agregar percentil
df_result['PERCENTILE'] = df_result['PRS_T2D'].rank(pct=True) * 100

# Guardar como TSV
out_tsv = f"{out_dir}/prs_diabetes_t2d_final.tsv"
df_result.to_csv(out_tsv, sep='\t', index=False, float_format='%.6f')
print(f"\n[{time.strftime('%H:%M:%S')}] Tabla final guardada: {out_tsv}")
print(f"  Muestras: {len(df_result):,}")
print(f"  Top 5 PRS más altos:")
print(df_result.head().to_string(index=False))
print(f"\n  Bottom 5 PRS más bajos:")
print(df_result.tail().to_string(index=False))

# También guardar el vector numpy por si se necesita downstream
np.save(f"{out_dir}/prs_t2d_total.npy", prs_total)

print(f"\n[{time.strftime('%H:%M:%S')}] PIPELINE PRS COMPLETADO.")
EOF

apptainer exec \
    --bind /mnt/cephfs:/mnt/cephfs \
    --bind "$TMPDIR:/tmp,$TMPDIR:/var/tmp" \
    "$SIF" \
    /opt/miniforge3/envs/vcf2zarr/bin/python "$PYTHON_SCRIPT"

rm -rf "$TMPDIR"
echo "==== MERGE PRS COMPLETADO ===="
