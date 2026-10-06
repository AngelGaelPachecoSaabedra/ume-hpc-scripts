#!/bin/bash
#SBATCH --job-name=check_x_zarr
#SBATCH --account=researchers
#SBATCH --qos=vip
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:15:00
#SBATCH --output=/mnt/cephfs/orgs/home/angel.pacheco/logs/check_x_%j.out
#SBATCH --error=/mnt/cephfs/orgs/home/angel.pacheco/logs/check_x_%j.err

set -euo pipefail

# Definir el contenedor Apptainer que usas
SIF="/mnt/cephfs/biocontainers/images/jupyter-biotools-1.2.sif"

TMPDIR="/mnt/cephfs/scratch/$USER/check_x_tmp/$SLURM_JOB_ID"
mkdir -p "$TMPDIR"
chmod 700 "$TMPDIR"
export TMPDIR TMP="$TMPDIR" TEMP="$TMPDIR"

PYTHON_SCRIPT="$TMPDIR/check_x.py"
cat << 'PYEOF' > "$PYTHON_SCRIPT"
import zarr
import numpy as np
import pandas as pd

# --- RUTAS ---
zarr_path = "/mnt/cephfs/hot_nvme/mcps/imputed-topmed/plink_to_zar/MCPS.HIS.TOPMED.b38.FR.2602.TARGET.imputed.plink.update2.zarr"
covar_file = "/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes/covar_ctrl_45.txt"

print("--- Iniciando Inspección de Cromosoma X ---")

# 1. Cargar Sexo
print("Cargando metadatos de sexo...")
df_covar = pd.read_csv(covar_file, sep='\t')
sex_map = dict(zip(df_covar['IID'].astype(str), df_covar['SEX']))

# 2. Abrir Zarr
print(f"Abriendo Zarr: {zarr_path}")
root = zarr.open_group(zarr_path, mode='r')

# 3. Mapear Muestras
sample_ids = root['sample_id'][:]
sample_ids = [s if isinstance(s, str) else s.decode() for s in sample_ids]

idx_males = []
idx_females = []

for i, sid in enumerate(sample_ids):
    s_sex = sex_map.get(sid, -1)
    if s_sex == 1 and len(idx_males) < 5:
        idx_males.append(i)
    elif s_sex == 0 and len(idx_females) < 5:
        idx_females.append(i)
    
    if len(idx_males) >= 5 and len(idx_females) >= 5:
        break

# 4. Localizar Cromosoma X
contigs = list(root['contig_id'][:])
if 'X' in contigs:
    x_idx = contigs.index('X')
elif '23' in contigs:
    x_idx = contigs.index('23')
else:
    raise ValueError("No se encontró cromosoma X o 23 en contig_id")

var_contig = root['variant_contig'][:]
start_x = np.searchsorted(var_contig, x_idx, side='left')

# Tomamos una variante un poco más adelante del inicio
test_idx = start_x + 1000 
pos = root['variant_position'][test_idx]
ref = root['variant_allele'][test_idx, 0]
alt = root['variant_allele'][test_idx, 1]

if isinstance(ref, bytes): ref = ref.decode()
if isinstance(alt, bytes): alt = alt.decode()

print(f"\n--- Inspeccionando Variante X:{pos} ({ref}/{alt}) ---")

gt_raw = root['call_genotype'][test_idx] 

print("\nResultados MUJERES (Esperado: 0/0, 0/1, 1/1):")
for idx in idx_females:
    raw = gt_raw[idx]
    suma = np.sum(np.maximum(raw, 0))
    print(f"  ID {sample_ids[idx]}: Raw={raw}, Suma={suma}")

print("\nResultados HOMBRES (CRÍTICO: ¿Es [1, -1] o [1, 1]?):")
for idx in idx_males:
    raw = gt_raw[idx]
    suma = np.sum(np.maximum(raw, 0))
    print(f"  ID {sample_ids[idx]}: Raw={raw}, Suma={suma}")

print("\n--- Fin del chequeo ---")
PYEOF

# Ejecutar dentro del contenedor
apptainer exec \
    --bind /mnt/cephfs:/mnt/cephfs \
    --bind "$TMPDIR:/tmp,$TMPDIR:/var/tmp" \
    "$SIF" \
    /opt/miniforge3/envs/vcf2zarr/bin/python "$PYTHON_SCRIPT"

rm -rf "$TMPDIR"
