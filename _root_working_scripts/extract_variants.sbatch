#!/bin/bash
#SBATCH --job-name=extract_vars
#SBATCH --account=researchers
#SBATCH --qos=vip
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=/mnt/cephfs/orgs/home/angel.pacheco/logs/extract_vars_%j.out
#SBATCH --error=/mnt/cephfs/orgs/home/angel.pacheco/logs/extract_vars_%j.err

set -euo pipefail

SIF="/mnt/cephfs/biocontainers/images/jupyter-biotools-1.2.sif"
CONDA_PYTHON="/opt/miniforge3/envs/vcf2zarr/bin/python"

ZARR_IN="/mnt/cephfs/hot_nvme/mcps/imputed-topmed/plink_to_zar/MCPS.HIS.TOPMED.b38.FR.2602.TARGET.imputed.plink.update2.zarr"
CSV_OUT="/mnt/cephfs/orgs/home/angel.pacheco/update2_variants_map.csv"

mkdir -p "/mnt/cephfs/orgs/home/angel.pacheco/logs"

# VOLVEMOS AL SCRATCH (Protegiendo el nodo)
BASE_TMP="/mnt/cephfs/scratch/$USER/plink_tmp"
TMPDIR="$BASE_TMP/$SLURM_JOB_ID"
mkdir -p "$TMPDIR"
chmod 700 "$TMPDIR"
export TMPDIR TMP="$TMPDIR" TEMP="$TMPDIR"
export APPTAINER_TMPDIR="$TMPDIR/apptainer-tmp"

echo "==== EXTRACCIÓN DE VARIANTES ZARR ===="
echo "Input: $ZARR_IN"
echo "Output: $CSV_OUT"
echo "TMP (Scratch): $TMPDIR"

PYTHON_SCRIPT="$TMPDIR/extract.py"
cat << 'EOF' > "$PYTHON_SCRIPT"
import sys
import zarr
import pandas as pd
import numpy as np
import time

zarr_path = sys.argv[1]
out_path = sys.argv[2]

print(f"[{time.strftime('%H:%M:%S')}] Abriendo Zarr: {zarr_path}")
try:
    root = zarr.open_group(zarr_path, mode='r')

    print(f"[{time.strftime('%H:%M:%S')}] Extrayendo vectores...")
    chrom_indices = root['variant_contig'][:]
    pos = root['variant_position'][:]
    
    # Extraemos los alelos
    alleles = root['variant_allele'][:]
    ref = alleles[:, 0]
    alt = alleles[:, 1]
    
    # ¡LA MAGIA DE LA TRADUCCIÓN!
    print(f"[{time.strftime('%H:%M:%S')}] Mapeando índices a nombres reales de cromosomas...")
    try:
        # Buscamos la lista maestra de nombres en el zarr
        contig_names = root['contig_id'][:]
        # Traducimos el arreglo de índices (0, 1, 2) a los nombres reales ('1', '2', 'X')
        chrom_real = np.array(contig_names)[chrom_indices]
    except KeyError:
        print("⚠️ Advertencia: No se encontró 'contig_id', usando índices crudos.")
        chrom_real = chrom_indices.astype(str)
    
    # Creamos el dataframe usando los nombres traducidos
    print(f"[{time.strftime('%H:%M:%S')}] Armando catálogo de variantes...")
    df = pd.DataFrame({
        'CHR': chrom_real,
        'POS': pos,
        'REF': ref.astype(str),
        'ALT': alt.astype(str)
    })
    
    print(f"[{time.strftime('%H:%M:%S')}] Guardando {len(df):,} variantes en {out_path}...")
    df.to_csv(out_path, index=False)
    print("¡Exito!")

except Exception as e:
    import traceback
    print("Error leyendo estructura Zarr:")
    traceback.print_exc()
    sys.exit(1)
EOF


# Ejecución (Con bind seguro a los temporales)
apptainer exec \
    --bind /mnt/cephfs:/mnt/cephfs \
    --bind "$TMPDIR:/tmp,$TMPDIR:/var/tmp" \
    "$SIF" \
    "$CONDA_PYTHON" "$PYTHON_SCRIPT" "$ZARR_IN" "$CSV_OUT"

rm -rf "$TMPDIR"
echo "Trabajo Terminado."
