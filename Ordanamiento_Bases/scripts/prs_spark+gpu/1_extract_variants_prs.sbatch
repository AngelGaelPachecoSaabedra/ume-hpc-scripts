#!/bin/bash
#SBATCH --job-name=extract_prs
#SBATCH --account=researchers
#SBATCH --qos=vip
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=/mnt/cephfs/orgs/home/angel.pacheco/logs/extract_prs_%j.out
#SBATCH --error=/mnt/cephfs/orgs/home/angel.pacheco/logs/extract_prs_%j.err

set -euo pipefail

SIF="/mnt/cephfs/biocontainers/images/jupyter-biotools-1.2.sif"
CONDA_PYTHON="/opt/miniforge3/envs/vcf2zarr/bin/python"

ZARR_IN="/mnt/cephfs/hot_nvme/mcps/imputed-topmed/plink_to_zar/MCPS.HIS.TOPMED.b38.FR.2602.TARGET.imputed.plink.update2.zarr"
CSV_OUT="/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes/update2_variants_map.csv"

# Archivo PGS Catalog (ya descargado y descomprimido)
PGS_RAW="/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes/PGS000014_hmPOS_GRCh38.txt"
PGS_CLEAN="/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes/PGS000014_clean_hg38.csv"

mkdir -p "/mnt/cephfs/orgs/home/angel.pacheco/logs"
mkdir -p "/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes"

# Temporales en Scratch
BASE_TMP="/mnt/cephfs/scratch/$USER/plink_tmp"
TMPDIR="$BASE_TMP/$SLURM_JOB_ID"
mkdir -p "$TMPDIR"
chmod 700 "$TMPDIR"
export TMPDIR TMP="$TMPDIR" TEMP="$TMPDIR"
export APPTAINER_TMPDIR="$TMPDIR/apptainer-tmp"

echo "==== ETAPA 1: EXTRACCIÓN DE VARIANTES ZARR + LIMPIEZA PGS ===="
echo "Zarr Input:  $ZARR_IN"
echo "CSV Output:  $CSV_OUT"
echo "PGS Input:   $PGS_RAW"
echo "PGS Output:  $PGS_CLEAN"
echo "TMP:         $TMPDIR"

# ============================================================
# PARTE A: Extraer catálogo de variantes del Zarr
# ============================================================
PYTHON_SCRIPT="$TMPDIR/extract.py"
cat << 'EOF' > "$PYTHON_SCRIPT"
import sys
import zarr
import pandas as pd
import numpy as np
import time

zarr_path = sys.argv[1]
out_path  = sys.argv[2]

print(f"[{time.strftime('%H:%M:%S')}] Abriendo Zarr: {zarr_path}")
try:
    root = zarr.open_group(zarr_path, mode='r')

    print(f"[{time.strftime('%H:%M:%S')}] Extrayendo vectores...")
    chrom_indices = root['variant_contig'][:]
    pos = root['variant_position'][:]

    alleles = root['variant_allele'][:]
    ref = alleles[:, 0]
    alt = alleles[:, 1]

    print(f"[{time.strftime('%H:%M:%S')}] Mapeando índices a nombres reales de cromosomas...")
    try:
        contig_names = root['contig_id'][:]
        chrom_real = np.array(contig_names)[chrom_indices]
    except KeyError:
        print("Advertencia: No se encontró 'contig_id', usando índices crudos.")
        chrom_real = chrom_indices.astype(str)

    print(f"[{time.strftime('%H:%M:%S')}] Armando catálogo de variantes...")
    df = pd.DataFrame({
        'CHR': chrom_real,
        'POS': pos,
        'REF': ref.astype(str),
        'ALT': alt.astype(str)
    })

    print(f"[{time.strftime('%H:%M:%S')}] Guardando {len(df):,} variantes en {out_path}...")
    df.to_csv(out_path, index=False)
    print(f"Catálogo Zarr extraído exitosamente.")

except Exception as e:
    import traceback
    print("Error leyendo estructura Zarr:")
    traceback.print_exc()
    sys.exit(1)
EOF

echo ">> Extrayendo catálogo de variantes del Zarr..."
apptainer exec \
    --bind /mnt/cephfs:/mnt/cephfs \
    --bind "$TMPDIR:/tmp,$TMPDIR:/var/tmp" \
    "$SIF" \
    "$CONDA_PYTHON" "$PYTHON_SCRIPT" "$ZARR_IN" "$CSV_OUT"

# ============================================================
# PARTE B: Limpiar archivo PGS Catalog
# Usamos las columnas armonizadas (hm_chr, hm_pos) que están
# en coordenadas GRCh38, y descartamos filas sin coordenadas.
# ============================================================
echo ""
echo ">> Limpiando archivo PGS Catalog..."

CLEAN_SCRIPT="$TMPDIR/clean_pgs.py"
cat << 'EOF' > "$CLEAN_SCRIPT"
import pandas as pd
import sys
import time

pgs_raw  = sys.argv[1]
pgs_out  = sys.argv[2]

print(f"[{time.strftime('%H:%M:%S')}] Leyendo PGS Catalog: {pgs_raw}")

# Leer saltando las líneas de comentario (#)
df = pd.read_csv(pgs_raw, sep='\t', comment='#')

print(f"[{time.strftime('%H:%M:%S')}] Variantes totales en PGS: {len(df):,}")

# Filtrar filas sin coordenadas armonizadas hg38
before = len(df)
df = df.dropna(subset=['hm_chr', 'hm_pos'])
df = df[df['hm_pos'] != '']
after = len(df)
print(f"[{time.strftime('%H:%M:%S')}] Descartadas {before - after:,} variantes sin coordenadas hg38")

# Renombrar columnas para que coincidan con el catálogo del Zarr
# hm_chr -> CHR, hm_pos -> POS
# effect_allele -> EA, other_allele -> OA (NO asumimos cuál es REF/ALT)
df_clean = pd.DataFrame({
    'CHR': df['hm_chr'].astype(str),
    'POS': df['hm_pos'].astype(int),
    'EA':  df['effect_allele'].astype(str),
    'OA':  df['other_allele'].astype(str),
    'BETA': df['effect_weight'].astype(float)
})

# Verificación: solo SNVs de una letra
mask_snv = (df_clean['EA'].str.len() == 1) & (df_clean['OA'].str.len() == 1)
n_non_snv = (~mask_snv).sum()
if n_non_snv > 0:
    print(f"[{time.strftime('%H:%M:%S')}] Descartando {n_non_snv:,} variantes no-SNV (indels)")
    df_clean = df_clean[mask_snv]

# Eliminar duplicados por posición + alelos
before_dedup = len(df_clean)
df_clean = df_clean.drop_duplicates(subset=['CHR', 'POS', 'EA', 'OA'])
after_dedup = len(df_clean)
if before_dedup != after_dedup:
    print(f"[{time.strftime('%H:%M:%S')}] Eliminados {before_dedup - after_dedup:,} duplicados")

print(f"[{time.strftime('%H:%M:%S')}] Guardando {len(df_clean):,} variantes limpias en {pgs_out}")
df_clean.to_csv(pgs_out, index=False)
print("Limpieza PGS completada.")
EOF

apptainer exec \
    --bind /mnt/cephfs:/mnt/cephfs \
    --bind "$TMPDIR:/tmp,$TMPDIR:/var/tmp" \
    "$SIF" \
    "$CONDA_PYTHON" "$CLEAN_SCRIPT" "$PGS_RAW" "$PGS_CLEAN"

rm -rf "$TMPDIR"
echo ""
echo "==== ETAPA 1 COMPLETADA ===="
echo "Archivos generados:"
echo "  1. $CSV_OUT        (catálogo de variantes del Zarr)"
echo "  2. $PGS_CLEAN      (pesos PGS limpios en hg38)"
