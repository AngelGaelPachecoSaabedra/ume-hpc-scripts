#!/bin/bash
#SBATCH --job-name=sbrc_download
#SBATCH --account=researchers
#SBATCH --qos=vip
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/mnt/cephfs/orgs/home/angel.pacheco/logs/sbrc_download_%j.out
#SBATCH --error=/mnt/cephfs/orgs/home/angel.pacheco/logs/sbrc_download_%j.err

set -euo pipefail

SIF="/mnt/cephfs/biocontainers/images/jupyter-biotools-1.4.sif"
OUT_BASE="/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes/sbayesrc"
mkdir -p "$OUT_BASE"

TMPDIR="/mnt/cephfs/scratch/$USER/sbrc_tmp/$SLURM_JOB_ID"
mkdir -p "$TMPDIR"
chmod 700 "$TMPDIR"
export TMPDIR TMP="$TMPDIR" TEMP="$TMPDIR"

echo "============================================================"
echo "  PIPELINE A — ETAPA 1: Formateo T2DGGI"
echo "  $(date)"
echo "============================================================"

# ================================================================
# 1C. Convertir a formato COJO (.ma) para GCTB/SBayesRC
# ================================================================
echo ""
echo "--- 1C: Convirtiendo a formato COJO (.ma) ---"

PYTHON_SCRIPT="$TMPDIR/format_t2dggi.py"
cat << 'PYEOF' > "$PYTHON_SCRIPT"
import pandas as pd
import numpy as np
import sys
import os
import gzip
import time

out_base = "/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes/sbayesrc"

t2dggi_path = os.path.join(out_base, "T2DGGI_all_ancestry_summary_stats.txt.gz")
print(f"[{time.strftime('%H:%M:%S')}] Leyendo {t2dggi_path}...")

# Leer header para identificar columnas
with gzip.open(t2dggi_path, 'rt') as f:
    header = f.readline().strip().split('\t')
print(f"  Columnas encontradas: {header}")

# Cargar datos (quitamos el dtype de CHR porque el nombre de la columna cambia)
df = pd.read_csv(t2dggi_path, sep='\t')
print(f"  Variantes leídas: {len(df):,}")

# Mapear columnas (Añadidas las variaciones de ortografía de T2DGGI)
col_map = {}
for c in df.columns:
    cl = c.lower()
    if cl in ('snp', 'rsid', 'markername', 'variant_id', 'id'): col_map['SNP'] = c
    elif cl in ('a1', 'ea', 'effect_allele', 'allele1', 'effectallele'): col_map['A1'] = c
    elif cl in ('a2', 'nea', 'other_allele', 'allele2', 'non_effect_allele', 'noneffectallele'): col_map['A2'] = c
    elif cl in ('freq', 'eaf', 'a1freq', 'frq', 'maf'): col_map['freq'] = c
    elif cl in ('beta', 'b', 'effect'): col_map['b'] = c
    elif cl in ('se', 'stderr', 'stderrbeta', 'betase'): col_map['se'] = c
    elif cl in ('p', 'pval', 'p_value', 'pvalue'): col_map['p'] = c
    elif cl in ('n', 'n_samples', 'neff', 'n_eff', 'nmiss'): col_map['N'] = c
    elif cl in ('chr', 'chrom', 'chromosome', 'chromsome'): col_map['CHR'] = c
    elif cl in ('pos', 'bp', 'position', 'genpos'): col_map['POS'] = c

print(f"  Mapeo de columnas: {col_map}")

# Si no hay SNP ID, crear uno con CHR:POS:A1:A2
if 'SNP' not in col_map:
    print("  Creando SNP ID desde CHR:POS:A1:A2...")
    df['SNP_ID'] = df[col_map['CHR']].astype(str) + ':' + df[col_map['POS']].astype(str) + ':' + \
                   df[col_map['A1']].astype(str) + ':' + df[col_map['A2']].astype(str)
    col_map['SNP'] = 'SNP_ID'

# Crear DataFrame COJO
df_cojo = pd.DataFrame({
    'SNP': df[col_map['SNP']],
    'A1':  df[col_map['A1']],
    'A2':  df[col_map['A2']],
    'freq': df[col_map['freq']].astype(float),
    'b':   df[col_map['b']].astype(float),
    'se':  df[col_map['se']].astype(float),
    'p':   df[col_map['p']].astype(float),
    'N':   df[col_map['N']].astype(float) if 'N' in col_map else 2535586
})

# QC básico
n_before = len(df_cojo)
df_cojo = df_cojo.dropna()
df_cojo = df_cojo[df_cojo['se'] > 0]
df_cojo = df_cojo[df_cojo['p'] > 0]
df_cojo = df_cojo[(df_cojo['freq'] > 0.001) & (df_cojo['freq'] < 0.999)]
n_after = len(df_cojo)
print(f"  QC: {n_before:,} → {n_after:,} variantes ({n_before - n_after:,} removidas)")

# Solo autosomas
if 'CHR' in col_map:
    df_cojo['_CHR'] = df[col_map['CHR']].astype(str)
    df_cojo = df_cojo[~df_cojo['_CHR'].isin(['X', 'Y', 'MT', 'chrX', 'chrY', '23', '24', '25'])]
    df_cojo = df_cojo.drop(columns=['_CHR'])
    print(f"  Solo autosomas: {len(df_cojo):,} variantes")

# Guardar
ma_file = os.path.join(out_base, "T2DGGI_ALL_autosomas.ma")
df_cojo.to_csv(ma_file, sep='\t', index=False)
print(f"\n[{time.strftime('%H:%M:%S')}] COJO guardado: {ma_file}")
PYEOF

apptainer exec --bind /mnt/cephfs:/mnt/cephfs --bind "$TMPDIR:/tmp,$TMPDIR:/var/tmp" \
    "$SIF" /opt/miniforge3/envs/vcf2zarr/bin/python "$PYTHON_SCRIPT"

rm -rf "$TMPDIR"
echo ""
echo "============================================================"
echo "  ETAPA 1 COMPLETADA."
echo "============================================================"
