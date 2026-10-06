#!/bin/bash
#SBATCH --job-name=check_alleles
#SBATCH --account=researchers
#SBATCH --qos=vip
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=00:30:00
#SBATCH --output=/mnt/cephfs/orgs/home/angel.pacheco/logs/check_alleles_%j.out
#SBATCH --error=/mnt/cephfs/orgs/home/angel.pacheco/logs/check_alleles_%j.err

set -euo pipefail

SIF="/mnt/cephfs/biocontainers/images/jupyter-biotools-1.2.sif"

TMPDIR="/mnt/cephfs/scratch/$USER/check_alleles/$SLURM_JOB_ID"
mkdir -p "$TMPDIR"
chmod 700 "$TMPDIR"
export TMPDIR TMP="$TMPDIR" TEMP="$TMPDIR"

echo "==== DIAGNÓSTICO: REPRESENTACIÓN DE ALELOS GWAS vs ZARR ===="

PYTHON_SCRIPT="$TMPDIR/check_alleles.py"
cat << 'PYEOF' > "$PYTHON_SCRIPT"
import zarr
import numpy as np
import gzip
import time

zarr_main_path = "/mnt/cephfs/hot_nvme/mcps/imputed-topmed/plink_to_zar/MCPS.HIS.TOPMED.b38.FR.2602.TARGET.imputed.plink.update2.zarr"
gwas_file = "/mnt/cephfs/orgs/home/angel.pacheco/gwas_miami/combined-gwas-results-all-ctrl-females.txt.gz"

print(f"[{time.strftime('%H:%M:%S')}] Abriendo Zarr...")
root = zarr.open_group(zarr_main_path, mode='r')

contig_ids = [c if isinstance(c, str) else c.decode() for c in root['contig_id'][:]]
z_chroms = root['variant_contig'][:]
z_pos = root['variant_position'][:]
z_alleles = root['variant_allele'][:]

print(f"  Zarr: {len(z_pos):,} variantes, contigs: {contig_ids}")

# Construir diccionario de lookup del Zarr: (chrom_str, pos) -> [(ref, alt), ...]
# Solo para un subconjunto (las variantes significativas del GWAS)
# Primero extraemos las variantes significativas (indels incluidos) del GWAS

print(f"\n[{time.strftime('%H:%M:%S')}] Extrayendo variantes significativas del GWAS (p < 5e-8, SIN filtro de longitud)...")

gwas_variants = []
gwas_snvs = 0
gwas_indels = 0
gwas_multi = 0

with gzip.open(gwas_file, 'rt') as f:
    header = f.readline().strip().split('\t')
    print(f"  Header GWAS: {header[:8]}...")
    
    for line in f:
        cols = line.strip().split('\t')
        p_val = float(cols[14])  # columna P
        if p_val >= 5e-8:
            continue
        
        chrom = cols[0]       # CHROM
        pos = int(cols[1])    # GENPOS
        snp_id = cols[2]      # ID
        allele0 = cols[3]     # ALLELE0 (REF/other)
        allele1 = cols[4]     # ALLELE1 (ALT/effect)
        beta = float(cols[9])
        
        len0, len1 = len(allele0), len(allele1)
        
        if len0 == 1 and len1 == 1:
            vtype = "SNV"
            gwas_snvs += 1
        elif len0 != len1:
            vtype = "INDEL"
            gwas_indels += 1
        else:
            vtype = "MNV"  # multi-nucleotide variant
            gwas_multi += 1
        
        gwas_variants.append({
            'chrom': chrom, 'pos': pos, 'id': snp_id,
            'a0': allele0, 'a1': allele1, 'beta': beta,
            'type': vtype, 'len0': len0, 'len1': len1
        })

print(f"\n  Total variantes significativas (p < 5e-8): {len(gwas_variants):,}")
print(f"  SNVs (1 vs 1): {gwas_snvs:,}")
print(f"  Indels: {gwas_indels:,}")
print(f"  MNVs: {gwas_multi:,}")

# ================================================================
# Mostrar ejemplos de indels del GWAS
# ================================================================
indels = [v for v in gwas_variants if v['type'] == 'INDEL']
print(f"\n{'='*70}")
print(f"EJEMPLOS DE INDELS EN EL GWAS (primeros 20):")
print(f"{'='*70}")
print(f"{'CHROM':<6} {'POS':<12} {'ALLELE0':<15} {'ALLELE1':<15} {'TIPO':<10} {'BETA':<10} {'ID'}")
for v in indels[:20]:
    print(f"{v['chrom']:<6} {v['pos']:<12} {v['a0']:<15} {v['a1']:<15} {v['type']:<10} {v['beta']:<10.4f} {v['id']}")

# ================================================================
# Ahora buscar CADA variante significativa en el Zarr
# ================================================================
print(f"\n[{time.strftime('%H:%M:%S')}] Buscando variantes en el Zarr...")

# Estadísticas de matching
match_exact = 0        # GWAS A0/A1 coincide exactamente con Zarr REF/ALT
match_flipped = 0      # GWAS A0/A1 es el inverso del Zarr (A0=ALT, A1=REF)
match_pos_diff_alleles = 0  # Misma posición pero alelos diferentes
no_match_pos = 0       # Posición no encontrada en el Zarr

# Detalles de los que NO coinciden
mismatches = []
matches_detail = []

for v in gwas_variants:
    chrom_str = v['chrom']
    # Manejar chrX: GWAS usa "23"
    if chrom_str == '23':
        continue  # Saltamos chrX por ahora, está en otro Zarr
    
    try:
        chrom_idx = contig_ids.index(chrom_str)
    except ValueError:
        no_match_pos += 1
        continue
    
    idx_start = np.searchsorted(z_chroms, chrom_idx, side='left')
    idx_end = np.searchsorted(z_chroms, chrom_idx, side='right')
    
    positions_chr = z_pos[idx_start:idx_end]
    pos_matches = np.where(positions_chr == v['pos'])[0]
    
    if len(pos_matches) == 0:
        no_match_pos += 1
        if v['type'] == 'INDEL':
            mismatches.append({**v, 'reason': 'POS_NOT_FOUND', 'zarr_ref': '-', 'zarr_alt': '-'})
        continue
    
    found = False
    for m in pos_matches:
        global_idx = idx_start + m
        zarr_ref = z_alleles[global_idx, 0]
        zarr_alt = z_alleles[global_idx, 1]
        if isinstance(zarr_ref, bytes):
            zarr_ref = zarr_ref.decode()
            zarr_alt = zarr_alt.decode()
        
        # Caso exacto: GWAS A0 == Zarr REF y GWAS A1 == Zarr ALT
        if v['a0'] == zarr_ref and v['a1'] == zarr_alt:
            match_exact += 1
            found = True
            if v['type'] == 'INDEL':
                matches_detail.append({**v, 'zarr_ref': zarr_ref, 'zarr_alt': zarr_alt, 'match': 'EXACT'})
            break
        
        # Caso invertido: GWAS A0 == Zarr ALT y GWAS A1 == Zarr REF
        if v['a0'] == zarr_alt and v['a1'] == zarr_ref:
            match_flipped += 1
            found = True
            if v['type'] == 'INDEL':
                matches_detail.append({**v, 'zarr_ref': zarr_ref, 'zarr_alt': zarr_alt, 'match': 'FLIPPED'})
            break
    
    if not found:
        match_pos_diff_alleles += 1
        # Guardar detalle del mismatch
        global_idx = idx_start + pos_matches[0]
        zarr_ref = z_alleles[global_idx, 0]
        zarr_alt = z_alleles[global_idx, 1]
        if isinstance(zarr_ref, bytes):
            zarr_ref = zarr_ref.decode()
            zarr_alt = zarr_alt.decode()
        
        if v['type'] == 'INDEL':
            mismatches.append({**v, 'reason': 'ALLELE_MISMATCH', 'zarr_ref': zarr_ref, 'zarr_alt': zarr_alt})

# ================================================================
# RESULTADOS
# ================================================================
total_checked = match_exact + match_flipped + match_pos_diff_alleles + no_match_pos
print(f"\n{'='*70}")
print(f"RESULTADOS DE ALINEACIÓN GWAS vs ZARR")
print(f"{'='*70}")
print(f"Total variantes checadas (autosomas): {total_checked:,}")
print(f"")
print(f"  Match exacto (A0=REF, A1=ALT):     {match_exact:>6,} ({match_exact/total_checked*100:.1f}%)")
print(f"  Match invertido (A0=ALT, A1=REF):   {match_flipped:>6,} ({match_flipped/total_checked*100:.1f}%)")
print(f"  Posición encontrada, alelos difieren:{match_pos_diff_alleles:>6,} ({match_pos_diff_alleles/total_checked*100:.1f}%)")
print(f"  Posición no encontrada en Zarr:      {no_match_pos:>6,} ({no_match_pos/total_checked*100:.1f}%)")

# Desglose por tipo
for vtype in ['SNV', 'INDEL', 'MNV']:
    sub = [v for v in gwas_variants if v['type'] == vtype and v['chrom'] != '23']
    if not sub:
        continue
    print(f"\n  --- {vtype}s ({len(sub):,} variantes) ---")

# ================================================================
# Detalle de indels que SÍ coinciden
# ================================================================
if matches_detail:
    print(f"\n{'='*70}")
    print(f"INDELS QUE SÍ COINCIDEN CON EL ZARR ({len(matches_detail)}):")
    print(f"{'='*70}")
    print(f"{'CHROM':<6} {'POS':<12} {'GWAS_A0':<15} {'GWAS_A1':<15} {'ZARR_REF':<15} {'ZARR_ALT':<15} {'MATCH'}")
    for v in matches_detail[:30]:
        print(f"{v['chrom']:<6} {v['pos']:<12} {v['a0']:<15} {v['a1']:<15} {v['zarr_ref']:<15} {v['zarr_alt']:<15} {v['match']}")

# ================================================================
# Detalle de indels que NO coinciden
# ================================================================
if mismatches:
    print(f"\n{'='*70}")
    print(f"INDELS QUE NO COINCIDEN ({len(mismatches)}):")
    print(f"{'='*70}")
    print(f"{'CHROM':<6} {'POS':<12} {'GWAS_A0':<15} {'GWAS_A1':<15} {'ZARR_REF':<15} {'ZARR_ALT':<15} {'RAZON'}")
    for v in mismatches[:30]:
        print(f"{v['chrom']:<6} {v['pos']:<12} {v['a0']:<15} {v['a1']:<15} {v['zarr_ref']:<15} {v['zarr_alt']:<15} {v['reason']}")

# ================================================================
# Distribución de longitudes de alelos en indels significativos
# ================================================================
print(f"\n{'='*70}")
print(f"DISTRIBUCIÓN DE LONGITUDES DE ALELOS (indels significativos)")
print(f"{'='*70}")
from collections import Counter
len_pairs_gwas = Counter((v['len0'], v['len1']) for v in gwas_variants if v['type'] == 'INDEL')
print("GWAS (ALLELE0_len, ALLELE1_len) -> count:")
for (l0, l1), cnt in sorted(len_pairs_gwas.items(), key=lambda x: -x[1])[:15]:
    print(f"  ({l0}, {l1}): {cnt}")

print(f"\n[{time.strftime('%H:%M:%S')}] Diagnóstico completado.")
PYEOF

apptainer exec \
    --bind /mnt/cephfs:/mnt/cephfs \
    --bind "$TMPDIR:/tmp,$TMPDIR:/var/tmp" \
    "$SIF" \
    /opt/miniforge3/envs/vcf2zarr/bin/python "$PYTHON_SCRIPT"

rm -rf "$TMPDIR"
echo "==== COMPLETADO ===="
