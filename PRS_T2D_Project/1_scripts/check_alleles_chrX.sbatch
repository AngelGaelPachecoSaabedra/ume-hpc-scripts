#!/bin/bash
#SBATCH --job-name=check_chrX
#SBATCH --account=researchers
#SBATCH --qos=vip
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=00:30:00
#SBATCH --output=/mnt/cephfs/orgs/home/angel.pacheco/logs/check_chrX_%j.out
#SBATCH --error=/mnt/cephfs/orgs/home/angel.pacheco/logs/check_chrX_%j.err

set -euo pipefail

SIF="/mnt/cephfs/biocontainers/images/jupyter-biotools-1.2.sif"

TMPDIR="/mnt/cephfs/scratch/$USER/check_chrX/$SLURM_JOB_ID"
mkdir -p "$TMPDIR"
chmod 700 "$TMPDIR"
export TMPDIR TMP="$TMPDIR" TEMP="$TMPDIR"

echo "==== DIAGNÓSTICO: GWAS vs ZARR chrX (VCF directo) ===="

PYTHON_SCRIPT="$TMPDIR/check_chrX.py"
cat << 'PYEOF' > "$PYTHON_SCRIPT"
import zarr
import numpy as np
import gzip
import time
from collections import Counter

zarr_x_path = "/mnt/cephfs/hot_nvme/mcps/imputed-topmed/zar_files/chrX.zarr"

# ================================================================
# 1. Abrir chrX Zarr y explorar su contenido
# ================================================================
print(f"[{time.strftime('%H:%M:%S')}] Abriendo Zarr chrX...")
root_x = zarr.open_group(zarr_x_path, mode='r')

print(f"\n  Arrays disponibles en chrX.zarr:")
for k in root_x:
    arr = root_x[k]
    print(f"    {k}: shape={arr.shape}, dtype={arr.dtype}")

contig_ids = [c if isinstance(c, str) else c.decode() for c in root_x['contig_id'][:]]
z_pos = root_x['variant_position'][:]
z_alleles = root_x['variant_allele'][:]
z_chroms = root_x['variant_contig'][:]

print(f"\n  Contigs: {contig_ids}")
print(f"  Total variantes chrX: {len(z_pos):,}")
print(f"  Total muestras: {root_x['sample_id'].shape[0]:,}")

# Distribución de longitudes de alelos en el Zarr chrX
print(f"\n  Distribución de longitudes de alelos en Zarr chrX:")
ref_lens = []
alt_lens = []
n_check = min(len(z_pos), 500000)  # Revisar primer medio millón
for i in range(n_check):
    ref = z_alleles[i, 0]
    alt = z_alleles[i, 1]
    if isinstance(ref, bytes):
        ref = ref.decode()
        alt = alt.decode()
    ref_lens.append(len(ref))
    alt_lens.append(len(alt))

ref_counter = Counter(ref_lens)
alt_counter = Counter(alt_lens)
pair_counter = Counter(zip(ref_lens, alt_lens))

print(f"  (Basado en primeras {n_check:,} variantes)")
print(f"\n  REF length distribution:")
for l, c in sorted(ref_counter.items()):
    print(f"    len={l}: {c:,} ({c/n_check*100:.1f}%)")

print(f"\n  ALT length distribution:")
for l, c in sorted(alt_counter.items()):
    print(f"    len={l}: {c:,} ({c/n_check*100:.1f}%)")

n_snvs_zarr = pair_counter.get((1,1), 0)
n_indels_zarr = sum(c for (r,a), c in pair_counter.items() if r != a)
print(f"\n  SNVs en Zarr (1,1): {n_snvs_zarr:,}")
print(f"  Indels en Zarr (r!=a): {n_indels_zarr:,}")

# ================================================================
# 2. Extraer variantes significativas del GWAS chrX
# ================================================================
print(f"\n[{time.strftime('%H:%M:%S')}] Extrayendo variantes significativas de chrX del GWAS...")

gwas_variants = []
for sex, gwas_file in [("females", "/mnt/cephfs/orgs/home/angel.pacheco/gwas_miami/combined-gwas-results-all-ctrl-females.txt.gz"),
                        ("males", "/mnt/cephfs/orgs/home/angel.pacheco/gwas_miami/combined-gwas-results-all-ctrl-males.txt.gz")]:
    n_sig = 0
    n_indels = 0
    with gzip.open(gwas_file, 'rt') as f:
        f.readline()
        for line in f:
            cols = line.strip().split('\t')
            if cols[0] != '23':  # Solo chrX
                continue
            p_val = float(cols[14])
            if p_val >= 5e-8:
                continue
            
            pos = int(cols[1])
            a0, a1 = cols[3], cols[4]
            beta = float(cols[9])
            vtype = "SNV" if len(a0) == 1 and len(a1) == 1 else "INDEL"
            
            gwas_variants.append({
                'sex': sex, 'pos': pos, 'id': cols[2],
                'a0': a0, 'a1': a1, 'beta': beta, 'type': vtype
            })
            n_sig += 1
            if vtype == 'INDEL':
                n_indels += 1
    
    print(f"  {sex}: {n_sig} variantes significativas en chrX ({n_indels} indels)")

print(f"  Total chrX significativas (ambos sexos): {len(gwas_variants)}")

# ================================================================
# 3. Buscar cada variante del GWAS en el Zarr chrX
# ================================================================
print(f"\n[{time.strftime('%H:%M:%S')}] Buscando variantes en Zarr chrX...")

match_exact = 0
match_flipped = 0
match_pos_diff_alleles = 0
no_match_pos = 0

indel_matches = []
indel_mismatches = []
snv_matches = []
snv_mismatches = []

for v in gwas_variants:
    pos_matches = np.where(z_pos == v['pos'])[0]
    
    if len(pos_matches) == 0:
        no_match_pos += 1
        if v['type'] == 'INDEL':
            indel_mismatches.append({**v, 'reason': 'POS_NOT_FOUND', 'zarr_ref': '-', 'zarr_alt': '-'})
        else:
            snv_mismatches.append({**v, 'reason': 'POS_NOT_FOUND', 'zarr_ref': '-', 'zarr_alt': '-'})
        continue
    
    found = False
    for m in pos_matches:
        zarr_ref = z_alleles[m, 0]
        zarr_alt = z_alleles[m, 1]
        if isinstance(zarr_ref, bytes):
            zarr_ref = zarr_ref.decode()
            zarr_alt = zarr_alt.decode()
        
        if v['a0'] == zarr_ref and v['a1'] == zarr_alt:
            match_exact += 1
            found = True
            if v['type'] == 'INDEL':
                indel_matches.append({**v, 'zarr_ref': zarr_ref, 'zarr_alt': zarr_alt, 'match': 'EXACT'})
            break
        
        if v['a0'] == zarr_alt and v['a1'] == zarr_ref:
            match_flipped += 1
            found = True
            if v['type'] == 'INDEL':
                indel_matches.append({**v, 'zarr_ref': zarr_ref, 'zarr_alt': zarr_alt, 'match': 'FLIPPED'})
            break
    
    if not found:
        match_pos_diff_alleles += 1
        zarr_ref = z_alleles[pos_matches[0], 0]
        zarr_alt = z_alleles[pos_matches[0], 1]
        if isinstance(zarr_ref, bytes):
            zarr_ref = zarr_ref.decode()
            zarr_alt = zarr_alt.decode()
        
        if v['type'] == 'INDEL':
            indel_mismatches.append({**v, 'reason': 'ALLELE_MISMATCH', 'zarr_ref': zarr_ref, 'zarr_alt': zarr_alt})
        else:
            snv_mismatches.append({**v, 'reason': 'ALLELE_MISMATCH', 'zarr_ref': zarr_ref, 'zarr_alt': zarr_alt})

# ================================================================
# RESULTADOS
# ================================================================
total = match_exact + match_flipped + match_pos_diff_alleles + no_match_pos

print(f"\n{'='*70}")
print(f"RESULTADOS: GWAS chrX vs ZARR chrX (VCF directo)")
print(f"{'='*70}")
print(f"Total variantes: {total}")
print(f"")
print(f"  Match exacto (A0=REF, A1=ALT):      {match_exact:>5} ({match_exact/total*100:.1f}%)")
print(f"  Match invertido (A0=ALT, A1=REF):    {match_flipped:>5} ({match_flipped/total*100:.1f}%)")
print(f"  Posición OK, alelos difieren:        {match_pos_diff_alleles:>5} ({match_pos_diff_alleles/total*100:.1f}%)")
print(f"  Posición no encontrada:              {no_match_pos:>5} ({no_match_pos/total*100:.1f}%)")

# Desglose SNV vs INDEL
snvs = [v for v in gwas_variants if v['type'] == 'SNV']
indels = [v for v in gwas_variants if v['type'] == 'INDEL']

print(f"\n  --- SNVs ({len(snvs)}) ---")
snv_found = len(snvs) - len(snv_mismatches)
print(f"    Encontrados: {snv_found} ({snv_found/max(len(snvs),1)*100:.1f}%)")
print(f"    No encontrados: {len(snv_mismatches)}")

print(f"\n  --- INDELs ({len(indels)}) ---")
print(f"    Encontrados (match): {len(indel_matches)} ({len(indel_matches)/max(len(indels),1)*100:.1f}%)")
print(f"    No encontrados: {len(indel_mismatches)}")

if indel_matches:
    print(f"\n{'='*70}")
    print(f"INDELS QUE SÍ COINCIDEN ({len(indel_matches)}):")
    print(f"{'='*70}")
    print(f"{'SEX':<8} {'POS':<12} {'GWAS_A0':<15} {'GWAS_A1':<15} {'ZARR_REF':<15} {'ZARR_ALT':<15} {'MATCH'}")
    for v in indel_matches:
        print(f"{v['sex']:<8} {v['pos']:<12} {v['a0']:<15} {v['a1']:<15} {v['zarr_ref']:<15} {v['zarr_alt']:<15} {v['match']}")

if indel_mismatches:
    print(f"\n{'='*70}")
    print(f"INDELS QUE NO COINCIDEN ({len(indel_mismatches)}):")
    print(f"{'='*70}")
    print(f"{'SEX':<8} {'POS':<12} {'GWAS_A0':<15} {'GWAS_A1':<15} {'ZARR_REF':<15} {'ZARR_ALT':<15} {'RAZON'}")
    for v in indel_mismatches[:20]:
        print(f"{v['sex']:<8} {v['pos']:<12} {v['a0']:<15} {v['a1']:<15} {v['zarr_ref']:<15} {v['zarr_alt']:<15} {v['reason']}")

if snv_mismatches:
    print(f"\n{'='*70}")
    print(f"SNVs QUE NO COINCIDEN ({len(snv_mismatches)}):")
    print(f"{'='*70}")
    print(f"{'SEX':<8} {'POS':<12} {'GWAS_A0':<15} {'GWAS_A1':<15} {'ZARR_REF':<15} {'ZARR_ALT':<15} {'RAZON'}")
    for v in snv_mismatches[:20]:
        print(f"{v['sex']:<8} {v['pos']:<12} {v['a0']:<15} {v['a1']:<15} {v['zarr_ref']:<15} {v['zarr_alt']:<15} {v['reason']}")

print(f"\n[{time.strftime('%H:%M:%S')}] Diagnóstico chrX completado.")
PYEOF

apptainer exec \
    --bind /mnt/cephfs:/mnt/cephfs \
    --bind "$TMPDIR:/tmp,$TMPDIR:/var/tmp" \
    "$SIF" \
    /opt/miniforge3/envs/vcf2zarr/bin/python "$PYTHON_SCRIPT"

rm -rf "$TMPDIR"
echo "==== COMPLETADO ===="
