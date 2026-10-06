#!/bin/bash
#SBATCH --job-name=prs_mcps
#SBATCH --account=researchers
#SBATCH --qos=vip
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=01:00:00
#SBATCH --output=/mnt/cephfs/orgs/home/angel.pacheco/logs/prs_mcps_%j.out
#SBATCH --error=/mnt/cephfs/orgs/home/angel.pacheco/logs/prs_mcps_%j.err

set -euo pipefail

SIF="/mnt/cephfs/biocontainers/images/jupyter-biotools-1.2.sif"

TMPDIR="/mnt/cephfs/scratch/$USER/prs_mcps_tmp/$SLURM_JOB_ID"
mkdir -p "$TMPDIR"
chmod 700 "$TMPDIR"
export TMPDIR TMP="$TMPDIR" TEMP="$TMPDIR"

echo "============================================================"
echo "PRS DIABETES T2 - BETAS GWAS MCPS (AUTOSOMAS + ChrX)"
echo "Cohorte Completa: 140K pacientes"
echo "============================================================"

PYTHON_SCRIPT="$TMPDIR/compute_prs_mcps.py"
cat << 'PYEOF' > "$PYTHON_SCRIPT"
import zarr
import numpy as np
import pandas as pd
import time
import gzip

print(f"[{time.strftime('%H:%M:%S')}] Iniciando PRS con betas MCPS...")

# ================================================================
# RUTAS
# ================================================================
zarr_main_path = "/mnt/cephfs/hot_nvme/mcps/imputed-topmed/plink_to_zar/MCPS.HIS.TOPMED.b38.FR.2602.TARGET.imputed.plink.update2.zarr"
zarr_x_path    = "/mnt/cephfs/hot_nvme/mcps/imputed-topmed/zar_files/chrX.zarr"
gwas_dir       = "/mnt/cephfs/orgs/home/angel.pacheco/gwas_miami"
ct_dir         = "/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes/ct_gwas"
out_dir        = "/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes"

# Archivos de datos clínicos y puente de IDs
covar_file     = "/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes/MCPS BASELINE.csv"
link_file      = "/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes/RGN_LINK_IID.csv"

# ================================================================
# 1. CARGAR VARIANTES C+T
# ================================================================
print(f"[{time.strftime('%H:%M:%S')}] Cargando variantes C+T...")
df_fem = pd.read_csv(f"{ct_dir}/clumped_females_all.tsv", sep='\t')
df_mal = pd.read_csv(f"{ct_dir}/clumped_males_all.tsv", sep='\t')

# ================================================================
# 2. QC: Filtrar por INFO y MAF del GWAS original
# ================================================================
print(f"[{time.strftime('%H:%M:%S')}] Aplicando QC (INFO >= 0.8, MAF >= 0.01)...")

def load_gwas_info(gwas_gz, variant_ids):
    ids_no_chr = {}
    for v in variant_ids:
        v_clean = v.replace("chrX:", "X:").replace("chr", "", 1) if v.startswith("chr") else v
        ids_no_chr[v_clean] = v
    
    info_dict = {}
    with gzip.open(gwas_gz, 'rt') as f:
        header = f.readline()
        for line in f:
            cols = line.strip().split('\t')
            snp_id = cols[2]
            if snp_id in ids_no_chr:
                info_score = float(cols[6])
                a1freq = float(cols[5])
                maf = min(a1freq, 1 - a1freq)
                info_dict[ids_no_chr[snp_id]] = {'INFO': info_score, 'MAF': maf}
    return info_dict

info_fem = load_gwas_info(f"{gwas_dir}/combined-gwas-results-all-ctrl-females.txt.gz", df_fem['SNP'].tolist())
df_fem['INFO'] = df_fem['SNP'].map(lambda x: info_fem.get(x, {}).get('INFO', 0))
df_fem['MAF']  = df_fem['SNP'].map(lambda x: info_fem.get(x, {}).get('MAF', 0))
df_fem = df_fem[(df_fem['INFO'] >= 0.8) & (df_fem['MAF'] >= 0.01)]

info_mal = load_gwas_info(f"{gwas_dir}/combined-gwas-results-all-ctrl-males.txt.gz", df_mal['SNP'].tolist())
df_mal['INFO'] = df_mal['SNP'].map(lambda x: info_mal.get(x, {}).get('INFO', 0))
df_mal['MAF']  = df_mal['SNP'].map(lambda x: info_mal.get(x, {}).get('MAF', 0))
df_mal = df_mal[(df_mal['INFO'] >= 0.8) & (df_mal['MAF'] >= 0.01)]

# ================================================================
# 3. CARGAR SEXO DE LOS PACIENTES Y MAPEO DE IDs
# ================================================================
print(f"[{time.strftime('%H:%M:%S')}] Cargando información de sexo y conectando IDs...")

# Leer archivo clínico (PATID -> MALE)
df_baseline = pd.read_csv(covar_file)
df_baseline = df_baseline.dropna(subset=['MALE']) # Quitar NaNs
patid_to_sex = dict(zip(df_baseline['PATID'].astype(str), df_baseline['MALE'].astype(int)))

# Leer archivo puente (IID -> PATID)
df_link = pd.read_csv(link_file)
iid_to_patid = dict(zip(df_link['IID'].astype(str), df_link['PATID'].astype(str)))

# Crear diccionario directo (IID -> Sexo)
sex_map = {}
for iid, patid in iid_to_patid.items():
    if patid in patid_to_sex:
        sex_map[iid] = patid_to_sex[patid]

print(f"  Pacientes clínicos cargados: {len(patid_to_sex):,}")
print(f"  Mapeos IID->PATID encontrados: {len(iid_to_patid):,}")
print(f"  IIDs con sexo confirmado listos para el PRS: {len(sex_map):,}")

# ================================================================
# 4. ABRIR AMBOS ZARRS Y PREPARAR ÍNDICES
# ================================================================
print(f"[{time.strftime('%H:%M:%S')}] Abriendo Zarrs y alineando muestras...")

root_main = zarr.open_group(zarr_main_path, mode='r')
root_x    = zarr.open_group(zarr_x_path, mode='r')

# Procesar muestras del Zarr principal como la referencia (140,831)
sample_ids_main = [s if isinstance(s, str) else s.decode() for s in root_main['sample_id'][:]]
n_samples = len(sample_ids_main)

# Determinar sexo para la matriz principal usando nuestro sex_map construido
sample_sex = np.array([sex_map.get(sid, -1) for sid in sample_ids_main])
idx_female  = np.where(sample_sex == 0)[0]
idx_male    = np.where(sample_sex == 1)[0]
idx_unknown = np.where(sample_sex == -1)[0]

print(f"  Mujeres identificadas en Zarr: {len(idx_female):,}")
print(f"  Hombres identificados en Zarr: {len(idx_male):,}")
print(f"  Muestras sin mapeo de sexo: {len(idx_unknown):,}")

# Mapear las muestras del ChrX para asegurar el mismo orden
sample_ids_x = [s if isinstance(s, str) else s.decode() for s in root_x['sample_id'][:]]
main_id_to_idx = {sid: i for i, sid in enumerate(sample_ids_main)}
x_to_main_order = np.array([main_id_to_idx.get(sid, -1) for sid in sample_ids_x])

# Metadata para localización
contigs_main = list(root_main['contig_id'][:])
z_chroms_main = root_main['variant_contig'][:]
z_pos_main    = root_main['variant_position'][:]
z_alleles_main= root_main['variant_allele'][:]

contigs_x = list(root_x['contig_id'][:])
z_chroms_x = root_x['variant_contig'][:]
z_pos_x    = root_x['variant_position'][:]
z_alleles_x= root_x['variant_allele'][:]

# ================================================================
# 5. FUNCIÓN PARA BUSCAR VARIANTE EN EL ZARR CORRECTO
# ================================================================
def extract_dosages(root_zarr, global_idx):
    """Extrae las dosis. Usa call_DS si existe, si no suma call_genotype."""
    if 'call_DS' in root_zarr:
        return root_zarr['call_DS'][global_idx, :]
    else:
        g_raw = root_zarr['call_genotype'][global_idx, :, :]
        return np.sum(np.maximum(g_raw, 0), axis=1).astype(np.float32)

def find_variant_genotypes(chrom_str, pos, a1_effect, a2_other):
    # Definir qué Zarr usar
    if chrom_str == 'X' or chrom_str == '23':
        root_target = root_x
        c_list = contigs_x
        z_chroms = z_chroms_x
        z_pos = z_pos_x
        z_alleles = z_alleles_x
        is_x = True
        try:
            chrom_idx = c_list.index('X')
        except ValueError:
            return None
    else:
        root_target = root_main
        c_list = contigs_main
        z_chroms = z_chroms_main
        z_pos = z_pos_main
        z_alleles = z_alleles_main
        is_x = False
        try:
            chrom_idx = c_list.index(chrom_str)
        except ValueError:
            return None

    idx_start = np.searchsorted(z_chroms, chrom_idx, side='left')
    idx_end   = np.searchsorted(z_chroms, chrom_idx, side='right')
    
    positions_chr = z_pos[idx_start:idx_end]
    matches = np.where(positions_chr == pos)[0]
    
    for m in matches:
        global_idx = idx_start + m
        ref = z_alleles[global_idx, 0]
        alt = z_alleles[global_idx, 1]
        if isinstance(ref, bytes):
            ref = ref.decode()
            alt = alt.decode()
        
        # Caso directo: effect_allele == ALT
        if a1_effect == alt and a2_other == ref:
            g_sum = extract_dosages(root_target, global_idx)
            # Reordenar si viene del ChrX
            if is_x:
                aligned_g = np.zeros(n_samples, dtype=np.float32)
                valid_mask = x_to_main_order >= 0
                aligned_g[x_to_main_order[valid_mask]] = g_sum[valid_mask]
                return aligned_g
            return g_sum
            
        # Caso invertido: effect_allele == REF
        if a1_effect == ref and a2_other == alt:
            g_sum = extract_dosages(root_target, global_idx)
            g_inverted = 2.0 - g_sum
            if is_x:
                aligned_g = np.zeros(n_samples, dtype=np.float32)
                valid_mask = x_to_main_order >= 0
                aligned_g[x_to_main_order[valid_mask]] = g_inverted[valid_mask]
                return aligned_g
            return g_inverted
            
    return None

# ================================================================
# 6. CALCULAR PRS
# ================================================================
prs_scores = np.zeros(n_samples, dtype=np.float32)

def calculate_for_group(df, idx_group, group_name):
    print(f"\n[{time.strftime('%H:%M:%S')}] Calculando PRS {group_name} ({len(df)} variantes)...")
    n_found = 0
    n_miss = 0
    for _, row in df.iterrows():
        chrom_raw = str(int(row['CHR']))
        chrom_str = 'X' if chrom_raw == '23' else chrom_raw
        pos = int(row['BP'])
        
        geno = find_variant_genotypes(chrom_str, pos, row['A1'], row['A2'])
        
        if geno is not None:
            prs_scores[idx_group] += geno[idx_group] * float(row['BETA'])
            n_found += 1
        else:
            n_miss += 1
    print(f"  Encontradas: {n_found}, No encontradas: {n_miss}")

calculate_for_group(df_fem, idx_female, "mujeres")
calculate_for_group(df_mal, idx_male, "hombres")

# ================================================================
# 7. CONSTRUIR TABLA FINAL Y Z-SCORE
# ================================================================
print(f"\n[{time.strftime('%H:%M:%S')}] Generando tabla final...")

# Al guardar, queremos conservar ambos IDs para que te sea fácil cruzarlo después
# Recuperar el PATID para el resultado final usando el mapa inverso (o None si no existe)
patid_list = [iid_to_patid.get(sid, 'UNKNOWN') for sid in sample_ids_main]

df_result = pd.DataFrame({
    'PATID': patid_list,
    'IID': sample_ids_main, 
    'SEX': sample_sex, 
    'PRS_T2D': prs_scores
})

# Filtramos los que tengan sexo conocido (los que se calcularon)
df_with_sex = df_result[df_result['SEX'] >= 0].copy()

for sex_val in [0, 1]:
    mask = df_with_sex['SEX'] == sex_val
    vals = df_with_sex.loc[mask, 'PRS_T2D']
    if len(vals) > 0 and vals.std() > 0:
        df_with_sex.loc[mask, 'PRS_Z'] = (vals - vals.mean()) / vals.std()
    else:
        df_with_sex.loc[mask, 'PRS_Z'] = 0.0
    df_with_sex.loc[mask, 'PERCENTILE'] = df_with_sex.loc[mask, 'PRS_T2D'].rank(pct=True) * 100

df_with_sex = df_with_sex.sort_values('PRS_T2D', ascending=False)
out_tsv = f"{out_dir}/prs_diabetes_mcps_gwas_final.tsv"
df_with_sex.to_csv(out_tsv, sep='\t', index=False, float_format='%.6f')

# ================================================================
# 8. RESUMEN
# ================================================================
print(f"\n{'='*60}")
print(f"RESUMEN PRS DIABETES T2 - GWAS MCPS")
print(f"{'='*60}")

for sex_val, sex_name in [(0, 'Mujeres'), (1, 'Hombres')]:
    sub = df_with_sex[df_with_sex['SEX'] == sex_val]
    if len(sub) > 0:
        print(f"\n{sex_name} (n={len(sub):,}):")
        print(f"  PRS: min={sub['PRS_T2D'].min():.6f}, max={sub['PRS_T2D'].max():.6f}, mean={sub['PRS_T2D'].mean():.6f}, std={sub['PRS_T2D'].std():.6f}")
        print(f"  Top 3 pacientes:")
        for _, r in sub.head(3).iterrows():
            print(f"    {r['PATID']} ({r['IID']}): PRS={r['PRS_T2D']:.6f} (z={r['PRS_Z']:.2f}, P{r['PERCENTILE']:.1f})")

print(f"\nArchivo guardado en: {out_tsv}")
print(f"\n[{time.strftime('%H:%M:%S')}] PIPELINE PRS MCPS COMPLETADO.")
PYEOF

apptainer exec \
    --bind /mnt/cephfs:/mnt/cephfs \
    --bind "$TMPDIR:/tmp,$TMPDIR:/var/tmp" \
    "$SIF" \
    /opt/miniforge3/envs/vcf2zarr/bin/python "$PYTHON_SCRIPT"

rm -rf "$TMPDIR"
echo "==== COMPLETADO ===="
