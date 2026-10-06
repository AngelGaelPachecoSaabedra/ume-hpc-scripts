#!/bin/bash
#SBATCH --job-name=prs_end2end
#SBATCH --account=researchers
#SBATCH --qos=vip
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=80G
#SBATCH --time=04:00:00
#SBATCH --output=/mnt/cephfs/orgs/home/angel.pacheco/logs/prs_end2end_%j.out
#SBATCH --error=/mnt/cephfs/orgs/home/angel.pacheco/logs/prs_end2end_%j.err

set -euo pipefail

SIF="/mnt/cephfs/biocontainers/images/jupyter-biotools-1.2.sif"

TMPDIR="/mnt/cephfs/scratch/$USER/prs_end2end_tmp/$SLURM_JOB_ID"
mkdir -p "$TMPDIR"
chmod 700 "$TMPDIR"
export TMPDIR TMP="$TMPDIR" TEMP="$TMPDIR"

# ================================================================
# RUTAS GLOBALES
# ================================================================
GWAS_DIR="/mnt/cephfs/orgs/home/angel.pacheco/gwas_miami"
PLINK_DIR="/mnt/cephfs/hot_nvme/mcps/imputed-topmed/plink_files/maximally_unrelated"
OUT_DIR="/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes"
CT_DIR="${OUT_DIR}/ct_gwas_v2"
mkdir -p "$CT_DIR"

echo "============================================================"
echo "PIPELINE END-TO-END: PRS DIABETES T2 (MCPS)"
echo "1. Clumping C+T (Panel Maximally Unrelated)"
echo "2. PRS Scoring (Cohorte Zarr 140K + Desglose por Cromosoma)"
echo "============================================================"

# ================================================================
# FASE 1: CLUMPING + THRESHOLDING (C+T)
# ================================================================
echo -e "\n[$(date +'%H:%M:%S')] INICIANDO FASE 1: CLUMPING (C+T)"

for SEX in females males; do
    echo "  >> Preparando stats de ${SEX}..."
    GWAS_FILE="${GWAS_DIR}/combined-gwas-results-all-ctrl-${SEX}.txt.gz"
    CLEAN_FILE="${CT_DIR}/gwas_clean_${SEX}.tsv"
    
    zcat "$GWAS_FILE" | awk -F'\t' '
        NR == 1 { print "SNP\tCHR\tBP\tA1\tA2\tBETA\tP"; next }
        $15 < 5e-8 && length($4) == 1 && length($5) == 1 {
            snp = "chr"$3
            print snp"\t"$1"\t"$2"\t"$5"\t"$4"\t"$10"\t"$15
        }
    ' > "$CLEAN_FILE"
done

for SEX in females males; do
    echo "  >> Ejecutando Clumping para ${SEX}..."
    CLEAN_FILE="${CT_DIR}/gwas_clean_${SEX}.tsv"
    CLUMP_MERGED="${CT_DIR}/clumped_${SEX}_all.tsv"
    
    echo -e "CHR\tBP\tSNP\tA1\tA2\tBETA\tP" > "$CLUMP_MERGED"
    
    for CHR in {1..22} X; do
        PLINK_PREFIX="${PLINK_DIR}/mcps-freeze150k_qcd_chr${CHR}_ivs"
        CHR_CODE=$([ "$CHR" == "X" ] && echo "23" || echo "$CHR")
        CLUMP_OUT="${TMPDIR}/clump_${SEX}_chr${CHR}"
        
        if [ ! -f "${PLINK_PREFIX}.bed" ]; then continue; fi
        
        N_CHR=$(awk -F'\t' -v c="$CHR_CODE" 'NR>1 && $2==c' "$CLEAN_FILE" | wc -l)
        if [ "$N_CHR" -eq 0 ]; then continue; fi
        
        apptainer exec --bind /mnt/cephfs:/mnt/cephfs --bind "$TMPDIR:/tmp,$TMPDIR:/var/tmp" "$SIF" \
            plink2 --bfile "$PLINK_PREFIX" --clump "$CLEAN_FILE" --clump-p1 5e-8 \
            --clump-r2 0.1 --clump-kb 250 --clump-snp-field SNP --clump-p-field P \
            --threads 8 --out "$CLUMP_OUT" 2>/dev/null || true
            
        if [ -f "${CLUMP_OUT}.clumps" ]; then
            awk 'NR>1 {print $3}' "${CLUMP_OUT}.clumps" > "${TMPDIR}/clump_ids_${SEX}_chr${CHR}.txt"
            grep -F -f "${TMPDIR}/clump_ids_${SEX}_chr${CHR}.txt" "$CLEAN_FILE" | \
                awk -F'\t' '{print $2"\t"$3"\t"$1"\t"$4"\t"$5"\t"$6"\t"$7}' >> "$CLUMP_MERGED"
        fi
    done
    echo "     TOTAL ${SEX}: $(tail -n +2 "$CLUMP_MERGED" | wc -l) variantes independientes."
done

# ================================================================
# FASE 2: PRS SCORING EN COHORTE COMPLETA (PYTHON)
# ================================================================
echo -e "\n[$(date +'%H:%M:%S')] INICIANDO FASE 2: CÁLCULO PRS (ZARR)"

PYTHON_SCRIPT="$TMPDIR/compute_prs_mcps.py"
cat << 'PYEOF' > "$PYTHON_SCRIPT"
import zarr
import numpy as np
import pandas as pd
import time
import gzip

# Rutas dentro de Python
zarr_main_path = "/mnt/cephfs/hot_nvme/mcps/imputed-topmed/plink_to_zar/MCPS.HIS.TOPMED.b38.FR.2602.TARGET.imputed.plink.update2.zarr"
zarr_x_path    = "/mnt/cephfs/hot_nvme/mcps/imputed-topmed/zar_files/chrX.zarr"
gwas_dir       = "/mnt/cephfs/orgs/home/angel.pacheco/gwas_miami"
ct_dir         = "/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes/ct_gwas_v2"
out_dir        = "/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes"
covar_file     = "/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes/MCPS BASELINE.csv"
link_file      = "/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes/RGN_LINK_IID.csv"

# 1. CARGAR VARIANTES C+T v2
print("  >> Cargando variantes curadas (v2)...")
df_fem = pd.read_csv(f"{ct_dir}/clumped_females_all.tsv", sep='\t')
df_mal = pd.read_csv(f"{ct_dir}/clumped_males_all.tsv", sep='\t')

# 2. QC FINAL (INFO > 0.8, MAF > 0.01)
print("  >> Aplicando QC final de imputación...")
def load_gwas_info(gwas_gz, variant_ids):
    ids_no_chr = {v.replace("chrX:", "X:").replace("chr", "", 1) if v.startswith("chr") else v: v for v in variant_ids}
    info_dict = {}
    with gzip.open(gwas_gz, 'rt') as f:
        f.readline()
        for line in f:
            cols = line.strip().split('\t')
            if cols[2] in ids_no_chr:
                a1freq = float(cols[5])
                info_dict[ids_no_chr[cols[2]]] = {'INFO': float(cols[6]), 'MAF': min(a1freq, 1 - a1freq)}
    return info_dict

for df, gwas_file in [(df_fem, "females"), (df_mal, "males")]:
    info = load_gwas_info(f"{gwas_dir}/combined-gwas-results-all-ctrl-{gwas_file}.txt.gz", df['SNP'].tolist())
    df['INFO'] = df['SNP'].map(lambda x: info.get(x, {}).get('INFO', 0))
    df['MAF']  = df['SNP'].map(lambda x: info.get(x, {}).get('MAF', 0))
    df.drop(df[(df['INFO'] < 0.8) | (df['MAF'] < 0.01)].index, inplace=True)

# 3. CARGA DE CLINICA E IDs
print("  >> Conectando IDs genómicos con clínicos...")
df_baseline = pd.read_csv(covar_file).dropna(subset=['MALE'])
patid_to_sex = dict(zip(df_baseline['PATID'].astype(str), df_baseline['MALE'].astype(int)))

df_link = pd.read_csv(link_file)
iid_to_patid = dict(zip(df_link['IID'].astype(str), df_link['PATID'].astype(str)))

sex_map = {iid: patid_to_sex[patid] for iid, patid in iid_to_patid.items() if patid in patid_to_sex}

# 4. ABRIR ZARRS Y ALINEAR
print("  >> Abriendo Zarrs y alineando 140K pacientes...")
root_main = zarr.open_group(zarr_main_path, mode='r')
root_x    = zarr.open_group(zarr_x_path, mode='r')

sample_ids_main = [s if isinstance(s, str) else s.decode() for s in root_main['sample_id'][:]]
n_samples = len(sample_ids_main)

sample_sex = np.array([sex_map.get(sid, -1) for sid in sample_ids_main])
idx_female, idx_male = np.where(sample_sex == 0)[0], np.where(sample_sex == 1)[0]

sample_ids_x = [s if isinstance(s, str) else s.decode() for s in root_x['sample_id'][:]]
main_id_to_idx = {sid: i for i, sid in enumerate(sample_ids_main)}
x_to_main_order = np.array([main_id_to_idx.get(sid, -1) for sid in sample_ids_x])

contigs_main, z_chroms_main, z_pos_main, z_alleles_main = list(root_main['contig_id'][:]), root_main['variant_contig'][:], root_main['variant_position'][:], root_main['variant_allele'][:]
contigs_x, z_chroms_x, z_pos_x, z_alleles_x = list(root_x['contig_id'][:]), root_x['variant_contig'][:], root_x['variant_position'][:], root_x['variant_allele'][:]

# 5. EXTRACCIÓN Y CÁLCULO PRS (Total y por Cromosoma)
def extract_dosages(root_zarr, global_idx):
    if 'call_DS' in root_zarr: 
        # APLASTAR A 1D PARA EVITAR SHAPE MISMATCH
        return np.ravel(root_zarr['call_DS'][global_idx]).astype(np.float32)
    
    g_raw = root_zarr['call_genotype'][global_idx, :, :]
    return np.sum(np.maximum(g_raw, 0), axis=1).astype(np.float32)

def find_variant_genotypes(chrom_str, pos, a1_effect, a2_other):
    is_x = chrom_str in ('X', '23')
    root_t, c_list, z_chroms, z_pos, z_alleles = (root_x, contigs_x, z_chroms_x, z_pos_x, z_alleles_x) if is_x else (root_main, contigs_main, z_chroms_main, z_pos_main, z_alleles_main)
    try: chrom_idx = c_list.index('X' if is_x else chrom_str)
    except ValueError: return None

    idx_start, idx_end = np.searchsorted(z_chroms, chrom_idx, side='left'), np.searchsorted(z_chroms, chrom_idx, side='right')
    matches = np.where(z_pos[idx_start:idx_end] == pos)[0]
    
    for m in matches:
        global_idx = idx_start + m
        ref, alt = z_alleles[global_idx, 0], z_alleles[global_idx, 1]
        if isinstance(ref, bytes): ref, alt = ref.decode(), alt.decode()
        
        if a1_effect == alt and a2_other == ref:
            g_sum = extract_dosages(root_t, global_idx)
        elif a1_effect == ref and a2_other == alt:
            g_sum = 2.0 - extract_dosages(root_t, global_idx)
        else: continue
            
        if is_x:
            aligned_g, valid_mask = np.zeros(n_samples, dtype=np.float32), x_to_main_order >= 0
            aligned_g[x_to_main_order[valid_mask]] = g_sum[valid_mask]
            return aligned_g
        return g_sum
    return None

# Inicializar diccionario para cromosomas y total
chromosomes = [str(i) for i in range(1, 23)] + ['X']
prs_by_chr = {chrom: np.zeros(n_samples, dtype=np.float32) for chrom in chromosomes}
prs_total = np.zeros(n_samples, dtype=np.float32)

for df, idx_group, g_name in [(df_fem, idx_female, "mujeres"), (df_mal, idx_male, "hombres")]:
    print(f"  >> Extrayendo dosis y calculando score por cromosoma ({g_name})...")
    for _, row in df.iterrows():
        chrom_raw = str(int(row['CHR'])) if row['CHR'] != 23 else 'X'
        geno = find_variant_genotypes(chrom_raw, int(row['BP']), row['A1'], row['A2'])
        if geno is not None:
            contribucion = geno[idx_group] * float(row['BETA'])
            prs_by_chr[chrom_raw][idx_group] += contribucion
            prs_total[idx_group] += contribucion

# 6. EXPORTAR RESULTADOS (Generando todas las columnas requeridas)
print("  >> Generando tabla final...")

# Construir el diccionario base con los IDs y sexo
df_data = {
    'PATID': [iid_to_patid.get(s, 'UNKNOWN') for s in sample_ids_main],
    'IID': sample_ids_main,
    'SEX': sample_sex
}

# Agregar las 23 columnas de los cromosomas (PRS_chr1 a PRS_chrX)
for chrom in chromosomes:
    df_data[f'PRS_chr{chrom}'] = prs_by_chr[chrom]

# Agregar el Total
df_data['PRS_TOTAL'] = prs_total

df_with_sex = pd.DataFrame(df_data)
df_with_sex = df_with_sex[df_with_sex['SEX'] >= 0].copy()

# Estandarizar el PRS_TOTAL (Z-score y Percentiles por sexo)
for sex_val in [0, 1]:
    mask = df_with_sex['SEX'] == sex_val
    vals = df_with_sex.loc[mask, 'PRS_TOTAL']
    df_with_sex.loc[mask, 'PRS_Z'] = (vals - vals.mean()) / vals.std() if len(vals) > 0 and vals.std() > 0 else 0.0
    df_with_sex.loc[mask, 'PERCENTILE'] = vals.rank(pct=True) * 100

# Ordenar por el riesgo total y guardar
df_with_sex.sort_values('PRS_TOTAL', ascending=False).to_csv(f"{out_dir}/prs_diabetes_mcps_gwas_final.tsv", sep='\t', index=False, float_format='%.6f')
print(f"  >> ¡PRS FINAL GUARDADO con {(len(chromosomes) + 6)} columnas en total!")
PYEOF

apptainer exec --bind /mnt/cephfs:/mnt/cephfs --bind "$TMPDIR:/tmp,$TMPDIR:/var/tmp" "$SIF" /opt/miniforge3/envs/vcf2zarr/bin/python "$PYTHON_SCRIPT"

rm -rf "$TMPDIR"
echo -e "\n[$(date +'%H:%M:%S')] ==== PIPELINE COMPLETADO EXITOSAMENTE ===="
