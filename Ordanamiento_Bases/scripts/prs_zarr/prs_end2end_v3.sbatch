#!/bin/bash
#SBATCH --job-name=prs_end2end_v4
#SBATCH --account=researchers
#SBATCH --qos=vip
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=80G
#SBATCH --time=04:00:00
#SBATCH --output=/mnt/cephfs/orgs/home/angel.pacheco/logs/prs_end2end_v4_%j.out
#SBATCH --error=/mnt/cephfs/orgs/home/angel.pacheco/logs/prs_end2end_v4_%j.err

set -euo pipefail

SIF="/mnt/cephfs/biocontainers/images/jupyter-biotools-1.2.sif"

TMPDIR="/mnt/cephfs/scratch/$USER/prs_end2end_v4_tmp/$SLURM_JOB_ID"
mkdir -p "$TMPDIR"
chmod 700 "$TMPDIR"
export TMPDIR TMP="$TMPDIR" TEMP="$TMPDIR"

# ================================================================
# RUTAS GLOBALES
# ================================================================
GWAS_DIR="/mnt/cephfs/orgs/home/angel.pacheco/gwas_miami"
PLINK_DIR="/mnt/cephfs/hot_nvme/mcps/imputed-topmed/plink_files/maximally_unrelated"
ZARR_DIR="/mnt/cephfs/hot_nvme/mcps/imputed-topmed/zar_files"
OUT_DIR="/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes"
CT_DIR="${OUT_DIR}/ct_gwas_v4"
mkdir -p "$CT_DIR"

echo "============================================================"
echo "PIPELINE END-TO-END v4: PRS DIABETES T2 (MCPS)"
echo "  - Fase 1: Clumping C+T Nativo (SNVs + Indels juntos)"
echo "  - Fase 2: PRS Scoring (Zarr VCF directo, call_DS puro)"
echo "============================================================"

# ================================================================
# FASE 1: CLUMPING + THRESHOLDING (C+T) NATIVO
# ================================================================
echo -e "\n[$(date +'%H:%M:%S')] INICIANDO FASE 1: CLUMPING (C+T)"

for SEX in females males; do
    echo "  >> Preparando stats de ${SEX} (P < 5e-8)..."
    GWAS_FILE="${GWAS_DIR}/combined-gwas-results-all-ctrl-${SEX}.txt.gz"
    CLEAN_FILE="${CT_DIR}/gwas_clean_${SEX}.tsv"

    # Filtramos por P-value sin discriminar Indels
    zcat "$GWAS_FILE" | awk -F'\t' '
        NR == 1 {
            print "SNP\tCHR\tBP\tA1\tA2\tBETA\tP"
            next
        }
        $15 < 5e-8 {
            snp = "chr"$3
            print snp"\t"$1"\t"$2"\t"$5"\t"$4"\t"$10"\t"$15
        }
    ' > "$CLEAN_FILE"

    N_VARS=$(tail -n +2 "$CLEAN_FILE" | wc -l)
    echo "     ${SEX}: ${N_VARS} variantes significativas (SNVs + Indels)"
done

for SEX in females males; do
    echo "  >> Ejecutando PLINK2 Clumping para ${SEX}..."
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

        # PLINK2 ahora procesará SNVs e Indels juntos calculando el LD perfecto
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

    N_TOTAL=$(tail -n +2 "$CLUMP_MERGED" | wc -l)
    echo "     TOTAL ${SEX}: ${N_TOTAL} variantes verdaderamente independientes."
done

# ================================================================
# FASE 2: PRS SCORING — ZARR POR CROMOSOMA (call_DS)
# ================================================================
echo -e "\n[$(date +'%H:%M:%S')] INICIANDO FASE 2: CÁLCULO PRS (ZARR POR CROMOSOMA)"

PYTHON_SCRIPT="$TMPDIR/compute_prs_mcps_v4.py"
cat << 'PYEOF' > "$PYTHON_SCRIPT"
import zarr
import numpy as np
import pandas as pd
import time
import gzip
import os

zarr_dir   = "/mnt/cephfs/hot_nvme/mcps/imputed-topmed/zar_files"
gwas_dir   = "/mnt/cephfs/orgs/home/angel.pacheco/gwas_miami"
ct_dir     = "/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes/ct_gwas_v4"
out_dir    = "/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes"
covar_file = "/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes/MCPS BASELINE.csv"
link_file  = "/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes/RGN_LINK_IID.csv"

chromosomes = [str(i) for i in range(1, 23)] + ['X']

print("  >> Cargando listas maestras C+T v4...")
df_fem = pd.read_csv(f"{ct_dir}/clumped_females_all.tsv", sep='\t')
df_mal = pd.read_csv(f"{ct_dir}/clumped_males_all.tsv", sep='\t')

print("  >> Aplicando QC final de imputación (INFO >= 0.8, MAF >= 0.01)...")
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

print("  >> Conectando IDs genómicos con clínicos...")
df_baseline = pd.read_csv(covar_file).dropna(subset=['MALE'])
patid_to_sex = dict(zip(df_baseline['PATID'].astype(str), df_baseline['MALE'].astype(int)))
df_link = pd.read_csv(link_file)
iid_to_patid = dict(zip(df_link['IID'].astype(str), df_link['PATID'].astype(str)))
sex_map = {iid: patid_to_sex[patid] for iid, patid in iid_to_patid.items() if patid in patid_to_sex}

print("  >> Abriendo Zarr por cromosoma y construyendo índices...")
ref_root = zarr.open_group(f"{zarr_dir}/chr1.zarr", mode='r')
sample_ids_ref = [s if isinstance(s, str) else s.decode() for s in ref_root['sample_id'][:]]
n_samples = len(sample_ids_ref)
ref_id_to_idx = {sid: i for i, sid in enumerate(sample_ids_ref)}

sample_sex = np.array([sex_map.get(sid, -1) for sid in sample_ids_ref])
idx_female = np.where(sample_sex == 0)[0]
idx_male   = np.where(sample_sex == 1)[0]

zarr_data = {}
for chrom in chromosomes:
    zarr_path = f"{zarr_dir}/chr{chrom}.zarr"
    if not os.path.exists(zarr_path): continue
    root = zarr.open_group(zarr_path, mode='r')
    sample_ids_this = [s if isinstance(s, str) else s.decode() for s in root['sample_id'][:]]
    reorder_map = np.array([ref_id_to_idx.get(sid, -1) for sid in sample_ids_this])
    
    zarr_data[chrom] = {
        'root': root, 'positions': root['variant_position'][:], 'alleles': root['variant_allele'][:],
        'reorder': reorder_map, 'same_order': np.array_equal(reorder_map, np.arange(n_samples))
    }

def extract_dosages(root_zarr, global_idx):
    if 'call_DS' in root_zarr:
        return np.ravel(root_zarr['call_DS'][global_idx]).astype(np.float32)
    return np.sum(np.maximum(root_zarr['call_genotype'][global_idx, :, :], 0), axis=1).astype(np.float32)

def find_variant_genotypes(chrom_str, pos, a1_effect, a2_other):
    chrom_key = 'X' if chrom_str in ('X', '23') else chrom_str
    if chrom_key not in zarr_data: return None
    zd = zarr_data[chrom_key]
    
    matches = np.where(zd['positions'] == pos)[0]
    for idx in matches:
        ref, alt = zd['alleles'][idx, 0], zd['alleles'][idx, 1]
        if isinstance(ref, bytes): ref, alt = ref.decode(), alt.decode()
        
        if a1_effect == alt and a2_other == ref: g_sum = extract_dosages(zd['root'], idx)
        elif a1_effect == ref and a2_other == alt: g_sum = 2.0 - extract_dosages(zd['root'], idx)
        else: continue
            
        if not zd['same_order']:
            aligned = np.zeros(n_samples, dtype=np.float32)
            valid = zd['reorder'] >= 0
            aligned[zd['reorder'][valid]] = g_sum[valid]
            return aligned
        return g_sum
    return None

prs_by_chr = {c: np.zeros(n_samples, dtype=np.float32) for c in chromosomes}
prs_total  = np.zeros(n_samples, dtype=np.float32)

for df, idx_group, g_name in [(df_fem, idx_female, "mujeres"), (df_mal, idx_male, "hombres")]:
    t0 = time.time()
    print(f"\n  >> Calculando PRS para {g_name} ({len(df)} variantes)...")
    for _, row in df.iterrows():
        chrom_raw = str(int(row['CHR'])) if row['CHR'] != 23 else 'X'
        geno = find_variant_genotypes(chrom_raw, int(row['BP']), str(row['A1']), str(row['A2']))
        if geno is not None:
            contribucion = geno[idx_group] * float(row['BETA'])
            prs_by_chr[chrom_raw][idx_group] += contribucion
            prs_total[idx_group] += contribucion
    print(f"     Completado en {time.time() - t0:.1f}s")

print("\n  >> Generando tabla y estandarizando Z-Score...")
df_data = {'PATID': [iid_to_patid.get(s, 'UNKNOWN') for s in sample_ids_ref], 'IID': sample_ids_ref, 'SEX': sample_sex}
for chrom in chromosomes: df_data[f'PRS_chr{chrom}'] = prs_by_chr[chrom]
df_data['PRS_TOTAL'] = prs_total

df_out = pd.DataFrame(df_data)
df_out = df_out[df_out['SEX'] >= 0].copy()

for sex_val in [0, 1]:
    mask = df_out['SEX'] == sex_val
    vals = df_out.loc[mask, 'PRS_TOTAL']
    if len(vals) > 0 and vals.std() > 0:
        df_out.loc[mask, 'PRS_Z'] = (vals - vals.mean()) / vals.std()
        df_out.loc[mask, 'PERCENTILE'] = vals.rank(pct=True) * 100

out_file = f"{out_dir}/prs_diabetes_mcps_gwas_v4_final.tsv"
df_out.sort_values('PRS_TOTAL', ascending=False).to_csv(out_file, sep='\t', index=False, float_format='%.6f')
print(f"  >> PRS v4 GUARDADO: {out_file}")
PYEOF

apptainer exec --bind /mnt/cephfs:/mnt/cephfs --bind "$TMPDIR:/tmp,$TMPDIR:/var/tmp" "$SIF" \
    /opt/miniforge3/envs/vcf2zarr/bin/python "$PYTHON_SCRIPT"

rm -rf "$TMPDIR"
echo -e "\n[$(date +'%H:%M:%S')] ==== PIPELINE v4 COMPLETADO EXITOSAMENTE ===="
