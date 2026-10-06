#!/bin/bash
#SBATCH --job-name=prs_cancer
#SBATCH --account=researchers
#SBATCH --qos=vip
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=02:30:00
#SBATCH --output=/mnt/cephfs/orgs/home/angel.pacheco/logs/prs_cancer_%j.out
#SBATCH --error=/mnt/cephfs/orgs/home/angel.pacheco/logs/prs_cancer_%j.err

set -euo pipefail

SIF="/mnt/cephfs/biocontainers/images/jupyter-biotools-1.5.sif"
PYTHON="/opt/miniforge3/envs/vcf2zarr/bin/python"

TMPDIR="/mnt/cephfs/scratch/$USER/prs_cancer_tmp/$SLURM_JOB_ID"
mkdir -p "$TMPDIR"
chmod 700 "$TMPDIR"
export TMPDIR TMP="$TMPDIR" TEMP="$TMPDIR"

PGS_DIR="/mnt/cephfs/orgs/home/angel.pacheco/Ordanamiento_Bases/pgs_scores"
OUT_DIR="/mnt/cephfs/orgs/home/angel.pacheco/prs_cancer"
mkdir -p "$OUT_DIR"

echo "================================================================"
echo "PIPELINE PRS CANCER ONCOLÓGICO (MCPS) — 6 SCORES"
echo "  PGS003766 (Próstata)  | PGS000028 (Mama)"
echo "  PGS000078 (Pulmón)    | PGS000073 (Cervical)"
echo "  PGS002265 (GI/CRC)    | PGS000874 (Hematológico/CLL)"
echo "Contenedor: jupyter-biotools-1.5.sif | Python: vcf2zarr env"
echo "================================================================"

PYTHON_SCRIPT="$TMPDIR/compute_prs_cancer.py"
cat << 'PYEOF' > "$PYTHON_SCRIPT"
import zarr
import numpy as np
import pandas as pd
import gzip
import sys
import time

# ================================================================
# RUTAS
# ================================================================
zarr_dir   = "/mnt/cephfs/hot_nvme/mcps/imputed-topmed/zar_files"
pgs_dir    = "/mnt/cephfs/orgs/home/angel.pacheco/Ordanamiento_Bases/pgs_scores"
baseline_f = "/mnt/cephfs/orgs/home/angel.pacheco/Ordanamiento_Bases/Bases_DNAnexus/MCPS BASELINE.csv"
out_dir    = "/mnt/cephfs/orgs/home/angel.pacheco/prs_cancer"

# ================================================================
# SNPs AMBIGUOS — excluidos siempre, ya sea por OA del scoring file
# o por OA inferida desde los alelos REF/ALT del zarr.
# ================================================================
AMBIGUOUS_PAIRS = {('A','T'), ('T','A'), ('C','G'), ('G','C')}

VALID_ALLELES = {'A', 'T', 'C', 'G'}

# ================================================================
# CONFIGURACIÓN: 6 SCORES
# sex: 'all' | 'male_only' | 'female_only'
# ================================================================
PGS_CONFIGS = {
    'PGS003766': {
        'file':  f"{pgs_dir}/PGS003766_hmPOS_GRCh38.txt.gz",
        'sex':   'male_only',
        'site':  'PROSTATE',
        'label': 'Próstata (Conti et al. 2021)',
    },
    'PGS000028': {
        'file':  f"{pgs_dir}/PGS000028_hmPOS_GRCh38.txt.gz",
        'sex':   'female_only',
        'site':  'BREAST',
        'label': 'Mama (Mavaddat et al. 2019)',
    },
    'PGS000078': {
        'file':  f"{pgs_dir}/PGS000078_hmPOS_GRCh38.txt.gz",
        'sex':   'all',
        'site':  'LUNG',
        'label': 'Pulmón (McKay et al. 2017)',
    },
    'PGS000073': {
        'file':  f"{pgs_dir}/PGS000073_hmPOS_GRCh38.txt.gz",
        'sex':   'female_only',
        'site':  'CERVICAL',
        'label': 'Cérvix (Chen et al. 2024)',
    },
    'PGS002265': {
        'file':  f"{pgs_dir}/PGS002265_hmPOS_GRCh38.txt.gz",
        'sex':   'all',
        'site':  'GI',
        'label': 'GI/Colorrectal (Huyghe et al. AJHG 2020)',
    },
    'PGS000874': {
        'file':  f"{pgs_dir}/PGS000874_hmPOS_GRCh38.txt.gz",
        'sex':   'all',
        'site':  'HEMATOLOGIC',
        'label': 'Hematológico/CLL (Berndt et al. Blood 2018)',
    },
}

# Contadores de variantes por PGS
counters = {
    pgs_id: {
        'total_in_file':      0,
        'ambiguous_file':     0,   # par (EA,OA) ambiguo desde el archivo
        'ambiguous_zarr':     0,   # par (EA,inferred_OA) ambiguo al matchear
        'no_hm_pos':          0,
        'usable':             0,
        'found_exact':        0,
        'allele_mismatch':    0,
        'not_in_zarr':        0,
    }
    for pgs_id in PGS_CONFIGS
}

# ================================================================
# FUNCIÓN: CARGAR SCORING FILE
# Maneja dos formatos:
#   (A) Con 'other_allele' explícito (PGS000028, PGS000073, etc.)
#   (B) Sin 'other_allele', solo 'hm_inferOtherAllele' (PGS002265, PGS000874)
#
# Para (B): OA se almacena como None cuando hm_inferOtherAllele es
# multi-valor o vacío. El filtro de ambigüedad se aplica en el matching
# usando los alelos REF/ALT del zarr.
# ================================================================
def load_pgs_file(path, pgs_id):
    rows = []
    cnt  = counters[pgs_id]

    with gzip.open(path, 'rt') as fh:
        for line in fh:
            if not line.startswith('#'):
                header = line.strip().split('\t')
                break

        # Detectar columna de posición harmonizada
        has_hm_chr = 'hm_chr'   in header
        has_hm_pos = 'hm_pos'   in header
        has_oa     = 'other_allele' in header
        has_infer  = 'hm_inferOtherAllele' in header

        idx_chr = header.index('hm_chr')   if has_hm_chr else header.index('chr_name')
        idx_pos = header.index('hm_pos')   if has_hm_pos else header.index('chr_position')
        idx_ea  = header.index('effect_allele')
        idx_bw  = header.index('effect_weight')
        idx_oa  = header.index('other_allele')        if has_oa    else None
        idx_inf = header.index('hm_inferOtherAllele') if has_infer else None

        for line in fh:
            cols = line.strip().split('\t')
            cnt['total_in_file'] += 1

            chr_val = cols[idx_chr].strip() if idx_chr < len(cols) else ''
            pos_val = cols[idx_pos].strip() if idx_pos < len(cols) else ''
            ea_val  = cols[idx_ea ].strip().upper() if idx_ea  < len(cols) else ''
            bw_val  = cols[idx_bw ].strip() if idx_bw  < len(cols) else ''

            # Extraer OA: primero desde other_allele, luego desde hm_inferOtherAllele
            oa_val = None
            if idx_oa is not None and idx_oa < len(cols):
                v = cols[idx_oa].strip().upper()
                if v in VALID_ALLELES:
                    oa_val = v
            elif idx_inf is not None and idx_inf < len(cols):
                v = cols[idx_inf].strip().upper()
                if v in VALID_ALLELES:   # solo si es un único nucleótido
                    oa_val = v
                # Si es multi-valor (A/T, G/C, …) → oa_val queda None

            # Descartar variantes sin posición harmonizada
            if not chr_val or not pos_val or chr_val in ('', 'NA', 'nan', '.'):
                cnt['no_hm_pos'] += 1
                continue

            # Filtrar EA inválidos
            if ea_val not in VALID_ALLELES:
                cnt['no_hm_pos'] += 1
                continue

            # Filtro de ambigüedad cuando OA es conocida desde el archivo
            if oa_val is not None and (ea_val, oa_val) in AMBIGUOUS_PAIRS:
                cnt['ambiguous_file'] += 1
                continue

            try:
                pos_int = int(float(pos_val))
                beta    = float(bw_val)
            except (ValueError, TypeError):
                cnt['no_hm_pos'] += 1
                continue

            chr_clean = chr_val.replace('chr', '').replace('CHR', '')
            cnt['usable'] += 1
            rows.append({
                'CHR':  chr_clean,
                'POS':  pos_int,
                'EA':   ea_val,
                'OA':   oa_val,      # puede ser None
                'BETA': beta,
            })

    return pd.DataFrame(rows)

# ================================================================
# CARGAR TODOS LOS SCORING FILES
# ================================================================
print("\n[1] Cargando scoring files del PGS Catalog...", flush=True)
for pgs_id, cfg in PGS_CONFIGS.items():
    print(f"  {pgs_id} ({cfg['label']}):", flush=True)
    cfg['df']     = load_pgs_file(cfg['file'], pgs_id)
    cfg['by_chr'] = {c: g for c, g in cfg['df'].groupby('CHR')}
    cnt = counters[pgs_id]
    pct_amb = (cnt['ambiguous_file'] / cnt['total_in_file'] * 100) if cnt['total_in_file'] else 0
    n_oa_none = (cfg['df']['OA'].isna()).sum()
    print(f"    total en archivo:        {cnt['total_in_file']}", flush=True)
    print(f"    sin posición hm:         {cnt['no_hm_pos']}", flush=True)
    print(f"    ambiguos (OA conocida):  {cnt['ambiguous_file']} ({pct_amb:.1f}%)", flush=True)
    print(f"    usables para PRS:        {cnt['usable']}  (OA desconocida: {n_oa_none})", flush=True)

# ================================================================
# CARGAR SEXO DESDE BASELINE
# ================================================================
print("\n[2] Cargando sexo biológico de BASELINE...", flush=True)
df_base = pd.read_csv(baseline_f, usecols=['PATID', 'MALE'])
df_base['PATID'] = df_base['PATID'].astype(str)
patid_to_male = dict(zip(df_base['PATID'], df_base['MALE'].astype(float)))

# ================================================================
# LISTA DE MUESTRAS DESDE zarr MASTER
# ================================================================
print("\n[3] Leyendo lista de muestras del zarr (chr1)...", flush=True)
root_master    = zarr.open_group(f"{zarr_dir}/chr1.zarr", mode='r')
sample_ids_raw = root_master['sample_id'][:]
sample_ids     = [s if isinstance(s, str) else s.decode() for s in sample_ids_raw]
n_samples      = len(sample_ids)
print(f"    Total muestras: {n_samples:,}", flush=True)

def iid_to_patid(iid):
    parts = iid.split('_')
    return parts[1] if len(parts) >= 3 else iid

sample_patids = [iid_to_patid(iid) for iid in sample_ids]
sample_sex    = np.array([patid_to_male.get(pid, np.nan) for pid in sample_patids],
                          dtype=np.float32)

idx_male   = np.where(sample_sex == 1.0)[0]
idx_female = np.where(sample_sex == 0.0)[0]
idx_unk    = np.where(np.isnan(sample_sex))[0]
print(f"    Hombres: {len(idx_male):,} | Mujeres: {len(idx_female):,} "
      f"| Sexo desconocido: {len(idx_unk):,}", flush=True)

# ================================================================
# INICIALIZAR ARRAYS PRS
# ================================================================
prs_arrays = {pgs_id: np.zeros(n_samples, dtype=np.float64) for pgs_id in PGS_CONFIGS}

# ================================================================
# LOOP PRINCIPAL — UN ZARR POR CROMOSOMA, TODOS LOS SCORES
# ================================================================
print("\n[4] Calculando PRS por cromosoma...", flush=True)
chromosomes = [str(c) for c in range(1, 23)] + ['X']

for chrom in chromosomes:
    scores_on_chrom = [pid for pid, cfg in PGS_CONFIGS.items() if chrom in cfg['by_chr']]
    if not scores_on_chrom:
        continue

    t0   = time.time()
    root = zarr.open_group(f"{zarr_dir}/chr{chrom}.zarr", mode='r')
    z_pos     = root['variant_position'][:]   # ordenado — verificado
    z_alleles = root['variant_allele'][:]     # (n_var, 2) = [REF, ALT]

    for pgs_id in scores_on_chrom:
        df_chr = PGS_CONFIGS[pgs_id]['by_chr'][chrom]
        cnt    = counters[pgs_id]

        for _, row in df_chr.iterrows():
            target_pos = int(row['POS'])
            ea         = str(row['EA']).upper()
            oa_known   = row['OA']        # str o None
            beta       = float(row['BETA'])

            # Búsqueda binaria (variant_position está ordenado)
            lo = np.searchsorted(z_pos, target_pos, side='left')
            hi = np.searchsorted(z_pos, target_pos, side='right')

            if lo == hi:
                cnt['not_in_zarr'] += 1
                continue

            found            = False
            any_pos_match    = False

            for vi in range(lo, hi):
                any_pos_match = True
                ref = z_alleles[vi, 0]
                alt = z_alleles[vi, 1]
                if isinstance(ref, bytes): ref = ref.decode()
                if isinstance(alt, bytes): alt = alt.decode()
                ref = ref.upper()
                alt = alt.upper()

                # Solo SNPs bialélicos simples
                if ref not in VALID_ALLELES or alt not in VALID_ALLELES:
                    continue

                # Intentar alinear EA con ALT o REF
                if ea == alt:
                    # EA = ALT → inferred_OA = REF
                    inferred_oa = ref
                    # Filtro de ambigüedad usando alelos del zarr
                    if (ea, inferred_oa) in AMBIGUOUS_PAIRS:
                        cnt['ambiguous_zarr'] += 1
                        continue
                    # Si el archivo proporcionó OA, verificar consistencia
                    if oa_known is not None and oa_known != ref:
                        cnt['allele_mismatch'] += 1
                        continue
                    ds = np.ravel(root['call_DS'][vi]).astype(np.float64)
                    prs_arrays[pgs_id] += ds * beta
                    cnt['found_exact'] += 1
                    found = True
                    break

                elif ea == ref:
                    # EA = REF → inferred_OA = ALT → flip dosage
                    inferred_oa = alt
                    if (ea, inferred_oa) in AMBIGUOUS_PAIRS:
                        cnt['ambiguous_zarr'] += 1
                        continue
                    if oa_known is not None and oa_known != alt:
                        cnt['allele_mismatch'] += 1
                        continue
                    ds = np.ravel(root['call_DS'][vi]).astype(np.float64)
                    prs_arrays[pgs_id] += (2.0 - ds) * beta
                    cnt['found_exact'] += 1
                    found = True
                    break

            if not found:
                if any_pos_match:
                    cnt['allele_mismatch'] += 1
                else:
                    cnt['not_in_zarr'] += 1

    elapsed = time.time() - t0
    print(f"  chr{chrom}: scores={scores_on_chrom} | {elapsed:.1f}s", flush=True)

# ================================================================
# REPORTE DE COBERTURA
# ================================================================
print("\n" + "="*65, flush=True)
print("REPORTE DE COBERTURA DE VARIANTES POR SCORE", flush=True)
print("="*65, flush=True)
for pgs_id, cfg in PGS_CONFIGS.items():
    cnt   = counters[pgs_id]
    total = cnt['total_in_file']
    usbl  = cnt['usable']
    found = cnt['found_exact']
    pct   = found / usbl * 100 if usbl > 0 else 0

    print(f"\n{pgs_id} — {cfg['label']}:", flush=True)
    print(f"  Total en scoring file:             {total:>5}", flush=True)
    print(f"  Sin posición hm / EA inválido:     {cnt['no_hm_pos']:>5}", flush=True)
    print(f"  Excluidos ambiguos (OA en archivo):{cnt['ambiguous_file']:>5}", flush=True)
    print(f"  Usables tras filtro inicial:       {usbl:>5}", flush=True)
    print(f"  Excluidos ambiguos (OA zarr):      {cnt['ambiguous_zarr']:>5}", flush=True)
    print(f"  Encontrados y usados en PRS:       {found:>5}  ({pct:.1f}%)", flush=True)
    print(f"  No en zarr (posición):             {cnt['not_in_zarr']:>5}", flush=True)
    print(f"  Mismatch de alelos:                {cnt['allele_mismatch']:>5}", flush=True)
print("="*65, flush=True)

# ================================================================
# FILTRO DE SEXO (NaN donde sexo es incompatible con el score)
# ================================================================
print("\n[5] Aplicando restricción de sexo...", flush=True)
for pgs_id, cfg in PGS_CONFIGS.items():
    sf = cfg['sex']
    if sf == 'male_only':
        prs_arrays[pgs_id][idx_female] = np.nan
        prs_arrays[pgs_id][idx_unk]    = np.nan
    elif sf == 'female_only':
        prs_arrays[pgs_id][idx_male] = np.nan
        prs_arrays[pgs_id][idx_unk]  = np.nan
    print(f"  {pgs_id} ({sf}): OK", flush=True)

# ================================================================
# ESTANDARIZAR (Z-SCORE Y PERCENTIL)
# ================================================================
print("\n[6] Estandarizando PRS...", flush=True)

sex_masks = {
    'male_only':   (sample_sex == 1.0),
    'female_only': (sample_sex == 0.0),
    'all':         (sample_sex == 0.0) | (sample_sex == 1.0),
}

df_out = pd.DataFrame({'IID': sample_ids, 'PATID': sample_patids, 'MALE': sample_sex})

for pgs_id, cfg in PGS_CONFIGS.items():
    raw_col = f"PRS_{pgs_id}_{cfg['site']}"
    z_col   = f"PRS_Z_{pgs_id}_{cfg['site']}"
    pct_col = f"PRS_PCT_{pgs_id}_{cfg['site']}"

    df_out[raw_col] = prs_arrays[pgs_id]

    mask  = sex_masks[cfg['sex']]
    vals  = pd.Series(prs_arrays[pgs_id])
    valid = mask & ~np.isnan(prs_arrays[pgs_id])

    df_out[z_col]   = np.nan
    df_out[pct_col] = np.nan

    if valid.sum() > 1:
        mu    = vals[valid].mean()
        sigma = vals[valid].std()
        if sigma > 0:
            df_out.loc[valid, z_col] = (vals[valid] - mu) / sigma
        df_out.loc[valid, pct_col] = vals[valid].rank(pct=True) * 100

    print(f"  {pgs_id}: {valid.sum():,} muestras con PRS válido", flush=True)

# ================================================================
# GUARDAR
# ================================================================
out_path = f"{out_dir}/prs_cancer_mcps_all_samples.tsv"
df_out.to_csv(out_path, sep='\t', index=False, float_format='%.6f', na_rep='NA')
print(f"\n[7] Tabla guardada: {out_path}", flush=True)
print(f"    Filas: {len(df_out):,} | Columnas: {len(df_out.columns)}", flush=True)
print(f"    Columnas: {list(df_out.columns)}", flush=True)
print("\n==== PRS CANCER SCORING COMPLETADO ====", flush=True)
PYEOF

echo "[$(date +'%H:%M:%S')] Iniciando cálculo PRS dentro del contenedor 1.5..."

apptainer exec \
    --bind /mnt/cephfs:/mnt/cephfs \
    --bind "$TMPDIR:/tmp,$TMPDIR:/var/tmp" \
    "$SIF" \
    "$PYTHON" "$PYTHON_SCRIPT"

rm -rf "$TMPDIR"
echo "[$(date +'%H:%M:%S')] ==== JOB COMPLETADO EXITOSAMENTE ===="
