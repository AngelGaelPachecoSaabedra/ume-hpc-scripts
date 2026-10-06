#!/bin/bash
#SBATCH --job-name=ld_check_indels
#SBATCH --account=researchers
#SBATCH --qos=vip
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=/mnt/cephfs/orgs/home/angel.pacheco/logs/ld_check_indels_%j.out
#SBATCH --error=/mnt/cephfs/orgs/home/angel.pacheco/logs/ld_check_indels_%j.err

set -euo pipefail

SIF="/mnt/cephfs/biocontainers/images/jupyter-biotools-1.2.sif"

TMPDIR="/mnt/cephfs/scratch/$USER/ld_check_tmp/$SLURM_JOB_ID"
mkdir -p "$TMPDIR"
chmod 700 "$TMPDIR"
export TMPDIR TMP="$TMPDIR" TEMP="$TMPDIR"

echo "============================================================"
echo "VERIFICACIÓN DE LD: INDELS vs SNVs CLUMPED"
echo "  Panel: Zarr VCF directo (140K, filtrado a 79K unrelated)"
echo "  Ventana: ±250 kb | Umbral: r² > 0.1"
echo "============================================================"

PYTHON_SCRIPT="$TMPDIR/ld_check_indels.py"
cat << 'PYEOF' > "$PYTHON_SCRIPT"
import zarr
import numpy as np
import pandas as pd
import os
import time

# ================================================================
# RUTAS
# ================================================================
zarr_dir   = "/mnt/cephfs/hot_nvme/mcps/imputed-topmed/zar_files"
ct_dir     = "/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes/ct_gwas_v3"
out_dir    = "/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes"
link_file  = "/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes/RGN_LINK_IID.csv"

# Archivo con los IIDs de los 79,612 maximally_unrelated
# Lo extraemos del .fam de PLINK
plink_dir  = "/mnt/cephfs/hot_nvme/mcps/imputed-topmed/plink_files/maximally_unrelated"

LD_WINDOW  = 250_000   # ±250 kb (mismo que clumping)
LD_THRESH  = 0.1       # r² umbral (mismo que --clump-r2)

# ================================================================
# 1. IDENTIFICAR LOS 79K MAXIMALLY UNRELATED EN EL ZARR
# ================================================================
print("=" * 60)
print("PASO 1: Identificando los 79,612 individuos en el Zarr")
print("=" * 60)

# Leer IIDs del .fam de PLINK (los 79K unrelated)
fam_file = f"{plink_dir}/mcps-freeze150k_qcd_chr1_ivs.fam"
fam_iids = set()
with open(fam_file) as f:
    for line in f:
        cols = line.strip().split()
        fam_iids.add(cols[1])  # Columna 2 del .fam = IID

print(f"  IIDs en .fam (unrelated): {len(fam_iids):,}")

# Abrir chr1 para obtener el orden de muestras del Zarr
ref_root = zarr.open_group(f"{zarr_dir}/chr1.zarr", mode='r')
sample_ids = [s if isinstance(s, str) else s.decode()
              for s in ref_root['sample_id'][:]]
n_total = len(sample_ids)

# Índices de los 79K dentro del Zarr
unrelated_idx = np.array([i for i, sid in enumerate(sample_ids) if sid in fam_iids])
print(f"  Muestras en Zarr: {n_total:,}")
print(f"  Unrelated encontrados en Zarr: {len(unrelated_idx):,}")

if len(unrelated_idx) < 70000:
    print("  [WARN] Menos de 70K encontrados. Verificar formato de IIDs.")
    print(f"  Ejemplo IID .fam: {list(fam_iids)[:3]}")
    print(f"  Ejemplo IID Zarr: {sample_ids[:3]}")

# ================================================================
# 2. CARGAR VARIANTES C+T v3 Y SEPARAR SNVs/INDELS
# ================================================================
print("\n" + "=" * 60)
print("PASO 2: Cargando variantes C+T v3")
print("=" * 60)

results_all = []

for sex in ["females", "males"]:
    print(f"\n  --- {sex.upper()} ---")
    df = pd.read_csv(f"{ct_dir}/clumped_{sex}_all.tsv", sep='\t')
    df['A1'] = df['A1'].astype(str)
    df['A2'] = df['A2'].astype(str)
    df['is_indel'] = (df['A1'].str.len() != 1) | (df['A2'].str.len() != 1)

    indels = df[df['is_indel']].copy()
    snvs   = df[~df['is_indel']].copy()

    print(f"  Total: {len(df)} variantes ({len(snvs)} SNVs + {len(indels)} indels)")

    if len(indels) == 0:
        print("  No hay indels para verificar, saltando.")
        continue

    # ================================================================
    # 3. PARA CADA INDEL, BUSCAR LD CON VARIANTES CERCANAS
    # ================================================================
    print(f"\n  Verificando LD de {len(indels)} indels contra variantes cercanas...")

    # Agrupar por cromosoma para abrir cada Zarr una sola vez
    for chrom_val in sorted(indels['CHR'].unique()):
        chrom_str = 'X' if chrom_val == 23 else str(int(chrom_val))
        zarr_path = f"{zarr_dir}/chr{chrom_str}.zarr"

        if not os.path.exists(zarr_path):
            print(f"    [WARN] {zarr_path} no existe, saltando")
            continue

        root = zarr.open_group(zarr_path, mode='r')
        positions = root['variant_position'][:]
        alleles   = root['variant_allele'][:]
        n_vars    = len(positions)

        # Índices de unrelated en ESTE Zarr (verificar orden de muestras)
        sample_ids_this = [s if isinstance(s, str) else s.decode()
                           for s in root['sample_id'][:]]
        if sample_ids_this == sample_ids:
            ur_idx = unrelated_idx  # Mismo orden
        else:
            sid_set = {sid: i for i, sid in enumerate(sample_ids_this)}
            ur_idx = np.array([sid_set[sid] for sid in sample_ids
                               if sid in sid_set and sid in fam_iids])

        indels_chr = indels[indels['CHR'] == chrom_val]
        snvs_chr   = snvs[snvs['CHR'] == chrom_val]

        # También buscar contra TODOS los vecinos del Zarr (no solo SNVs clumped)
        # Esto detecta LD con variantes que no pasaron el clumping original

        print(f"\n    chr{chrom_str}: {len(indels_chr)} indels a verificar, "
              f"{len(snvs_chr)} SNVs clumped en este cromosoma")

        for _, indel_row in indels_chr.iterrows():
            indel_pos = int(indel_row['BP'])
            indel_a1  = indel_row['A1']
            indel_a2  = indel_row['A2']
            indel_snp = indel_row['SNP']

            # Buscar el indel en el Zarr
            indel_matches = np.where(positions == indel_pos)[0]
            indel_zarr_idx = None
            indel_flip = False

            for m in indel_matches:
                ref = alleles[m, 0]
                alt = alleles[m, 1]
                if isinstance(ref, bytes):
                    ref, alt = ref.decode(), alt.decode()
                if indel_a1 == alt and indel_a2 == ref:
                    indel_zarr_idx = m
                    indel_flip = False
                    break
                elif indel_a1 == ref and indel_a2 == alt:
                    indel_zarr_idx = m
                    indel_flip = True
                    break

            if indel_zarr_idx is None:
                print(f"      [MISS] {indel_snp} no encontrado en Zarr")
                results_all.append({
                    'sex': sex, 'indel_snp': indel_snp, 'indel_chr': chrom_str,
                    'indel_pos': indel_pos, 'indel_a1': indel_a1, 'indel_a2': indel_a2,
                    'indel_beta': indel_row['BETA'], 'indel_p': indel_row['P'],
                    'status': 'NOT_IN_ZARR', 'max_r2': np.nan,
                    'max_r2_variant': '', 'max_r2_pos': np.nan,
                    'n_neighbors_checked': 0, 'n_neighbors_ld_high': 0
                })
                continue

            # Extraer dosificación del indel (solo los 79K unrelated)
            ds_indel = np.ravel(root['call_DS'][indel_zarr_idx])[ur_idx].astype(np.float64)
            if indel_flip:
                ds_indel = 2.0 - ds_indel

            # Varianza del indel
            var_indel = np.var(ds_indel)
            if var_indel < 1e-10:
                print(f"      [MONO] {indel_snp} es monomórfico en los 79K")
                results_all.append({
                    'sex': sex, 'indel_snp': indel_snp, 'indel_chr': chrom_str,
                    'indel_pos': indel_pos, 'indel_a1': indel_a1, 'indel_a2': indel_a2,
                    'indel_beta': indel_row['BETA'], 'indel_p': indel_row['P'],
                    'status': 'MONOMORPHIC', 'max_r2': 0.0,
                    'max_r2_variant': '', 'max_r2_pos': np.nan,
                    'n_neighbors_checked': 0, 'n_neighbors_ld_high': 0
                })
                continue

            # Encontrar TODAS las variantes del Zarr dentro de ±250kb
            window_mask = (positions >= indel_pos - LD_WINDOW) & \
                          (positions <= indel_pos + LD_WINDOW) & \
                          (np.arange(n_vars) != indel_zarr_idx)
            neighbor_indices = np.where(window_mask)[0]

            # ---- Verificar contra SNVs clumped del mismo cromosoma ----
            snv_positions_set = set(snvs_chr['BP'].astype(int).tolist())

            max_r2 = 0.0
            max_r2_var = ""
            max_r2_pos = 0
            n_checked = 0
            n_high_ld = 0
            high_ld_details = []

            for ni in neighbor_indices:
                n_ref = alleles[ni, 0]
                n_alt = alleles[ni, 1]
                if isinstance(n_ref, bytes):
                    n_ref, n_alt = n_ref.decode(), n_alt.decode()

                ds_neighbor = np.ravel(root['call_DS'][ni])[ur_idx].astype(np.float64)

                var_neighbor = np.var(ds_neighbor)
                if var_neighbor < 1e-10:
                    continue

                n_checked += 1

                # Calcular r²
                cov = np.mean(ds_indel * ds_neighbor) - np.mean(ds_indel) * np.mean(ds_neighbor)
                r2 = (cov ** 2) / (var_indel * var_neighbor)

                neighbor_pos = int(positions[ni])
                neighbor_id = f"chr{chrom_str}:{neighbor_pos}:{n_ref}:{n_alt}"
                is_clumped_snv = neighbor_pos in snv_positions_set

                if r2 > max_r2:
                    max_r2 = r2
                    max_r2_var = neighbor_id
                    max_r2_pos = neighbor_pos

                if r2 > LD_THRESH:
                    n_high_ld += 1
                    tag = "CLUMPED_SNV" if is_clumped_snv else "OTHER"
                    high_ld_details.append({
                        'variant': neighbor_id,
                        'pos': neighbor_pos,
                        'r2': r2,
                        'type': tag,
                        'is_indel': len(n_ref) != 1 or len(n_alt) != 1
                    })

            # Clasificar resultado
            if max_r2 > LD_THRESH:
                # ¿Está en LD con un SNV ya clumped?
                clumped_ld = [d for d in high_ld_details if d['type'] == 'CLUMPED_SNV']
                if clumped_ld:
                    status = "REDUNDANT_WITH_CLUMPED"
                    top_clumped = max(clumped_ld, key=lambda x: x['r2'])
                    print(f"      [REDUNDANTE] {indel_snp} r²={top_clumped['r2']:.4f} "
                          f"con SNV clumped {top_clumped['variant']}")
                else:
                    status = "LD_WITH_OTHER"
                    print(f"      [LD_OTRO] {indel_snp} max_r²={max_r2:.4f} con {max_r2_var} "
                          f"(no es SNV clumped)")
            else:
                status = "INDEPENDENT"
                print(f"      [OK] {indel_snp} independiente (max_r²={max_r2:.4f}, "
                      f"{n_checked} vecinos revisados)")

            results_all.append({
                'sex': sex, 'indel_snp': indel_snp, 'indel_chr': chrom_str,
                'indel_pos': indel_pos, 'indel_a1': indel_a1, 'indel_a2': indel_a2,
                'indel_beta': indel_row['BETA'], 'indel_p': indel_row['P'],
                'status': status, 'max_r2': max_r2,
                'max_r2_variant': max_r2_var, 'max_r2_pos': max_r2_pos,
                'n_neighbors_checked': n_checked, 'n_neighbors_ld_high': n_high_ld
            })

# ================================================================
# 4. RESUMEN Y EXPORTAR
# ================================================================
print("\n" + "=" * 60)
print("RESUMEN DE VERIFICACIÓN DE LD")
print("=" * 60)

df_results = pd.DataFrame(results_all)

if len(df_results) == 0:
    print("  No se encontraron indels para verificar.")
else:
    for sex in df_results['sex'].unique():
        df_sex = df_results[df_results['sex'] == sex]
        print(f"\n  --- {sex.upper()} ({len(df_sex)} indels) ---")
        for status in ['INDEPENDENT', 'REDUNDANT_WITH_CLUMPED', 'LD_WITH_OTHER',
                        'NOT_IN_ZARR', 'MONOMORPHIC']:
            n = sum(df_sex['status'] == status)
            if n > 0:
                print(f"    {status}: {n}")

        # Detalle de redundantes
        redundant = df_sex[df_sex['status'] == 'REDUNDANT_WITH_CLUMPED']
        if len(redundant) > 0:
            print(f"\n    Indels REDUNDANTES con SNVs ya clumped (deberían eliminarse):")
            for _, r in redundant.iterrows():
                print(f"      {r['indel_snp']} (beta={r['indel_beta']:.4f}, p={r['indel_p']:.2e}) "
                      f"→ r²={r['max_r2']:.4f} con {r['max_r2_variant']}")

    # Guardar reporte completo
    report_file = f"{out_dir}/ld_check_indels_v3_report.tsv"
    df_results.to_csv(report_file, sep='\t', index=False, float_format='%.6f')
    print(f"\n  Reporte completo: {report_file}")

    # Generar lista de indels que SÍ son independientes (para usar en el PRS)
    independent = df_results[df_results['status'] == 'INDEPENDENT']
    print(f"\n  CONCLUSIÓN: {len(independent)} de {len(df_results)} indels son "
          f"genuinamente independientes (r² < {LD_THRESH} con todos los vecinos).")

    redundant_total = df_results[df_results['status'] == 'REDUNDANT_WITH_CLUMPED']
    if len(redundant_total) > 0:
        print(f"  ⚠ {len(redundant_total)} indels son REDUNDANTES con SNVs ya clumped "
              f"y deberían eliminarse del PRS para evitar doble conteo.")

print("\n==== VERIFICACIÓN COMPLETADA ====")
PYEOF

apptainer exec --bind /mnt/cephfs:/mnt/cephfs --bind "$TMPDIR:/tmp,$TMPDIR:/var/tmp" "$SIF" \
    /opt/miniforge3/envs/vcf2zarr/bin/python "$PYTHON_SCRIPT"

rm -rf "$TMPDIR"
echo -e "\n[$(date +'%H:%M:%S')] ==== LD CHECK COMPLETADO ===="
