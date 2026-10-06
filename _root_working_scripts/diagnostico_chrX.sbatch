#!/bin/bash
#SBATCH --job-name=diag2_chrX
#SBATCH --account=researchers
#SBATCH --qos=vip
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:15:00
#SBATCH --output=/mnt/cephfs/orgs/home/angel.pacheco/logs/diag2_chrX_%j.out
#SBATCH --error=/mnt/cephfs/orgs/home/angel.pacheco/logs/diag2_chrX_%j.err

set -euo pipefail

SIF="/mnt/cephfs/biocontainers/images/jupyter-biotools-1.2.sif"
TMPDIR="/mnt/cephfs/scratch/$USER/diag2_tmp/$SLURM_JOB_ID"
mkdir -p "$TMPDIR"
chmod 700 "$TMPDIR"
export TMPDIR TMP="$TMPDIR" TEMP="$TMPDIR"

CHRX_BED="/mnt/cephfs/hot_nvme/mcps/imputed-topmed/plink_files/maximally_unrelated/mcps-freeze150k_qcd_chrX_ivs"

echo "============================================================"
echo "  DIAGNÓSTICO 2: chrX codificación + software"
echo "  $(date)"
echo "============================================================"

# ================================================================
# 1. Encontrar PLINK2 en r-bio
# ================================================================
echo ""
echo "--- 1. Buscando PLINK2 en entornos conda ---"
apptainer exec --bind /mnt/cephfs:/mnt/cephfs "$SIF" \
    bash -c '
    for env in /opt/miniforge3/envs/*/bin/plink2; do
        if [ -x "$env" ]; then
            echo "✅ PLINK2 encontrado: $env"
            "$env" --version 2>&1 | head -1
        fi
    done
    # También buscar en PATH general
    find /opt/miniforge3 -name "plink2" -type f 2>/dev/null
    '

# ================================================================
# 2. Usar PLINK2 de r-bio para geno-counts en chrX
# ================================================================
echo ""
echo "--- 2. Conteo de genotipos chrX (geno-counts) ---"
apptainer exec --bind /mnt/cephfs:/mnt/cephfs --bind "$TMPDIR:/tmp,$TMPDIR:/var/tmp" "$SIF" \
    bash -c "
    PLINK2=/opt/miniforge3/envs/r-bio/bin/plink2
    if [ -x \"\$PLINK2\" ]; then
        \$PLINK2 --bfile ${CHRX_BED} \\
            --chr 23 \\
            --thin-count 500 \\
            --geno-counts \\
            --out $TMPDIR/chrX_geno \\
            --threads 4 \\
            --memory 12000

        echo ''
        echo 'Primeras 15 líneas de geno-counts:'
        head -16 $TMPDIR/chrX_geno.gcount

        echo ''
        echo 'Columnas con HAP (hemicigoto):'
        head -1 $TMPDIR/chrX_geno.gcount

        # Contar si hay columnas HAP_REF_CT / HAP_ALT_CT
        if head -1 $TMPDIR/chrX_geno.gcount | grep -q 'HAP'; then
            echo '✅ PLINK2 reconoce hemicigotos — columnas HAP presentes'
            echo 'Ejemplo de conteos HAP (primeras 10 variantes):'
            awk 'NR<=11' $TMPDIR/chrX_geno.gcount | column -t
        else
            echo '⚠️ Sin columnas HAP — todo codificado como diploide'
            awk 'NR<=11' $TMPDIR/chrX_geno.gcount | column -t
        fi
    else
        echo '❌ PLINK2 no encontrado en r-bio'
    fi
    "

# ================================================================
# 3. Verificar .nosex y .fam más detallado
# ================================================================
echo ""
echo "--- 3. Contenido del .nosex ---"
wc -l "${CHRX_BED}.nosex"
echo "Primeras 5 líneas:"
head -5 "${CHRX_BED}.nosex"
echo ""
echo "Primeras 5 líneas del .fam:"
head -5 "${CHRX_BED}.fam"
echo ""
echo "Sexo en .fam — col 5 (1=male, 2=female, 0=unknown):"
awk '{print $5}' "${CHRX_BED}.fam" | sort | uniq -c | sort -rn

# ================================================================
# 4. Log original de la creación del bfile chrX
# ================================================================
echo ""
echo "--- 4. Log de creación de chrX bfiles ---"
cat "${CHRX_BED}.log" 2>/dev/null || echo "Log no disponible"

# ================================================================
# 5. Verificar Zarr chrX — muestreo de hombres vs mujeres
# ================================================================
echo ""
echo "--- 5. Zarr chrX: call_DS en hombres vs mujeres ---"
cat << 'PYEOF' > "$TMPDIR/check_ds_sex.py"
import zarr
import numpy as np
import pandas as pd
import time

print(f"[{time.strftime('%H:%M:%S')}] Analizando call_DS de chrX por sexo...")

# Abrir Zarr chrX
root = zarr.open_group("/mnt/cephfs/hot_nvme/mcps/imputed-topmed/zar_files/chrX.zarr", mode='r')
sample_ids = [s if isinstance(s, str) else s.decode() for s in root['sample_id'][:]]

# Leer sexo del .fam de autosomas (chr1) para tener sexo de las 140K
# El .fam del panel 79K tiene sexo, pero necesitamos el de todas las 140K
# Usar el chrX.zarr y los genotipos para inferir sexo:
# En chrX, hombres hemicigotos deberían tener call_genotype con solo 0/0 o 1/1 (no 0/1)

# Tomar 1000 variantes aleatorias y 1000 muestras
np.random.seed(42)
n_vars = min(500, root['call_DS'].shape[0])
n_samp = min(2000, root['call_DS'].shape[1])
var_idx = np.sort(np.random.choice(root['call_DS'].shape[0], n_vars, replace=False))
samp_idx = np.arange(n_samp)

# Leer genotipos discretos para inferir sexo
print(f"[{time.strftime('%H:%M:%S')}] Leyendo genotipos de {n_vars} variantes x {n_samp} muestras...")
gt = root['call_genotype'][var_idx[:, None], samp_idx[None, :], :]  # (vars, samps, 2)
gt_sum = np.sum(np.maximum(gt, 0), axis=2)  # (vars, samps)

# Contar heterocigotos por muestra
het_count = np.sum(gt_sum == 1, axis=0)  # por muestra
het_rate = het_count / n_vars

# Hombres: het_rate ≈ 0 (hemicigotos, no deberían tener het)
# Mujeres: het_rate > 0 (diploides, sí tienen het)
print(f"\n  Tasa de heterocigosidad por muestra (primeras {n_samp}):")
print(f"    min={het_rate.min():.4f}, max={het_rate.max():.4f}, mean={het_rate.mean():.4f}")
print(f"    mediana={np.median(het_rate):.4f}")

# Clasificar
threshold = 0.02  # < 2% het = probablemente hombre
inferred_male = het_rate < threshold
inferred_female = ~inferred_male
print(f"\n  Inferencia de sexo (umbral het < {threshold}):")
print(f"    Probables hombres: {inferred_male.sum()}")
print(f"    Probables mujeres: {inferred_female.sum()}")
print(f"    Ratio F/M: {inferred_female.sum() / max(inferred_male.sum(), 1):.2f}")

# Ahora verificar call_DS para hombres vs mujeres
print(f"\n[{time.strftime('%H:%M:%S')}] Comparando call_DS por sexo inferido...")
ds = root['call_DS'][var_idx[:, None], samp_idx[None, :], 0]  # (vars, samps)

ds_male = ds[:, inferred_male]
ds_female = ds[:, inferred_female]

print(f"\n  call_DS HOMBRES (n={inferred_male.sum()}):")
print(f"    min={ds_male.min():.4f}, max={ds_male.max():.4f}, mean={ds_male.mean():.4f}")
vals_m = np.round(ds_male.flatten(), 1)
unique_m, counts_m = np.unique(vals_m, return_counts=True)
print(f"    Top 10 valores (redondeado):")
for v, c in sorted(zip(unique_m, counts_m), key=lambda x: -x[1])[:10]:
    print(f"      {v:.1f}: {c:,} ({100*c/len(vals_m):.1f}%)")

# ¿Cuántos valores son intermedios (0.1-0.9 o 1.1-1.9)?
intermediate_m = ((vals_m > 0.05) & (vals_m < 0.95)) | ((vals_m > 1.05) & (vals_m < 1.95))
print(f"    Valores intermedios (no 0/1/2): {intermediate_m.sum():,} ({100*intermediate_m.sum()/len(vals_m):.1f}%)")

print(f"\n  call_DS MUJERES (n={inferred_female.sum()}):")
print(f"    min={ds_female.min():.4f}, max={ds_female.max():.4f}, mean={ds_female.mean():.4f}")
vals_f = np.round(ds_female.flatten(), 1)
unique_f, counts_f = np.unique(vals_f, return_counts=True)
print(f"    Top 10 valores (redondeado):")
for v, c in sorted(zip(unique_f, counts_f), key=lambda x: -x[1])[:10]:
    print(f"      {v:.1f}: {c:,} ({100*c/len(vals_f):.1f}%)")

intermediate_f = ((vals_f > 0.05) & (vals_f < 0.95)) | ((vals_f > 1.05) & (vals_f < 1.95))
print(f"    Valores intermedios (no 0/1/2): {intermediate_f.sum():,} ({100*intermediate_f.sum()/len(vals_f):.1f}%)")

print(f"\n[{time.strftime('%H:%M:%S')}] Diagnóstico completado.")
PYEOF

apptainer exec \
    --bind /mnt/cephfs:/mnt/cephfs \
    --bind "$TMPDIR:/tmp,$TMPDIR:/var/tmp" \
    "$SIF" \
    /opt/miniforge3/envs/vcf2zarr/bin/python "$TMPDIR/check_ds_sex.py"

# ================================================================
# 6. Software: R packages y GCTB
# ================================================================
echo ""
echo "--- 6. R packages disponibles en r-bio ---"
apptainer exec --bind /mnt/cephfs:/mnt/cephfs "$SIF" \
    /opt/miniforge3/envs/r-bio/bin/Rscript -e '
pkgs <- c("bigsnpr", "bigstatsr", "bigsparser", "data.table", "R.utils")
for (p in pkgs) {
    if (requireNamespace(p, quietly=TRUE)) {
        cat(sprintf("✅ %s: %s\n", p, as.character(packageVersion(p))))
    } else {
        cat(sprintf("❌ %s: NO instalado\n", p))
    }
}
cat(sprintf("\nR version: %s\n", R.version.string))
' 2>/dev/null

echo ""
echo "--- 7. Buscar GCTB en todo el sistema ---"
find /mnt/cephfs -name "gctb" -type f 2>/dev/null | head -5
find /opt -name "gctb" -type f 2>/dev/null | head -5
which gctb 2>/dev/null || echo "GCTB no en PATH"

echo ""
echo "--- 8. Contenedor 1.3 en Podman ---"
# Verificar si podman está accesible
which podman 2>/dev/null && podman ps --format "{{.ID}} {{.Image}} {{.Status}}" 2>/dev/null | head -5 || echo "Podman no accesible desde nodo de cómputo"

rm -rf "$TMPDIR"
echo ""
echo "============================================================"
echo "  DIAGNÓSTICO 2 COMPLETADO: $(date)"
echo "============================================================"
