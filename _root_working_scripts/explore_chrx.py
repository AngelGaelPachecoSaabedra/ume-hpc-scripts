import zarr
import numpy as np
import pandas as pd

zarr_x_path = "/mnt/cephfs/hot_nvme/mcps/imputed-topmed/zar_files/chrX.zarr"
covar_file = "/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes/covar_ctrl_45.txt"

print(f"Abriendo Zarr: {zarr_x_path}")
root = zarr.open_group(zarr_x_path, mode='r')

print("\n--- Metadatos del Zarr ---")
print("Arrays disponibles:", list(root.keys()))
if 'contig_id' in root:
    print("Contigs presentes en este archivo:", list(root['contig_id'][:]))
else:
    print("No hay array de 'contig_id'.")

# Cargar covariables para ubicar hombres y mujeres
print("\nBuscando hombres y mujeres en la covariable...")
df_covar = pd.read_csv(covar_file, sep='\t')
sex_map = dict(zip(df_covar['IID'].astype(str), df_covar['SEX']))

sample_ids = root['sample_id'][:]
sample_ids = [s if isinstance(s, str) else s.decode() for s in sample_ids]

idx_males = []
idx_females = []
for i, sid in enumerate(sample_ids):
    s_sex = sex_map.get(sid, -1)
    if s_sex == 1 and len(idx_males) < 5:
        idx_males.append(i)
    elif s_sex == 0 and len(idx_females) < 5:
        idx_females.append(i)
    if len(idx_males) == 5 and len(idx_females) == 5:
        break

# Inspeccionar los genotipos de una variante cualquiera (índice 1000)
pos = root['variant_position'][1000]
ref = root['variant_allele'][1000, 0]
alt = root['variant_allele'][1000, 1]
if isinstance(ref, bytes): ref = ref.decode()
if isinstance(alt, bytes): alt = alt.decode()

print(f"\n--- Inspeccionando Variante en Posición {pos} ({ref}/{alt}) ---")
gt_raw = root['call_genotype'][1000]

print("\nResultados MUJERES:")
for idx in idx_females:
    raw = gt_raw[idx]
    suma = np.sum(np.maximum(raw, 0))
    print(f"  ID {sample_ids[idx]}: Raw={raw}, Suma dosificación={suma}")

print("\nResultados HOMBRES:")
for idx in idx_males:
    raw = gt_raw[idx]
    suma = np.sum(np.maximum(raw, 0))
    print(f"  ID {sample_ids[idx]}: Raw={raw}, Suma dosificación={suma}")

