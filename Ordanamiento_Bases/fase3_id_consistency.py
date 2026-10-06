import pandas as pd

BASE_DIR   = "/mnt/cephfs/orgs/home/angel.pacheco/Ordanamiento_Bases/"
BASAL      = BASE_DIR + "Bases_DNAnexus/MCPS BASELINE.csv"
MORT       = BASE_DIR + "Bases_DNAnexus/MCPS MORTALITY (DNAnexus version 2.1).csv"
DENOM_PATH = BASE_DIR + "denominador_genomico_F145K.csv"

df_b = pd.read_csv(BASAL,  usecols=["PATID"])
df_m = pd.read_csv(MORT,   usecols=["PATID","STATUS","D040"])
df_d = pd.read_csv(DENOM_PATH, usecols=["PATID"], dtype={"PATID": str})

set_b = set(df_b["PATID"])
set_m = set(df_m["PATID"])
set_d = set(df_d["PATID"])

print("=" * 65)
print("VALIDACIÓN CONSISTENCIA DE ID")
print("=" * 65)
print(f"\nBASELINE  → N PATID únicos:   {len(set_b):>10,}")
print(f"MORTALITY → N PATID únicos:   {len(set_m):>10,}")
print(f"DENOMINADOR → N PATID únicos: {len(set_d):>10,}")

print(f"\nBaseline ∩ Mortality:         {len(set_b & set_m):>10,}  (esperado: 159517)")
print(f"Baseline ∩ Denominador:       {len(set_b & set_d):>10,}")
print(f"Mortality ∩ Denominador:      {len(set_m & set_d):>10,}")

print(f"\nEn denominador, NO en baseline: {len(set_d - set_b):>7,}")
print(f"En denominador, NO en mortality: {len(set_d - set_m):>6,}")

print()
print("=== Prefijos PATID en denominador ===")
mcpsds_in_denom  = df_d["PATID"].str.startswith("MCPSDS").sum()
mcpsrgn_in_denom = df_d["PATID"].str.startswith("MCPSRGN").sum()
print(f"  MCPSDS  en denominador: {mcpsds_in_denom:>8,}")
print(f"  MCPSRGN en denominador: {mcpsrgn_in_denom:>8,}")

print()
print("=== PATIDs MCPSDS en baseline que están en denominador ===")
ds_ids_baseline = set(df_b[df_b["PATID"].str.startswith("MCPSDS")]["PATID"])
print(f"  MCPSDS en baseline:     {len(ds_ids_baseline):>8,}")
overlap_ds_denom = ds_ids_baseline & set_d
print(f"  MCPSDS en denominador:  {len(overlap_ds_denom):>8,}  ← si 0: MCPSDS sin genómica")

print()
print("=== Pérdida por inner join Baseline ∩ Denominador ===")
print(f"  Baseline total:         {len(set_b):>8,}")
print(f"  Post inner join:        {len(set_b & set_d):>8,}")
print(f"  Pérdida:                {len(set_b) - len(set_b & set_d):>8,}  ({(len(set_b)-len(set_b&set_d))/len(set_b)*100:.1f}%)")

print()
print("=== D040 NaN asociado a STATUS='U' (confirmar) ===")
df_m_full = pd.read_csv(MORT, usecols=["PATID","STATUS","D040"])
nan_not_u = df_m_full[(df_m_full["D040"].isna()) & (df_m_full["STATUS"] != "U")]
print(f"  NaN en D040 con STATUS ≠ 'U': {len(nan_not_u)}  (esperado: 0)")
print(f"  NaN en D040 con STATUS == 'U': {df_m_full[(df_m_full['D040'].isna()) & (df_m_full['STATUS']=='U')].shape[0]}")
