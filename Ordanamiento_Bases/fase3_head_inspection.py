import pandas as pd
import numpy as np

BASAL = "/mnt/cephfs/orgs/home/angel.pacheco/Ordanamiento_Bases/Bases_DNAnexus/MCPS BASELINE.csv"
MORT  = "/mnt/cephfs/orgs/home/angel.pacheco/Ordanamiento_Bases/Bases_DNAnexus/MCPS MORTALITY (DNAnexus version 2.1).csv"

CANCER_VARS = ["BASE_LUNGCANCER","BASE_OTHCANCER","BASE_PROSTATECANCER",
               "BASE_CERVCANCER","BASE_BREASTCANCER","BASE_STOMCANCER","BASE_ORALCANCER"]
COVAR_VARS  = ["AGE","MALE","BASE_DIABETES","EVER_SMOK","CURR_SMOK"]
MORT_VARS   = ["PATID","STATUS","D040","ICD10_UNDERLYING"]

# ── BASAL ────────────────────────────────────────────────────────
LOAD_COLS = ["PATID"] + CANCER_VARS + COVAR_VARS
df_b = pd.read_csv(BASAL, usecols=LOAD_COLS)
print("=" * 65)
print("BASAL — head(8)")
print("=" * 65)
print(df_b.head(8).to_string())
print()
print("Dtypes:")
print(df_b.dtypes.to_string())
print()
print("Missing values:")
print(df_b.isnull().sum().to_string())
print()
print("Estadísticos descriptivos (variables cáncer):")
print(df_b[CANCER_VARS].describe().to_string())
print()
print("Prevalencias basales (% positivo):")
for v in CANCER_VARS:
    pct = df_b[v].mean()*100
    n   = df_b[v].sum()
    print(f"  {v:30s}: n={int(n):>5,}  ({pct:.2f}%)")
print()
print("Valores únicos clave:")
print(f"  AGE  rango: {df_b['AGE'].min()} - {df_b['AGE'].max()}")
print(f"  MALE valores únicos: {df_b['MALE'].unique()}")
print(f"  BASE_DIABETES valores únicos: {df_b['BASE_DIABETES'].unique()}")
print(f"  EVER_SMOK valores únicos: {df_b['EVER_SMOK'].unique()}")

# ── MORTALIDAD ───────────────────────────────────────────────────
df_m = pd.read_csv(MORT, usecols=MORT_VARS)
print()
print("=" * 65)
print("MORTALIDAD — head(8)")
print("=" * 65)
print(df_m.head(8).to_string())
print()
print("Dtypes:")
print(df_m.dtypes.to_string())
print()
print("Missing values:")
print(df_m.isnull().sum().to_string())
print()
print("STATUS distribución:")
print(df_m["STATUS"].value_counts(dropna=False).to_string())
print()
print("D040 distribución (todos):")
print(df_m["D040"].value_counts(dropna=False).to_string())
print()
print("D040 distribución (STATUS='A' vivos):")
alive = df_m[df_m["STATUS"]=="A"]
print(alive["D040"].value_counts(dropna=False).to_string())
print()
print(f"PATID únicos en MORTALITY: {df_m['PATID'].nunique()}")
print(f"PATID únicos en BASELINE:  {df_b['PATID'].nunique()}")
