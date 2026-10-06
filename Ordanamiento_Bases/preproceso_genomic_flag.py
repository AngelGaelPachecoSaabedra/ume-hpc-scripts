import pandas as pd

BASE_DIR = "/mnt/cephfs/orgs/home/angel.pacheco/Ordanamiento_Bases/"
PATH = BASE_DIR + "denominador_genomico_F145K.csv"

df = pd.read_csv(PATH)

print("Columnas antes:", df.columns.tolist())
print("Shape antes:", df.shape)
print()
print("EXOME_MCPS_FLAG — distribución:")
print(df["EXOME_MCPS_FLAG"].value_counts(dropna=False).to_string())
print()

# ── Regla GENOMIC_QC_FLAG ─────────────────────────────────────────
# linkage_risk = 1 si FLAG == '6' OR contiene '6' en su valor string
# NaN (sin flag) → 0, flags sin '6' → 0
def flag_linkage_risk(val):
    if pd.isna(val):
        return 0
    return 1 if "6" in str(val) else 0

df["GENOMIC_QC_FLAG"] = df["EXOME_MCPS_FLAG"].apply(flag_linkage_risk)

print("GENOMIC_QC_FLAG — distribución:")
print(df["GENOMIC_QC_FLAG"].value_counts(dropna=False).to_string())
print()

# Verificar: filas con flag '6', '5|6', '6|7' → deben ser 1
check = df[df["EXOME_MCPS_FLAG"].isin(["6","5|6","6|7"])]["GENOMIC_QC_FLAG"].unique()
print(f"Valores GENOMIC_QC_FLAG para flags con '6': {check}  ← esperado: [1]")
check2 = df[df["EXOME_MCPS_FLAG"].isin(["4","5","7"])]["GENOMIC_QC_FLAG"].unique()
print(f"Valores GENOMIC_QC_FLAG para flags sin '6': {check2}  ← esperado: [0]")
print()

# Guardar
col_order = ["PATID","IID_EXOME","IID_WGS","HAS_EXOME","HAS_WGS","EXOME_MCPS_FLAG","GENOMIC_QC_FLAG"]
df = df[col_order]
df.to_csv(PATH, index=False)
print("Archivo actualizado guardado:", PATH)
print()
print("head(6):")
print(df.head(6).to_string())
print()
print("Columnas finales:", df.columns.tolist())
print("Shape final:", df.shape)
