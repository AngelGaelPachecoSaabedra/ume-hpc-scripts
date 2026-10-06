import pandas as pd

BASAL = "/mnt/cephfs/orgs/home/angel.pacheco/Ordanamiento_Bases/Bases_DNAnexus/MCPS BASELINE.csv"

print("=" * 65)
print("FASE 3 — VALIDACIÓN ARCHIVO BASAL")
print("=" * 65)
print(f"Archivo: {BASAL.split('/')[-1]}\n")

# Cargar sólo cabecera primero
df_cols = pd.read_csv(BASAL, nrows=0)
all_cols = df_cols.columns.tolist()
print(f"Total columnas: {len(all_cols)}")
print(f"Todas las columnas:\n{all_cols}\n")

# ── Variables oncológicas obligatorias (FASE 1) ──────────────────
CANCER_VARS = [
    "BASE_LUNGCANCER",
    "BASE_OTHCANCER",
    "BASE_PROSTATECANCER",
    "BASE_CERVCANCER",
    "BASE_BREASTCANCER",
    "BASE_STOMCANCER",
    "BASE_ORALCANCER"
]
print("=== Verificación variables oncológicas ===")
for v in CANCER_VARS:
    status = "PRESENTE ✓" if v in all_cols else "AUSENTE ✗"
    print(f"  {v:30s} → {status}")

# ── Identificar covariables ──────────────────────────────────────
print("\n=== Búsqueda de covariables ===")
# Edad
edad_cols = [c for c in all_cols if "age" in c.lower() or "edad" in c.lower()]
print(f"  EDAD candidatas: {edad_cols}")

# Sexo
sexo_cols = [c for c in all_cols if "sex" in c.lower() or "male" in c.lower() or "genero" in c.lower()]
print(f"  SEXO candidatas: {sexo_cols}")

# Diabetes
dm_cols = [c for c in all_cols if "diab" in c.lower()]
print(f"  DIABETES candidatas: {dm_cols}")

# Tabaquismo
smok_cols = [c for c in all_cols if "smok" in c.lower() or "tabac" in c.lower() or "cigar" in c.lower()]
print(f"  TABAQUISMO candidatas: {smok_cols}")
