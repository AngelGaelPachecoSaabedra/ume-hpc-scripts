"""
FASE 4: CONSTRUCCIÓN DEL FENOTIPO ONCOLÓGICO
MCPS — Estudio Prospectivo Ciudad de México
Script reproducible — pandas
"""

import pandas as pd
import numpy as np

# ─── RUTAS ────────────────────────────────────────────────────────────────────
BASE_DIR  = "/mnt/cephfs/orgs/home/angel.pacheco/Ordanamiento_Bases/"
PATH_DENOM = BASE_DIR + "denominador_genomico_F145K.csv"
PATH_BASAL = BASE_DIR + "Bases_DNAnexus/MCPS BASELINE.csv"
PATH_MORT  = BASE_DIR + "Bases_DNAnexus/MCPS MORTALITY (DNAnexus version 2.1).csv"

CANCER_VARS = [
    "BASE_LUNGCANCER", "BASE_OTHCANCER", "BASE_PROSTATECANCER",
    "BASE_CERVCANCER", "BASE_BREASTCANCER", "BASE_STOMCANCER", "BASE_ORALCANCER"
]

BASAL_COLS = ["PATID", "AGE", "MALE", "BASE_DIABETES", "EVER_SMOK"] + CANCER_VARS

MORT_COLS  = ["PATID", "STATUS", "D040", "ICD10_UNDERLYING",
              "CAUSA1", "CAUSA2", "CAUSA3", "CAUSA4", "CAUSA5", "CAUSA6"]

FINAL_COLS = [
    "PATID", "IID_EXOME", "IID_WGS", "GENOMIC_QC_FLAG",
    "AGE", "MALE", "BASE_DIABETES", "EVER_SMOK",
    "EVIDENCIA_CASO_BASAL", "EVIDENCIA_INCIDENCIA_MORTAL",
    "CANCER_SOURCE", "PHENOTYPE_CANCER_STATUS"
]

SEP = "=" * 68

# ══════════════════════════════════════════════════════════════════════════════
# PASO 1 — CARGA DE DATOS
# ══════════════════════════════════════════════════════════════════════════════
print(SEP)
print("PASO 1 — CARGA DE DATOS")
print(SEP)

df_denom = pd.read_csv(PATH_DENOM, dtype={"PATID": str, "EXOME_MCPS_FLAG": str})
df_basal = pd.read_csv(PATH_BASAL, usecols=BASAL_COLS, dtype={"PATID": str})
df_mort  = pd.read_csv(PATH_MORT,  usecols=MORT_COLS,  dtype={"PATID": str})

print(f"  Denominador genómico:  {len(df_denom):>8,} filas x {df_denom.shape[1]} cols")
print(f"  Basal:                 {len(df_basal):>8,} filas x {df_basal.shape[1]} cols")
print(f"  Mortalidad:            {len(df_mort):>8,} filas x {df_mort.shape[1]} cols")

# ══════════════════════════════════════════════════════════════════════════════
# PASO 2 — INTEGRACIÓN (INNER JOINS) + EXCLUSIÓN STATUS='U'
# ══════════════════════════════════════════════════════════════════════════════
print(f"\n{SEP}")
print("PASO 2 — INTEGRACIÓN Y FILTROS")
print(SEP)

# Inner join 1: denominador ∩ basal
df = df_denom.merge(df_basal, on="PATID", how="inner")
n_post_join1 = len(df)
print(f"  Post inner join (denominador ∩ basal):       {n_post_join1:>8,}")

# Inner join 2: resultado ∩ mortalidad
df = df.merge(df_mort, on="PATID", how="inner")
n_post_join2 = len(df)
print(f"  Post inner join (∩ mortalidad):              {n_post_join2:>8,}")

# Verificar integridad: no debe haber pérdida entre join1 y join2
assert n_post_join1 == n_post_join2, \
    f"ERROR: pérdida inesperada en join2 ({n_post_join1} → {n_post_join2})"
print(f"  Sin pérdida entre joins ✓")

# Impacto STATUS antes de excluir
n_status_u = (df["STATUS"] == "U").sum()
n_status_a = (df["STATUS"] == "A").sum()
n_status_d = (df["STATUS"] == "D").sum()
print(f"\n  STATUS distribución (pre-exclusión):")
print(f"    A (vivo):      {n_status_a:>8,}")
print(f"    D (fallecido): {n_status_d:>8,}")
print(f"    U (incierto):  {n_status_u:>8,}  ← a excluir")

# Excluir STATUS='U'
df = df[df["STATUS"] != "U"].copy()
n_post_excl_u = len(df)
print(f"\n  Post exclusión STATUS='U':                   {n_post_excl_u:>8,}")
print(f"  Participantes excluidos por STATUS='U':      {n_post_join2 - n_post_excl_u:>8,}")

# ══════════════════════════════════════════════════════════════════════════════
# PASO 3 — ESTANDARIZACIÓN PREVIA
# ══════════════════════════════════════════════════════════════════════════════
print(f"\n{SEP}")
print("PASO 3 — ESTANDARIZACIÓN")
print(SEP)

# D040 → binaria entera (1.0 → 1, cualquier otro → 0, NaN no debería existir post-exclusión U)
nan_d040_post_excl = df["D040"].isna().sum()
print(f"  NaN en D040 post-exclusión STATUS='U': {nan_d040_post_excl}  (esperado: 0)")
df["D040"] = (df["D040"] == 1).astype(int)
print(f"  D040 convertido a int (0/1). Distribución: {df['D040'].value_counts().to_dict()}")

# EVER_SMOK — reporte de NA, sin imputar
nan_ever_smok = df["EVER_SMOK"].isna().sum()
print(f"\n  NaN en EVER_SMOK: {nan_ever_smok}  ({nan_ever_smok/len(df)*100:.3f}% del dataset post-exclusión)")

# GENOMIC_QC_FLAG — confirmar presencia
print(f"\n  GENOMIC_QC_FLAG presente: {'GENOMIC_QC_FLAG' in df.columns}")
print(f"  GENOMIC_QC_FLAG distribución: {df['GENOMIC_QC_FLAG'].value_counts().to_dict()}")

# ══════════════════════════════════════════════════════════════════════════════
# PASO 4 — DERIVACIÓN DE VARIABLES CATEGÓRICAS
# ══════════════════════════════════════════════════════════════════════════════
print(f"\n{SEP}")
print("PASO 4 — DERIVACIÓN DEL FENOTIPO")
print(SEP)

# ── EVIDENCIA_CASO_BASAL ─────────────────────────────────────────────────────
# 1 si cualquier variable basal == 1; 0 si todas son 0 o NaN
df["EVIDENCIA_CASO_BASAL"] = df[CANCER_VARS].max(axis=1).clip(upper=1).fillna(0).astype(int)

n_basal_casos = df["EVIDENCIA_CASO_BASAL"].sum()
print(f"  EVIDENCIA_CASO_BASAL=1 (autorreporte positivo): {n_basal_casos:>6,}")
print(f"  Distribución: {df['EVIDENCIA_CASO_BASAL'].value_counts().to_dict()}")

# Desglose por tipo de cáncer basal
print("\n  Desglose por variable basal (pueden solaparse):")
for v in CANCER_VARS:
    n = int(df[v].sum())
    print(f"    {v:30s}: {n:>5,}")

# ── EVIDENCIA_INCIDENCIA_MORTAL ───────────────────────────────────────────────
df["EVIDENCIA_INCIDENCIA_MORTAL"] = df["D040"].astype(int)

n_mort_casos = df["EVIDENCIA_INCIDENCIA_MORTAL"].sum()
print(f"\n  EVIDENCIA_INCIDENCIA_MORTAL=1 (D040, mortalidad oncológica): {n_mort_casos:>5,}")
print(f"  Distribución: {df['EVIDENCIA_INCIDENCIA_MORTAL'].value_counts().to_dict()}")

# ── PHENOTYPE_CANCER_STATUS ───────────────────────────────────────────────────
caso_mask    = (df["EVIDENCIA_CASO_BASAL"] == 1) | (df["EVIDENCIA_INCIDENCIA_MORTAL"] == 1)
control_mask = (df["EVIDENCIA_CASO_BASAL"] == 0) & (df["EVIDENCIA_INCIDENCIA_MORTAL"] == 0)

df["PHENOTYPE_CANCER_STATUS"] = "INDEFINIDO"  # Fallback de seguridad
df.loc[control_mask, "PHENOTYPE_CANCER_STATUS"] = "CONTROL"
df.loc[caso_mask,    "PHENOTYPE_CANCER_STATUS"] = "CASO"

# Verificar ausencia de INDEFINIDO
n_indefinido = (df["PHENOTYPE_CANCER_STATUS"] == "INDEFINIDO").sum()
assert n_indefinido == 0, f"ERROR: {n_indefinido} registros con fenotipo INDEFINIDO"
print(f"\n  PHENOTYPE_CANCER_STATUS — sin valores INDEFINIDOS ✓")

# ── CANCER_SOURCE ─────────────────────────────────────────────────────────────
conditions = [
    (df["EVIDENCIA_CASO_BASAL"] == 1) & (df["EVIDENCIA_INCIDENCIA_MORTAL"] == 0),
    (df["EVIDENCIA_CASO_BASAL"] == 0) & (df["EVIDENCIA_INCIDENCIA_MORTAL"] == 1),
    (df["EVIDENCIA_CASO_BASAL"] == 1) & (df["EVIDENCIA_INCIDENCIA_MORTAL"] == 1),
    (df["EVIDENCIA_CASO_BASAL"] == 0) & (df["EVIDENCIA_INCIDENCIA_MORTAL"] == 0),
]
choices = ["SOLO_BASAL", "SOLO_MORTALIDAD", "AMBOS", "NINGUNO"]

df["CANCER_SOURCE"] = np.select(conditions, choices, default="ERROR")

n_error_source = (df["CANCER_SOURCE"] == "ERROR").sum()
assert n_error_source == 0, f"ERROR: {n_error_source} registros sin CANCER_SOURCE asignado"
print(f"  CANCER_SOURCE — sin valores ERROR ✓")

# ══════════════════════════════════════════════════════════════════════════════
# PASO 5 — SELECCIÓN DE VARIABLES FINALES
# ══════════════════════════════════════════════════════════════════════════════
print(f"\n{SEP}")
print("PASO 5 — SELECCIÓN DE VARIABLES FINALES")
print(SEP)

# Verificar que todas las columnas finales existen
missing_cols = [c for c in FINAL_COLS if c not in df.columns]
if missing_cols:
    print(f"  ALERTA — columnas faltantes: {missing_cols}")
else:
    print(f"  Todas las columnas finales presentes ✓")

df_final = df[FINAL_COLS].copy()
print(f"  Shape final: {df_final.shape[0]:,} filas x {df_final.shape[1]} columnas")
print(f"  Columnas: {df_final.columns.tolist()}")

# ══════════════════════════════════════════════════════════════════════════════
# PASO 6 — CONTROL DE CALIDAD (REPORTE COMPLETO)
# ══════════════════════════════════════════════════════════════════════════════
print(f"\n{SEP}")
print("PASO 6 — REPORTE DE CONTROL DE CALIDAD")
print(SEP)

print(f"\n{'─'*50}")
print("6.1  FLUJO DE REGISTROS")
print(f"{'─'*50}")
print(f"  Denominador genómico (input):       {len(df_denom):>8,}")
print(f"  Post inner join (∩ basal):          {n_post_join1:>8,}  (Δ={len(df_denom)-n_post_join1:,})")
print(f"  Post inner join (∩ mortalidad):     {n_post_join2:>8,}  (Δ=0)")
print(f"  Excluidos STATUS='U':               {n_status_u:>8,}")
print(f"  Dataset analítico final:            {n_post_excl_u:>8,}")

print(f"\n{'─'*50}")
print("6.2  MISSING DATA")
print(f"{'─'*50}")
for col in ["BASE_DIABETES", "EVER_SMOK"]:
    n_na = df_final[col].isna().sum()
    pct  = n_na / len(df_final) * 100
    print(f"  {col:20s}  NaN={n_na:>5,}  ({pct:.3f}%)")

print(f"\n{'─'*50}")
print("6.3  PHENOTYPE_CANCER_STATUS")
print(f"{'─'*50}")
pheno_dist = df_final["PHENOTYPE_CANCER_STATUS"].value_counts()
total_final = len(df_final)
for label, n in pheno_dist.items():
    pct = n / total_final * 100
    print(f"  {label:10s}  n={n:>7,}  ({pct:.2f}%)")
print(f"  {'TOTAL':10s}  n={total_final:>7,}")

print(f"\n{'─'*50}")
print("6.4  DESGLOSE DE CASOS POR FUENTE (CANCER_SOURCE)")
print(f"{'─'*50}")
source_dist = df_final["CANCER_SOURCE"].value_counts()
n_casos_total = (df_final["PHENOTYPE_CANCER_STATUS"] == "CASO").sum()
for label, n in source_dist.items():
    pct_all   = n / total_final * 100
    pct_casos = n / n_casos_total * 100 if n_casos_total > 0 and label != "NINGUNO" else 0
    if label == "NINGUNO":
        print(f"  {label:20s}  n={n:>7,}  ({pct_all:.2f}% del total)")
    else:
        print(f"  {label:20s}  n={n:>7,}  ({pct_all:.2f}% del total | {pct_casos:.1f}% de casos)")

# Verificar coherencia: NINGUNO == CONTROL
n_ninguno  = (df_final["CANCER_SOURCE"] == "NINGUNO").sum()
n_controles = (df_final["PHENOTYPE_CANCER_STATUS"] == "CONTROL").sum()
match_ctrl = n_ninguno == n_controles
print(f"\n  NINGUNO ({n_ninguno:,}) == CONTROL ({n_controles:,}): {match_ctrl} ✓" if match_ctrl
      else f"\n  ALERTA: NINGUNO ({n_ninguno}) ≠ CONTROL ({n_controles})")

print(f"\n{'─'*50}")
print("6.5  IMPACTO STATUS='U' SOBRE CASOS")
print(f"{'─'*50}")
# Recuperar registros U del join para verificar cuántos serían casos
df_u_check = df_denom.merge(df_basal, on="PATID", how="inner") \
                      .merge(df_mort,  on="PATID", how="inner")
df_u_check = df_u_check[df_u_check["STATUS"] == "U"].copy()
df_u_check["_basal_flag"] = df_u_check[CANCER_VARS].max(axis=1).clip(upper=1).fillna(0).astype(int)
# D040 es NaN para STATUS=U por definición MCPS
n_u_basal_casos = df_u_check["_basal_flag"].sum()
print(f"  Participantes STATUS='U' excluidos:       {n_status_u:>6,}")
print(f"  De esos, con EVIDENCIA_CASO_BASAL=1:      {n_u_basal_casos:>6,}  ← no clasificables como caso/control")
print(f"  (D040 es NaN para todos los STATUS='U' — confirmado en FASE 3)")

print(f"\n{'─'*50}")
print("6.6  CONTEO CÓDIGOS ICD-10 INESPECÍFICOS (D038)")
print(f"{'─'*50}")
# Recuperar ICD10_UNDERLYING del dataset pre-selección (df contiene CAUSA1-6 + ICD10)
# Usar df_mort completo filtrado a no-U con D040=1
df_mort_casos = df_mort[
    (df_mort["STATUS"] != "U") &
    (df_mort["D040"] == 1.0)
].copy()

all_causa_cols = ["ICD10_UNDERLYING","CAUSA1","CAUSA2","CAUSA3","CAUSA4","CAUSA5","CAUSA6"]
for prefix, label in [("C76","C76x (ill-defined sites)"),
                       ("C80","C80x (site unspecified)"),
                       ("C97","C97X (multiple primary)")]:
    # Contar en ICD10_UNDERLYING (causa adjudicada, la más relevante)
    n_under = df_mort_casos["ICD10_UNDERLYING"].fillna("").str.startswith(prefix).sum()
    # Contar apariciones en CUALQUIER campo CAUSA
    n_any = 0
    for col in all_causa_cols:
        if col in df_mort_casos.columns:
            n_any += df_mort_casos[col].fillna("").str.startswith(prefix).sum()
    print(f"  {label}:")
    print(f"    En ICD10_UNDERLYING (causa adjudicada): {n_under:>5,}")
    print(f"    En cualquier campo CAUSA (CAUSA1-6):    {n_any:>5,}")

# Verificación cruzada D040 con denominador
print(f"\n{'─'*50}")
print("6.7  VERIFICACIÓN CRUZADA")
print(f"{'─'*50}")
print(f"  Suma EVIDENCIA_CASO_BASAL:          {df_final['EVIDENCIA_CASO_BASAL'].sum():>6,}")
print(f"  Suma EVIDENCIA_INCIDENCIA_MORTAL:   {df_final['EVIDENCIA_INCIDENCIA_MORTAL'].sum():>6,}")
print(f"  Casos (basal OR mortalidad):        {(df_final['PHENOTYPE_CANCER_STATUS']=='CASO').sum():>6,}")
print(f"  Controles:                          {(df_final['PHENOTYPE_CANCER_STATUS']=='CONTROL').sum():>6,}")
# Casos solo basal
n_sb = (df_final["CANCER_SOURCE"] == "SOLO_BASAL").sum()
n_sm = (df_final["CANCER_SOURCE"] == "SOLO_MORTALIDAD").sum()
n_ab = (df_final["CANCER_SOURCE"] == "AMBOS").sum()
assert n_sb + n_sm + n_ab == n_casos_total, "ERROR suma de fuentes ≠ total casos"
print(f"  SOLO_BASAL + SOLO_MORTALIDAD + AMBOS = {n_sb}+{n_sm}+{n_ab} = {n_sb+n_sm+n_ab} ✓")

# GENOMIC_QC_FLAG en casos vs controles
print(f"\n{'─'*50}")
print("6.8  GENOMIC_QC_FLAG POR FENOTIPO")
print(f"{'─'*50}")
flag_cross = df_final.groupby("PHENOTYPE_CANCER_STATUS")["GENOMIC_QC_FLAG"].value_counts().unstack(fill_value=0)
print(flag_cross.to_string())

print(f"\n{SEP}")
print("FASE 4 COMPLETADA — DATASET LISTO PARA FASE 5")
print(SEP)
print(f"\n  df_final shape: {df_final.shape}")
print(f"  Columnas: {df_final.columns.tolist()}")

# Exportar referencia al namespace global para depuración si se importa
# (no exportar CSV aún — eso es FASE 5)
print("\n  [CSV NO exportado — pendiente autorización FASE 5]")
