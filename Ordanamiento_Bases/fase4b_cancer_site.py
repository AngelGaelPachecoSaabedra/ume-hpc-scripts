"""
FASE 4b: ESTRATIFICACIÓN POR SITIO DE CÁNCER
MCPS — Estudio Prospectivo Ciudad de México
Script reproducible — pandas
Trabaja sobre la base construida en FASE 4 (re-ejecuta pipeline interno)
"""

import pandas as pd
import numpy as np

# ─── RUTAS ────────────────────────────────────────────────────────────────────
BASE_DIR   = "/mnt/cephfs/orgs/home/angel.pacheco/Ordanamiento_Bases/"
PATH_DENOM = BASE_DIR + "denominador_genomico_F145K.csv"
PATH_BASAL = BASE_DIR + "Bases_DNAnexus/MCPS BASELINE.csv"
PATH_MORT  = BASE_DIR + "Bases_DNAnexus/MCPS MORTALITY (DNAnexus version 2.1).csv"

CANCER_VARS = [
    "BASE_LUNGCANCER", "BASE_OTHCANCER", "BASE_PROSTATECANCER",
    "BASE_CERVCANCER", "BASE_BREASTCANCER", "BASE_STOMCANCER", "BASE_ORALCANCER"
]

# Columnas de basal — incluye vars individuales de cáncer para clasificación de sitio
BASAL_COLS = ["PATID", "AGE", "MALE", "BASE_DIABETES", "EVER_SMOK"] + CANCER_VARS

# Mortalidad — añadir ICD10_UNDERLYING para clasificación de sitio
MORT_COLS  = ["PATID", "STATUS", "D040", "ICD10_UNDERLYING"]

SEP  = "=" * 68
LINE = "─" * 50

# ══════════════════════════════════════════════════════════════════════════════
# BLOQUE A — REPRODUCIR PIPELINE FASE 4 (sin modificar PHENOTYPE_CANCER_STATUS)
# ══════════════════════════════════════════════════════════════════════════════
print(SEP)
print("BLOQUE A — REPRODUCCIÓN PIPELINE FASE 4")
print(SEP)

df_denom = pd.read_csv(PATH_DENOM, dtype={"PATID": str, "EXOME_MCPS_FLAG": str})
df_basal = pd.read_csv(PATH_BASAL, usecols=BASAL_COLS, dtype={"PATID": str})
df_mort  = pd.read_csv(PATH_MORT,  usecols=MORT_COLS,  dtype={"PATID": str})

# Inner joins
df = df_denom.merge(df_basal, on="PATID", how="inner")
df = df.merge(df_mort, on="PATID", how="inner")

# Excluir STATUS='U'
df = df[df["STATUS"] != "U"].copy()

# D040 → binaria entera
df["D040"] = (df["D040"] == 1).astype(int)

# Variables fenotípicas (FASE 4)
df["EVIDENCIA_CASO_BASAL"]       = df[CANCER_VARS].max(axis=1).clip(upper=1).fillna(0).astype(int)
df["EVIDENCIA_INCIDENCIA_MORTAL"] = df["D040"].astype(int)

caso_mask    = (df["EVIDENCIA_CASO_BASAL"] == 1) | (df["EVIDENCIA_INCIDENCIA_MORTAL"] == 1)
control_mask = (df["EVIDENCIA_CASO_BASAL"] == 0) & (df["EVIDENCIA_INCIDENCIA_MORTAL"] == 0)

df["PHENOTYPE_CANCER_STATUS"] = np.where(caso_mask, "CASO", "CONTROL")

conditions_src = [
    (df["EVIDENCIA_CASO_BASAL"] == 1) & (df["EVIDENCIA_INCIDENCIA_MORTAL"] == 0),
    (df["EVIDENCIA_CASO_BASAL"] == 0) & (df["EVIDENCIA_INCIDENCIA_MORTAL"] == 1),
    (df["EVIDENCIA_CASO_BASAL"] == 1) & (df["EVIDENCIA_INCIDENCIA_MORTAL"] == 1),
    (df["EVIDENCIA_CASO_BASAL"] == 0) & (df["EVIDENCIA_INCIDENCIA_MORTAL"] == 0),
]
df["CANCER_SOURCE"] = np.select(conditions_src,
                                ["SOLO_BASAL","SOLO_MORTALIDAD","AMBOS","NINGUNO"],
                                default="ERROR")

print(f"  Dataset post-FASE4: {len(df):,} filas")
print(f"  CASOS:    {(df['PHENOTYPE_CANCER_STATUS']=='CASO').sum():,}")
print(f"  CONTROLES:{(df['PHENOTYPE_CANCER_STATUS']=='CONTROL').sum():,}")

# ──────────────────────────────────────────────────────────────────────────────
# Trabajar sobre COPIA del DataFrame — PHENOTYPE_CANCER_STATUS no se modifica
# ──────────────────────────────────────────────────────────────────────────────
df_work = df.copy()

# ══════════════════════════════════════════════════════════════════════════════
# BLOQUE B — CLASIFICACIÓN CANCER_SITE
# ══════════════════════════════════════════════════════════════════════════════
print(f"\n{SEP}")
print("BLOQUE B — DERIVACIÓN DE CANCER_SITE")
print(SEP)

# ─── B.1 Inicializar CANCER_SITE en "NONE" para todos ─────────────────────
df_work["CANCER_SITE"] = "NONE"

# ─── B.2 CLASIFICACIÓN BASAL (prioridad: BREAST > LUNG > PROSTATE >
#         CERVICAL > GI > ORAL > OTHER) ──────────────────────────────────────
# Se aplica de menor a mayor prioridad (el último en escribir "gana")
# Usando np.where en cascada — orden: OTHER primero, BREAST al final
print("  Aplicando clasificación basal (prioridad decreciente últim→primer)...")

basal_rules = [
    ("BASE_OTHCANCER",    "OTHER"),
    ("BASE_ORALCANCER",   "ORAL"),
    ("BASE_STOMCANCER",   "GI"),
    ("BASE_CERVCANCER",   "CERVICAL"),
    ("BASE_PROSTATECANCER","PROSTATE"),
    ("BASE_LUNGCANCER",   "LUNG"),
    ("BASE_BREASTCANCER", "BREAST"),
]

for col, site_label in basal_rules:
    mask = df_work[col] == 1
    df_work.loc[mask, "CANCER_SITE"] = site_label

n_basal_classified = (
    df_work[df_work["PHENOTYPE_CANCER_STATUS"] == "CASO"]["CANCER_SITE"] != "NONE"
).sum()
n_casos_total = (df_work["PHENOTYPE_CANCER_STATUS"] == "CASO").sum()
print(f"  Casos clasificados por basal: {n_basal_classified:,} / {n_casos_total:,}")

# ─── B.3 CLASIFICACIÓN ICD10 para CASOS sin sitio basal ──────────────────────
print("\n  Aplicando clasificación por ICD10_UNDERLYING...")

# Candidatos: CASO + CANCER_SITE == "NONE" + EVIDENCIA_INCIDENCIA_MORTAL == 1
mortal_only_mask = (
    (df_work["PHENOTYPE_CANCER_STATUS"] == "CASO") &
    (df_work["CANCER_SITE"] == "NONE") &
    (df_work["EVIDENCIA_INCIDENCIA_MORTAL"] == 1)
)
print(f"  Candidatos a clasificar por ICD10: {mortal_only_mask.sum():,}")

def classify_icd10(icd_code):
    """Clasificar ICD-10 a sitio oncológico según reglas FASE 4b."""
    if pd.isna(icd_code) or str(icd_code).strip() == "":
        return "UNKNOWN"
    code = str(icd_code).strip().upper()
    # Extraer prefijo numérico para comparaciones de rango
    # ICD-10 formato: letra + 2-3 dígitos (ej. C169, C50X)
    if len(code) < 3:
        return "OTHER"
    letter = code[0]
    if letter != "C":
        # Códigos D3xx (benignos/borderline in situ) — clasificar como OTHER
        return "OTHER"
    try:
        num = int(code[1:3])  # Primeros 2 dígitos numéricos
    except ValueError:
        return "OTHER"

    # ── Reglas de clasificación (orden: más específico primero) ──
    # LUNG: C33–C34
    if num in (33, 34):
        return "LUNG"
    # BREAST: C50
    if num == 50:
        return "BREAST"
    # PROSTATE: C61
    if num == 61:
        return "PROSTATE"
    # CERVICAL: C53
    if num == 53:
        return "CERVICAL"
    # GI: C15–C26 (esófago, estómago, intestino delgado, colon, recto,
    #              ano, hígado, vías biliares, páncreas, otros digestivos)
    if 15 <= num <= 26:
        return "GI"
    # HEMATOLOGIC: C81–C96 (linfomas, leucemias, mieloma)
    if 81 <= num <= 96:
        return "HEMATOLOGIC"
    # UNKNOWN (inespecíficos): C76, C80, C97
    if num in (76, 80, 97):
        return "UNKNOWN"
    # OTHER: todo lo demás C00–C99 no clasificado arriba
    if 0 <= num <= 99:
        return "OTHER"
    return "OTHER"

# Aplicar solo a los candidatos
df_work.loc[mortal_only_mask, "CANCER_SITE"] = \
    df_work.loc[mortal_only_mask, "ICD10_UNDERLYING"].apply(classify_icd10)

n_post_icd = (
    (df_work["PHENOTYPE_CANCER_STATUS"] == "CASO") &
    (df_work["CANCER_SITE"] != "NONE")
).sum()
print(f"  Casos clasificados (basal + ICD10): {n_post_icd:,} / {n_casos_total:,}")

# ──────────────────────────────────────────────────────────────────────────────
# BLOQUE C — VALIDACIONES OBLIGATORIAS
# ──────────────────────────────────────────────────────────────────────────────
print(f"\n{SEP}")
print("BLOQUE C — VALIDACIONES")
print(SEP)

# V1: Todo CASO debe tener CANCER_SITE != "NONE"
casos_sin_sitio = (
    (df_work["PHENOTYPE_CANCER_STATUS"] == "CASO") &
    (df_work["CANCER_SITE"] == "NONE")
).sum()
v1_ok = casos_sin_sitio == 0
print(f"  [V1] Casos con CANCER_SITE='NONE':    {casos_sin_sitio}")
print(f"       {'PASS ✓' if v1_ok else 'FAIL ✗ — revisar pipeline'}")

if not v1_ok:
    # Diagnóstico de casos sin sitio
    df_sin_sitio = df_work[
        (df_work["PHENOTYPE_CANCER_STATUS"] == "CASO") &
        (df_work["CANCER_SITE"] == "NONE")
    ][["PATID","EVIDENCIA_CASO_BASAL","EVIDENCIA_INCIDENCIA_MORTAL",
       "ICD10_UNDERLYING","CANCER_SOURCE"]]
    print("       DIAGNÓSTICO — primeras filas sin CANCER_SITE:")
    print(df_sin_sitio.head(10).to_string())

# V2: Todo CONTROL debe tener CANCER_SITE == "NONE"
controles_con_sitio = (
    (df_work["PHENOTYPE_CANCER_STATUS"] == "CONTROL") &
    (df_work["CANCER_SITE"] != "NONE")
).sum()
v2_ok = controles_con_sitio == 0
print(f"\n  [V2] Controles con CANCER_SITE!='NONE': {controles_con_sitio}")
print(f"       {'PASS ✓' if v2_ok else 'FAIL ✗ — revisar pipeline'}")

# ──────────────────────────────────────────────────────────────────────────────
# BLOQUE D — REPORTE DE DISTRIBUCIONES
# ──────────────────────────────────────────────────────────────────────────────
print(f"\n{SEP}")
print("BLOQUE D — REPORTE DE DISTRIBUCIONES")
print(SEP)

# D.1 Distribución total de CANCER_SITE
print(f"\n{LINE}")
print("D.1  DISTRIBUCIÓN TOTAL DE CANCER_SITE (todos los participantes)")
print(LINE)
site_total = df_work["CANCER_SITE"].value_counts(dropna=False)
for site, n in site_total.items():
    pct = n / len(df_work) * 100
    print(f"  {site:15s}  n={n:>8,}  ({pct:.3f}%)")

# D.2 Distribución solo en CASOS
print(f"\n{LINE}")
print("D.2  DISTRIBUCIÓN DE CANCER_SITE — SÓLO CASOS (n={:,})".format(n_casos_total))
print(LINE)
df_casos = df_work[df_work["PHENOTYPE_CANCER_STATUS"] == "CASO"]
site_casos = df_casos["CANCER_SITE"].value_counts(dropna=False)
for site, n in site_casos.items():
    pct_casos = n / n_casos_total * 100
    pct_total = n / len(df_work) * 100
    print(f"  {site:15s}  n={n:>5,}  ({pct_casos:5.1f}% de casos | {pct_total:.2f}% del total)")

# D.3 Tabla cruzada CANCER_SITE vs CANCER_SOURCE
print(f"\n{LINE}")
print("D.3  TABLA CRUZADA: CANCER_SITE × CANCER_SOURCE (sólo CASOS)")
print(LINE)
cross = pd.crosstab(
    df_casos["CANCER_SITE"],
    df_casos["CANCER_SOURCE"],
    margins=True,
    margins_name="TOTAL"
)
# Ordenar filas por total descendente (excluyendo fila TOTAL)
site_order = (cross.drop("TOTAL")
                   .sort_values("TOTAL", ascending=False)
                   .index.tolist()) + ["TOTAL"]
cross = cross.loc[site_order]
print(cross.to_string())

# D.4 Conteo UNKNOWN y OTHER
print(f"\n{LINE}")
print("D.4  DETALLE UNKNOWN y OTHER en CASOS")
print(LINE)
for target in ["UNKNOWN", "OTHER"]:
    sub = df_casos[df_casos["CANCER_SITE"] == target]
    print(f"\n  {target} — n={len(sub):,}")
    if len(sub) > 0:
        # Distribución de ICD10_UNDERLYING dentro del grupo
        icd_dist = sub["ICD10_UNDERLYING"].value_counts(dropna=False).head(15)
        print(f"  ICD10_UNDERLYING más frecuentes:")
        for code, cnt in icd_dist.items():
            print(f"    {str(code):10s}  n={cnt:>4,}")
        # Fuente (basal o mortalidad)
        print(f"  Fuente de clasificación:")
        print(f"    SOLO_BASAL:       {(sub['CANCER_SOURCE']=='SOLO_BASAL').sum():>5,}")
        print(f"    SOLO_MORTALIDAD:  {(sub['CANCER_SOURCE']=='SOLO_MORTALIDAD').sum():>5,}")
        print(f"    AMBOS:            {(sub['CANCER_SOURCE']=='AMBOS').sum():>5,}")

# D.5 Resumen ejecutivo
print(f"\n{SEP}")
print("RESUMEN EJECUTIVO — FASE 4b")
print(SEP)
print(f"\n  Dataset de trabajo (copia FASE 4): {len(df_work):>8,} filas")
print(f"  CASOS total:                       {n_casos_total:>8,}")
print(f"  CONTROLES total:                   {(df_work['PHENOTYPE_CANCER_STATUS']=='CONTROL').sum():>8,}")
print(f"\n  Clasificados por BASAL:            {n_basal_classified:>8,}  ({n_basal_classified/n_casos_total*100:.1f}% de casos)")
mortal_class = n_post_icd - n_basal_classified
print(f"  Clasificados por ICD10_UNDERLYING: {mortal_class:>8,}  ({mortal_class/n_casos_total*100:.1f}% de casos)")
print(f"\n  Validaciones:")
print(f"    V1 (casos sin sitio = 0):  {'PASS ✓' if v1_ok else 'FAIL ✗'}")
print(f"    V2 (ctrl con sitio = 0):   {'PASS ✓' if v2_ok else 'FAIL ✗'}")
print(f"\n  CANCER_SITE disponible. CSV NO exportado — pendiente autorización FASE 5.")

# Columnas del DataFrame de trabajo disponibles para FASE 5
COLS_DISPONIBLES = [
    "PATID","IID_EXOME","IID_WGS","GENOMIC_QC_FLAG",
    "AGE","MALE","BASE_DIABETES","EVER_SMOK",
    "EVIDENCIA_CASO_BASAL","EVIDENCIA_INCIDENCIA_MORTAL",
    "CANCER_SOURCE","PHENOTYPE_CANCER_STATUS","CANCER_SITE"
]
missing_cols_check = [c for c in COLS_DISPONIBLES if c not in df_work.columns]
print(f"\n  Columnas para FASE 5: {COLS_DISPONIBLES}")
if missing_cols_check:
    print(f"  ALERTA columnas faltantes: {missing_cols_check}")
else:
    print(f"  Todas las columnas requeridas presentes ✓")
