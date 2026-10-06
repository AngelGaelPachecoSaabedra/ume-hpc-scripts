"""
FASE 5: EXPORTACIÓN Y REPORTE FINAL
MCPS — Estudio Prospectivo Ciudad de México
Pipeline completo reproducible — pandas
"""

import pandas as pd
import numpy as np
import os
from datetime import datetime

# ─── RUTAS ────────────────────────────────────────────────────────────────────
BASE_DIR   = "/mnt/cephfs/orgs/home/angel.pacheco/Ordanamiento_Bases/"
PATH_DENOM = BASE_DIR + "denominador_genomico_F145K.csv"
PATH_BASAL = BASE_DIR + "Bases_DNAnexus/MCPS BASELINE.csv"
PATH_MORT  = BASE_DIR + "Bases_DNAnexus/MCPS MORTALITY (DNAnexus version 2.1).csv"
OUT_CSV    = BASE_DIR + "MCPS_Cohorte_Casos_Controles_Oncologicos_F145K.csv"
OUT_META   = BASE_DIR + "MCPS_Cohorte_Metadata.txt"

CANCER_VARS = [
    "BASE_LUNGCANCER", "BASE_OTHCANCER", "BASE_PROSTATECANCER",
    "BASE_CERVCANCER", "BASE_BREASTCANCER", "BASE_STOMCANCER", "BASE_ORALCANCER"
]
BASAL_COLS = ["PATID", "AGE", "MALE", "BASE_DIABETES", "EVER_SMOK"] + CANCER_VARS
MORT_COLS  = ["PATID", "STATUS", "D040", "ICD10_UNDERLYING"]

VALID_SITES   = {"NONE","BREAST","LUNG","PROSTATE","CERVICAL",
                 "GI","ORAL","HEMATOLOGIC","OTHER","UNKNOWN"}
VALID_PHENO   = {"CASO", "CONTROL"}
VALID_BINARY  = {0, 1}

FINAL_COLS = [
    "PATID", "IID_EXOME", "IID_WGS", "GENOMIC_QC_FLAG",
    "AGE", "MALE", "BASE_DIABETES", "EVER_SMOK",
    "EVIDENCIA_CASO_BASAL", "EVIDENCIA_INCIDENCIA_MORTAL",
    "CANCER_SOURCE", "CANCER_SITE", "PHENOTYPE_CANCER_STATUS"
]

SEP  = "=" * 68
LINE = "─" * 50
ERRORS = []

# ══════════════════════════════════════════════════════════════════════════════
# BLOQUE 0 — REPRODUCIR PIPELINE FASE 4 + 4b
# ══════════════════════════════════════════════════════════════════════════════
print(SEP)
print("BLOQUE 0 — REPRODUCCIÓN PIPELINE FASE 4 + 4b")
print(SEP)

df_denom = pd.read_csv(PATH_DENOM, dtype={"PATID": str, "EXOME_MCPS_FLAG": str})
df_basal = pd.read_csv(PATH_BASAL, usecols=BASAL_COLS, dtype={"PATID": str})
df_mort  = pd.read_csv(PATH_MORT,  usecols=MORT_COLS,  dtype={"PATID": str})

df = df_denom.merge(df_basal, on="PATID", how="inner")
df = df.merge(df_mort, on="PATID", how="inner")
df = df[df["STATUS"] != "U"].copy()

df["D040"] = (df["D040"] == 1).astype(int)
df["EVIDENCIA_CASO_BASAL"]        = df[CANCER_VARS].max(axis=1).clip(upper=1).fillna(0).astype(int)
df["EVIDENCIA_INCIDENCIA_MORTAL"] = df["D040"].astype(int)

caso_mask    = (df["EVIDENCIA_CASO_BASAL"] == 1) | (df["EVIDENCIA_INCIDENCIA_MORTAL"] == 1)
df["PHENOTYPE_CANCER_STATUS"] = np.where(caso_mask, "CASO", "CONTROL")

conds_src = [
    (df["EVIDENCIA_CASO_BASAL"] == 1) & (df["EVIDENCIA_INCIDENCIA_MORTAL"] == 0),
    (df["EVIDENCIA_CASO_BASAL"] == 0) & (df["EVIDENCIA_INCIDENCIA_MORTAL"] == 1),
    (df["EVIDENCIA_CASO_BASAL"] == 1) & (df["EVIDENCIA_INCIDENCIA_MORTAL"] == 1),
    (df["EVIDENCIA_CASO_BASAL"] == 0) & (df["EVIDENCIA_INCIDENCIA_MORTAL"] == 0),
]
df["CANCER_SOURCE"] = np.select(conds_src,
                                ["SOLO_BASAL","SOLO_MORTALIDAD","AMBOS","NINGUNO"],
                                default="ERROR")

# CANCER_SITE (FASE 4b)
df["CANCER_SITE"] = "NONE"
for col, label in [
    ("BASE_OTHCANCER",    "OTHER"),
    ("BASE_ORALCANCER",   "ORAL"),
    ("BASE_STOMCANCER",   "GI"),
    ("BASE_CERVCANCER",   "CERVICAL"),
    ("BASE_PROSTATECANCER","PROSTATE"),
    ("BASE_LUNGCANCER",   "LUNG"),
    ("BASE_BREASTCANCER", "BREAST"),
]:
    df.loc[df[col] == 1, "CANCER_SITE"] = label

def classify_icd10(code):
    if pd.isna(code) or str(code).strip() == "":
        return "UNKNOWN"
    c = str(code).strip().upper()
    if len(c) < 3 or c[0] != "C":
        return "OTHER"
    try:
        num = int(c[1:3])
    except ValueError:
        return "OTHER"
    if num in (33, 34):   return "LUNG"
    if num == 50:         return "BREAST"
    if num == 61:         return "PROSTATE"
    if num == 53:         return "CERVICAL"
    if 15 <= num <= 26:   return "GI"
    if 81 <= num <= 96:   return "HEMATOLOGIC"
    if num in (76, 80, 97): return "UNKNOWN"
    return "OTHER"

mortal_only = (
    (df["PHENOTYPE_CANCER_STATUS"] == "CASO") &
    (df["CANCER_SITE"] == "NONE") &
    (df["EVIDENCIA_INCIDENCIA_MORTAL"] == 1)
)
df.loc[mortal_only, "CANCER_SITE"] = df.loc[mortal_only, "ICD10_UNDERLYING"].apply(classify_icd10)

print(f"  Pipeline reproducido. Shape: {df.shape[0]:,} filas x {df.shape[1]} cols")

# Copiar para trabajo
df_final = df.copy()

# ══════════════════════════════════════════════════════════════════════════════
# BLOQUE 1 — VERIFICACIONES PREVIAS
# ══════════════════════════════════════════════════════════════════════════════
print(f"\n{SEP}")
print("BLOQUE 1 — VERIFICACIONES PREVIAS")
print(SEP)

# V1: PATID único
n_patid_dup = df_final["PATID"].duplicated().sum()
v1_ok = n_patid_dup == 0
print(f"\n  [V1] PATID únicos — duplicados encontrados: {n_patid_dup}")
print(f"       {'PASS ✓' if v1_ok else 'FAIL ✗'}")
if not v1_ok:
    ERRORS.append(f"V1: {n_patid_dup} PATIDs duplicados")

# V2: sin filas duplicadas completas
n_row_dup = df_final.duplicated().sum()
v2_ok = n_row_dup == 0
print(f"\n  [V2] Filas completamente duplicadas: {n_row_dup}")
print(f"       {'PASS ✓' if v2_ok else 'FAIL ✗'}")
if not v2_ok:
    ERRORS.append(f"V2: {n_row_dup} filas duplicadas")

# V3: variables binarias solo contienen 0/1
for varbin in ["GENOMIC_QC_FLAG", "EVIDENCIA_CASO_BASAL", "EVIDENCIA_INCIDENCIA_MORTAL"]:
    vals = set(df_final[varbin].dropna().unique())
    ok   = vals <= VALID_BINARY
    print(f"\n  [V3] {varbin} — valores únicos: {vals}")
    print(f"       {'PASS ✓' if ok else 'FAIL ✗ — valores inesperados: ' + str(vals - VALID_BINARY)}")
    if not ok:
        ERRORS.append(f"V3: {varbin} contiene valores no binarios: {vals - VALID_BINARY}")

# V4: PHENOTYPE_CANCER_STATUS solo CASO/CONTROL
pheno_vals = set(df_final["PHENOTYPE_CANCER_STATUS"].unique())
v4_ok = pheno_vals <= VALID_PHENO
print(f"\n  [V4] PHENOTYPE_CANCER_STATUS — valores únicos: {pheno_vals}")
print(f"       {'PASS ✓' if v4_ok else 'FAIL ✗ — valores inesperados: ' + str(pheno_vals - VALID_PHENO)}")
if not v4_ok:
    ERRORS.append(f"V4: PHENOTYPE_CANCER_STATUS valores inesperados: {pheno_vals - VALID_PHENO}")

# V5: CANCER_SITE solo valores válidos
site_vals = set(df_final["CANCER_SITE"].unique())
v5_ok = site_vals <= VALID_SITES
print(f"\n  [V5] CANCER_SITE — valores únicos: {site_vals}")
print(f"       {'PASS ✓' if v5_ok else 'FAIL ✗ — valores inesperados: ' + str(site_vals - VALID_SITES)}")
if not v5_ok:
    ERRORS.append(f"V5: CANCER_SITE valores inesperados: {site_vals - VALID_SITES}")

# V6: dimensiones esperadas
v6_rows = len(df_final) == 138688
v6_ok   = v6_rows
print(f"\n  [V6] Filas esperadas 138,688 — filas actuales: {len(df_final):,}")
print(f"       {'PASS ✓' if v6_ok else 'FAIL ✗'}")
if not v6_ok:
    ERRORS.append(f"V6: filas esperadas 138688, encontradas {len(df_final)}")

# ── ABORTAR si hay errores ────────────────────────────────────────────────────
if ERRORS:
    print(f"\n{'!'*68}")
    print("PIPELINE ABORTADO — ERRORES DETECTADOS:")
    for e in ERRORS:
        print(f"  ✗ {e}")
    print(f"{'!'*68}")
    raise SystemExit(1)

print(f"\n  Todas las verificaciones superadas. Procediendo a exportación.")

# ══════════════════════════════════════════════════════════════════════════════
# BLOQUE 2 — SELECCIÓN ESTRICTA DE COLUMNAS (13, en orden exacto)
# ══════════════════════════════════════════════════════════════════════════════
print(f"\n{SEP}")
print("BLOQUE 2 — SELECCIÓN ESTRICTA DE COLUMNAS")
print(SEP)

missing_final = [c for c in FINAL_COLS if c not in df_final.columns]
if missing_final:
    raise ValueError(f"Columnas faltantes en DataFrame: {missing_final}")

df_export = df_final[FINAL_COLS].copy()

print(f"  Columnas seleccionadas ({len(FINAL_COLS)}):")
for i, col in enumerate(FINAL_COLS, 1):
    dtype = df_export[col].dtype
    n_na  = df_export[col].isna().sum()
    print(f"    {i:2d}. {col:35s} dtype={str(dtype):10s}  NaN={n_na:,}")

assert list(df_export.columns) == FINAL_COLS, "ERROR: orden de columnas incorrecto"
print(f"\n  Orden de columnas verificado ✓")
print(f"  Shape final: {df_export.shape[0]:,} filas x {df_export.shape[1]} columnas")

# ══════════════════════════════════════════════════════════════════════════════
# BLOQUE 3 — EXPORTACIÓN CSV
# ══════════════════════════════════════════════════════════════════════════════
print(f"\n{SEP}")
print("BLOQUE 3 — EXPORTACIÓN CSV")
print(SEP)

# EVER_SMOK: NaN se conserva como vacío en CSV (comportamiento default de pandas)
df_export.to_csv(OUT_CSV, index=False, na_rep="")

csv_size_mb = os.path.getsize(OUT_CSV) / (1024 * 1024)
print(f"  CSV escrito: {OUT_CSV}")
print(f"  Tamaño en disco: {csv_size_mb:.2f} MB")

# Verificación post-escritura: recargar y comprobar dimensiones
df_verify = pd.read_csv(OUT_CSV, dtype={"PATID": str}, nrows=5)
cols_verify = pd.read_csv(OUT_CSV, nrows=0).columns.tolist()
n_rows_verify = sum(1 for _ in open(OUT_CSV)) - 1  # contar líneas menos header
print(f"  Post-escritura — columnas: {len(cols_verify)}, filas (conteo líneas): {n_rows_verify:,}")
assert len(cols_verify) == 13, f"ERROR columnas post-escritura: {len(cols_verify)}"
assert n_rows_verify == 138688, f"ERROR filas post-escritura: {n_rows_verify}"
print(f"  Verificación post-escritura: PASS ✓")

# ══════════════════════════════════════════════════════════════════════════════
# BLOQUE 4 — ARCHIVO METADATA
# ══════════════════════════════════════════════════════════════════════════════
print(f"\n{SEP}")
print("BLOQUE 4 — GENERACIÓN ARCHIVO METADATA")
print(SEP)

ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
n_casos     = (df_export["PHENOTYPE_CANCER_STATUS"] == "CASO").sum()
n_controles = (df_export["PHENOTYPE_CANCER_STATUS"] == "CONTROL").sum()
n_solo_basal = (df_export["CANCER_SOURCE"] == "SOLO_BASAL").sum()
n_solo_mort  = (df_export["CANCER_SOURCE"] == "SOLO_MORTALIDAD").sum()
n_ambos      = (df_export["CANCER_SOURCE"] == "AMBOS").sum()
n_flag6      = (df_export["GENOMIC_QC_FLAG"] == 1).sum()
site_dist    = df_export[df_export["PHENOTYPE_CANCER_STATUS"]=="CASO"]["CANCER_SITE"].value_counts()

meta_lines = f"""================================================================================
MCPS — METADATA COHORTE ONCOLÓGICA
Archivo: MCPS_Cohorte_Casos_Controles_Oncologicos_F145K.csv
Generado: {ts}
================================================================================

1. FUENTES DE DATOS
   - Denominador genómico: denominador_genomico_F145K.csv
     (Exoma QC Pass N=141,046 + WGS QC Pass N=9,950; WGS ⊂ Exoma)
   - Base basal: MCPS BASELINE.csv (DNAnexus)
   - Mortalidad: MCPS MORTALITY (DNAnexus version 2.1).csv
   - Seguimiento: hasta 31 de diciembre de 2020

2. DEFINICIÓN DE CASO vs CONTROL
   CASO:
     EVIDENCIA_CASO_BASAL = 1   (autorreporte de cáncer al basal)
     OR
     EVIDENCIA_INCIDENCIA_MORTAL = 1   (muerte oncológica D040=1)

   CONTROL:
     EVIDENCIA_CASO_BASAL = 0
     AND
     EVIDENCIA_INCIDENCIA_MORTAL = 0

3. EVIDENCIA_CASO_BASAL
   Variable = 1 si cualquiera de las siguientes es 1:
     BASE_LUNGCANCER     (cáncer de pulmón, autorreporte)
     BASE_OTHCANCER      (otro cáncer, autorreporte)
     BASE_PROSTATECANCER (cáncer de próstata, autorreporte)
     BASE_CERVCANCER     (cáncer cervical, autorreporte)
     BASE_BREASTCANCER   (cáncer de mama, autorreporte)
     BASE_STOMCANCER     (cáncer digestivo, autorreporte)
     BASE_ORALCANCER     (cáncer oral/vías altas, autorreporte)
   Codificación basal: 1=sí, 0=no. Fuente: data dictionary MCPS baseline.

4. EVIDENCIA_INCIDENCIA_MORTAL
   Variable = 1 si D040 == 1.
   D040 = ANY NEOPLASTIC DEATH (suma de D025–D039), variable derivada MCPS.
   Basado en ICD10_UNDERLYING (causa adjudicada única de muerte).
   Incluye códigos inespecíficos: C76x (n=41), C80x (n=101), C97X (n=2)
   en causa adjudicada — decisión ontológica de incluirlos (D038/D039 MCPS).

5. EXCLUSIÓN STATUS='U'
   Participantes con STATUS='U' (vital status incierto al 31-dic-2020)
   excluidos del análisis prospectivo, conforme a recomendación MCPS.
   Excluidos: 2,358 participantes (23 con autorreporte basal positivo,
   no clasificables sin D040).

6. DEFINICIÓN CANCER_SITE
   Clasificación por sitio anatómico. Prioridad para basal:
     BREAST > LUNG > PROSTATE > CERVICAL > GI > ORAL > OTHER
   Para SOLO_MORTALIDAD (sin sitio basal), se usa ICD10_UNDERLYING:
     LUNG:        C33–C34
     BREAST:      C50
     PROSTATE:    C61
     CERVICAL:    C53
     GI:          C15–C26
     HEMATOLOGIC: C81–C96
     UNKNOWN:     C76, C80, C97 (ill-defined/unspecified/multiple)
     OTHER:       todos los demás C00–C99 no clasificados arriba
     NONE:        controles (PHENOTYPE_CANCER_STATUS = CONTROL)

7. GENOMIC_QC_FLAG
   Flag de riesgo de linkage genómico.
   = 1 si EXOME_MCPS_FLAG contiene '6' (Linkage may be unreliable).
   Recomendación: análisis de sensibilidad excluyendo GENOMIC_QC_FLAG=1.
   N con flag=1: {n_flag6:,} ({n_flag6/len(df_export)*100:.2f}% del total)

8. TAMAÑO FINAL DE LA COHORTE ANALÍTICA
   Total participantes:   {len(df_export):>8,}
   CASOS:                 {n_casos:>8,}  ({n_casos/len(df_export)*100:.2f}%)
   CONTROLES:             {n_controles:>8,}  ({n_controles/len(df_export)*100:.2f}%)

9. DESGLOSE DE CASOS
   SOLO_BASAL:            {n_solo_basal:>8,}  ({n_solo_basal/n_casos*100:.1f}% de casos)
   SOLO_MORTALIDAD:       {n_solo_mort:>8,}  ({n_solo_mort/n_casos*100:.1f}% de casos)
   AMBOS:                 {n_ambos:>8,}  ({n_ambos/n_casos*100:.1f}% de casos)

10. DISTRIBUCIÓN CANCER_SITE (en CASOS)
{chr(10).join(f"   {s:15s}  n={n:>5,}  ({n/n_casos*100:.1f}%)" for s, n in site_dist.items())}

11. COVARIABLES INCLUIDAS
   AGE:           edad autorreportada al reclutamiento (años; rango: 35–112)
   MALE:          sexo biológico (1=hombre, 0=mujer)
   BASE_DIABETES: diabetes al basal, autorreporte (1=sí, 0=no; NaN=0)
   EVER_SMOK:     alguna vez fumador (1=sí, 0=no; NaN=43 participantes)

12. NOTAS METODOLÓGICAS
   - Estudio MCPS: cohorte prospectiva de adultos ≥35 años reclutados
     en Ciudad de México (Coyoacán e Iztapalapa), ~1998–2004.
   - El denominador genómico F145K (exoma QC Pass) es subconjunto
     de los 159,517 participantes MCPS totales (88.4% de cobertura).
   - 18,471 participantes sin datos genómicos excluidos (4,064 MCPSDS
     sin renombramiento genómico + 14,407 MCPSRGN sin QC pass).
   - Los 9,949 participantes con WGS (IID_WGS no vacío) forman un
     subconjunto del denominador exómico.
   - La restricción genómica no introduce sesgo diferencial sobre
     GENOMIC_QC_FLAG entre casos y controles (0.55% vs 0.59%).
   - Historia de enfermedad: autorreporte únicamente (basal).
     No hay datos de incidencia clínica no letal en este pipeline.

================================================================================
Fin del archivo de metadata.
================================================================================
"""

with open(OUT_META, "w", encoding="utf-8") as fh:
    fh.write(meta_lines)

meta_size_kb = os.path.getsize(OUT_META) / 1024
print(f"  Metadata escrita: {OUT_META}")
print(f"  Tamaño metadata: {meta_size_kb:.1f} KB")

# ══════════════════════════════════════════════════════════════════════════════
# BLOQUE 5 — VALIDACIÓN FINAL Y MENSAJE DE ÉXITO
# ══════════════════════════════════════════════════════════════════════════════
print(f"\n{SEP}")
print("BLOQUE 5 — VALIDACIÓN FINAL Y CONFIRMACIÓN DE ÉXITO")
print(SEP)

print(f"""
  ╔══════════════════════════════════════════════════════════════╗
  ║           EXPORTACIÓN COMPLETADA EXITOSAMENTE               ║
  ╠══════════════════════════════════════════════════════════════╣
  ║                                                              ║
  ║  CSV exportado:                                              ║
  ║    {OUT_CSV}
  ║                                                              ║
  ║  Metadata exportada:                                         ║
  ║    {OUT_META}
  ║                                                              ║
  ╠══════════════════════════════════════════════════════════════╣
  ║  Filas exportadas:        {len(df_export):>10,}                       ║
  ║  Columnas exportadas:     {len(df_export.columns):>10,}                       ║
  ║  Tamaño CSV en disco:     {csv_size_mb:>9.2f} MB                      ║
  ╠══════════════════════════════════════════════════════════════╣
  ║  CASOS:                   {n_casos:>10,}                       ║
  ║  CONTROLES:               {n_controles:>10,}                       ║
  ╠══════════════════════════════════════════════════════════════╣
  ║  Timestamp:  {ts}                       ║
  ╚══════════════════════════════════════════════════════════════╝
""")

print(f"  Verificaciones finales:")
print(f"    Filas == 138,688:  {len(df_export) == 138688} ✓")
print(f"    Columnas == 13:    {len(df_export.columns) == 13} ✓")
print(f"    PATID únicos:      {not df_export['PATID'].duplicated().any()} ✓")
print(f"    CSV existe en disco: {os.path.isfile(OUT_CSV)} ✓")
print(f"    Metadata existe:   {os.path.isfile(OUT_META)} ✓")

print(f"\n{SEP}")
print("PIPELINE MCPS ONCOLÓGICO F145K — COMPLETADO")
print(SEP)
