import pandas as pd

BASE_DIR = "/mnt/cephfs/orgs/home/angel.pacheco/Ordanamiento_Bases/"

# Cargar datos ya validados
ex_qc = pd.read_excel(
    BASE_DIR + "2-11-21 - MCPS Exome Sequencing Update (F145K).xlsx",
    sheet_name="QC Pass (N=141,046)",
    usecols=["Sample Name", "MCPS Flag (N=3,349)"]
).rename(columns={"Sample Name": "IID", "MCPS Flag (N=3,349)": "EXOME_MCPS_FLAG"})
ex_qc = ex_qc[ex_qc["IID"].notna()]
# Convertir flag a string para uniformidad
ex_qc["EXOME_MCPS_FLAG"] = ex_qc["EXOME_MCPS_FLAG"].astype(str).replace("nan", pd.NA)

wgs_qc = pd.read_excel(
    BASE_DIR + "6-15-21 - MCPS Whole Genome Sequencing Update (Post-renaming).xlsx",
    sheet_name="QC Pass (N=9,950)",
    usecols=["Sample Name"]
).rename(columns={"Sample Name": "IID"})
wgs_qc = wgs_qc[wgs_qc["IID"].notna()]

link = pd.read_csv(BASE_DIR + "Bases_DNAnexus/RGN_LINK_IID.csv", usecols=["PATID", "IID"])

# Linkage
ex_linked  = ex_qc.merge(link, on="IID", how="left")
wgs_linked = wgs_qc.merge(link, on="IID", how="left")

# Sets de PATIDs únicos
set_exome = set(ex_linked[ex_linked['PATID'].notna()]['PATID'])
set_wgs   = set(wgs_linked[wgs_linked['PATID'].notna()]['PATID'])

# ── Deduplicar exoma: por PATID, preferir sin flag (flag is NA) ──────────
ex_linked['_flag_null'] = ex_linked['EXOME_MCPS_FLAG'].isna().astype(int)
ex_dedup = (ex_linked[ex_linked['PATID'].notna()]
            .sort_values('_flag_null', ascending=False)   # sin-flag (1) primero
            .drop_duplicates(subset='PATID', keep='first')
            [['PATID', 'IID', 'EXOME_MCPS_FLAG']]
            .rename(columns={'IID': 'IID_EXOME'}))

# ── Deduplicar WGS ───────────────────────────────────────────────────────
wgs_dedup = (wgs_linked[wgs_linked['PATID'].notna()]
             .drop_duplicates(subset='PATID', keep='first')
             [['PATID', 'IID']]
             .rename(columns={'IID': 'IID_WGS'}))

# ── Denominador: outer join unión ────────────────────────────────────────
denom = pd.merge(ex_dedup, wgs_dedup, on='PATID', how='outer')
denom['HAS_EXOME'] = denom['IID_EXOME'].notna().astype(int)
denom['HAS_WGS']   = denom['IID_WGS'].notna().astype(int)

print("=" * 65)
print("RESUMEN FASE 2 — DENOMINADOR GENÓMICO")
print("=" * 65)
print(f"\nExoma QC Pass (PATID únicos linked): {len(set_exome):>10,}")
print(f"WGS   QC Pass (PATID únicos linked): {len(set_wgs):>10,}")
print(f"  WGS sin link PATID:                {wgs_linked['PATID'].isna().sum():>10,}")
print(f"Intersección (exoma ∩ WGS):          {len(set_exome & set_wgs):>10,}")
print(f"Unión (exoma ∪ WGS):                 {len(set_exome | set_wgs):>10,}")
print(f"Solo exoma (no WGS):                 {len(set_exome - set_wgs):>10,}")
print(f"Solo WGS   (no exoma):               {len(set_wgs - set_exome):>10,}")
print()
print(f"Denominador final (tabla, unión):    {len(denom):>10,}")
print(f"  Con exoma:                         {denom['HAS_EXOME'].sum():>10,}")
print(f"  Con WGS:                           {denom['HAS_WGS'].sum():>10,}")
print(f"  Con ambos (exoma + WGS):           {((denom['HAS_EXOME']==1)&(denom['HAS_WGS']==1)).sum():>10,}")
print()
print("Distribución flags exoma en denominador final:")
print(denom['EXOME_MCPS_FLAG'].value_counts(dropna=False).to_string())
print()
print("Primeras 5 filas:")
print(denom.head().to_string())

# Guardar
out_path = BASE_DIR + "denominador_genomico_F145K.csv"
denom.to_csv(out_path, index=False)
print(f"\nGuardado: {out_path}")
