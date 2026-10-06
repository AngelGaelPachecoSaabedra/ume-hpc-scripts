import pandas as pd

BASE_DIR = "/mnt/cephfs/orgs/home/angel.pacheco/Ordanamiento_Bases/"

print("=" * 65)
print("FASE 2: CONSTRUCCIÓN DEL DENOMINADOR GENÓMICO")
print("=" * 65)

# ─────────────────────────────────────────────────────────────
# 1. CARGAR QC PASS — EXOMA
# ─────────────────────────────────────────────────────────────
print("\n[1] Cargando Exoma QC Pass...")
ex_qc = pd.read_excel(
    BASE_DIR + "2-11-21 - MCPS Exome Sequencing Update (F145K).xlsx",
    sheet_name="QC Pass (N=141,046)",
    usecols=["Sample Name", "MCPS Flag (N=3,349)"]
)
print(f"    Filas brutas: {len(ex_qc)}")

# IID = Sample Name
ex_qc = ex_qc.rename(columns={
    "Sample Name": "IID",
    "MCPS Flag (N=3,349)": "EXOME_MCPS_FLAG"
})

# Limpiar filas vacías de IID
ex_qc = ex_qc[ex_qc["IID"].notna() & (ex_qc["IID"] != "")]
print(f"    Filas con IID válido: {len(ex_qc)}")
print(f"    Con MCPS Flag (cualquiera): {ex_qc['EXOME_MCPS_FLAG'].notna().sum()}")
print(f"    Distribución flags: {ex_qc['EXOME_MCPS_FLAG'].value_counts(dropna=False).to_dict()}")

# Muestra IID
print(f"    Muestra IID: {ex_qc['IID'].head(3).tolist()}")

# ─────────────────────────────────────────────────────────────
# 2. CARGAR QC PASS — WGS
# ─────────────────────────────────────────────────────────────
print("\n[2] Cargando WGS QC Pass...")
wgs_qc = pd.read_excel(
    BASE_DIR + "6-15-21 - MCPS Whole Genome Sequencing Update (Post-renaming).xlsx",
    sheet_name="QC Pass (N=9,950)",
    usecols=["Sample Name"]
)
print(f"    Filas brutas: {len(wgs_qc)}")

wgs_qc = wgs_qc.rename(columns={"Sample Name": "IID"})
wgs_qc = wgs_qc[wgs_qc["IID"].notna() & (wgs_qc["IID"] != "")]
print(f"    Filas con IID válido: {len(wgs_qc)}")
print(f"    Muestra IID: {wgs_qc['IID'].head(3).tolist()}")

# ─────────────────────────────────────────────────────────────
# 3. CARGAR RGN_LINK_IID (IID → PATID)
# ─────────────────────────────────────────────────────────────
print("\n[3] Cargando RGN_LINK_IID...")
link = pd.read_csv(BASE_DIR + "Bases_DNAnexus/RGN_LINK_IID.csv",
                   usecols=["PATID", "IID"])
print(f"    Total filas: {len(link)}")
print(f"    IIDs únicos: {link['IID'].nunique()}")
print(f"    PATIDs únicos: {link['PATID'].nunique()}")
print(f"    Muestra IID en link: {link['IID'].head(3).tolist()}")

# Verificar formato IID coincide
print(f"\n    Muestra formato IID Exome: {ex_qc['IID'].iloc[0]}")
print(f"    Muestra formato IID RGN_LINK: {link['IID'].iloc[0]}")

# ─────────────────────────────────────────────────────────────
# 4. LINKAGE IID → PATID para EXOMA
# ─────────────────────────────────────────────────────────────
print("\n[4] Linkage Exoma → PATID...")
ex_linked = ex_qc.merge(link, on="IID", how="left")
print(f"    Exoma QC Pass total: {len(ex_linked)}")
print(f"    Con PATID (linked): {ex_linked['PATID'].notna().sum()}")
print(f"    Sin PATID (no link): {ex_linked['PATID'].isna().sum()}")

# IIDs duplicados en exoma (mismo IID, diferente entrada)
dup_iid_ex = ex_linked['IID'].duplicated().sum()
print(f"    IIDs duplicados en exoma: {dup_iid_ex}")

# PATIDs duplicados (mismo PATID, múltiples muestras)
dup_patid_ex = ex_linked[ex_linked['PATID'].notna()]['PATID'].duplicated().sum()
print(f"    PATIDs duplicados en exoma (múltiples muestras): {dup_patid_ex}")

# ─────────────────────────────────────────────────────────────
# 5. LINKAGE IID → PATID para WGS
# ─────────────────────────────────────────────────────────────
print("\n[5] Linkage WGS → PATID...")
wgs_linked = wgs_qc.merge(link, on="IID", how="left")
print(f"    WGS QC Pass total: {len(wgs_linked)}")
print(f"    Con PATID (linked): {wgs_linked['PATID'].notna().sum()}")
print(f"    Sin PATID (no link): {wgs_linked['PATID'].isna().sum()}")

dup_patid_wgs = wgs_linked[wgs_linked['PATID'].notna()]['PATID'].duplicated().sum()
print(f"    PATIDs duplicados en WGS: {dup_patid_wgs}")

# ─────────────────────────────────────────────────────────────
# 6. CONJUNTOS DE PATID ÚNICOS
# ─────────────────────────────────────────────────────────────
print("\n[6] Conjuntos de PATIDs únicos QC-viables...")

set_exome = set(ex_linked[ex_linked['PATID'].notna()]['PATID'])
set_wgs   = set(wgs_linked[wgs_linked['PATID'].notna()]['PATID'])

intersection = set_exome & set_wgs
union_set    = set_exome | set_wgs

print(f"\n    N Exoma (PATIDs únicos, linked): {len(set_exome)}")
print(f"    N WGS   (PATIDs únicos, linked): {len(set_wgs)}")
print(f"    Intersección (ambas plataformas): {len(intersection)}")
print(f"    Unión (al menos una plataforma): {len(union_set)}")
print(f"    Solo Exoma (no WGS): {len(set_exome - set_wgs)}")
print(f"    Solo WGS (no Exoma): {len(set_wgs - set_exome)}")

# ─────────────────────────────────────────────────────────────
# 7. GUARDAR DENOMINADOR FINAL
# ─────────────────────────────────────────────────────────────
print("\n[7] Construyendo tabla denominador final (unión)...")

# Para Exoma: tomar un PATID por participante (el más limpio = sin flag si posible)
ex_dedup = (ex_linked[ex_linked['PATID'].notna()]
            .sort_values('EXOME_MCPS_FLAG', na_position='first')  # NaN (sin flag) primero
            .drop_duplicates(subset='PATID', keep='first')
            [['PATID', 'IID', 'EXOME_MCPS_FLAG']]
            .rename(columns={'IID': 'IID_EXOME'}))

# Para WGS
wgs_dedup = (wgs_linked[wgs_linked['PATID'].notna()]
             .drop_duplicates(subset='PATID', keep='first')
             [['PATID', 'IID']]
             .rename(columns={'IID': 'IID_WGS'}))

# Denominador: outer join para unión
denom = pd.merge(ex_dedup, wgs_dedup, on='PATID', how='outer')
denom['HAS_EXOME'] = denom['IID_EXOME'].notna().astype(int)
denom['HAS_WGS']   = denom['IID_WGS'].notna().astype(int)

print(f"    Total denominador (unión): {len(denom)}")
print(f"    Con exoma:       {denom['HAS_EXOME'].sum()}")
print(f"    Con WGS:         {denom['HAS_WGS'].sum()}")
print(f"    Con ambos:       {((denom['HAS_EXOME']==1) & (denom['HAS_WGS']==1)).sum()}")

# Guardar
out_path = BASE_DIR + "denominador_genomico_F145K.csv"
denom.to_csv(out_path, index=False)
print(f"\n    Guardado: {out_path}")
print("\nPrimeras 5 filas del denominador:")
print(denom.head().to_string())
