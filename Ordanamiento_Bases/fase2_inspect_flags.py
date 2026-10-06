import pandas as pd

BASE_DIR = "/mnt/cephfs/orgs/home/angel.pacheco/Ordanamiento_Bases/"

# ─── EXOME: hoja MCPS Flags ──────────────────────────────
print("=" * 60)
print("EXOME — hoja 'MCPS Flags (N=4,535)'")
print("=" * 60)
ex_flags = pd.read_excel(
    BASE_DIR + "2-11-21 - MCPS Exome Sequencing Update (F145K).xlsx",
    sheet_name="MCPS Flags (N=4,535)",
    nrows=10
)
print(f"Columnas: {ex_flags.columns.tolist()}")
print(ex_flags.head(10).to_string())

# ─── EXOME: valor MCPS Flag en QC Pass ───────────────────
print("\n" + "=" * 60)
print("EXOME QC Pass — distribución columna 'MCPS Flag (N=3,349)'")
print("=" * 60)
ex_qc = pd.read_excel(
    BASE_DIR + "2-11-21 - MCPS Exome Sequencing Update (F145K).xlsx",
    sheet_name="QC Pass (N=141,046)"
)
print(f"Total filas cargadas: {len(ex_qc)}")
print(f"Columnas: {ex_qc.columns.tolist()}")
print("\nDistribución MCPS Flag:")
print(ex_qc['MCPS Flag (N=3,349)'].value_counts(dropna=False).head(10))
print("\nNo-NaN (flagged en QC Pass):", ex_qc['MCPS Flag (N=3,349)'].notna().sum())

# ─── WGS: hoja MCPS QC Flags ─────────────────────────────
print("\n" + "=" * 60)
print("WGS — hoja 'MCPS QC Flags (N=177, 23 Fail)'")
print("=" * 60)
wgs_flags = pd.read_excel(
    BASE_DIR + "6-15-21 - MCPS Whole Genome Sequencing Update (Post-renaming).xlsx",
    sheet_name="MCPS QC Flags (N=177, 23 Fail)",
    nrows=10
)
print(f"Columnas: {wgs_flags.columns.tolist()}")
print(wgs_flags.head(10).to_string())
