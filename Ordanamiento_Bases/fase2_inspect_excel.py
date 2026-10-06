import pandas as pd

BASE_DIR = "/mnt/cephfs/orgs/home/angel.pacheco/Ordanamiento_Bases/"

# ─── EXOME ───────────────────────────────────────────────
print("=" * 60)
print("EXOME — hoja 'QC Pass (N=141,046)'")
print("=" * 60)

ex_qc = pd.read_excel(
    BASE_DIR + "2-11-21 - MCPS Exome Sequencing Update (F145K).xlsx",
    sheet_name="QC Pass (N=141,046)",
    nrows=5
)
print(f"Columnas ({len(ex_qc.columns)}): {ex_qc.columns.tolist()}")
print("\nPrimeras 5 filas:")
print(ex_qc.to_string())

# ─── EXOME Summary ───────────────────────────────────────
print("\n" + "=" * 60)
print("EXOME — hoja 'Summary'")
print("=" * 60)
ex_sum = pd.read_excel(
    BASE_DIR + "2-11-21 - MCPS Exome Sequencing Update (F145K).xlsx",
    sheet_name="Summary",
    nrows=30
)
print(ex_sum.to_string())

# ─── WGS ─────────────────────────────────────────────────
print("\n" + "=" * 60)
print("WGS — hoja 'QC Pass (N=9,950)'")
print("=" * 60)

wgs_qc = pd.read_excel(
    BASE_DIR + "6-15-21 - MCPS Whole Genome Sequencing Update (Post-renaming).xlsx",
    sheet_name="QC Pass (N=9,950)",
    nrows=5
)
print(f"Columnas ({len(wgs_qc.columns)}): {wgs_qc.columns.tolist()}")
print("\nPrimeras 5 filas:")
print(wgs_qc.to_string())

# ─── WGS Summary ─────────────────────────────────────────
print("\n" + "=" * 60)
print("WGS — hoja 'Summary'")
print("=" * 60)
wgs_sum = pd.read_excel(
    BASE_DIR + "6-15-21 - MCPS Whole Genome Sequencing Update (Post-renaming).xlsx",
    sheet_name="Summary",
    nrows=30
)
print(wgs_sum.to_string())
