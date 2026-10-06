"""
Validación PRS pipeline: Spark+CuPy (nuevo) vs baseline (previo)
PGS000363

Python: estadísticas + exporta tabla para R
Gráficos: generados por graficar_comparativa_prs_pgs000363.R
"""

import glob
import numpy as np
import pandas as pd
from scipy import stats
import warnings
warnings.filterwarnings("ignore")

# ─── RUTAS ────────────────────────────────────────────────────────────────────
DIR_NEW  = "/mnt/cephfs/orgs/home/angel.pacheco/Ordanamiento_Bases/scripts/resultados/PGS000363"
FILE_OLD = "/mnt/cephfs/hot_nvme/pgscatalog/scores/PGS000363/PGS000363_PRS_total.tsv"
OUT_DIR  = "/mnt/cephfs/orgs/home/angel.pacheco/Ordanamiento_Bases/scripts"
OUT_CSV  = f"{OUT_DIR}/prs_comparativa_pgs000363.csv"

# ─── PASO 1: PRS TOTAL NUEVO ──────────────────────────────────────────────────
print("=" * 60)
print("PASO 1 — Construyendo PRS total nuevo")
print("=" * 60)

chr_files = sorted(glob.glob(f"{DIR_NEW}/PGS000363_chr*_scores.tsv"))
print(f"  Archivos encontrados: {len(chr_files)}")

frames = []
for f in chr_files:
    df = pd.read_csv(f, sep="\t")
    frames.append(df)

df_all = pd.concat(frames, ignore_index=True)
print(f"  Filas totales (todos los cromosomas): {len(df_all):,}")

df_new = (
    df_all
    .dropna(subset=["sample_id", "PRS"])
    .groupby("sample_id", as_index=False)["PRS"]
    .sum()
    .rename(columns={"PRS": "PRS_NEW"})
)
print(f"  Muestras únicas (nuevo): {len(df_new):,}")

# ─── PASO 2: CARGAR BASELINE Y UNIR ──────────────────────────────────────────
print()
print("=" * 60)
print("PASO 2 — Cargando baseline y haciendo inner join")
print("=" * 60)

df_old = pd.read_csv(FILE_OLD, sep="\t")[["sample_id", "PRS_total"]].rename(
    columns={"PRS_total": "PRS_OLD"}
)
df_old = df_old.dropna(subset=["sample_id", "PRS_OLD"])
print(f"  Muestras baseline: {len(df_old):,}")

df = df_new.merge(df_old, on="sample_id", how="inner")
df = df.dropna(subset=["PRS_NEW", "PRS_OLD"])
print(f"  Muestras en inner join: {len(df):,}")
print(f"  Solo en nuevo (no en baseline): {len(df_new) - len(df):,}")
print(f"  Solo en baseline (no en nuevo): {len(df_old) - len(df):,}")

# ─── PASO 3: VALIDACIÓN NUMÉRICA ─────────────────────────────────────────────
print()
print("=" * 60)
print("PASO 3 — Validación numérica")
print("=" * 60)

r_pearson, p_pearson   = stats.pearsonr(df["PRS_NEW"], df["PRS_OLD"])
r_spearman, p_spearman = stats.spearmanr(df["PRS_NEW"], df["PRS_OLD"])

abs_diff = (df["PRS_NEW"] - df["PRS_OLD"]).abs()
mae      = float(abs_diff.mean())
max_diff = float(abs_diff.max())

p50 = float(np.percentile(abs_diff, 50))
p90 = float(np.percentile(abs_diff, 90))
p99 = float(np.percentile(abs_diff, 99))

print(f"\n  Correlaciones:")
print(f"    Pearson  r  = {r_pearson:.8f}   (p = {p_pearson:.2e})")
print(f"    Spearman rho= {r_spearman:.8f}   (p = {p_spearman:.2e})")

print(f"\n  Error absoluto:")
print(f"    MAE              = {mae:.6e}")
print(f"    Max abs diff     = {max_diff:.6e}")

print(f"\n  Percentiles del error absoluto:")
print(f"    p50 = {p50:.6e}")
print(f"    p90 = {p90:.6e}")
print(f"    p99 = {p99:.6e}")

# ─── PASO 4: VALIDACIÓN POR PERCENTILES ──────────────────────────────────────
print()
print("=" * 60)
print("PASO 4 — Validación por percentiles (1–100)")
print("=" * 60)

def assign_percentile(series):
    ranks = series.rank(pct=True)
    pctls = (ranks * 100).clip(1, 100).round().astype(int)
    return pctls

df["pct_NEW"] = assign_percentile(df["PRS_NEW"])
df["pct_OLD"] = assign_percentile(df["PRS_OLD"])
df["pct_diff"] = (df["pct_NEW"] - df["pct_OLD"]).abs()

n = len(df)
exact   = int((df["pct_diff"] == 0).sum())
within1 = int((df["pct_diff"] <= 1).sum())
within2 = int((df["pct_diff"] <= 2).sum())
within5 = int((df["pct_diff"] <= 5).sum())

print(f"\n  Concordancia de percentiles (n={n:,}):")
print(f"    Match exacto (±0)  : {exact/n*100:.2f}%  ({exact:,})")
print(f"    Dentro de ±1 pctil : {within1/n*100:.2f}%  ({within1:,})")
print(f"    Dentro de ±2 pctil : {within2/n*100:.2f}%  ({within2:,})")
print(f"    Dentro de ±5 pctil : {within5/n*100:.2f}%  ({within5:,})")

pct_err_p50 = float(np.percentile(df["pct_diff"], 50))
pct_err_p90 = float(np.percentile(df["pct_diff"], 90))
pct_err_p99 = float(np.percentile(df["pct_diff"], 99))
print(f"\n  Distribución del error en percentiles:")
print(f"    p50 = {pct_err_p50:.1f}")
print(f"    p90 = {pct_err_p90:.1f}")
print(f"    p99 = {pct_err_p99:.1f}")

# ─── EXPORTAR CSV PARA R ──────────────────────────────────────────────────────
print()
print("=" * 60)
print("PASO 5 — Exportando tabla para R")
print("=" * 60)

df[["sample_id", "PRS_NEW", "PRS_OLD", "pct_NEW", "pct_OLD", "pct_diff"]].to_csv(
    OUT_CSV, index=False
)
print(f"  Tabla guardada: {OUT_CSV}")
print(f"  Filas: {len(df):,}  |  Columnas: sample_id, PRS_NEW, PRS_OLD, pct_NEW, pct_OLD, pct_diff")

print()
print("=" * 60)
print("ANÁLISIS PYTHON COMPLETADO — ejecutar R para gráficos")
print("=" * 60)
