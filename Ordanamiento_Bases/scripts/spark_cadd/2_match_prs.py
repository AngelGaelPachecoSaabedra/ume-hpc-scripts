from pyspark.sql import SparkSession
import pyspark.sql.functions as F
from pyspark.sql.functions import broadcast

spark = SparkSession.builder.appName("Match_PRS_Diabetes").getOrCreate()

# ----------------------------------------------------------------
# 1. Cargar catálogo de variantes del Zarr (CHR, POS, REF, ALT)
# ----------------------------------------------------------------
print("1. Cargando catálogo de variantes del Zarr...")
df_vars = spark.read.csv(
    "/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes/update2_variants_map.csv",
    header=True, inferSchema=True
)
df_vars = df_vars.withColumn("POS", F.col("POS").cast("integer"))

# ----------------------------------------------------------------
# 2. Cargar pesos PGS limpios (CHR, POS, EA, OA, BETA)
# ----------------------------------------------------------------
print("2. Cargando pesos PGS (PGS000014 - Diabetes T2)...")
df_pgs = spark.read.csv(
    "/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes/PGS000014_clean_hg38.csv",
    header=True, inferSchema=True
)
df_pgs = df_pgs.withColumn("POS", F.col("POS").cast("integer"))

# ----------------------------------------------------------------
# 3. Broadcast Join con alineación de alelos
#
# El PGS tiene effect_allele (EA) y other_allele (OA).
# En el Zarr tenemos REF y ALT.
# Hay dos casos válidos:
#   Caso A: EA == ALT y OA == REF  →  BETA se usa directo
#   Caso B: EA == REF y OA == ALT  →  BETA se invierte (flip = true)
# ----------------------------------------------------------------
print("3. Ejecutando Broadcast Join con alineación de alelos...")

# Caso A: effect_allele coincide con ALT (caso directo)
match_direct = df_vars.join(
    broadcast(df_pgs),
    (df_vars["CHR"] == df_pgs["CHR"]) &
    (df_vars["POS"] == df_pgs["POS"]) &
    (df_vars["ALT"] == df_pgs["EA"]) &
    (df_vars["REF"] == df_pgs["OA"]),
    "inner"
).select(
    df_vars["CHR"],
    df_vars["POS"],
    df_vars["REF"],
    df_vars["ALT"],
    df_pgs["BETA"],
    F.lit(False).alias("FLIP")
)

# Caso B: effect_allele coincide con REF (caso invertido)
match_flip = df_vars.join(
    broadcast(df_pgs),
    (df_vars["CHR"] == df_pgs["CHR"]) &
    (df_vars["POS"] == df_pgs["POS"]) &
    (df_vars["ALT"] == df_pgs["OA"]) &
    (df_vars["REF"] == df_pgs["EA"]),
    "inner"
).select(
    df_vars["CHR"],
    df_vars["POS"],
    df_vars["REF"],
    df_vars["ALT"],
    df_pgs["BETA"],
    F.lit(True).alias("FLIP")
)

# Unir ambos casos
matched_df = match_direct.unionByName(match_flip)

# Estadísticas
n_direct = match_direct.count()
n_flip   = match_flip.count()
n_total  = n_direct + n_flip
print(f"   Matches directos (EA=ALT): {n_direct:,}")
print(f"   Matches invertidos (EA=REF): {n_flip:,}")
print(f"   Total matches: {n_total:,}")

# ----------------------------------------------------------------
# 4. Guardar como Parquet
# ----------------------------------------------------------------
print("4. Guardando pesos alineados en Parquet...")
out_path = "/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes/matched_prs_weights.parquet"
matched_df.write.mode("overwrite").parquet(out_path)

print(f"Cruce PRS exitoso! Archivo: {out_path}")
spark.stop()
