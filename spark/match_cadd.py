from pyspark.sql import SparkSession
import pyspark.sql.functions as F
from pyspark.sql.functions import broadcast

spark = SparkSession.builder.appName("Match_CADD").getOrCreate()

print("1. Cargando catálogo de variantes del paciente...")
df_vars = spark.read.csv("/mnt/cephfs/orgs/home/angel.pacheco/update2_variants_map.csv", header=True, inferSchema=True)
df_vars = df_vars.withColumn("POS", F.col("POS").cast("integer"))

print("2. Cargando base de datos masiva CADD...")
# CADD tiene un encabezado que empieza con #, lo saltamos y nombramos las columnas
df_cadd = spark.read.csv("/mnt/cephfs/orgs/shared/CADD/whole_genome_SNVs.tsv.gz", sep='\t', comment='#', header=False)

# CADD columns: 0=Chrom, 1=Pos, 2=Ref, 3=Alt, 4=RawScore, 5=PHRED
df_cadd = df_cadd.select(
    F.col("_c0").alias("CHR"),
    F.col("_c1").cast("integer").alias("POS"),
    F.col("_c2").alias("REF"),
    F.col("_c3").alias("ALT"),
    F.col("_c5").cast("float").alias("PHRED")  # Este es el score de patogenicidad real
)

print("3. Ejecutando Broadcast Inner Join...")
# Forzamos a que el df_vars se copie a la RAM de todos los workers para cruzar volando
matched_df = df_cadd.join(broadcast(df_vars), ["CHR", "POS", "REF", "ALT"], "inner")

print("4. Guardando matriz de pesos en formato Parquet...")
out_path = "/mnt/cephfs/orgs/home/angel.pacheco/matched_weights.parquet"
matched_df.write.mode("overwrite").parquet(out_path)

print("¡Cruce exitoso!")
spark.stop()
