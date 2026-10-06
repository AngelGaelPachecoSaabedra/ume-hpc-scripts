#!/bin/bash
#SBATCH --job-name=spark_prs
#SBATCH --account=researchers
#SBATCH --qos=vip
#SBATCH --partition=normal
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=/mnt/cephfs/orgs/home/angel.pacheco/spark/spark_prs_%j.out
#SBATCH --error=/mnt/cephfs/orgs/home/angel.pacheco/spark/spark_prs_%j.err
#SBATCH -D /mnt/cephfs/orgs/home/angel.pacheco/spark

echo "=============================================="
echo "SPARK DISTRIBUIDO - CRUCE PRS DIABETES (PGS000014)"
echo "Job ID:        $SLURM_JOB_ID"
echo "Nodos:         $SLURM_NODELIST"
echo "=============================================="

export SPARK_HOME=/mnt/cephfs/biocontainers/spark/spark
export JAVA_HOME=/usr/lib/jvm/java-17-openjdk-amd64
export PATH="$JAVA_HOME/bin:$SPARK_HOME/bin:$SPARK_HOME/sbin:$PATH"
export PYSPARK_PYTHON=/usr/bin/python3

JOB_DIR=/mnt/cephfs/orgs/home/angel.pacheco/spark/job_${SLURM_JOB_ID}
export SPARK_LOG_DIR=${JOB_DIR}/logs
export SPARK_WORKER_DIR=${JOB_DIR}/work
export SPARK_LOCAL_DIRS=${JOB_DIR}/tmp

mkdir -p "$SPARK_LOG_DIR" "$SPARK_WORKER_DIR" "$SPARK_LOCAL_DIRS"

MASTER_HOST=$(scontrol show hostnames "$SLURM_NODELIST" | head -n 1)
MASTER_IP=$(getent hosts "$MASTER_HOST" | awk '{print $1}')
MASTER_URL="spark://$MASTER_IP:7077"

# ========================================
# INICIAR MASTER
# ========================================
echo ">> Iniciando Master en $MASTER_HOST..."
srun --overlap -N1 -w "$MASTER_HOST" --job-name=spark-master bash -c "
    export SPARK_HOME='$SPARK_HOME'
    export JAVA_HOME='$JAVA_HOME'
    export PATH='$PATH'
    export SPARK_LOG_DIR='$SPARK_LOG_DIR'
    cd \$SPARK_HOME
    nohup ./sbin/start-master.sh --host $MASTER_IP --port 7077 >/dev/null 2>&1
    for i in {1..30}; do
        if ss -tlnp 2>/dev/null | grep -q ':7077'; then
            while ps aux | grep -q '[o]rg.apache.spark.deploy.master.Master'; do sleep 5; done
            exit 0
        fi
        sleep 2
    done
    exit 1
" &
MASTER_SRUN_PID=$!
sleep 10

# ========================================
# INICIAR WORKERS
# ========================================
echo ">> Iniciando Workers..."
NODES=($(scontrol show hostnames "$SLURM_NODELIST"))
for NODE in "${NODES[@]}"; do
    srun --overlap -N1 -w "$NODE" --job-name=spark-worker bash -c "
        export SPARK_HOME='$SPARK_HOME'
        export JAVA_HOME='$JAVA_HOME'
        export PATH='$PATH'
        export SPARK_LOG_DIR='$SPARK_LOG_DIR'
        export SPARK_WORKER_DIR='$SPARK_WORKER_DIR'
        HOSTNAME=\$(hostname)
        cd \$SPARK_HOME
        nohup ./sbin/start-worker.sh '$MASTER_URL' \
            --host \$HOSTNAME \
            --cores $SLURM_CPUS_PER_TASK \
            --memory 50G >/dev/null 2>&1
        sleep 5
        while ps aux | grep -q '[o]rg.apache.spark.deploy.worker.Worker'; do sleep 5; done
    " &
done
sleep 20

# ========================================
# GENERAR SCRIPT PYSPARK
# ========================================
cat << 'EOF' > match_prs.py
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
EOF

# ========================================
# EJECUTAR APLICACIÓN PYSPARK
# ========================================
echo ">> Ejecutando Cruce PRS..."
$SPARK_HOME/bin/spark-submit \
    --master "$MASTER_URL" \
    --deploy-mode client \
    --driver-memory 8g \
    --executor-cores $SLURM_CPUS_PER_TASK \
    --executor-memory 50g \
    --conf spark.sql.autoBroadcastJoinThreshold=1048576000 \
    match_prs.py

STATUS=$?

# ========================================
# DETENER CLUSTER
# ========================================
echo ">> Deteniendo cluster..."
for NODE in "${NODES[@]}"; do
    ssh "$NODE" "cd '$SPARK_HOME' && ./sbin/stop-worker.sh; pkill -f 'org.apache.spark.deploy.worker'" 2>/dev/null &
done
ssh "$MASTER_HOST" "cd '$SPARK_HOME' && ./sbin/stop-master.sh; pkill -f 'org.apache.spark.deploy.master'" 2>/dev/null
wait
pkill -P $$ 2>/dev/null

echo "Fin: $(date)"
exit $STATUS
