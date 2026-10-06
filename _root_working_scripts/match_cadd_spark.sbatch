#!/bin/bash
#SBATCH --job-name=spark_cadd
#SBATCH --account=researchers
#SBATCH --qos=vip
#SBATCH --partition=normal
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=64G              # Aumentamos la RAM porque CADD es masivo
#SBATCH --time=04:00:00
#SBATCH --output=/mnt/cephfs/orgs/home/angel.pacheco/spark/spark_cadd_%j.out
#SBATCH --error=/mnt/cephfs/orgs/home/angel.pacheco/spark/spark_cadd_%j.err
#SBATCH -D /mnt/cephfs/orgs/home/angel.pacheco/spark

echo "=============================================="
echo "SPARK DISTRIBUIDO - CROSS MATCH CADD"
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
# GENERAR SCRIPT PYSPARK AL VUELO
# ========================================
cat << 'EOF' > match_cadd.py
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
EOF

# ========================================
# EJECUTAR APLICACIÓN PYSPARK
# ========================================
echo ">> Ejecutando Cruce Masivo..."
$SPARK_HOME/bin/spark-submit \
    --master "$MASTER_URL" \
    --deploy-mode client \
    --driver-memory 8g \
    --executor-cores $SLURM_CPUS_PER_TASK \
    --executor-memory 50g \
    --conf spark.sql.autoBroadcastJoinThreshold=1048576000 \
    match_cadd.py

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
