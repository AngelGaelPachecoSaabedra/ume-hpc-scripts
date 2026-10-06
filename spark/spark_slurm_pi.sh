#!/bin/bash
#SBATCH --job-name=spark_pi
#SBATCH --account=researchers
#SBATCH --qos=vip
#SBATCH --partition=normal
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/mnt/cephfs/orgs/home/angel.pacheco/spark/spark_pi_%j.out
#SBATCH --error=/mnt/cephfs/orgs/home/angel.pacheco/spark/spark_pi_%j.err
#SBATCH -D /mnt/cephfs/orgs/home/angel.pacheco/spark

echo "=============================================="
echo "SPARK DISTRIBUIDO - $SLURM_JOB_NUM_NODES nodos"
echo "Job ID:        $SLURM_JOB_ID"
echo "Nodos:         $SLURM_NODELIST"
echo "Inicio:        $(date)"
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

echo ""
echo "Master: $MASTER_HOST ($MASTER_IP)"
echo ""

# ========================================
# INICIAR MASTER EN BACKGROUND PERSISTENTE
# ========================================
echo ">> Iniciando Master..."

srun -N1 -w "$MASTER_HOST" --job-name=spark-master bash -c "
    export SPARK_HOME='$SPARK_HOME'
    export JAVA_HOME='$JAVA_HOME'
    export PATH='$PATH'
    export SPARK_LOG_DIR='$SPARK_LOG_DIR'
    
    cd \$SPARK_HOME
    
    # Iniciar master y mantenerlo vivo
    nohup ./sbin/start-master.sh --host $MASTER_IP --port 7077 >/dev/null 2>&1
    
    # Esperar a que esté listo
    for i in {1..30}; do
        if ss -tlnp 2>/dev/null | grep -q ':7077'; then
            echo '✓ Master listo'
            # MANTENER PROCESO VIVO
            while ps aux | grep -q '[o]rg.apache.spark.deploy.master.Master'; do
                sleep 5
            done
            exit 0
        fi
        sleep 2
    done
    echo '✗ Master falló'
    exit 1
" &

MASTER_SRUN_PID=$!

# Esperar confirmación de que master está listo
sleep 10

# Verificar que master está corriendo
ssh "$MASTER_HOST" "ss -tlnp | grep ':7077'" > /dev/null
if [ $? -ne 0 ]; then
    echo "ERROR: Master no está escuchando en 7077"
    kill $MASTER_SRUN_PID 2>/dev/null
    exit 1
fi

echo "✓ Master confirmado en $MASTER_IP:7077"
echo ""

# ========================================
# INICIAR WORKERS EN CADA NODO
# ========================================
echo ">> Iniciando Workers..."

NODES=($(scontrol show hostnames "$SLURM_NODELIST"))

for NODE in "${NODES[@]}"; do
    echo "   - Worker en $NODE"
    
    srun -N1 -w "$NODE" --job-name=spark-worker bash -c "
        export SPARK_HOME='$SPARK_HOME'
        export JAVA_HOME='$JAVA_HOME'
        export PATH='$PATH'
        export SPARK_LOG_DIR='$SPARK_LOG_DIR'
        export SPARK_WORKER_DIR='$SPARK_WORKER_DIR'
        
        HOSTNAME=\$(hostname)
        echo \"Starting worker on \$HOSTNAME\"
        
        cd \$SPARK_HOME
        
        # Iniciar worker
        nohup ./sbin/start-worker.sh '$MASTER_URL' \
            --host \$HOSTNAME \
            --cores $SLURM_CPUS_PER_TASK \
            --memory 14G >/dev/null 2>&1
        
        # MANTENER PROCESO VIVO
        sleep 5
        while ps aux | grep -q '[o]rg.apache.spark.deploy.worker.Worker'; do
            sleep 5
        done
    " &
done

sleep 20

# ========================================
# VERIFICAR CLUSTER
# ========================================
echo ""
echo ">> Estado del cluster:"

ssh "$MASTER_HOST" "ps aux | grep -E 'Master|Worker' | grep -v grep | wc -l | xargs echo '  Procesos Spark en master:'"
ssh "$MASTER_HOST" "ss -tlnp | grep ':7077' && echo '  ✓ Puerto 7077 activo' || echo '  ✗ Puerto 7077 no activo'"

echo ""
echo "  Workers por nodo:"
for NODE in "${NODES[@]}"; do
    COUNT=$(ssh "$NODE" "ps aux | grep 'org.apache.spark.deploy.worker.Worker' | grep -v grep | wc -l")
    echo "    $NODE: $COUNT"
done

echo ""

# ========================================
# EJECUTAR APLICACIÓN
# ========================================
echo ">> Ejecutando SparkPi..."
echo ""

$SPARK_HOME/bin/spark-submit \
    --master "$MASTER_URL" \
    --deploy-mode client \
    --driver-memory 2g \
    --executor-cores 4 \
    --executor-memory 8g \
    --total-executor-cores 16 \
    --conf spark.dynamicAllocation.enabled=false \
    --conf spark.network.timeout=800s \
    --class org.apache.spark.examples.SparkPi \
    $SPARK_HOME/examples/jars/spark-examples*.jar \
    1000

STATUS=$?

echo ""
echo ">> Resultado: $STATUS"
echo ""

# ========================================
# DETENER TODO
# ========================================
echo ">> Deteniendo cluster..."

for NODE in "${NODES[@]}"; do
    ssh "$NODE" "cd '$SPARK_HOME' && ./sbin/stop-worker.sh; pkill -f 'org.apache.spark.deploy.worker'" 2>/dev/null &
done

ssh "$MASTER_HOST" "cd '$SPARK_HOME' && ./sbin/stop-master.sh; pkill -f 'org.apache.spark.deploy.master'" 2>/dev/null

wait

# Matar los srun que mantienen procesos vivos
pkill -P $$ 2>/dev/null

echo "✓ Cluster detenido"
echo ""
echo "=============================================="
echo "Fin: $(date)"
echo "=============================================="

exit $STATUS