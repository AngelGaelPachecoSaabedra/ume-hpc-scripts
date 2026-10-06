#!/bin/bash
# =============================================================================
# run_prs_spark_gpu.sbatch
# =============================================================================
# PRS Spark+GPU: matching de variantes con PySpark local + scoring CuPy GPU.
#
# Arquitectura de integración Spark ↔ Apptainer:
# ──────────────────────────────────────────────
#   • Spark: usa la instalación real en CephFS (sin pip, sin runtime alterno)
#       SPARK_HOME=/mnt/cephfs/biocontainers/spark/spark
#   • PySpark: importado dentro del contenedor via PYTHONPATH → $SPARK_HOME/python
#   • Java:    expuesto al contenedor vía --bind /usr/lib/jvm:/usr/lib/jvm
#   • Modo:    local[N] (N = CPUs Slurm) — nodo único, sin levantar master/workers
#   • Python:  /opt/miniforge3/envs/jax-gpu/bin/python (tiene CuPy + zarr)
#              ejecutado DENTRO del contenedor 1.5.sif
#
# Uso:
#   sbatch --export=PGS_ID=PGS000001,CHROM=1 run_prs_spark_gpu.sbatch
#   sbatch --export=PGS_ID=PGS000001,CHROM=X,WINDOW_SIZE=20000 run_prs_spark_gpu.sbatch
#   sbatch --export=PGS_ID=PGS000363,CHROM=1,GPU_BATCH_MATCHES=16000 run_prs_spark_gpu.sbatch
#
# Variables obligatorias:
#   PGS_ID        PGS Catalog ID (ej. PGS000001)
#   CHROM         Cromosoma: 1-22 o X
#
# Variables opcionales (con defaults):
#   WINDOW_SIZE        Variantes zarr por ventana      (default: 10000)
#   MISSING_STRATEGY   mean | zero | skip              (default: mean)
#   DRIVER_MEM         Memoria driver Spark            (default: 60g)
#   GPU_BATCH_MATCHES  Variantes matched por bloque GPU (default: 4000)
# =============================================================================

#SBATCH --job-name=prs_spark_gpu
#SBATCH --account=researchers
#SBATCH --qos=vip
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --gres=gpu:1
#SBATCH --time=06:00:00
#SBATCH --output=/mnt/cephfs/orgs/home/angel.pacheco/logs/prs_spark_gpu_%j.out
#SBATCH --error=/mnt/cephfs/orgs/home/angel.pacheco/logs/prs_spark_gpu_%j.err

set -euo pipefail

# =============================================================================
# Validación de variables de entrada
# =============================================================================
: "${PGS_ID:?[ERROR] PGS_ID no definido. Uso: sbatch --export=PGS_ID=PGS000001,CHROM=1 ...}"
: "${CHROM:?[ERROR]  CHROM no definido. Uso: CHROM=1 o CHROM=X}"

# =============================================================================
# Rutas — contenedor y script
# =============================================================================
SIF="/mnt/cephfs/biocontainers/images/jupyter-biotools-1.5.sif"
SCRIPT_DIR="/mnt/cephfs/orgs/home/angel.pacheco/Ordanamiento_Bases/scripts"

# =============================================================================
# Spark — instalación real en CephFS
# Nota: /mnt/cephfs/biocontainers/spark/spark puede ser un symlink a la
# versión concreta spark-3.5.3-bin-hadoop3. Ambas rutas son válidas;
# esta es la que usan los scripts existentes (match_cadd_spark, 2_match_prs_spark).
# Si se quiere apuntar directamente a la versión:
#   SPARK_HOME="/mnt/cephfs/biocontainers/spark/spark/spark-3.5.3-bin-hadoop3"
# =============================================================================
SPARK_HOME="/mnt/cephfs/biocontainers/spark/spark"

# Java del host — misma ruta que usan todos tus scripts Spark actuales.
# Se expone al contenedor vía --bind /usr/lib/jvm:/usr/lib/jvm (read-only).
JAVA_HOME="/usr/lib/jvm/java-17-openjdk-amd64"

# Python dentro del contenedor que tiene CuPy + zarr + numpy
# (entorno jax-gpu del contenedor jupyter-biotools-1.5.sif)
CUPY_PYTHON="/opt/miniforge3/envs/jax-gpu/bin/python"

# =============================================================================
# Parámetros de ejecución con defaults
# =============================================================================
WINDOW_SIZE="${WINDOW_SIZE:-10000}"
MISSING_STRATEGY="${MISSING_STRATEGY:-mean}"
DRIVER_MEM="${DRIVER_MEM:-60g}"
GPU_BATCH_MATCHES="${GPU_BATCH_MATCHES:-4000}"

# Spark local mode: usa todos los CPUs asignados por Slurm
SPARK_MASTER="local[${SLURM_CPUS_PER_TASK}]"

# =============================================================================
# Scratch aislado por job (evita colisiones entre jobs simultáneos)
# =============================================================================
BASE_TMP="/mnt/cephfs/scratch/${USER}/prs_spark_gpu"
TMPDIR="${BASE_TMP}/${SLURM_JOB_ID}"
mkdir -p "${TMPDIR}"
chmod 700 "${TMPDIR}"

# Directorios internos de Spark — todos en scratch para evitar I/O en /tmp
SPARK_LOCAL_DIRS="${TMPDIR}/spark-local"
SPARK_LOG_DIR="${TMPDIR}/spark-logs"
SPARK_WORKER_DIR="${TMPDIR}/spark-worker"
mkdir -p "${SPARK_LOCAL_DIRS}" "${SPARK_LOG_DIR}" "${SPARK_WORKER_DIR}"

# Apptainer: sus propios temporales también en scratch
APPTAINER_TMP="${TMPDIR}/apptainer-tmp"
APPTAINER_CACHEDIR="${TMPDIR}/apptainer-cache"
mkdir -p "${APPTAINER_TMP}" "${APPTAINER_CACHEDIR}"

export TMPDIR TMP="${TMPDIR}" TEMP="${TMPDIR}"
export APPTAINER_TMPDIR="${APPTAINER_TMP}"
export APPTAINER_CACHEDIR

mkdir -p "/mnt/cephfs/orgs/home/angel.pacheco/logs"

# =============================================================================
# Cabecera informativa
# =============================================================================
echo "================================================================"
echo "  PRS Spark+GPU Pipeline"
echo "  PGS_ID:      ${PGS_ID}"
echo "  CHROM:       ${CHROM}"
echo "  WINDOW:           ${WINDOW_SIZE}"
echo "  MISSING:          ${MISSING_STRATEGY}"
echo "  DRIVER_MEM:       ${DRIVER_MEM}"
echo "  GPU_BATCH_MATCHES:${GPU_BATCH_MATCHES}"
echo "  SPARK_MASTER:${SPARK_MASTER}"
echo "----------------------------------------------------------------"
echo "  Job ID:      ${SLURM_JOB_ID}"
echo "  Nodo:        ${SLURMD_NODENAME}"
echo "  CPUs:        ${SLURM_CPUS_PER_TASK}"
echo "  RAM:         ${SLURM_MEM_PER_NODE:-??} MB"
echo "  Scratch:     ${TMPDIR}"
echo "----------------------------------------------------------------"
echo "  SPARK_HOME:  ${SPARK_HOME}"
echo "  JAVA_HOME:   ${JAVA_HOME}"
echo "  Python:      ${CUPY_PYTHON}"
echo "  Contenedor:  ${SIF}"
echo "================================================================"
echo ""

# =============================================================================
# Ejecución dentro del contenedor
# =============================================================================
# --nv                → inyecta librerías NVIDIA del host (necesario para CuPy)
# --bind /mnt/cephfs  → acceso a zarr (solo lectura) + Spark instalado en CephFS
# --bind /usr/lib/jvm → expone Java del host al contenedor (necesario para Spark JVM)
#                       /usr/lib/jvm es una ruta estándar del OS, no datos de usuario
# bash -lc "…"        → login shell para cargar variables del entorno del contenedor
#
# Dentro del contenedor:
#   1. Se construye PYTHONPATH con $SPARK_HOME/python y py4j (sin pip install)
#   2. Se configura JAVA_HOME y PATH para que el JVM sea accesible
#   3. Se configura LD_LIBRARY_PATH con librerías CUDA/cuBLAS para CuPy
#   4. Se lanza compute_prs_spark_gpu.py con el Python del entorno jax-gpu
# =============================================================================
apptainer exec \
    --nv \
    --bind /mnt/cephfs:/mnt/cephfs \
    --bind /usr/lib/jvm:/usr/lib/jvm \
    --bind "${TMPDIR}:/tmp,${TMPDIR}:/var/tmp" \
    "${SIF}" \
    bash -lc "

    # ── PySpark desde instalación real en CephFS (sin pip) ───────────────────
    # SPARK_HOME está montado en el contenedor vía --bind /mnt/cephfs
    export SPARK_HOME='${SPARK_HOME}'

    # py4j: buscar el zip exacto (Spark 3.5.x usa py4j-0.10.9.7)
    PY4J_ZIP=\$(ls \"\${SPARK_HOME}/python/lib/\"/py4j-*.zip 2>/dev/null | head -1)
    if [ -z \"\${PY4J_ZIP}\" ]; then
        echo '[ERROR] py4j zip no encontrado en \${SPARK_HOME}/python/lib/' >&2
        exit 1
    fi

    # PYTHONPATH: pyspark desde SPARK_HOME + py4j
    export PYTHONPATH=\"\${SPARK_HOME}/python:\${PY4J_ZIP}\${PYTHONPATH:+:\${PYTHONPATH}}\"

    # ── Java del host (expuesto vía --bind /usr/lib/jvm) ────────────────────
    export JAVA_HOME='${JAVA_HOME}'
    export PATH=\"\${JAVA_HOME}/bin:\${SPARK_HOME}/bin:\${PATH}\"

    # ── Spark: directorios de trabajo en scratch ──────────────────────────────
    export SPARK_LOCAL_DIRS='${SPARK_LOCAL_DIRS}'
    export SPARK_LOG_DIR='${SPARK_LOG_DIR}'
    export SPARK_WORKER_DIR='${SPARK_WORKER_DIR}'

    # PYSPARK_PYTHON: el Python dentro del contenedor que usará Spark para
    # ejecutar el código del usuario (driver + workers en modo local)
    export PYSPARK_PYTHON='${CUPY_PYTHON}'
    export PYSPARK_DRIVER_PYTHON='${CUPY_PYTHON}'

    # ── CUDA / cuBLAS para CuPy ──────────────────────────────────────────────
    JAX_LIBS=\$(find /opt/miniforge3/envs/jax-gpu/lib/python3.10/site-packages/nvidia/ \
        -type d -name lib 2>/dev/null | tr '\n' ':')
    export LD_LIBRARY_PATH=\"/.singularity.d/libs:/opt/miniforge3/envs/jax-gpu/lib:\${JAX_LIBS}\${LD_LIBRARY_PATH:-}\"

    # ── Temporales ───────────────────────────────────────────────────────────
    export TMPDIR='${TMPDIR}'
    export TMP='${TMPDIR}'
    export TEMP='${TMPDIR}'

    # ── Diagnóstico pre-ejecución ────────────────────────────────────────────
    echo '[INFO] Java:'
    java -version 2>&1 | head -1

    echo '[INFO] PySpark disponible:'
    '${CUPY_PYTHON}' -c 'import pyspark; print(\"  pyspark\", pyspark.__version__)'

    echo '[INFO] CuPy disponible:'
    '${CUPY_PYTHON}' -c 'import cupy; print(\"  cupy\", cupy.__version__)'

    echo '[INFO] Zarr disponible:'
    '${CUPY_PYTHON}' -c 'import zarr; print(\"  zarr\", zarr.__version__)'

    echo '[INFO] py4j: '\${PY4J_ZIP}
    echo '[INFO] PYTHONPATH (primeros 200): '\${PYTHONPATH:0:200}
    echo ''

    # ── Pipeline ─────────────────────────────────────────────────────────────
    '${CUPY_PYTHON}' '${SCRIPT_DIR}/compute_prs_spark_gpu.py' \
        --pgs-id                '${PGS_ID}'              \
        --chrom                 '${CHROM}'                \
        --window-size           '${WINDOW_SIZE}'          \
        --missing-strategy      '${MISSING_STRATEGY}'     \
        --spark-master          '${SPARK_MASTER}'         \
        --spark-driver-memory   '${DRIVER_MEM}'           \
        --spark-executor-memory '16g'                     \
        --gpu-batch-matches     '${GPU_BATCH_MATCHES}'
    "

STATUS=$?

# =============================================================================
# Limpieza del scratch
# =============================================================================
rm -rf "${TMPDIR}"

# =============================================================================
# Resultado final
# =============================================================================
echo ""
if [ "${STATUS}" -eq 0 ]; then
    echo "================================================================"
    echo "  ÉXITO: PRS chr${CHROM} (${PGS_ID}) completado."
    echo "  Salida en: /mnt/cephfs/hot_nvme/pgscatalog/scores/${PGS_ID}/"
    echo "    ${PGS_ID}_chr${CHROM}_scores.tsv"
    echo "    ${PGS_ID}_chr${CHROM}_metadata.json"
    echo "================================================================"
else
    echo "================================================================" >&2
    echo "  ERROR: PRS chr${CHROM} (${PGS_ID}) — código ${STATUS}"        >&2
    echo "  Log: /mnt/cephfs/orgs/home/angel.pacheco/logs/prs_spark_gpu_${SLURM_JOB_ID}.out" >&2
    echo "================================================================" >&2
fi

exit "${STATUS}"
