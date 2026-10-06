#!/bin/bash
#SBATCH --job-name=match_espiro
#SBATCH --account=researchers
#SBATCH --qos=vip
#SBATCH --nodes=4
#SBATCH --ntasks=4
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=48G
#SBATCH --time=04:00:00
#SBATCH --output=/mnt/cephfs/orgs/home/angel.pacheco/espirometrias/match_%j.out
#SBATCH --error=/mnt/cephfs/orgs/home/angel.pacheco/espirometrias/match_%j.err
#SBATCH -D /mnt/cephfs/orgs/home/angel.pacheco/espirometrias

echo "========================================"
echo "MATCHING DE ESPIROMETRÍAS - JOB SLURM"
echo "========================================"

WORKDIR="/mnt/cephfs/orgs/home/angel.pacheco/espirometrias"
cd "$WORKDIR" || exit 1

echo "Directorio de trabajo: $WORKDIR"
echo ""

# ❌ NADA de module load python/3.10 ni gcc/11.2.0
# module purge
# module load python/3.10
# module load gcc/11.2.0

# Activar el entorno correcto (.venv)
if [ ! -d ".venv" ]; then
    echo "ERROR: No existe .venv. Créalo en el login con 'uv venv .venv' e instala pandas/numpy."
    exit 1
fi

echo "Activando entorno .venv..."
source .venv/bin/activate
echo "Python usado:"
which python
python --version
echo ""

echo "========================================"
echo "INICIANDO PROCESAMIENTO DISTRIBUIDO"
echo "========================================"


# Número total de tasks
TOTAL_TASKS=4

# Lanzar tasks en paralelo usando srun
# Cada task se ejecuta en un nodo diferente con sus 24 CPUs
for task_id in $(seq 0 $((TOTAL_TASKS-1))); do
    echo "Lanzando Task $((task_id+1))/$TOTAL_TASKS..."

    srun --nodes=1 --ntasks=1 --cpus-per-task=$SLURM_CPUS_PER_TASK \
         --exclusive \
         python matching_distribuido.py \
         --task-id=$task_id \
         --total-tasks=$TOTAL_TASKS \
         --cpus=$SLURM_CPUS_PER_TASK \
         &
done

echo ""
echo "Todas las tasks lanzadas. Esperando a que terminen..."
echo ""

# Esperar a que todas las tasks terminen
wait

echo ""
echo "========================================"
echo "TODAS LAS TASKS COMPLETADAS"
echo "========================================"
echo ""

# Verificar que se generaron todos los archivos temporales
TEMP_DIR="$WORKDIR/temp_chunks"
echo "Verificando archivos temporales en: $TEMP_DIR"

if [ -d "$TEMP_DIR" ]; then
    num_files=$(ls -1 $TEMP_DIR/match_task_*.csv 2>/dev/null | wc -l)
    echo "Archivos temporales encontrados: $num_files"

    if [ $num_files -eq $TOTAL_TASKS ]; then
        echo "✓ Todos los archivos temporales generados correctamente"
    else
        echo "⚠ ADVERTENCIA: Se esperaban $TOTAL_TASKS archivos, se encontraron $num_files"
    fi
else
    echo "⚠ ADVERTENCIA: Directorio temporal no encontrado"
fi

echo ""

# Combinar resultados (esto ya lo hace la última task automáticamente)
# Pero verificamos que el archivo final exista
FINAL_FILE="$WORKDIR/matching_4fuentes_completo.csv"

if [ -f "$FINAL_FILE" ]; then
    echo "✓ Archivo final generado: $FINAL_FILE"

    # Estadísticas del archivo
    num_lines=$(wc -l < "$FINAL_FILE")
    file_size=$(du -h "$FINAL_FILE" | cut -f1)

    echo ""
    echo "Estadísticas del archivo final:"
    echo "  - Líneas: $num_lines"
    echo "  - Tamaño: $file_size"
    echo ""

    # Mostrar primeras líneas del header
    echo "Primeras columnas del archivo:"
    head -n 1 "$FINAL_FILE" | tr ',' '\n' | head -n 20
    echo "  ..."

else
    echo "⚠ ERROR: No se encontró el archivo final"
    echo "Intentando combinar manualmente..."

    # Combinar manualmente si es necesario
    python -c "
import pandas as pd
from pathlib import Path
import os

temp_dir = Path('$TEMP_DIR')
files = sorted(temp_dir.glob('match_task_*.csv'))

if files:
    print(f'Combinando {len(files)} archivos...')
    dfs = [pd.read_csv(f) for f in files]
    df_final = pd.concat(dfs, ignore_index=True)
    df_final.to_csv('$FINAL_FILE', index=False)
    print(f'✓ Archivo final creado: {len(df_final):,} registros')
else:
    print('ERROR: No se encontraron archivos temporales')
"
fi

# Limpiar archivos temporales
if [ -d "$TEMP_DIR" ]; then
    echo ""
    echo "Limpiando archivos temporales..."
    rm -rf "$TEMP_DIR"
    echo "✓ Archivos temporales eliminados"
fi

echo ""
echo "========================================"
echo "JOB COMPLETADO CON ÉXITO"
echo "========================================"
echo ""
echo "Job ID: $SLURM_JOB_ID"
echo "Tiempo de ejecución:"
sacct -j $SLURM_JOB_ID --format=JobID,JobName,Elapsed,State,MaxRSS

echo ""
echo "Para ver el archivo de resultados:"
echo "  less $FINAL_FILE"
echo ""
echo "Para ver estadísticas:"
echo "  python -c \"import pandas as pd; df=pd.read_csv('$FINAL_FILE'); print(df['MATCH'].value_counts()); print(df['match_confianza'].value_counts())\""
echo ""

exit 0
