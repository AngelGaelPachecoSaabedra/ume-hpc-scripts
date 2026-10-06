#!/bin/bash
#SBATCH --job-name=vcf2zarr_chr10
#SBATCH --account=researchers
#SBATCH --qos=large
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=150G
#SBATCH --time=12:00:00
#SBATCH --output=/mnt/cephfs/orgs/home/angel.pacheco/vcf2zarr/logs/vcf2zarr_chr10_%j.out
#SBATCH --error=/mnt/cephfs/orgs/home/angel.pacheco/vcf2zarr/logs/vcf2zarr_chr10_%j.err
#SBATCH -D /mnt/cephfs/orgs/home/angel.pacheco/vcf2zarr

echo "======================================"
echo " VCF → Zarr conversión cromosoma 10"
echo " Job ID:     $SLURM_JOB_ID"
echo " Nodo:       $SLURM_NODELIST"
echo " CPUs:       $SLURM_CPUS_PER_TASK"
echo " Memoria:    $SLURM_MEM_PER_NODE MB"
echo " Inicio:     $(date)"
echo "======================================"
echo ""

# Variables
VENV_PATH=/mnt/cephfs/orgs/home/angel.pacheco/uv-envs/vcf2zarr
VCF_FILE="/mnt/cephfs/hot_nvme/mcps/imputed-topmed/vcf_files/mcps-freeze150k_qcd_chr10.vcf.gz"
ZARR_OUTPUT="/mnt/cephfs/orgs/home/angel.pacheco/vcf2zarr/zarr_output/mcps-freeze150k_qcd_chr10.zarr"

# Verificar que el venv existe
if [ ! -d "$VENV_PATH" ]; then
    echo "ERROR: Virtual environment no existe: $VENV_PATH"
    echo ""
    echo "Créalo con:"
    echo "  python3 -m venv $VENV_PATH"
    echo "  source $VENV_PATH/bin/activate"
    echo "  pip install scikit-allel zarr numpy numcodecs"
    exit 1
fi

# Activar venv
echo "Activando ambiente: $VENV_PATH"
source ${VENV_PATH}/bin/activate

# Verificar Python
echo "Python: $(which python)"
echo "Versión: $(python --version)"
echo ""

# Verificar paquetes
echo "Verificando dependencias..."
python -c "import allel; print('  ✓ scikit-allel:', allel.__version__)" || {
    echo "  ✗ ERROR: scikit-allel no disponible"
    exit 1
}
python -c "import zarr; print('  ✓ zarr:', zarr.__version__)" || exit 1
python -c "import numpy; print('  ✓ numpy:', numpy.__version__)" || exit 1
python -c "import numcodecs; print('  ✓ numcodecs:', numcodecs.__version__)" || exit 1
echo ""

# Verificar VCF
if [ ! -f "$VCF_FILE" ]; then
    echo "ERROR: VCF no encontrado: $VCF_FILE"
    exit 1
fi

echo "Archivos:"
echo "  Input:  $VCF_FILE"
echo "  Output: $ZARR_OUTPUT"
echo ""

# Crear directorio de salida
mkdir -p "$(dirname $ZARR_OUTPUT)"

# Conversión
echo "Iniciando conversión..."
echo "======================================"

python << 'EOF'
import allel
import sys
from datetime import datetime

vcf_file = "/mnt/cephfs/hot_nvme/mcps/imputed-topmed/vcf_files/mcps-freeze150k_qcd_chr10.vcf.gz"
zarr_output = "/mnt/cephfs/orgs/home/angel.pacheco/vcf2zarr/zarr_output/mcps-freeze150k_qcd_chr10.zarr"

print(f"Timestamp inicio: {datetime.now()}")
print(f"VCF: {vcf_file}")
print(f"Zarr: {zarr_output}")
print()

try:
    allel.vcf_to_zarr(
        vcf_file,
        zarr_output,
        fields="*",
        alt_number=1,
        overwrite=True,
        log=sys.stdout,
        chunk_length=10000,
        chunk_width=1000
    )
    
    print()
    print("=" * 70)
    print("✓ CONVERSIÓN EXITOSA")
    print("=" * 70)
    
    # Información del resultado
    import zarr
    z = zarr.open(zarr_output, mode='r')
    
    print()
    print("Información del Zarr:")
    print(f"  Ubicación: {zarr_output}")
    
    if 'variants/POS' in z:
        n_variants = z['variants/POS'].shape[0]
        print(f"  Variantes: {n_variants:,}")
    
    if 'samples' in z:
        n_samples = len(z['samples'])
        print(f"  Muestras: {n_samples:,}")
    
    print(f"  Grupos: {list(z.group_keys())}")
    
    # Tamaño
    import subprocess
    try:
        size_output = subprocess.check_output(['du', '-sh', zarr_output]).decode()
        size = size_output.split()[0]
        print(f"  Tamaño: {size}")
    except:
        pass
    
    print()
    print(f"Timestamp fin: {datetime.now()}")
    print("=" * 70)
    
except Exception as e:
    print()
    print("=" * 70)
    print(f"✗ ERROR: {e}")
    print("=" * 70)
    import traceback
    traceback.print_exc()
    sys.exit(1)
EOF

STATUS=$?

echo ""
echo "======================================"
echo " Finalizado: $(date)"
echo " Exit code: $STATUS"
echo "======================================"

deactivate
exit $STATUS