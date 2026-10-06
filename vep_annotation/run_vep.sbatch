#!/bin/bash
#SBATCH --job-name=vep_annot
#SBATCH --account=researchers
#SBATCH --qos=vip
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=/mnt/cephfs/orgs/home/angel.pacheco/vep_annotation/results/vep_run_%j.out
#SBATCH --error=/mnt/cephfs/orgs/home/angel.pacheco/vep_annotation/results/vep_run_%j.err
#SBATCH --mail-user=angelpachecosaavedra@gmail.com
#SBATCH --mail-type=BEGIN,END,FAIL

set -euo pipefail

# ================================================================
# VEP 115 con Plugins (REVEL, CADD, SpliceAI, PrimateAI)
#
# Input: Archivo VCF
# Output: Anotación en TSV
# ================================================================

VEP_SIF="/mnt/cephfs/orgs/home/angel.pacheco/vep_annotation/ensembl-vep-115.sif"
WORK="/mnt/cephfs/orgs/home/angel.pacheco/vep_annotation"
PLUGINS_DIR="${WORK}/plugins"

# Variables de I/O
INPUT_VCF="${WORK}/input/test_variants.vcf"
RESULTS_DIR="${WORK}/results"
OUTPUT_TSV="${RESULTS_DIR}/test_4plugins.tsv"

mkdir -p "$RESULTS_DIR"

# ================================================================
# Manejo de temporales
# ================================================================
TMPDIR="/mnt/cephfs/scratch/$USER/vep_tmp/$SLURM_JOB_ID"
mkdir -p "$TMPDIR"; chmod 700 "$TMPDIR"
export TMPDIR TMP="$TMPDIR" TEMP="$TMPDIR"

echo "============================================================"
echo "  PIPELINE VEP — Anotación con Plugins"
echo "  Nodo: $SLURMD_NODENAME"
echo "  RAM: $(free -g | awk 'NR==2{print $2}')G"
echo "  CPUs: $SLURM_CPUS_PER_TASK"
echo "  $(date)"
echo "============================================================"

# ================================================================
# Pre-flight checks
# ================================================================
echo ""
echo "--- Pre-flight ---"

if [[ ! -f "$INPUT_VCF" ]]; then
  echo "ERROR: No existe el input $INPUT_VCF"
  exit 1
fi

if [[ ! -f "$VEP_SIF" ]]; then
  echo "ERROR: No existe el contenedor VEP en $VEP_SIF"
  exit 1
fi

echo "   Input VCF     : $INPUT_VCF"
echo "   VEP SIF       : $VEP_SIF"
echo "   Plugins dir   : $PLUGINS_DIR"
echo "   Output TSV    : $OUTPUT_TSV"

# ================================================================
# PASO 1: Ejecución de VEP 115
# ================================================================
echo ""
echo "============================================================"
echo "  Ejecutando Ensembl VEP (Modo Offline)"
echo "============================================================"
echo "  Inicio: $(date)"

apptainer exec \
    --bind /mnt/cephfs:/mnt/cephfs \
    --bind "$TMPDIR:/tmp,$TMPDIR:/var/tmp" \
    --env VEP_DIR_PLUGINS="$PLUGINS_DIR" \
    "$VEP_SIF" \
    vep \
        --input_file "$INPUT_VCF" \
        --output_file "$OUTPUT_TSV" \
        --cache --offline --assembly GRCh38 \
        --dir_cache /mnt/cephfs/hot_nvme/vep \
        --dir_plugins "$PLUGINS_DIR" \
        --fasta /mnt/cephfs/hot_nvme/vep/homo_sapiens/115_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz \
        --plugin REVEL,"${PLUGINS_DIR}/revel_with_transcript_ids.GRCh38.header.tsv.gz" \
        --plugin CADD,snv=/mnt/cephfs/hot_nvme/cadd/whole_genome_SNVs.sorted.tsv.bgz,indels=/mnt/cephfs/hot_nvme/cadd/gnomad.genomes.r4.0.indel.tsv.gz \
        --plugin SpliceAI,snv=/mnt/cephfs/hot_nvme/genome_scores_v1.3/spliceai_scores.raw.snv.hg38.vcf.gz,indel=/mnt/cephfs/hot_nvme/genome_scores_v1.3/spliceai_scores.raw.indel.hg38.vcf.gz,cutoff=0.5 \
        --plugin PrimateAI,/mnt/cephfs/hot_nvme/primateai/PrimateAI_scores_v0.2_GRCh38_sorted.tsv.bgz \
        --tab --force_overwrite \
        --sift b --polyphen b \
        --symbol --canonical --numbers \
        --af --af_gnomade \
        --fork $SLURM_CPUS_PER_TASK \
        --verbose

echo "  Fin: $(date)"

# Limpieza segura de temporales
rm -rf "$TMPDIR"

echo ""
echo "============================================================"
echo "  PIPELINE VEP COMPLETADO — $(date)"
echo ""
echo "  Outputs en $RESULTS_DIR/:"
echo "    $(basename "$OUTPUT_TSV")"
echo "============================================================"
