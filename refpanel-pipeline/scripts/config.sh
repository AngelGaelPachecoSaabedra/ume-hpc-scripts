#!/bin/bash
# ============================================================
# config.sh — Configuración global (chr1-22 + chrX)
# ============================================================

# --- Containers ---
export BIO_SIF="/mnt/cephfs/biocontainers/images/jupyter-biotools-1.2.sif"
export IMP_SIF="/mnt/cephfs/biocontainers/images/imputation-tools-1.1.sif"
export BIND="--bind /mnt/cephfs:/mnt/cephfs"

# --- Input ---
export PVCF_DIR="/mnt/cephfs/hot_nvme/mcps/whole_genome_sequencing/pVCF/NF"

# --- Output final ---
export PANEL_DIR="/mnt/cephfs/panel_referencia"

# --- Scratch (temporales) ---
export SCRATCH="/mnt/cephfs/scratch/angel.pacheco/refpanel_build"

# --- Logs ---
export LOG_DIR="/mnt/cephfs/orgs/home/angel.pacheco/logs"

# --- Reference genome ---
export REF_FASTA="/mnt/cephfs/hot_nvme/hg38/ref/hg38.fa"

# --- Genetic maps ---
export GMAP_DIR="/mnt/cephfs/panel_referencia/maps"

# --- Sex file para chrX (SHAPEIT5 format: 1=male, 2=female) ---
export SEX_FILE="/mnt/cephfs/panel_referencia/wgs_sex.txt"

# --- QC ---
export MIN_GQ=20
export MIN_DP=10
export MAX_MISS=0.05
export HWE_THRESH="1e-10"

# --- Función: convertir SLURM_ARRAY_TASK_ID a nombre de cromosoma ---
# 1-22 → chr1-chr22, 23 → chrX
get_chr_name() {
    local task_id=$1
    if [ "$task_id" -eq 23 ]; then
        echo "X"
    else
        echo "$task_id"
    fi
}

# --- Chromosome → pVCF chunks (incluyendo chrX) ---
declare -A CHR_CHUNKS
CHR_CHUNKS=(
    [1]="MCPS_WGS_Freeze_Two.01_chr1_1_123400000.vcf.gz|MCPS_WGS_Freeze_Two.02_chr1_123400001_248956422.vcf.gz"
    [2]="MCPS_WGS_Freeze_Two.03_chr2_1_93900000.vcf.gz|MCPS_WGS_Freeze_Two.04_chr2_93900001_242193529.vcf.gz"
    [3]="MCPS_WGS_Freeze_Two.05_chr3_1_90900000.vcf.gz|MCPS_WGS_Freeze_Two.06_chr3_90900001_198295559.vcf.gz"
    [4]="MCPS_WGS_Freeze_Two.07_chr4_1_50000000.vcf.gz|MCPS_WGS_Freeze_Two.08_chr4_50000001_190214555.vcf.gz"
    [5]="MCPS_WGS_Freeze_Two.09_chr5_1_48800000.vcf.gz|MCPS_WGS_Freeze_Two.10_chr5_48800001_181538259.vcf.gz"
    [6]="MCPS_WGS_Freeze_Two.11_chr6_1_59800000.vcf.gz|MCPS_WGS_Freeze_Two.12_chr6_59800001_170805979.vcf.gz"
    [7]="MCPS_WGS_Freeze_Two.13_chr7_1_60100000.vcf.gz|MCPS_WGS_Freeze_Two.14_chr7_60100001_159345973.vcf.gz"
    [8]="MCPS_WGS_Freeze_Two.15_chr8_1_45200000.vcf.gz|MCPS_WGS_Freeze_Two.16_chr8_45200001_145138636.vcf.gz"
    [9]="MCPS_WGS_Freeze_Two.17_chr9_1_43000000.vcf.gz|MCPS_WGS_Freeze_Two.18_chr9_43000001_138394717.vcf.gz"
    [10]="MCPS_WGS_Freeze_Two.19_chr10_1_39800000.vcf.gz|MCPS_WGS_Freeze_Two.20_chr10_39800001_133797422.vcf.gz"
    [11]="MCPS_WGS_Freeze_Two.21_chr11_1_53400000.vcf.gz|MCPS_WGS_Freeze_Two.22_chr11_53400001_135086622.vcf.gz"
    [12]="MCPS_WGS_Freeze_Two.23_chr12_1_35500000.vcf.gz|MCPS_WGS_Freeze_Two.24_chr12_35500001_133275309.vcf.gz"
    [13]="MCPS_WGS_Freeze_Two.25_chr13_1_17700000.vcf.gz|MCPS_WGS_Freeze_Two.26_chr13_17700001_114364328.vcf.gz"
    [14]="MCPS_WGS_Freeze_Two.27_chr14_1_17200000.vcf.gz|MCPS_WGS_Freeze_Two.28_chr14_17200001_107043718.vcf.gz"
    [15]="MCPS_WGS_Freeze_Two.29_chr15_1_19000000.vcf.gz|MCPS_WGS_Freeze_Two.30_chr15_19000001_101991189.vcf.gz"
    [16]="MCPS_WGS_Freeze_Two.31_chr16_1_36800000.vcf.gz|MCPS_WGS_Freeze_Two.32_chr16_36800001_90338345.vcf.gz"
    [17]="MCPS_WGS_Freeze_Two.33_chr17_1_25100000.vcf.gz|MCPS_WGS_Freeze_Two.34_chr17_25100001_83257441.vcf.gz"
    [18]="MCPS_WGS_Freeze_Two.35_chr18_1_18500000.vcf.gz|MCPS_WGS_Freeze_Two.36_chr18_18500001_80373285.vcf.gz"
    [19]="MCPS_WGS_Freeze_Two.37_chr19_1_26200000.vcf.gz|MCPS_WGS_Freeze_Two.38_chr19_26200001_58617616.vcf.gz"
    [20]="MCPS_WGS_Freeze_Two.39_chr20_1_28100000.vcf.gz|MCPS_WGS_Freeze_Two.40_chr20_28100001_64444167.vcf.gz"
    [21]="MCPS_WGS_Freeze_Two.41_chr21_1_12000000.vcf.gz|MCPS_WGS_Freeze_Two.42_chr21_12000001_46709983.vcf.gz"
    [22]="MCPS_WGS_Freeze_Two.43_chr22_1_15000000.vcf.gz|MCPS_WGS_Freeze_Two.44_chr22_15000001_50818468.vcf.gz"
    [X]="MCPS_WGS_Freeze_Two.45_chrX_1_61000000.vcf.gz|MCPS_WGS_Freeze_Two.46_chrX_61000001_156040895.vcf.gz"
)
export CHR_CHUNKS
