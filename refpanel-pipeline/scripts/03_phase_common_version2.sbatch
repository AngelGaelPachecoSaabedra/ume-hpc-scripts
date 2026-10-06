#!/bin/bash
#SBATCH --job-name=refpanel_03_phase_common
#SBATCH --account=researchers
#SBATCH --qos=vip
#SBATCH --partition=normal
#SBATCH --array=1-23
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=150G
#SBATCH --time=UNLIMITED
#SBATCH --output=/mnt/cephfs/orgs/home/angel.pacheco/logs/refpanel_03_phase_common_%A_%a.out
#SBATCH --error=/mnt/cephfs/orgs/home/angel.pacheco/logs/refpanel_03_phase_common_%A_%a.err
#SBATCH --mail-user=angelpachecosaavedra@gmail.com
#SBATCH --mail-type=BEGIN,END,FAIL
set -euo pipefail
source /mnt/cephfs/orgs/home/angel.pacheco/refpanel-pipeline/scripts/config.sh

CHR=$(get_chr_name ${SLURM_ARRAY_TASK_ID})
IN_BCF="${SCRATCH}/02_filter/chr${CHR}.filtered.bcf"
OUTDIR="${SCRATCH}/03_phase"
OUT_BCF="${OUTDIR}/chr${CHR}.common.phased.bcf"
GMAP="${GMAP_DIR}/chr${CHR}.b38.gmap.gz"

echo "================================================================"
echo "  Paso 3: Phase Common — chr${CHR}"
echo "  Job: ${SLURM_JOB_ID}, Node: ${SLURMD_NODENAME}"
echo "  Start: $(date)"
echo "================================================================"

mkdir -p "${OUTDIR}"
[ -f "${IN_BCF}" ] || { echo "ERROR: No existe ${IN_BCF}"; exit 1; }

if [ -f "${OUT_BCF}" ] && [ -f "${OUT_BCF}.csi" ]; then
    echo "SKIP: ya existe"
    exit 0
fi

echo -n "  Variantes en input: "
apptainer exec ${BIND} "${BIO_SIF}" bcftools view -H "${IN_BCF}" | wc -l

# ====================================================================
if [ "$CHR" != "X" ]; then
    # ===== AUTOSOMAS: cromosoma completo =====
    T=${SLURM_CPUS_PER_TASK}
    echo ""
    echo "  shapeit5_phase_common --filter-maf 0.001 --region chr${CHR}"
    echo ""
    apptainer exec ${BIND} "${IMP_SIF}" \
        shapeit5_phase_common \
            --input "${IN_BCF}" \
            --filter-maf 0.001 \
            --map "${GMAP}" \
            --output "${OUT_BCF}" \
            --region "chr${CHR}" \
            --thread "${T}" \
            --log "${OUTDIR}/chr${CHR}.common.log"
    apptainer exec ${BIND} "${BIO_SIF}" bcftools index -f "${OUT_BCF}"

else
    # ===== chrX: chunks definidos en cM + ligate =====
    T=${SLURM_CPUS_PER_TASK}          # con chunks pequeños caben los 16 cores
    CHRLEN=156040895                  # longitud chrX GRCh38
    CHUNK_CM=30                       # tamaño de chunk en cM
    OVL_CM=8                          # solape en cM (2x ventana HMM de 4 cM)

    CHUNKDIR="${OUTDIR}/chrX_chunks"
    mkdir -p "${CHUNKDIR}"

    # --- lista de haploides (hombres): SEX_FILE = "IID\tSEX", 1=male ---
    HAPLOID_FILE="${OUTDIR}/chrX_haploid_males.txt"
    awk -F'\t' 'NR>1 && $2==1 {print $1}' "${SEX_FILE}" > "${HAPLOID_FILE}"
    echo "  chrX: $(wc -l < ${HAPLOID_FILE}) hombres (haploides)"

    # --- generar coordenadas de chunks: cM -> bp, desde el gmap ---
    # gmap (sin header): col1=pos(bp)  col2=chr  col3=cM
    COORDS="${CHUNKDIR}/chrX_chunks.coords.txt"
    zcat "${GMAP}" | awk -v CHUNK="${CHUNK_CM}" -v OVL="${OVL_CM}" -v CHRLEN="${CHRLEN}" '
        { pos[NR]=$1+0; cm[NR]=$3+0; n=NR }
        END {
            maxcm = cm[n]; step = CHUNK - OVL; first = 1
            for (s = 0; s < maxcm; s += step) {
                lo = s; hi = s + CHUNK; if (hi > maxcm) hi = maxcm
                bplo = pos[1]
                for (i = 1; i <= n; i++) if (cm[i] >= lo) { bplo = pos[i]; break }
                bphi = pos[n]
                for (i = n; i >= 1; i--) if (cm[i] <= hi) { bphi = pos[i]; break }
                if (first) { bplo = 1; first = 0 }
                if (hi >= maxcm) bphi = CHRLEN
                printf "%d-%d\n", bplo, bphi
                if (hi >= maxcm) break
            }
        }' > "${COORDS}"

    N_CHUNKS=$(wc -l < "${COORDS}")
    echo "  chrX: ${N_CHUNKS} chunks (${CHUNK_CM} cM / solape ${OVL_CM} cM)"
    echo "  --- coordenadas (bp) ---"
    cat -n "${COORDS}"
    echo "  ------------------------"

    # --- detectar binario de ligate dentro del contenedor ---
    LIGATE_BIN=$(apptainer exec ${BIND} "${IMP_SIF}" bash -lc \
        'for b in shapeit5_ligate ligate SHAPEIT5_ligate; do if command -v "$b" >/dev/null 2>&1; then command -v "$b"; break; fi; done' || true)
    [ -n "${LIGATE_BIN}" ] || { echo "ERROR: no encuentro binario ligate en ${IMP_SIF}"; exit 1; }
    echo "  ligate: ${LIGATE_BIN}"

    # --- fasear cada chunk (resumible vía .csi) ---
    FILELIST="${CHUNKDIR}/chunks.files.txt"
    : > "${FILELIST}"
    i=0
    while IFS= read -r RANGE; do
        CK=$(printf "%02d" "$i")
        CHUNK_BCF="${CHUNKDIR}/chrX.chunk${CK}.common.bcf"
        if [ -f "${CHUNK_BCF}.csi" ]; then
            echo "  chrX chunk${CK} (chrX:${RANGE}): SKIP (ya existe)"
        else
            echo "  chrX chunk${CK}: faseando chrX:${RANGE}"
            apptainer exec ${BIND} "${IMP_SIF}" \
                shapeit5_phase_common \
                    --input "${IN_BCF}" \
                    --filter-maf 0.001 \
                    --map "${GMAP}" \
                    --output "${CHUNK_BCF}" \
                    --region "chrX:${RANGE}" \
                    --haploids "${HAPLOID_FILE}" \
                    --thread "${T}" \
                    --log "${CHUNKDIR}/chrX.chunk${CK}.log"
            apptainer exec ${BIND} "${BIO_SIF}" bcftools index -f "${CHUNK_BCF}"
        fi
        echo "${CHUNK_BCF}" >> "${FILELIST}"
        i=$((i+1))
    done < "${COORDS}"

    # --- ligar los chunks en un chrX completo ---
    echo ""
    echo "  Ligando ${i} chunks -> ${OUT_BCF}"
    apptainer exec ${BIND} "${IMP_SIF}" \
        "${LIGATE_BIN}" \
            --input "${FILELIST}" \
            --output "${OUT_BCF}" \
            --thread "${T}" \
            --index
fi
# ====================================================================

echo ""
echo -n "  Variantes faseadas: "
apptainer exec ${BIND} "${BIO_SIF}" bcftools view -H "${OUT_BCF}" | wc -l
echo "  Output: ${OUT_BCF} ($(ls -lh ${OUT_BCF} | awk '{print $5}'))"
echo "  Done: $(date)"
