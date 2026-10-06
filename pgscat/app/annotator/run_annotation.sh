#!/usr/bin/env bash
# ═══════════════════════════════════════════════════════════════════════════════
#  run_annotation.sh  –  Apptainer wrapper for the Variant Annotator pipeline
#  Version 1.6  (gnomAD 4.1.1 + MCPS population frequencies + af_effective)
# ═══════════════════════════════════════════════════════════════════════════════
#
#  Usage:
#    ./annotator/run_annotation.sh <PGS_ID> [OPTIONS]
#
#  Required environment variables (or edit defaults below):
#    DATA_DIR          Directory containing {PGS_ID}/ subfolders with betamap
#    ANNOTATIONS_DIR   Output root; {PGS_ID}/ subfolder is created automatically
#    GFF3_PATH         Full path to GENCODE GFF3 file
#    SIF_PATH          Full path to the built variant_annotator.sif image
#
#  Optional environment variables:
#    FASTA_PATH        Full path to GRCh38 reference FASTA (.fa / .fa.gz)
#                      Enables coding-consequence annotation (missense, synonymous,
#                      stop-gained, frameshift) and splice donor/acceptor via GT/AG.
#    DBNSFP_PATH       Full path to dbNSFP5 tabix file (dbNSFP5.0a_grch38.gz)
#                      Enables CADD, REVEL, SIFT, PolyPhen2, ClinVar scores.
#    DBSNP_FREQ_PATH   Full path to dbSNP population frequency VCF (freq.vcf.gz)
#                      Enables rsid, af_global, af_max_population, af_population_summary,
#                      and rarity_class annotation. Complementary to DBNSFP_PATH.
#                      Default: /mnt/cephfs/hot_nvme/dbsnp/population_frequency/freq.vcf.gz
#    GNOMAD_BASE       Base directory of gnomAD 4.1.1 data (must contain variants/exome/
#                      and variants/genome/ with chrN.vcf.bgz + .tbi files).
#                      Default: /mnt/cephfs/hot_nvme/gnomad_4.1.1
#                      Set to empty string ("") to disable gnomAD annotation.
#    MCPS_SOURCE_DIR   Directory containing MCPS per-chromosome TSV.GZ files
#                      (chrN.freq.tsv.gz).  A per-run parquet index is built in
#                      {ANNOTATIONS_DIR}/{PGS_ID}/mcps_idx/ on first use.
#                      Default: /mnt/cephfs/hot_nvme/mcps/mcps-variant-browser-afs
#                      Set to empty string ("") to disable MCPS annotation.
#    REGULATORY_BED    Comma- or space-separated list of BED file paths for
#                      regulatory element annotation.
#                      Defaults to derived files in:
#                        /mnt/cephfs/hot/pgscat/work/regulatory_derived/
#                        encodeCcreCombined_GRCh38.bed.gz  (ENCODE cCRE V3, highest priority)
#                        encode_tfbs_GRCh38.bed.gz         (ENCODE TFBS clusters, GRCh38)
#                        encode_dhs_GRCh38.bed.gz          (DNase DHS open chromatin, GRCh38)
#                        oreganno_GRCh38.bed.gz            (ORegAnno curated elements, GRCh38)
#                      Build with:
#                        python scripts/build_regulatory_beds.py --all
#                      Set to empty string ("") to disable regulatory annotation.
#
#  Example (full v1.6 with FASTA + dbNSFP5 + dbSNP + gnomAD + MCPS):
#    DATA_DIR=/mnt/cephfs/hot/pgscat/data \
#    ANNOTATIONS_DIR=/mnt/cephfs/hot/pgscat/annotations \
#    GFF3_PATH=/mnt/cephfs/hot_nvme/gencode/GRCh38.14/gencode.v49.basic.annotation.gff3 \
#    SIF_PATH=/mnt/cephfs/hot/pgscat/apptainer/variant_annotator.sif \
#    FASTA_PATH=/mnt/cephfs/hot_nvme/hg38/ref/hg38.fa \
#    DBNSFP_PATH=/mnt/cephfs/hot_nvme/dbNSFP5/dbNSFP5.0a_grch38.gz \
#    DBSNP_FREQ_PATH=/mnt/cephfs/hot_nvme/dbsnp/population_frequency/freq.vcf.gz \
#    GNOMAD_BASE=/mnt/cephfs/hot_nvme/gnomad_4.1.1 \
#    MCPS_SOURCE_DIR=/mnt/cephfs/hot_nvme/mcps/mcps-variant-browser-afs \
#    ./annotator/run_annotation.sh PGS000001
#
#  Example (minimal — GFF3 only):
#    DATA_DIR=/mnt/cephfs/hot/pgscat/data \
#    ANNOTATIONS_DIR=/mnt/cephfs/hot/pgscat/annotations \
#    GFF3_PATH=/mnt/cephfs/hot_nvme/gencode/GRCh38.14/gencode.v49.basic.annotation.gff3 \
#    SIF_PATH=/mnt/cephfs/hot/pgscat/apptainer/variant_annotator.sif \
#    ./annotator/run_annotation.sh PGS000001
#
#  Slurm (optional):
#    sbatch --partition=highmem --mem=32G --cpus-per-task=4 \
#           --job-name=annotate_PGS000001 \
#           ./annotator/run_annotation.sh PGS000001
# ═══════════════════════════════════════════════════════════════════════════════
set -euo pipefail

# ── Defaults (override via env) ───────────────────────────────────────────────
DATA_DIR="${DATA_DIR:-/mnt/cephfs/hot/pgscat/data}"
ANNOTATIONS_DIR="${ANNOTATIONS_DIR:-/mnt/cephfs/hot/pgscat/annotations}"
GFF3_PATH="${GFF3_PATH:-/mnt/cephfs/hot_nvme/gencode/GRCh38.14/gencode.v49.basic.annotation.gff3}"
SIF_PATH="${SIF_PATH:-/mnt/cephfs/hot/pgscat/apptainer/variant_annotator.sif}"
ANNOTATOR_SRC="${ANNOTATOR_SRC:-/mnt/cephfs/hot/pgscat/app/annotator}"

# Optional inputs (empty string = not used)
FASTA_PATH="${FASTA_PATH:-}"
DBNSFP_PATH="${DBNSFP_PATH:-}"         # tabix-indexed dbNSFP5 file
DBSNP_FREQ_PATH="${DBSNP_FREQ_PATH:-/mnt/cephfs/hot_nvme/dbsnp/population_frequency/freq.vcf.gz}"
GNOMAD_BASE="${GNOMAD_BASE:-/mnt/cephfs/hot_nvme/gnomad_4.1.1}"
# Pre-built parquet index (build once with gnomad_index.py — eliminates CephFS seek latency).
# Default: ${GNOMAD_BASE}/parquet  Set to "" to force tabix fallback.
GNOMAD_PARQUET_DIR="${GNOMAD_PARQUET_DIR:-${GNOMAD_BASE}/parquet}"
MCPS_SOURCE_DIR="${MCPS_SOURCE_DIR:-/mnt/cephfs/hot_nvme/mcps/mcps-variant-browser-afs}"

# Regulatory BED: default to derived files if they exist, allow override or disable with REGULATORY_BED=""
# Priority order: cCRE V3 > TFBS > DHS > ORegAnno (matches _SOURCE_PRIORITY in regulatory.py)
_REG_DERIVED_DIR="/mnt/cephfs/hot/pgscat/work/regulatory_derived"
_REG_CCRE="${_REG_DERIVED_DIR}/encodeCcreCombined_GRCh38.bed.gz"
_REG_TFBS="${_REG_DERIVED_DIR}/encode_tfbs_GRCh38.bed.gz"
_REG_DHS="${_REG_DERIVED_DIR}/encode_dhs_GRCh38.bed.gz"
_REG_OREG="${_REG_DERIVED_DIR}/oreganno_GRCh38.bed.gz"
if [[ -z "${REGULATORY_BED+x}" ]]; then
    # REGULATORY_BED not set at all — auto-detect derived files
    _AUTO_REG=()
    [[ -f "${_REG_CCRE}" ]] && _AUTO_REG+=("${_REG_CCRE}")  # highest priority
    [[ -f "${_REG_TFBS}" ]] && _AUTO_REG+=("${_REG_TFBS}")
    [[ -f "${_REG_DHS}"  ]] && _AUTO_REG+=("${_REG_DHS}")
    [[ -f "${_REG_OREG}" ]] && _AUTO_REG+=("${_REG_OREG}")
    REGULATORY_BED="${_AUTO_REG[*]}"  # space-separated (may be empty if not built)
else
    # REGULATORY_BED was explicitly set (possibly to "" to disable)
    : # use as-is
fi

# ── Parse arguments ───────────────────────────────────────────────────────────
if [[ $# -lt 1 ]]; then
    echo "Usage: $0 <PGS_ID> [--no-parquet] [--verbose]" >&2
    exit 1
fi

PGS_ID="${1^^}"   # uppercase
shift
EXTRA_ARGS=("$@")

# ── Resolve paths ─────────────────────────────────────────────────────────────
BETAMAP_PATH="${DATA_DIR}/${PGS_ID}/${PGS_ID}_hmPOS_GRCh38.betamap.tsv.gz"
OUT_DIR="${ANNOTATIONS_DIR}/${PGS_ID}"
LOG_FILE="${OUT_DIR}/annotation.log"

# ── Pre-flight checks ─────────────────────────────────────────────────────────
if [[ ! -f "${SIF_PATH}" ]]; then
    echo "ERROR: Apptainer image not found: ${SIF_PATH}" >&2
    echo "Build it with:" >&2
    echo "  apptainer build ${SIF_PATH} /mnt/cephfs/hot/pgscat/app/apptainer/variant_annotator.def" >&2
    exit 1
fi

if [[ ! -f "${BETAMAP_PATH}" ]]; then
    echo "ERROR: Betamap not found: ${BETAMAP_PATH}" >&2
    echo "Run score preparation first via the web platform." >&2
    exit 1
fi

if [[ ! -f "${GFF3_PATH}" ]]; then
    echo "ERROR: GFF3 not found: ${GFF3_PATH}" >&2
    exit 1
fi

if [[ -n "${FASTA_PATH}" && ! -f "${FASTA_PATH}" ]]; then
    echo "ERROR: FASTA not found: ${FASTA_PATH}" >&2
    exit 1
fi

if [[ -n "${DBNSFP_PATH}" ]]; then
    if [[ ! -f "${DBNSFP_PATH}" ]]; then
        echo "ERROR: dbNSFP5 not found: ${DBNSFP_PATH}" >&2
        exit 1
    fi
    if [[ ! -f "${DBNSFP_PATH}.tbi" ]]; then
        echo "ERROR: dbNSFP5 tabix index not found: ${DBNSFP_PATH}.tbi" >&2
        exit 1
    fi
fi

if [[ -n "${DBSNP_FREQ_PATH}" ]]; then
    if [[ ! -f "${DBSNP_FREQ_PATH}" ]]; then
        echo "WARNING: dbSNP freq VCF not found: ${DBSNP_FREQ_PATH} — skipping population frequencies" >&2
        DBSNP_FREQ_PATH=""
    elif [[ ! -f "${DBSNP_FREQ_PATH}.tbi" && ! -f "${DBSNP_FREQ_PATH}.csi" ]]; then
        echo "WARNING: dbSNP freq VCF tabix index not found: ${DBSNP_FREQ_PATH}.tbi — skipping" >&2
        DBSNP_FREQ_PATH=""
    fi
fi

# ── Prepare output directory ──────────────────────────────────────────────────
mkdir -p "${OUT_DIR}"
echo "=== Variant Annotator v1.5  PGS_ID=${PGS_ID}  $(date -Iseconds) ===" | tee "${LOG_FILE}"

# ── Build bind mounts and Python args ─────────────────────────────────────────
BETAMAP_IN_CONTAINER="/data/${PGS_ID}/${PGS_ID}_hmPOS_GRCh38.betamap.tsv.gz"
GFF3_DIR="$(dirname "${GFF3_PATH}")"
GFF3_FILENAME="$(basename "${GFF3_PATH}")"
GFF3_IN_CONTAINER="/ref/${GFF3_FILENAME}"
OUT_IN_CONTAINER="/annotations/${PGS_ID}"

BIND_OPTS=(
    "--bind" "${DATA_DIR}:/data:ro"
    "--bind" "${GFF3_DIR}:/ref:ro"
    "--bind" "${ANNOTATIONS_DIR}:/annotations:rw"
    "--bind" "${ANNOTATOR_SRC}:/app/annotator:ro"
)

PYTHON_ARGS=(
    "--betamap" "${BETAMAP_IN_CONTAINER}"
    "--gff3"    "${GFF3_IN_CONTAINER}"
    "--outdir"  "${OUT_IN_CONTAINER}"
    "--pgs-id"  "${PGS_ID}"
)

# ── Optional FASTA ────────────────────────────────────────────────────────────
if [[ -n "${FASTA_PATH}" ]]; then
    FASTA_DIR="$(dirname "${FASTA_PATH}")"
    FASTA_FILENAME="$(basename "${FASTA_PATH}")"
    BIND_OPTS+=("--bind" "${FASTA_DIR}:/fasta:ro")
    PYTHON_ARGS+=("--fasta" "/fasta/${FASTA_FILENAME}")
    echo "  FASTA:   ${FASTA_PATH}" | tee -a "${LOG_FILE}"
fi

# ── Optional dbNSFP5 ──────────────────────────────────────────────────────────
if [[ -n "${DBNSFP_PATH}" ]]; then
    DBNSFP_DIR="$(dirname "${DBNSFP_PATH}")"
    DBNSFP_FILENAME="$(basename "${DBNSFP_PATH}")"
    BIND_OPTS+=("--bind" "${DBNSFP_DIR}:/dbnsfp:ro")
    PYTHON_ARGS+=("--dbnsfp" "/dbnsfp/${DBNSFP_FILENAME}")
    echo "  dbNSFP5: ${DBNSFP_PATH}" | tee -a "${LOG_FILE}"
fi

# ── Optional dbSNP population frequency VCF (v1.3) ───────────────────────────
if [[ -n "${DBSNP_FREQ_PATH}" ]]; then
    DBSNP_FREQ_DIR="$(dirname "${DBSNP_FREQ_PATH}")"
    DBSNP_FREQ_FILENAME="$(basename "${DBSNP_FREQ_PATH}")"
    BIND_OPTS+=("--bind" "${DBSNP_FREQ_DIR}:/dbsnp_freq:ro")
    PYTHON_ARGS+=("--dbsnp" "/dbsnp_freq/${DBSNP_FREQ_FILENAME}")
    echo "  dbSNP:   ${DBSNP_FREQ_PATH}" | tee -a "${LOG_FILE}"
fi

# ── Optional regulatory BED files ─────────────────────────────────────────────
REG_BED_ARGS=()
if [[ -n "${REGULATORY_BED}" ]]; then
    # REGULATORY_BED_DIR is a single shared parent directory for all BED files.
    # All provided BED paths must reside under the same directory tree.
    # For multi-directory setups, extend this section to add multiple binds.
    #
    # Simple single-directory approach: bind the directory of the first BED file.
    # More complex setups: set REGULATORY_BED_DIR explicitly.
    REGULATORY_BED_DIR="${REGULATORY_BED_DIR:-}"

    # shellcheck disable=SC2206
    REG_BED_ARRAY=(${REGULATORY_BED})

    if [[ -z "${REGULATORY_BED_DIR}" && ${#REG_BED_ARRAY[@]} -gt 0 ]]; then
        REGULATORY_BED_DIR="$(dirname "${REG_BED_ARRAY[0]}")"
    fi

    if [[ -n "${REGULATORY_BED_DIR}" ]]; then
        BIND_OPTS+=("--bind" "${REGULATORY_BED_DIR}:/regulatory:ro")
    fi

    for BED_PATH in "${REG_BED_ARRAY[@]}"; do
        BED_FILENAME="$(basename "${BED_PATH}")"
        REG_BED_ARGS+=("/regulatory/${BED_FILENAME}")
        echo "  RegBED:  ${BED_PATH}" | tee -a "${LOG_FILE}"
    done

    if [[ ${#REG_BED_ARGS[@]} -gt 0 ]]; then
        PYTHON_ARGS+=("--regulatory-bed" "${REG_BED_ARGS[@]}")
    fi
fi

# ── Optional gnomAD 4.1.1 (v1.4+) ───────────────────────────────────────────
# Fast path: if GNOMAD_PARQUET_DIR exists, bind it and pass --gnomad-parquet-dir.
#   Build once: python gnomad_index.py --gnomad-dir ${GNOMAD_BASE}/variants \
#                                      --out-dir ${GNOMAD_BASE}/parquet
# Slow path: tabix range queries (per-variant CephFS seeks — slow for large PGS).
if [[ -n "${GNOMAD_BASE}" ]]; then
    if [[ -d "${GNOMAD_BASE}/variants" ]]; then
        BIND_OPTS+=("--bind" "${GNOMAD_BASE}:/gnomad:ro")
        PYTHON_ARGS+=("--gnomad-dir" "/gnomad/variants")
        # Parquet index fast path
        if [[ -n "${GNOMAD_PARQUET_DIR}" && -d "${GNOMAD_PARQUET_DIR}" ]]; then
            BIND_OPTS+=("--bind" "${GNOMAD_PARQUET_DIR}:/gnomad_parquet:ro")
            PYTHON_ARGS+=("--gnomad-parquet-dir" "/gnomad_parquet")
            echo "  gnomAD:  ${GNOMAD_BASE} [active, parquet fast path: ${GNOMAD_PARQUET_DIR}]" | tee -a "${LOG_FILE}"
        else
            echo "  gnomAD:  ${GNOMAD_BASE} [active, tabix fallback — run gnomad_index.py to build parquet]" | tee -a "${LOG_FILE}"
        fi
    else
        echo "  gnomAD:  ${GNOMAD_BASE}/variants not found — skipped" | tee -a "${LOG_FILE}"
    fi
else
    echo "  gnomAD:  GNOMAD_BASE unset — skipped" | tee -a "${LOG_FILE}"
fi

# ── Optional MCPS population frequencies (v1.6) ───────────────────────────────
# Mounts the source TSV.GZ directory.  The annotator builds a per-run parquet
# index in {ANNOTATIONS_DIR}/{PGS_ID}/mcps_idx/ (writable via /annotations).
if [[ -n "${MCPS_SOURCE_DIR}" ]]; then
    if [[ -d "${MCPS_SOURCE_DIR}" ]]; then
        BIND_OPTS+=("--bind" "${MCPS_SOURCE_DIR}:/mcps:ro")
        PYTHON_ARGS+=("--mcps-dir" "/mcps")
        echo "  MCPS:    ${MCPS_SOURCE_DIR} [active]" | tee -a "${LOG_FILE}"
    else
        echo "  MCPS:    ${MCPS_SOURCE_DIR} not found — skipped" | tee -a "${LOG_FILE}"
    fi
else
    echo "  MCPS:    MCPS_SOURCE_DIR unset — skipped" | tee -a "${LOG_FILE}"
fi

# ── Log run parameters ────────────────────────────────────────────────────────
echo "Running Apptainer:" | tee -a "${LOG_FILE}"
echo "  SIF:     ${SIF_PATH}" | tee -a "${LOG_FILE}"
echo "  Betamap: ${BETAMAP_PATH}" | tee -a "${LOG_FILE}"
echo "  GFF3:    ${GFF3_PATH}" | tee -a "${LOG_FILE}"
echo "  Output:  ${OUT_DIR}" | tee -a "${LOG_FILE}"
echo "" | tee -a "${LOG_FILE}"

# ── Run inside Apptainer ──────────────────────────────────────────────────────
# Bind mount layout inside container:
#   /data        → betamap directory (read-only)
#   /ref         → GFF3 directory (read-only)
#   /fasta       → FASTA directory (read-only, if --fasta set)
#   /dbnsfp      → dbNSFP5 directory (read-only, if --dbnsfp set)
#   /dbsnp_freq  → dbSNP population freq VCF directory (read-only, if DBSNP_FREQ_PATH set)
#   /gnomad      → gnomAD 4.1.1 base directory (read-only, if GNOMAD_BASE set)
#   /mcps        → MCPS TSV.GZ source directory (read-only, if MCPS_SOURCE_DIR set)
#   /regulatory  → regulatory BED directory (read-only, if --regulatory-bed set)
#   /annotations → output root (read-write; mcps_idx/ written here)
#   /app/annotator → annotator source (read-only)

apptainer exec \
    "${BIND_OPTS[@]}" \
    "${SIF_PATH}" \
    python /app/annotator/annotate_variants.py \
        "${PYTHON_ARGS[@]}" \
        "${EXTRA_ARGS[@]}" \
    2>&1 | tee -a "${LOG_FILE}"

EXIT_CODE="${PIPESTATUS[0]}"

if [[ ${EXIT_CODE} -eq 0 ]]; then
    echo "" | tee -a "${LOG_FILE}"
    echo "=== SUCCESS  $(date -Iseconds) ===" | tee -a "${LOG_FILE}"
    echo "Outputs:"
    ls -lh "${OUT_DIR}"/*.tsv.gz "${OUT_DIR}"/*.parquet "${OUT_DIR}"/*.json 2>/dev/null | \
        awk '{print "  " $5 "  " $9}' | tee -a "${LOG_FILE}"
else
    echo "" | tee -a "${LOG_FILE}"
    echo "=== FAILED (exit ${EXIT_CODE})  $(date -Iseconds) ===" | tee -a "${LOG_FILE}"
fi

exit "${EXIT_CODE}"
