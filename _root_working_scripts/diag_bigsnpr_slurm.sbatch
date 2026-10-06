#!/bin/bash
#SBATCH --job-name=diag_bigsnpr
#SBATCH --account=researchers
#SBATCH --qos=vip
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=20G
#SBATCH --time=00:30:00
#SBATCH --output=/mnt/cephfs/orgs/home/angel.pacheco/logs/diag_bigsnpr_%j.out
#SBATCH --error=/mnt/cephfs/orgs/home/angel.pacheco/logs/diag_bigsnpr_%j.err

set -uo pipefail   # SIN -e: queremos que continúe aunque algo falle

# ================================================================
# DIAGNÓSTICO EMPÍRICO: qué vars SLURM ve bigsnpr dentro del cont.
# Lanzar separadamente en cada nodo a investigar:
#   sbatch --nodelist=MCPS-HPC02 diag_bigsnpr_slurm.sbatch
#   sbatch --nodelist=MCPS-HPC04 diag_bigsnpr_slurm.sbatch
# Comparar outputs.
# ================================================================

SIF="/mnt/cephfs/biocontainers/images/jupyter-biotools-1.5.sif"
PLINK_DIR="/mnt/cephfs/hot_nvme/mcps/imputed-topmed/plink_files/maximally_unrelated"

TMPDIR="/mnt/cephfs/scratch/$USER/diag_bigsnpr_$SLURM_JOB_ID"
mkdir -p "$TMPDIR"; chmod 700 "$TMPDIR"
export TMPDIR TMP="$TMPDIR" TEMP="$TMPDIR"

echo "============================================================"
echo "  DIAG bigsnpr — Nodo: $SLURMD_NODENAME"
echo "  Job: $SLURM_JOB_ID   $(date)"
echo "============================================================"

# ──────────────────────────────────────────────────────────────────
# [A] Variables SLURM_* en el HOST (antes de entrar al contenedor)
# ──────────────────────────────────────────────────────────────────
echo ""
echo "═══════════ [A] SLURM_* vars en el HOST ═══════════"
env | grep -E '^(SLURM|SLURMD)' | sort
echo ""
echo "uname -a: $(uname -a)"
echo "nproc:    $(nproc)"
echo "Memoria libre: $(free -g | awk 'NR==2{print $7}')G de $(free -g | awk 'NR==2{print $2}')G"
echo "Carga: $(uptime)"

# ──────────────────────────────────────────────────────────────────
# [B] Verificar cgroups v2 (lo que parallelly consulta también)
# ──────────────────────────────────────────────────────────────────
echo ""
echo "═══════════ [B] Cgroups v2 (CPU affinity) ═══════════"
for cg in /sys/fs/cgroup/cpuset.cpus /sys/fs/cgroup/cpuset.cpus.effective; do
  if [ -r "$cg" ]; then
    echo "$cg: $(cat $cg)"
  else
    echo "$cg: no legible"
  fi
done

# ──────────────────────────────────────────────────────────────────
# [C] Comparar variables DENTRO vs FUERA del contenedor
# ──────────────────────────────────────────────────────────────────
R_DIAG="$TMPDIR/diag.R"
cat << 'REOF' > "$R_DIAG"
cat("══════════ [C] DENTRO del contenedor: vars SLURM_* visibles ══════════\n")
slurm_vars <- Sys.getenv()
slurm_only <- slurm_vars[grepl("^(SLURM|R_PARALLELLY)", names(slurm_vars))]
for (n in sort(names(slurm_only))) {
  cat(sprintf("  %s = %s\n", n, slurm_only[[n]]))
}

cat("\n══════════ [D] parallelly view ══════════\n")
suppressMessages({
  library(parallelly)
  library(bigparallelr)
})
cat(sprintf("  parallelly version: %s\n", as.character(packageVersion("parallelly"))))
cat(sprintf("  bigparallelr version: %s\n", as.character(packageVersion("bigparallelr"))))

# Probar cada método de availableCores por separado (sin que falle el script)
methods <- c("system", "nproc", "mc.cores", "cgroups.cpuset",
             "cgroups.cpuquota", "Slurm", "SGE", "PBS", "LSF",
             "fallback", "custom")
for (m in methods) {
  r <- tryCatch(parallelly::availableCores(methods=m),
                error=function(e) sprintf("ERROR: %s", conditionMessage(e)),
                warning=function(w) sprintf("WARN: %s", conditionMessage(w)))
  cat(sprintf("  availableCores(methods='%s'): %s\n", m, as.character(r[1])))
}

cat(sprintf("\n  ► availableCores() default: %d\n", parallelly::availableCores()))
cat(sprintf("  ► availableCores(logical=FALSE): %d\n", parallelly::availableCores(logical=FALSE)))
cat(sprintf("  ► bigparallelr::nb_cores(): %d\n", bigparallelr::nb_cores()))

cat("\n══════════ [E] bigsnpr versión y settings ══════════\n")
suppressMessages(library(bigsnpr))
cat(sprintf("  bigsnpr version: %s\n", as.character(packageVersion("bigsnpr"))))
cat(sprintf("  bigsparser version: %s\n", as.character(packageVersion("bigsparser"))))
cat(sprintf("  bigstatsr version: %s\n", as.character(packageVersion("bigstatsr"))))
cat(sprintf("  OMP_NUM_THREADS = %s\n", Sys.getenv("OMP_NUM_THREADS")))
cat(sprintf("  OPENBLAS_NUM_THREADS = %s\n", Sys.getenv("OPENBLAS_NUM_THREADS")))

cat("\n══════════ [F] Test pequeño: snp_cor con 1000 variantes ══════════\n")
# Usar el .bed de chr22 (más pequeño) para test rápido
plink_dir <- Sys.getenv("PLINK_DIR")
tmp_dir <- Sys.getenv("TMPDIR")
bed_file <- file.path(plink_dir, "mcps-freeze150k_qcd_chr22_ivs.bed")

if (!file.exists(bed_file)) {
  cat(sprintf("  No existe %s — saltando test snp_cor\n", bed_file))
} else {
  cat(sprintf("  Leyendo %s...\n", bed_file))
  t0 <- Sys.time()
  bk_file <- file.path(tmp_dir, "diag_chr22")
  snp_readBed(bed_file, backingfile=bk_file)
  obj <- snp_attach(paste0(bk_file, ".rds"))
  G <- obj$genotypes
  cat(sprintf("  G: %d × %d, lectura: %.1fs\n", nrow(G), ncol(G),
              as.numeric(difftime(Sys.time(), t0, units="secs"))))

  # Tomar 1000 variantes random del medio del chr para test
  set.seed(42)
  m_total <- ncol(G)
  ind_test <- sort(sample(seq(2000, m_total-2000), 1000))
  pos_test <- obj$map$physical.pos[ind_test] / 1000

  cat("\n  Ejecutando snp_cor con 1000 variantes, ncores=4, thr_r2=0.005...\n")
  t0 <- Sys.time()
  corr <- tryCatch({
    snp_cor(Gna=G, ind.col=ind_test, infos.pos=pos_test,
            size=3000, ncores=4L, thr_r2=0.005)
  }, error=function(e) {
    cat(sprintf("  ✗ ERROR: %s\n", conditionMessage(e)))
    NULL
  })
  t_cor <- as.numeric(difftime(Sys.time(), t0, units="secs"))

  if (!is.null(corr)) {
    cat(sprintf("  ✓ snp_cor OK: %.1fs, dim=%dx%d, nnz=%s\n",
                t_cor, nrow(corr), ncol(corr),
                format(length(corr@x), big.mark=",")))
  }
}

cat("\n══════════ DIAG FIN ══════════\n")
REOF

echo ""
echo "═══════════ Ejecutando dentro del contenedor ═══════════"
apptainer exec \
  --bind /mnt/cephfs:/mnt/cephfs \
  --bind "$TMPDIR:/tmp,$TMPDIR:/var/tmp" \
  --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 --env MKL_NUM_THREADS=1 \
  --env BLIS_NUM_THREADS=1 --env NUMEXPR_NUM_THREADS=1 \
  --env PLINK_DIR="$PLINK_DIR" \
  "$SIF" /opt/miniforge3/envs/r-bio/bin/Rscript "$R_DIAG"

EXIT_CODE=$?
echo ""
echo "Exit code de R: $EXIT_CODE"

rm -rf "$TMPDIR"

echo ""
echo "============================================================"
echo "  DIAG FIN — $(date)"
echo "============================================================"
