#!/bin/bash
#SBATCH --job-name=sbrc_run
#SBATCH --account=researchers
#SBATCH --qos=vip
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=150G
#SBATCH --time=UNLIMITED
#SBATCH --output=/mnt/cephfs/orgs/home/angel.pacheco/logs/sbrc_run_%j.out
#SBATCH --error=/mnt/cephfs/orgs/home/angel.pacheco/logs/sbrc_run_%j.err
#SBATCH --mail-user=angelpachecosaavedra@gmail.com
#SBATCH --mail-type=BEGIN,END,FAIL

set -euo pipefail

SBRC_SIF="/mnt/cephfs/biocontainers/images/sbayesrc-0.2.6.sif"
SIF14="/mnt/cephfs/biocontainers/images/jupyter-biotools-1.4.sif"

PLINK_DIR="/mnt/cephfs/hot_nvme/mcps/imputed-topmed/plink_files/maximally_unrelated"
OUT_BASE="/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes/sbayesrc"
LD_DIR="${OUT_BASE}/ldm_79k"

MA_FILE="${OUT_BASE}/T2DGGI_ALL_hg38_matched.ma"

mkdir -p "${OUT_BASE}/results"

TMPDIR="/mnt/cephfs/scratch/$USER/sbrc_run_tmp/$SLURM_JOB_ID"
mkdir -p "$TMPDIR"
chmod 700 "$TMPDIR"
export TMPDIR TMP="$TMPDIR" TEMP="$TMPDIR"

echo "============================================================"
echo "  PIPELINE A — ETAPA 2: SBayesR (LD + MCMC)"
echo "  Nodo: $SLURMD_NODENAME"
echo "  RAM: $(free -g | awk 'NR==2{print $2}')G"
echo "  $(date)"
echo "============================================================"

if [ ! -f "$MA_FILE" ]; then
    echo " No se encontró $MA_FILE"
    exit 1
fi

echo "  Summary stats: $MA_FILE ($(wc -l < "$MA_FILE") líneas)"

R_SCRIPT="$TMPDIR/run_sbayesr.R"
cat << 'REOF' > "$R_SCRIPT"
library(SBayesRC)

ma_file   <- "/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes/sbayesrc/T2DGGI_ALL_hg38_matched.ma"
geno_dir  <- "/mnt/cephfs/hot_nvme/mcps/imputed-topmed/plink_files/maximally_unrelated/mcps-freeze150k_qcd_chr{CHR}_ivs"
out_dir   <- "/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes/sbayesrc/results"
ld_dir    <- "/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes/sbayesrc/ldm_79k"
threads   <- 24

# ================================================================
# PASO 1: Generando bloques LD
# ================================================================
cat("\n=== PASO 1: Generando bloques LD ===\n")

ld_script <- file.path(ld_dir, "ld.sh")

# Solo generar si no existe (LDstep1 ya corrió exitosamente antes)
if (!file.exists(ld_script)) {
    tryCatch({
        SBayesRC::LDstep1(
            mafile = ma_file, genoPrefix = geno_dir,
            outDir = ld_dir, genoCHR = "1-22", log2file = TRUE
        )
        cat(" Paso 1 completado\n")
    }, error = function(e) {
        cat(sprintf(" Error en LDstep1: %s\n", e$message))
        stop(e)
    })
} else {
    cat("  ld.sh ya existe, saltando LDstep1\n")
}

# ================================================================
# PASO 2: Generando matrices LD por bloque (PARALELO INTELIGENTE)
# ================================================================
cat("\n=== PASO 2: Generando matrices LD por bloque ===\n")

# 1. Leer el archivo original a la memoria RAM de R
ld_lines <- readLines(ld_script)
n_blocks <- length(ld_lines)
cat(sprintf("  Bloques totales en script: %d\n", n_blocks))

# 2. Reemplazar los hilos usando R puro (evita sed -i en CephFS/Apptainer)
ld_lines <- gsub("--thread [0-9]+ ", "", ld_lines)
ld_lines <- gsub("gctb ", "gctb --thread 4 ", ld_lines)

# 3. Detectar qué bloques ya terminaron para NO sobreescribirlos
completed_bins <- list.files(ld_dir, pattern = "\\.ldm\\.full\\.bin$")
completed <- length(completed_bins)
cat(sprintf("  Bloques ya completados en disco: %d\n", completed))

if (completed >= n_blocks * 0.95) {
    cat(sprintf("  %d/%d bloques completados (>95%%), saltando Paso 2 bulk\n", completed, n_blocks))
} else {
    lines_to_run <- ld_lines

    # Si ya hay bloques, los filtramos de la lista
    if (completed > 0) {
        completed_prefixes <- gsub("\\.ldm\\.full\\.bin$", "", completed_bins)
        for (pref in completed_prefixes) {
            pattern <- paste0("/", pref, "$|/", pref, " ")
            lines_to_run <- lines_to_run[!grepl(pattern, lines_to_run)]
        }
    }

    cat(sprintf("  Bloques pendientes a ejecutar: %d\n", length(lines_to_run)))

    # 4. Guardar SOLO los pendientes en un script temporal
    ld_pending <- file.path(ld_dir, "ld_pending.sh")
    writeLines(lines_to_run, ld_pending)

    # 5. Ejecutar paralelamente
    cat("  Lanzando ejecución paralela (6 bloques × 4 hilos = 24 cores)...\n")
    t0 <- Sys.time()

    system(sprintf("cat %s | xargs -d'\\n' -P 6 -I {} bash -c '{}'", ld_pending))

    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "hours"))
    completed_after <- length(list.files(ld_dir, pattern = "\\.ldm\\.full\\.bin$"))
    cat(sprintf(" Paso 2 completado en %.1f horas (%d bloques totales)\n", elapsed, completed_after))
}

# ================================================================
# PASO 2b: Re-intentar bloques fallidos (segfault) sin threading
# ================================================================
failed_blocks <- c("b248", "b307", "b406", "b529")
for (blk in failed_blocks) {
    bin_file <- file.path(ld_dir, paste0(blk, ".ldm.full.bin"))
    if (!file.exists(bin_file)) {
        # Buscar el comando original en ld.sh
        cmd_match <- grep(paste0("/", blk, " "), ld_lines, value = TRUE)
        if (length(cmd_match) > 0) {
            # Ejecutar sin --thread (single-threaded, más estable)
            cmd_clean <- gsub("--thread [0-9]+ ", "", cmd_match[1])
            cat(sprintf("  Re-intentando %s (single-thread)...\n", blk))
            system(cmd_clean)
            if (file.exists(bin_file)) {
                cat(sprintf("   %s completado\n", blk))
            } else {
                cat(sprintf("  %s falló de nuevo, continuando sin él\n", blk))
            }
        }
    }
}

# ================================================================
# PASO 3: Eigen-decomposición
# ================================================================
# PASO 3: Eigen-decomposición (un bloque a la vez)
cat("\n=== PASO 3: Eigen-decomposición ===\n")
bin_files <- list.files(ld_dir, pattern = "\\.ldm\\.full\\.bin$")
block_ids <- as.integer(gsub("b(\\d+)\\.ldm\\.full\\.bin", "\\1", bin_files))
block_ids <- sort(block_ids[!is.na(block_ids)])
cat(sprintf("  Descomponiendo %d bloques...\n", length(block_ids)))

n_done <- 0
n_fail <- 0
for (bi in block_ids) {
    tryCatch({
        SBayesRC::LDstep3(outDir = ld_dir, blockIndex = bi)
        n_done <- n_done + 1
        if (n_done %% 50 == 0) cat(sprintf("    Progreso: %d/%d\n", n_done, length(block_ids)))
    }, error = function(e) {
        cat(sprintf(" loque %d falló: %s\n", bi, e$message))
        n_fail <<- n_fail + 1
    })
}
cat(sprintf("  Eigen completado: %d ok, %d fallidos\n", n_done, n_fail))
# ================================================================
# PASO 4: SBayesR MCMC (sin anotaciones)
# ================================================================
cat("\n=== PASO 4: SBayesR MCMC ===\n")
tryCatch({
    # Intentar con threads primero
    SBayesRC::sbayesrc(
        mafile = ma_file, LDdir = ld_dir,
        outDir = out_dir, annot = NULL, threads = threads
    )
    cat("\n SBayesR MCMC completado\n")
}, error = function(e) {
    if (grepl("unused argument.*threads", e$message)) {
        cat("  threads no soportado, re-intentando sin threads...\n")
        SBayesRC::sbayesrc(
            mafile = ma_file, LDdir = ld_dir,
            outDir = out_dir, annot = NULL
        )
        cat("\n SBayesR MCMC completado (sin threads)\n")
    } else {
        cat(sprintf("\n Error en SBayesR: %s\n", e$message))
        stop(e)
    }
})

# ================================================================
# PASO 5: Exportando betas para scoring
# ================================================================
cat("\n=== PASO 5: Exportando betas para scoring ===\n")

snp_res <- file.path(out_dir, "sbayesrc.snpRes")
if (!file.exists(snp_res)) {
    snp_files <- list.files(out_dir, pattern = "\\.snpRes$", full.names = TRUE)
    if (length(snp_files) > 0) snp_res <- snp_files[1]
}

if (file.exists(snp_res)) {
    df_snp <- read.table(snp_res, header = TRUE)

    cat(sprintf("  Variantes totales: %d\n", nrow(df_snp)))
    cat(sprintf("  Variantes con efecto != 0: %d\n", sum(df_snp$A1Effect != 0)))

    out_weights <- file.path(dirname(out_dir), "sbayesrc_weights_for_scoring.tsv")
    write.table(
        df_snp[, c("Id", "A1", "A2", "A1Effect", "A1Frq")],
        out_weights, sep = "\t", row.names = FALSE, quote = FALSE,
        col.names = c("SNP", "A1", "A2", "BETA", "FREQ")
    )
    cat(sprintf("\n Pesos para scoring: %s\n", out_weights))
} else {
    cat(" No se encontró archivo .snpRes\n")
}

cat(sprintf("\n=== SBayesR completado: %s ===\n", Sys.time()))
REOF

if [ -f "$SBRC_SIF" ]; then
    apptainer exec --bind /mnt/cephfs:/mnt/cephfs --bind "$TMPDIR:/tmp,$TMPDIR:/var/tmp" \
        "$SBRC_SIF" Rscript "$R_SCRIPT"
else
    apptainer exec --bind /mnt/cephfs:/mnt/cephfs --bind "$TMPDIR:/tmp,$TMPDIR:/var/tmp" \
        "$SIF14" /opt/miniforge3/envs/r-bio/bin/Rscript "$R_SCRIPT"
fi

rm -rf "$TMPDIR"

echo ""
echo "============================================================"
echo "  PIPELINE A — ETAPA 2 COMPLETADA"
echo "  Siguiente: sbrc_03_scoring_gpu.sbatch"
echo "============================================================"
