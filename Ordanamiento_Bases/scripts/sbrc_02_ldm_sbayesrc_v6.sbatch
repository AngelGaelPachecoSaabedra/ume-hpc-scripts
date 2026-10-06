#!/bin/bash
#SBATCH --job-name=sbrc_run
#SBATCH --account=researchers
#SBATCH --qos=vip
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=150G
#SBATCH --time=UNLIMITED
#SBATCH --output=/mnt/cephfs/orgs/home/angel.pacheco/logs/sbrc_run_%j.out
#SBATCH --error=/mnt/cephfs/orgs/home/angel.pacheco/logs/sbrc_run_%j.err
#SBATCH --mail-user=angelpachecosaavedra@gmail.com
#SBATCH --mail-type=BEGIN,END,FAIL

set -euo pipefail

# ================================================================
# PIPELINE A — SBayesRC v6 (DEFINITIVO)
#
# Todas las lecciones aprendidas de v1-v5:
#
# API SBayesRC 0.2.6:
#   - LDstep3(outDir, blockIndex, thresh) — NO acepta threads
#   - sbayesrc(mafile, LDdir, outPrefix, annot="") — NO outDir, NO threads
#   - LDstep3 procesa UN bloque por llamada (blockIndex obligatorio)
#   - LDstep3 NO genera snp.info — hay que construirlo manualmente
#   - sbayesrc() lee snp.info y requiere columna "Block"
#
# Bloques problemáticos excluidos:
#   - bNA: variantes sin cromosoma (no se genera .ldm.full.bin)
#   - b114 (58,831 vars): nMarker*nMarker overflow int32 en LDstep3
#   - b371 (46,633 vars): mismo overflow
#   Total excluido: ~105K variantes (~1% de 10M)
#
# Formato eigen .bin (float32 todo):
#   m(4) + k(4) + lambdaSum(4) + thresh(4) + eigenvals[k](4) + eigenvecs[k*m](4)
#
# Flujo:
#   1. LDstep1 (generar ld.sh + snplists)
#   2. GCTB LD matrices (paralelo 6×4 threads)
#   2b. Re-intentar bloques fallidos en single-thread
#   3. LDstep3 eigen por bloque (skipea existentes)
#   3b. Construir snp.info y ldm.info consistentes
#   4. sbayesrc MCMC
#   5. Exportar betas
# ================================================================

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
echo "  PIPELINE A — SBayesRC v6 (DEFINITIVO)"
echo "  Nodo: $SLURMD_NODENAME"
echo "  RAM: $(free -g | awk 'NR==2{print $2}')G"
echo "  $(date)"
echo "============================================================"

if [ ! -f "$MA_FILE" ]; then
    echo "No se encontró $MA_FILE"
    exit 1
fi

echo "  Summary stats: $MA_FILE ($(wc -l < "$MA_FILE") líneas)"

R_SCRIPT="$TMPDIR/run_sbayesr_v6.R"
cat << 'REOF' > "$R_SCRIPT"
library(SBayesRC)
library(data.table)

ma_file   <- "/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes/sbayesrc/T2DGGI_ALL_hg38_matched.ma"
geno_dir  <- "/mnt/cephfs/hot_nvme/mcps/imputed-topmed/plink_files/maximally_unrelated/mcps-freeze150k_qcd_chr{CHR}_ivs"
out_dir   <- "/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes/sbayesrc/results"
ld_dir    <- "/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes/sbayesrc/ldm_79k"

# Bloques que causan integer overflow en LDstep3 (nMarker > 46340 → nMarker^2 > 2^31)
# y el bloque NA (variantes sin cromosoma asignado)
SKIP_BLOCKS <- c("NA", "114", "371")

# ================================================================
# PASO 1: LDstep1 — generar bloques LD
# ================================================================
cat("\n=== PASO 1: Generando bloques LD ===\n")

ld_script <- file.path(ld_dir, "ld.sh")

if (!file.exists(ld_script)) {
    SBayesRC::LDstep1(
        mafile = ma_file, genoPrefix = geno_dir,
        outDir = ld_dir, genoCHR = "1-22", log2file = TRUE
    )
    cat(" LDstep1 completado\n")
} else {
    cat("  ld.sh ya existe, saltando LDstep1\n")
}

# ================================================================
# PASO 2: GCTB LD matrices por bloque (paralelo)
# ================================================================
cat("\n=== PASO 2: Generando matrices LD por bloque ===\n")

ld_lines <- readLines(ld_script)
n_blocks <- length(ld_lines)
cat(sprintf("  Bloques totales en script: %d\n", n_blocks))

# Ajustar threads: 6 bloques × 4 hilos = 24 cores
ld_lines <- gsub("--thread [0-9]+ ", "", ld_lines)
ld_lines <- gsub("gctb ", "gctb --thread 4 ", ld_lines)

completed_bins <- list.files(ld_dir, pattern = "\\.ldm\\.full\\.bin$")
completed <- length(completed_bins)
cat(sprintf("  Bloques completados: %d\n", completed))

if (completed >= n_blocks * 0.95) {
    cat(sprintf("  %d/%d (>95%%), saltando bulk\n", completed, n_blocks))
} else {
    # Filtrar bloques ya completados
    lines_to_run <- ld_lines
    if (completed > 0) {
        completed_prefixes <- gsub("\\.ldm\\.full\\.bin$", "", completed_bins)
        for (pref in completed_prefixes) {
            pattern <- paste0("/", pref, "$|/", pref, " ")
            lines_to_run <- lines_to_run[!grepl(pattern, lines_to_run)]
        }
    }
    cat(sprintf("  Pendientes: %d\n", length(lines_to_run)))

    ld_pending <- file.path(ld_dir, "ld_pending.sh")
    writeLines(lines_to_run, ld_pending)

    cat("  Ejecutando (6 bloques × 4 hilos)...\n")
    t0 <- Sys.time()
    system(sprintf("cat %s | xargs -d'\\n' -P 6 -I {} bash -c '{}'", ld_pending))
    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "hours"))
    completed_after <- length(list.files(ld_dir, pattern = "\\.ldm\\.full\\.bin$"))
    cat(sprintf(" %.1f horas, %d bloques\n", elapsed, completed_after))
}

# ================================================================
# PASO 2b: Re-intentar bloques con segfault (single-thread)
# ================================================================
cat("\n=== PASO 2b: Re-intentando bloques fallidos ===\n")

known_segfault <- c("b248", "b307", "b406", "b529")
for (blk in known_segfault) {
    bin_file <- file.path(ld_dir, paste0(blk, ".ldm.full.bin"))
    if (!file.exists(bin_file)) {
        cmd_match <- grep(paste0("/", blk, " "), ld_lines, value = TRUE)
        if (length(cmd_match) > 0) {
            cmd_clean <- gsub("--thread [0-9]+ ", "", cmd_match[1])
            cat(sprintf("  Re-intentando %s (single-thread)...\n", blk))
            system(cmd_clean)
            if (file.exists(bin_file)) {
                cat(sprintf("  %s completado\n", blk))
            } else {
                cat(sprintf("  ⚠️%s falló, continuando\n", blk))
            }
        }
    }
}

# ================================================================
# PASO 3: Eigen-decomposición (bloque por bloque)
# ================================================================
cat("\n=== PASO 3: Eigen-decomposición ===\n")

# Leer ldm.info ORIGINAL para obtener la lista completa de bloques
ldm_info_orig <- file.path(ld_dir, "ldm.info")
if (!file.exists(ldm_info_orig)) {
    stop("No se encontró ldm.info — LDstep1 no se ejecutó correctamente")
}

dt_ldm <- fread(ldm_info_orig)
all_block_ids <- as.character(dt_ldm$Block)
cat(sprintf("  Bloques en ldm.info: %d\n", nrow(dt_ldm)))

# Determinar bloques válidos: tienen .ldm.full.bin, no están en SKIP_BLOCKS,
# y nMarker < 46340 (evitar integer overflow en LDstep3)
valid_blocks <- c()
skip_reasons <- list()

for (i in seq_len(nrow(dt_ldm))) {
    blk <- as.character(dt_ldm$Block[i])
    n_snps <- dt_ldm$NumSnps[i]
    bin_file <- file.path(ld_dir, paste0("b", blk, ".ldm.full.bin"))

    if (blk %in% SKIP_BLOCKS) {
        skip_reasons[[blk]] <- "en lista SKIP_BLOCKS"
        next
    }
    if (!file.exists(bin_file)) {
        skip_reasons[[blk]] <- "sin .ldm.full.bin"
        next
    }
    if (n_snps > 46340) {
        # nMarker * nMarker overflow int32 (sqrt(2^31) ≈ 46340)
        skip_reasons[[blk]] <- sprintf("nMarker=%d > 46340 (overflow)", n_snps)
        next
    }
    valid_blocks <- c(valid_blocks, blk)
}

cat(sprintf("  Bloques válidos: %d\n", length(valid_blocks)))
cat(sprintf("  Bloques excluidos: %d\n", length(skip_reasons)))
for (blk in names(skip_reasons)) {
    cat(sprintf("    b%s: %s\n", blk, skip_reasons[[blk]]))
}

# Generar eigen para bloques válidos que aún no lo tienen
# LDstep3 usa blockIndex = posición en ldm.info (1-indexed)
n_eigen_generated <- 0
n_eigen_skipped <- 0

for (blk in valid_blocks) {
    eigen_file <- file.path(ld_dir, paste0("block", blk, ".eigen.bin"))

    if (file.exists(eigen_file)) {
        n_eigen_skipped <- n_eigen_skipped + 1
        next
    }

    # Encontrar blockIndex (posición en ldm.info)
    blk_idx <- which(dt_ldm$Block == as.integer(blk))
    if (length(blk_idx) == 0) {
        blk_idx <- which(as.character(dt_ldm$Block) == blk)
    }
    if (length(blk_idx) == 0) {
        cat(sprintf(" b%s no encontrado en ldm.info\n", blk))
        next
    }

    cat(sprintf("  Generando eigen para b%s (blockIndex=%d)...\n", blk, blk_idx))
    tryCatch({
        SBayesRC::LDstep3(outDir = ld_dir, blockIndex = blk_idx)
        n_eigen_generated <- n_eigen_generated + 1
        cat(sprintf(" block%s.eigen.bin\n", blk))
    }, error = function(e) {
        cat(sprintf("b%s falló: %s\n", blk, e$message))
        # Agregar a skip si falló
        skip_reasons[[blk]] <<- sprintf("LDstep3 error: %s", e$message)
        valid_blocks <<- valid_blocks[valid_blocks != blk]
    })
}

cat(sprintf("\n  Eigen: %d existentes, %d nuevos, %d total\n",
    n_eigen_skipped, n_eigen_generated,
    n_eigen_skipped + n_eigen_generated))

# Filtrar valid_blocks a solo los que REALMENTE tienen eigen
final_blocks <- c()
for (blk in valid_blocks) {
    eigen_file <- file.path(ld_dir, paste0("block", blk, ".eigen.bin"))
    if (file.exists(eigen_file)) {
        final_blocks <- c(final_blocks, blk)
    }
}
cat(sprintf("  Bloques con eigen confirmados: %d\n", length(final_blocks)))

# ================================================================
# PASO 3b: Construir ldm.info y snp.info consistentes
# ================================================================
cat("\n=== PASO 3b: Construyendo snp.info y ldm.info ===\n")

# --- Reconstruir ldm.info solo con bloques finales ---
dt_final <- dt_ldm[as.character(dt_ldm$Block) %in% final_blocks, ]

# Recalcular StartSnpIdx y EndSnpIdx
idx <- 0L
new_start <- integer(nrow(dt_final))
new_end <- integer(nrow(dt_final))
for (i in seq_len(nrow(dt_final))) {
    new_start[i] <- idx
    new_end[i] <- idx + dt_final$NumSnps[i] - 1L
    idx <- idx + dt_final$NumSnps[i]
}
dt_final$StartSnpIdx <- new_start
dt_final$EndSnpIdx <- new_end

fwrite(dt_final, file.path(ld_dir, "ldm.info"), sep = "\t")
cat(sprintf("  ldm.info: %d bloques, %d variantes totales\n",
    nrow(dt_final), sum(dt_final$NumSnps)))

# --- Construir snp.info concatenando b{N}.ldm.full.info ---
# Necesita: todas las columnas originales + columna Block + Index global
cat("  Construyendo snp.info...\n")

snp_info_file <- file.path(ld_dir, "snp.info")

# Leer el header de cualquier .info
sample_info <- file.path(ld_dir, "b1.ldm.full.info")
header_line <- readLines(sample_info, n = 1)
# Convertir espacios a tabs y agregar columna Block
header_fields <- strsplit(trimws(header_line), "\\s+")[[1]]
header_out <- paste(c(header_fields, "Block"), collapse = "\t")

# Escribir header
writeLines(header_out, snp_info_file)

# Concatenar datos bloque por bloque en el orden de dt_final
global_idx <- 0L
for (i in seq_len(nrow(dt_final))) {
    blk <- as.character(dt_final$Block[i])
    info_file <- file.path(ld_dir, paste0("b", blk, ".ldm.full.info"))

    if (!file.exists(info_file)) {
        cat(sprintf("   Falta b%s.ldm.full.info\n", blk))
        next
    }

    # Leer datos (sin header)
    dt_block <- fread(info_file, skip = 1, header = FALSE)
    n_vars <- nrow(dt_block)

    # Reemplazar columna Index (V8) con índice global
    dt_block$V8 <- seq(global_idx, global_idx + n_vars - 1L)

    # Agregar columna Block
    dt_block$Block <- as.integer(blk)

    # Append al archivo
    fwrite(dt_block, snp_info_file, sep = "\t",
           col.names = FALSE, append = TRUE)

    global_idx <- global_idx + n_vars

    if (i %% 100 == 0) {
        cat(sprintf("    Progreso: %d/%d bloques (%d variantes)\n",
            i, nrow(dt_final), global_idx))
    }
}

cat(sprintf("  snp.info: %d variantes, %d bloques\n",
    global_idx, nrow(dt_final)))

# Verificación de consistencia
snp_count <- as.integer(system(sprintf("wc -l < %s", snp_info_file), intern = TRUE)) - 1L
ldm_count <- sum(dt_final$NumSnps)
if (snp_count != ldm_count) {
    cat(sprintf("  ICONSISTENCIA: snp.info=%d vs ldm.info=%d\n",
        snp_count, ldm_count))
} else {
    cat(sprintf(" Consistencia verificada: %d variantes\n", snp_count))
}

# Verificar un bloque específico (b529)
b529_row <- dt_final[dt_final$Block == 529, ]
if (nrow(b529_row) > 0) {
    b529_start <- b529_row$StartSnpIdx
    # Leer la línea correspondiente de snp.info (start + 2 por header)
    test_line <- system(sprintf("awk 'NR==%d' %s", b529_start + 2, snp_info_file),
                        intern = TRUE)
    cat(sprintf("  Verificación b529 (idx=%d): %s\n",
        b529_start, substr(test_line, 1, 80)))
}

# ================================================================
# PASO 4: SBayesRC MCMC
# ================================================================
cat("\n=== PASO 4: SBayesR MCMC ===\n")

out_prefix <- file.path(out_dir, "sbayesrc")
cat(sprintf("  mafile:    %s\n", ma_file))
cat(sprintf("  LDdir:     %s\n", ld_dir))
cat(sprintf("  outPrefix: %s\n", out_prefix))
cat(sprintf("  annot:     (ninguna)\n"))
cat(sprintf("  Bloques:   %d (de 591 originales)\n", length(final_blocks)))

tryCatch({
    SBayesRC::sbayesrc(
        mafile    = ma_file,
        LDdir     = ld_dir,
        outPrefix = out_prefix,
        annot     = "",
        thresh    = 0.9,
        bTune     = FALSE
    )
    cat("\n SBayesR MCMC completado\n")
}, error = function(e) {
    cat(sprintf("\n Error en SBayesR: %s\n", e$message))
    stop(e)
})

# ================================================================
# PASO 5: Exportar betas para scoring
# ================================================================
cat("\n=== PASO 5: Exportando betas para scoring ===\n")

snp_res <- paste0(out_prefix, ".snpRes")
if (!file.exists(snp_res)) {
    snp_files <- list.files(out_dir, pattern = "\\.snpRes$", full.names = TRUE)
    if (length(snp_files) > 0) {
        snp_res <- snp_files[1]
        cat(sprintf("  Usando: %s\n", basename(snp_res)))
    }
}

if (file.exists(snp_res)) {
    df_snp <- fread(snp_res)

    cat(sprintf("  Variantes totales: %d\n", nrow(df_snp)))
    cat(sprintf("  Variantes con efecto != 0: %d\n", sum(df_snp$A1Effect != 0)))

    out_weights <- file.path(dirname(out_dir), "sbayesrc_weights_for_scoring.tsv")
    fwrite(
        df_snp[, .(SNP = Id, A1, A2, BETA = A1Effect, FREQ = A1Frq)],
        out_weights, sep = "\t"
    )
    cat(sprintf("  Pesos: %s\n", out_weights))

    # Estadísticas
    nonzero <- sum(df_snp$A1Effect != 0)
    cat(sprintf("  Non-zero: %d (%.1f%%)\n", nonzero, nonzero/nrow(df_snp)*100))
    cat(sprintf("  |beta| range: [%.2e, %.2e]\n",
        min(abs(df_snp$A1Effect[df_snp$A1Effect != 0])),
        max(abs(df_snp$A1Effect))))

    par_res <- paste0(out_prefix, ".parRes")
    if (file.exists(par_res)) {
        cat("\n  Parámetros estimados:\n")
        print(fread(par_res))
    }
} else {
    cat("No se encontró .snpRes\n")
    cat(sprintf("  Archivos en %s:\n", out_dir))
    cat(paste("   ", list.files(out_dir), collapse = "\n"), "\n")
}

cat(sprintf("\n=== Pipeline A completado: %s ===\n", Sys.time()))
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
echo "  PIPELINE A — SBayesRC v6 COMPLETADO"
echo "  Siguiente: sbrc_03_scoring_gpu.sbatch"
echo "============================================================"
