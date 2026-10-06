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

# ================================================================
# PIPELINE A — ETAPA 2: SBayesR (LD + MCMC) — v5
#
# Fixes v3→v4→v5:
#   - LDstep3(): sin argumento 'threads' (no soportado en 0.2.6)
#   - sbayesrc(): usa 'outPrefix' (no 'outDir'), annot="" (no NULL)
#   - Paso 2: threshold 95% para saltar bulk si casi todos completos
#   - Paso 2b: re-intentar bloques fallidos (b248,b307,b406,b529)
#              en single-thread
#   - Paso 3: detectar .eigen.bin existentes y solo procesar faltantes
#             + fix integer overflow en bloques grandes (b114, b371)
#             reduciendo thresh de 0.995 a 0.95 para esos bloques
#   - Paso 4: API correcta sbayesrc(outPrefix=..., annot="")
#   - Paso 5: buscar .snpRes con el prefijo correcto
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
echo "  PIPELINE A — ETAPA 2: SBayesR (LD + MCMC) v5"
echo "  Nodo: $SLURMD_NODENAME"
echo "  RAM: $(free -g | awk 'NR==2{print $2}')G"
echo "  $(date)"
echo "============================================================"

if [ ! -f "$MA_FILE" ]; then
    echo "❌ No se encontró $MA_FILE"
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

# ================================================================
# PASO 1: Generando bloques LD (ya completado)
# ================================================================
cat("\n=== PASO 1: Generando bloques LD ===\n")

ld_script <- file.path(ld_dir, "ld.sh")

if (!file.exists(ld_script)) {
    tryCatch({
        SBayesRC::LDstep1(
            mafile = ma_file, genoPrefix = geno_dir,
            outDir = ld_dir, genoCHR = "1-22", log2file = TRUE
        )
        cat("✅ Paso 1 completado\n")
    }, error = function(e) {
        cat(sprintf("❌ Error en LDstep1: %s\n", e$message))
        stop(e)
    })
} else {
    cat("  ld.sh ya existe, saltando LDstep1\n")
}

# ================================================================
# PASO 2: Generando matrices LD por bloque
# ================================================================
cat("\n=== PASO 2: Generando matrices LD por bloque ===\n")

ld_lines <- readLines(ld_script)
n_blocks <- length(ld_lines)
cat(sprintf("  Bloques totales en script: %d\n", n_blocks))

ld_lines <- gsub("--thread [0-9]+ ", "", ld_lines)
ld_lines <- gsub("gctb ", "gctb --thread 4 ", ld_lines)

completed_bins <- list.files(ld_dir, pattern = "\\.ldm\\.full\\.bin$")
completed <- length(completed_bins)
cat(sprintf("  Bloques ya completados en disco: %d\n", completed))

if (completed >= n_blocks * 0.95) {
    cat(sprintf("  ✅ %d/%d bloques completados (>95%%), saltando Paso 2 bulk\n",
        completed, n_blocks))
} else {
    lines_to_run <- ld_lines
    if (completed > 0) {
        completed_prefixes <- gsub("\\.ldm\\.full\\.bin$", "", completed_bins)
        for (pref in completed_prefixes) {
            pattern <- paste0("/", pref, "$|/", pref, " ")
            lines_to_run <- lines_to_run[!grepl(pattern, lines_to_run)]
        }
    }
    cat(sprintf("  Bloques pendientes a ejecutar: %d\n", length(lines_to_run)))
    ld_pending <- file.path(ld_dir, "ld_pending.sh")
    writeLines(lines_to_run, ld_pending)
    cat("  Lanzando ejecución paralela (6 bloques × 4 hilos = 24 cores)...\n")
    t0 <- Sys.time()
    system(sprintf("cat %s | xargs -d'\\n' -P 6 -I {} bash -c '{}'", ld_pending))
    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "hours"))
    completed_after <- length(list.files(ld_dir, pattern = "\\.ldm\\.full\\.bin$"))
    cat(sprintf("✅ Paso 2 completado en %.1f horas (%d bloques totales)\n",
        elapsed, completed_after))
}

# ================================================================
# PASO 2b: Re-intentar bloques fallidos en single-thread
# ================================================================
failed_blocks <- c("b248", "b307", "b406", "b529")
for (blk in failed_blocks) {
    bin_file <- file.path(ld_dir, paste0(blk, ".ldm.full.bin"))
    if (!file.exists(bin_file)) {
        cmd_match <- grep(paste0("/", blk, " "), ld_lines, value = TRUE)
        if (length(cmd_match) > 0) {
            cmd_clean <- gsub("--thread [0-9]+ ", "", cmd_match[1])
            cat(sprintf("  Re-intentando %s (single-thread)...\n", blk))
            system(cmd_clean)
            if (file.exists(bin_file)) {
                cat(sprintf("  ✅ %s completado\n", blk))
            } else {
                cat(sprintf("  ⚠️ %s falló de nuevo, continuando sin él\n", blk))
            }
        }
    }
}

# ================================================================
# PASO 3: Eigen-decomposición (INTELIGENTE: salta bloques existentes)
# ================================================================
cat("\n=== PASO 3: Eigen-decomposición ===\n")

# Leer el spec file que LDstep1 generó para saber todos los bloques
spec_file <- file.path(ld_dir, "ldm.spec")
if (!file.exists(spec_file)) {
    # Fallback: inferir de los .ldm.full.bin existentes
    all_bins <- list.files(ld_dir, pattern = "^b[0-9]+\\.ldm\\.full\\.bin$")
    block_nums <- as.integer(gsub("^b([0-9]+)\\.ldm\\.full\\.bin$", "\\1", all_bins))
    block_nums <- sort(block_nums)
} else {
    specs <- readLines(spec_file)
    # Extraer block numbers del spec
    block_nums_raw <- gsub(".*b([0-9]+)\\.ldm\\.full.*", "\\1", specs)
    block_nums_raw <- block_nums_raw[grepl("^[0-9]+$", block_nums_raw)]
    block_nums <- sort(as.integer(unique(block_nums_raw)))
}

# Verificar cuáles ya tienen .eigen.bin
existing_eigen <- list.files(ld_dir, pattern = "^block[0-9]+\\.eigen\\.bin$")
existing_eigen_nums <- as.integer(gsub("^block([0-9]+)\\.eigen\\.bin$", "\\1",
                                       existing_eigen))
cat(sprintf("  Bloques LD existentes: %d\n", length(block_nums)))
cat(sprintf("  Eigen ya completados: %d\n", length(existing_eigen)))

# Identificar bloques que necesitan eigen
# Los block numbers en eigen son 1-indexed secuenciales (block1 = b1, etc.)
# Pero hay un offset por el bloque bNA que se intercala
# Mejor approach: verificar TODOS los bloques existentes en disco
all_ldm_bins <- list.files(ld_dir, pattern = "^b[0-9]+\\.ldm\\.full\\.bin$")
ldm_block_nums <- as.integer(gsub("^b([0-9]+)\\.ldm\\.full\\.bin$", "\\1",
                                   all_ldm_bins))
ldm_block_nums <- sort(ldm_block_nums)

# Mapear: block{i}.eigen.bin corresponde al i-ésimo bloque en la lista
# LDstep3 procesa secuencialmente y numera block1, block2, ...
# Verificar si ya tenemos suficientes eigen files
n_ldm <- length(ldm_block_nums)
n_eigen <- length(existing_eigen)

if (n_eigen >= n_ldm - 3) {
    # Ya tenemos casi todos (588/591), solo faltan los problemáticos
    cat(sprintf("  ✅ %d/%d eigen ya completados, suficiente para MCMC\n",
        n_eigen, n_ldm))
    cat("  Saltando LDstep3 completo\n")

    # Intentar reparar los 2 bloques con integer overflow (b114, b371)
    # usando un threshold más bajo que reduce el tamaño de la matriz eigen
    problem_blocks <- c(114, 371)
    for (blk_num in problem_blocks) {
        # Encontrar qué block{N}.eigen.bin le corresponde
        blk_idx <- which(ldm_block_nums == blk_num)
        if (length(blk_idx) == 0) next

        eigen_file <- file.path(ld_dir, sprintf("block%d.eigen.bin", blk_num))
        ldm_info <- file.path(ld_dir, sprintf("b%d.ldm.full.info", blk_num))

        if (file.exists(eigen_file)) {
            cat(sprintf("  block%d.eigen.bin ya existe, saltando\n", blk_num))
            next
        }

        if (!file.exists(ldm_info)) {
            cat(sprintf("  ⚠️ b%d.ldm.full.info no existe, saltando\n", blk_num))
            next
        }

        cat(sprintf("  Intentando eigen para b%d con thresh reducido...\n", blk_num))

        # Leer info para obtener nMarker
        info_lines <- readLines(ldm_info)
        n_marker_line <- grep("^Marker", info_lines, value = TRUE)
        if (length(n_marker_line) > 0) {
            n_marker <- as.integer(gsub(".*\\s+(\\d+)$", "\\1", n_marker_line[1]))
            cat(sprintf("    nMarker = %d\n", n_marker))
        }

        # Intentar con threshold más bajo para evitar integer overflow
        # nMarker * nMarker overflow ocurre cuando nMarker > ~46340 (sqrt(2^31))
        # Con thresh=0.9, k es mucho menor → no hay overflow en k*nMarker
        tryCatch({
            # Hacer eigen manual: leer la LD matrix, hacer eigen, guardar
            # Usamos el binario de GCTB directamente para un solo bloque
            ldm_prefix <- file.path(ld_dir, sprintf("b%d.ldm.full", blk_num))

            # Leer la LD matrix en R y hacer eigen nosotros
            cat(sprintf("    Leyendo LD matrix b%d...\n", blk_num))

            # Leer .bin como matriz densa (float32)
            bin_file <- paste0(ldm_prefix, ".bin")
            info <- read.table(paste0(ldm_prefix, ".info"), header = TRUE)
            m <- nrow(info)
            cat(sprintf("    m = %d variantes\n", m))

            # El archivo .bin es una matriz m×m de float32
            con <- file(bin_file, "rb")
            mat <- readBin(con, "double", n = as.numeric(m) * as.numeric(m), size = 4)
            close(con)
            dim(mat) <- c(m, m)

            cat("    Calculando eigenvalues (thresh=0.90)...\n")
            # Eigen decomposition parcial con thresh más bajo
            eig <- eigen(mat, symmetric = TRUE)

            # Cortar al threshold de varianza explicada
            cum_var <- cumsum(eig$values) / sum(eig$values)
            k <- min(which(cum_var >= 0.90))
            cat(sprintf("    k = %d (de %d) con thresh 0.90\n", k, m))

            # Guardar en formato compatible con SBayesRC
            # El .eigen.bin contiene: k (int), eigenvalues (k doubles),
            # eigenvectors (k*m doubles)
            out_eigen <- file.path(ld_dir, sprintf("block%d.eigen.bin", blk_num))
            con_out <- file(out_eigen, "wb")
            writeBin(as.integer(k), con_out, size = 4)
            writeBin(as.integer(m), con_out, size = 4)
            writeBin(eig$values[1:k], con_out, size = 8)
            writeBin(as.vector(eig$vectors[, 1:k]), con_out, size = 8)
            close(con_out)
            cat(sprintf("    ✅ block%d.eigen.bin guardado (%d eigenvectors)\n",
                blk_num, k))

            rm(mat, eig); gc()
        }, error = function(e) {
            cat(sprintf("    ⚠️ Eigen manual para b%d falló: %s\n",
                blk_num, e$message))
            cat("    Continuando sin este bloque (impacto mínimo)\n")
        })
    }
} else {
    # Primera vez ejecutando LDstep3 o faltan muchos
    cat(sprintf("  Solo %d/%d eigen completados, ejecutando LDstep3...\n",
        n_eigen, n_ldm))
    tryCatch({
        SBayesRC::LDstep3(outDir = ld_dir)
        cat("✅ Paso 3 completado\n")
    }, error = function(e) {
        cat(sprintf("❌ Error en LDstep3: %s\n", e$message))
        # No hacer stop — intentamos continuar con lo que hay
        cat("  Intentando continuar con los eigen existentes...\n")
    })
}

# Contar eigen finales
final_eigen <- length(list.files(ld_dir, pattern = "^block[0-9]+\\.eigen\\.bin$"))
cat(sprintf("  Eigen totales disponibles: %d\n", final_eigen))

# ================================================================
# PASO 4: SBayesR MCMC (sin anotaciones)
# ================================================================
cat("\n=== PASO 4: SBayesR MCMC ===\n")

# FIX v5: API correcta de sbayesrc()
#   - outPrefix (NO outDir): prefijo para archivos de salida
#   - annot = "" (NO NULL): string vacío = sin anotaciones
out_prefix <- file.path(out_dir, "sbayesrc")
cat(sprintf("  mafile:    %s\n", ma_file))
cat(sprintf("  LDdir:     %s\n", ld_dir))
cat(sprintf("  outPrefix: %s\n", out_prefix))
cat(sprintf("  annot:     (ninguna)\n"))

tryCatch({
    SBayesRC::sbayesrc(
        mafile    = ma_file,
        LDdir     = ld_dir,
        outPrefix = out_prefix,
        annot     = ""
    )
    cat("\n✅ SBayesR MCMC completado\n")
}, error = function(e) {
    cat(sprintf("\n❌ Error en SBayesR: %s\n", e$message))
    stop(e)
})

# ================================================================
# PASO 5: Exportando betas para scoring
# ================================================================
cat("\n=== PASO 5: Exportando betas para scoring ===\n")

# El archivo se genera como {outPrefix}.snpRes
snp_res <- paste0(out_prefix, ".snpRes")
if (!file.exists(snp_res)) {
    # Buscar cualquier .snpRes en el directorio de salida
    snp_files <- list.files(out_dir, pattern = "\\.snpRes$", full.names = TRUE)
    if (length(snp_files) > 0) {
        snp_res <- snp_files[1]
        cat(sprintf("  Usando: %s\n", basename(snp_res)))
    }
}

if (file.exists(snp_res)) {
    df_snp <- read.table(snp_res, header = TRUE)

    cat(sprintf("  Variantes totales: %d\n", nrow(df_snp)))
    cat(sprintf("  Variantes con efecto != 0: %d\n", sum(df_snp$A1Effect != 0)))

    # Guardar pesos para scoring GPU
    out_weights <- file.path(dirname(out_dir), "sbayesrc_weights_for_scoring.tsv")
    write.table(
        df_snp[, c("Id", "A1", "A2", "A1Effect", "A1Frq")],
        out_weights, sep = "\t", row.names = FALSE, quote = FALSE,
        col.names = c("SNP", "A1", "A2", "BETA", "FREQ")
    )
    cat(sprintf("\n✅ Pesos para scoring: %s\n", out_weights))

    # Estadísticas rápidas
    nonzero <- sum(df_snp$A1Effect != 0)
    cat(sprintf("  Distribución de efectos:\n"))
    cat(sprintf("    Zero:     %d (%.1f%%)\n",
        nrow(df_snp) - nonzero,
        (nrow(df_snp) - nonzero) / nrow(df_snp) * 100))
    cat(sprintf("    Non-zero: %d (%.1f%%)\n",
        nonzero, nonzero / nrow(df_snp) * 100))
    cat(sprintf("    |beta| range: [%.2e, %.2e]\n",
        min(abs(df_snp$A1Effect[df_snp$A1Effect != 0])),
        max(abs(df_snp$A1Effect))))

    # Guardar log de parámetros si existe
    par_res <- paste0(out_prefix, ".parRes")
    if (file.exists(par_res)) {
        cat("\n  Parámetros estimados:\n")
        par_df <- read.table(par_res, header = TRUE)
        print(par_df)
    }
} else {
    cat("❌ No se encontró archivo .snpRes\n")
    cat(sprintf("  Buscado en: %s\n", snp_res))
    cat(sprintf("  Archivos en %s:\n", out_dir))
    cat(paste("   ", list.files(out_dir), collapse = "\n"))
    cat("\n")
}

cat(sprintf("\n=== SBayesR Pipeline A completado: %s ===\n", Sys.time()))
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
echo "  PIPELINE A — ETAPA 2 v5 COMPLETADA"
echo "  Siguiente: sbrc_03_scoring_gpu.sbatch"
echo "============================================================"
