#!/bin/bash
# Diagnóstico bigsparser 0.7.3 — v3: evita crash en externalptr
# Uso: bash test_sfbm_reattach.sh

SIF="/mnt/cephfs/biocontainers/images/jupyter-biotools-1.4.sif"
PREV_TMPDIR="/mnt/cephfs/scratch/angel.pacheco/ldpred2_hm3_tmp/2931"
TEST_SBK="$PREV_TMPDIR/corr_chr22.sbk"

R_SCRIPT=$(mktemp /tmp/test_reattach_XXXX.R)
cat << 'REOF' > "$R_SCRIPT"
suppressMessages({ library(Matrix); library(bigsparser) })
cat(sprintf("bigsparser: %s\n\n", as.character(packageVersion("bigsparser"))))

safe_print_field <- function(f, val) {
  tryCatch({
    if (inherits(val, "externalptr")) {
      cat(sprintf("    %s: <externalptr>\n", f))
    } else if (is.function(val) || is.environment(val)) {
      cat(sprintf("    %s: <function/env>\n", f))
    } else if (length(val) == 0) {
      cat(sprintf("    %s: (empty, class=%s)\n", f, class(val)[1]))
    } else if (length(val) > 10) {
      cat(sprintf("    %s (len=%d, class=%s): [%s ...]\n",
          f, length(val), class(val)[1], paste(head(val, 6), collapse=",")))
    } else {
      cat(sprintf("    %s (class=%s): %s\n", f, class(val)[1], paste(val, collapse=" ")))
    }
  }, error = function(e) cat(sprintf("    %s: ERROR(%s)\n", f, conditionMessage(e))))
}

# ── 1) Crear SFBM compact de prueba (50x50) ──────────────────────────────────
cat("=== Crear SFBM_compact de prueba ===\n")
set.seed(42)
n <- 50
M <- matrix(rnorm(n*n), n, n); M <- M %*% t(M); diag(M) <- diag(M) + 1
M[abs(M) < 4] <- 0
spm <- as(Matrix(M, sparse=TRUE), "symmetricMatrix")
tmp_bk <- tempfile("test_sfbm")
sfbm <- as_SFBM(spm, backingfile=tmp_bk, compact=TRUE)
cat(sprintf("  clase: %s\n", paste(class(sfbm), collapse="/")))

# Campos del objeto creado
cat("  campos:\n")
fld_names <- names(sfbm$getRefClass()$fields())
for (f in fld_names) {
  val <- tryCatch(sfbm[[f]], error=function(e) paste("ERR:", conditionMessage(e)))
  safe_print_field(f, val)
}

sbk_path  <- paste0(tmp_bk, ".sbk")
sbk_size  <- file.info(sbk_path)$size
cat(sprintf("  .sbk size: %d bytes\n", sbk_size))
cat(sprintf("  nrow=%d, ncol=%d, p len=%d\n", sfbm$nrow, sfbm$ncol, length(sfbm$p)))
cat(sprintf("  p[1:6]: %s\n", paste(head(sfbm$p, 6), collapse=",")))

# ── 2) Mismo con SFBM_corr_compact_RC directamente ───────────────────────────
cat("\n=== Crear via SFBM_corr_compact_RC ===\n")
gen_cc <- get("SFBM_corr_compact_RC", envir=asNamespace("bigsparser"))
sfbm_cc <- tryCatch(
  gen_cc$new(spm=spm, backingfile=tempfile("test_cc")),
  error=function(e) { cat(sprintf("  FAIL gen_cc$new: %s\n", conditionMessage(e))); NULL }
)
if (!is.null(sfbm_cc)) {
  cat(sprintf("  clase: %s\n", paste(class(sfbm_cc), collapse="/")))
  cat("  campos:\n")
  for (f in names(sfbm_cc$getRefClass()$fields())) {
    val <- tryCatch(sfbm_cc[[f]], error=function(e) paste("ERR:", conditionMessage(e)))
    safe_print_field(f, val)
  }
  cat(sprintf("  first_i len=%d, head: %s\n",
      length(sfbm_cc$first_i), paste(head(sfbm_cc$first_i, 8), collapse=",")))
}

# ── 3) Tamaño del .sbk real (chr22) ──────────────────────────────────────────
real_sbk <- commandArgs(trailingOnly=TRUE)[1]
cat(sprintf("\n=== .sbk real: %s ===\n", basename(real_sbk)))
if (file.exists(real_sbk)) {
  real_size <- file.info(real_sbk)$size
  ncol_chr22 <- 14360L
  # p tiene ncol+1 enteros de 8 bytes al inicio en algunos formatos
  # Estimamos: si fuera dense-stripe con avg_width doubles por col
  avg_dbl <- (real_size / 8) / ncol_chr22
  cat(sprintf("  size: %d bytes = %d doubles\n", real_size, real_size %/% 8))
  cat(sprintf("  si ncol=%d: avg doubles/col = %.1f\n", ncol_chr22, avg_dbl))
  # nnz reportado en log: 9,933,050 → si .sbk solo guarda values:
  # 9933050 * 8 = 79,464,400 bytes (~79 MB)
  # El archivo es 43,591,550,600 → ~43 GB → incluye indices
  # Estimado con indices int32: values(8) + col_idx(4) = 12 bytes/nnz → ~119 MB
  # 43 GB >> 119 MB → no son solo los nonzeros del LD original
  # bigsparser compact: almacena ld_matrix ENTERA (dense stripes por SNP)
  cat(sprintf("  43 GB sugiere almacenamiento denso por stripe (no solo nnz)\n"))
}

# ── 4) Field-swap: ¿se puede modificar $backingfile y $p? ────────────────────
cat("\n=== Test field-swap ===\n")
# Crear dummy 3x3 para no contaminar sfbm de prueba
dummy_spm <- as(Matrix::Diagonal(3), "symmetricMatrix")
gen_cc2 <- get("SFBM_corr_compact_RC", envir=asNamespace("bigsparser"))
sfbm_d <- tryCatch(
  gen_cc2$new(spm=dummy_spm, backingfile=tempfile("dummy_cc")),
  error=function(e) { cat(sprintf("  FAIL dummy: %s\n", conditionMessage(e))); NULL }
)

if (!is.null(sfbm_d)) {
  cat(sprintf("  dummy: clase=%s, nrow=%d\n", class(sfbm_d)[1], sfbm_d$nrow))

  # Test 1: ¿podemos leer p del dummy?
  cat(sprintf("  dummy$p: %s\n", paste(sfbm_d$p, collapse=",")))
  cat(sprintf("  dummy$first_i: %s\n", paste(sfbm_d$first_i, collapse=",")))

  # Test 2: field-swap backingfile
  tryCatch({
    sfbm_d$backingfile <- real_sbk  # sin .sbk
    cat(sprintf("  backingfile swap: OK → %s\n", sfbm_d$backingfile))
  }, error=function(e) cat(sprintf("  backingfile swap: FAIL: %s\n", conditionMessage(e))))

  # Test 3: field swap nrow/ncol
  tryCatch({
    sfbm_d$nrow <- 14360L
    cat(sprintf("  nrow swap: OK → %d\n", sfbm_d$nrow))
  }, error=function(e) cat(sprintf("  nrow swap: FAIL: %s\n", conditionMessage(e))))

  # Test 4: ¿address activo sin crash? (recalcula extptr del .sbk)
  tryCatch({
    addr <- sfbm_d$address
    cat(sprintf("  address recalc: OK (class=%s)\n", class(addr)[1]))
  }, error=function(e) cat(sprintf("  address recalc: FAIL: %s\n", conditionMessage(e))))
}

# ── 5) ¿Qué hace $initialize de SFBM_compact? ────────────────────────────────
cat("\n=== Source de initialize (SFBM_compact_RC) ===\n")
gen3 <- get("SFBM_compact_RC", envir=asNamespace("bigsparser"))
tryCatch(print(body(gen3$def@refMethods$initialize)),
         error=function(e) cat(sprintf("  FAIL: %s\n", conditionMessage(e))))

cat("\nDiagnóstico v3 completado.\n")
REOF

apptainer exec --bind /mnt/cephfs:/mnt/cephfs \
    "$SIF" /opt/miniforge3/envs/r-bio/bin/Rscript "$R_SCRIPT" "$TEST_SBK"

rm -f "$R_SCRIPT"
