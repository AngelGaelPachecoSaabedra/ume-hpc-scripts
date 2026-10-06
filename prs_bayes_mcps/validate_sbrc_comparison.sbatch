#!/bin/bash
#SBATCH --job-name=val_sbrc_cmp
#SBATCH --account=researchers
#SBATCH --qos=vip
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=80G
#SBATCH --time=06:00:00
#SBATCH --exclude=mcps-hpc-01
#SBATCH --output=/mnt/cephfs/orgs/home/angel.pacheco/logs/val_sbrc_cmp_%j.out
#SBATCH --error=/mnt/cephfs/orgs/home/angel.pacheco/logs/val_sbrc_cmp_%j.err

set -euo pipefail

SIF="/mnt/cephfs/biocontainers/images/jupyter-biotools-1.4.sif"
TMPDIR="/mnt/cephfs/scratch/$USER/val_sbrc_cmp_tmp/$SLURM_JOB_ID"
mkdir -p "$TMPDIR"
chmod 700 "$TMPDIR"
export TMPDIR TMP="$TMPDIR" TEMP="$TMPDIR"
export APPTAINER_TMPDIR="$TMPDIR/apptainer-tmp"

echo "============================================================"
echo "  VALIDACIÓN CLÍNICA COMPARATIVA — SBayesRC"
echo "  Pipeline A (T2DGGI) vs Pipeline A2 (MCPS GWAS)"
echo "  Fenotipo: HbA1c >= 6.5% + Diagnóstico + Meds"
echo "  Nodo: $SLURMD_NODENAME"
echo "  $(date)"
echo "============================================================"

# ================================================================
# RUTAS
# ================================================================
BASE_DIR="/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes"
SBRC_DIR="${BASE_DIR}/sbayesrc"
PLINK_DIR="/mnt/cephfs/hot_nvme/mcps/imputed-topmed/plink_files/maximally_unrelated"
CLIN_FILE="${BASE_DIR}/MCPS BASELINE.csv"
LINK_FILE="${BASE_DIR}/RGN_LINK_IID.csv"
FIG_DIR="${SBRC_DIR}/figures_validation"
SCORE_DIR="${SBRC_DIR}/scoring"

mkdir -p "$FIG_DIR" "$SCORE_DIR"

# Weights
W_T2DGGI="${SBRC_DIR}/sbayesrc_weights_for_scoring.tsv"
W_MCPS="${SBRC_DIR}/results_mcps/sbayesrc_mcps_weights_for_scoring.tsv"

echo ""
echo "--- Pre-flight ---"
for F in "$W_T2DGGI" "$W_MCPS" "$CLIN_FILE" "$LINK_FILE"; do
    if [ -f "$F" ]; then
        echo "  ✓ $(basename "$F")"
    else
        echo "  ✗ FALTA: $F"
        exit 1
    fi
done

# ================================================================
# FASE 1: plink2 --score PARA AMBOS SETS DE WEIGHTS
# ================================================================
echo ""
echo "============================================================"
echo "  FASE 1: Scoring con plink2 --score (chr1..22)"
echo "============================================================"

# Inspeccionar columnas de los weights
echo ""
echo "--- Estructura de weights ---"
echo "  T2DGGI:"
head -1 "$W_T2DGGI" | tr '\t' '\n' | head -5
echo "  Variantes: $(tail -n +2 "$W_T2DGGI" | wc -l)"
echo ""
echo "  MCPS:"
head -1 "$W_MCPS" | tr '\t' '\n' | head -5
echo "  Variantes: $(tail -n +2 "$W_MCPS" | wc -l)"

# Crear script de scoring por cromosoma
SCORE_SCRIPT="$TMPDIR/score_both.sh"
cat << 'SCOREEOF' > "$SCORE_SCRIPT"
#!/bin/bash
set -euo pipefail

PLINK_DIR="$1"
W_FILE="$2"
OUT_PREFIX="$3"
LABEL="$4"
SIF="$5"

echo ""
echo "  >> Scoring $LABEL..."

# Detectar columnas del weights file
HEADER=$(head -1 "$W_FILE")
# Buscar nombre de columna de ID (SNP, Id, Name)
ID_COL=""
for cand in SNP Id Name ID; do
    if echo "$HEADER" | grep -qw "$cand"; then
        ID_COL="$cand"
        break
    fi
done

# Buscar columna de alelo (A1)
A1_COL=""
for cand in A1 EffectAllele; do
    if echo "$HEADER" | grep -qw "$cand"; then
        A1_COL="$cand"
        break
    fi
done

# Buscar columna de efecto (BETA, A1Effect)
BETA_COL=""
for cand in BETA A1Effect beta; do
    if echo "$HEADER" | grep -qw "$cand"; then
        BETA_COL="$cand"
        break
    fi
done

echo "     Columnas detectadas: ID=$ID_COL, A1=$A1_COL, BETA=$BETA_COL"

if [ -z "$ID_COL" ] || [ -z "$A1_COL" ] || [ -z "$BETA_COL" ]; then
    echo "ERROR: No se detectaron las columnas necesarias en $W_FILE"
    echo "  Header: $HEADER"
    exit 1
fi

# Scoring por cromosoma (1-22)
for CHR in $(seq 1 22); do
    PLINK_PREFIX="${PLINK_DIR}/mcps-freeze150k_qcd_chr${CHR}_ivs"
    if [ ! -f "${PLINK_PREFIX}.bed" ]; then
        echo "     chr${CHR}: SKIP (no .bed)"
        continue
    fi

    apptainer exec --bind /mnt/cephfs:/mnt/cephfs "$SIF" \
        plink2 --bfile "$PLINK_PREFIX" \
               --score "$W_FILE" header-read \
               cols=+scoresums \
               --score-col-nums 3 \
               --threads 4 \
               --memory 4000 \
               --out "${OUT_PREFIX}_chr${CHR}" \
               2>/dev/null || {
        echo "     chr${CHR}: plink2 falló, intentando con cols explícitos..."
        # Fallback: especificar columnas por número
        # Col 1 = ID, Col 2 = A1, Col 3 = BETA
        apptainer exec --bind /mnt/cephfs:/mnt/cephfs "$SIF" \
            plink2 --bfile "$PLINK_PREFIX" \
                   --score "$W_FILE" 1 2 3 header-read \
                   cols=+scoresums \
                   --threads 4 \
                   --memory 4000 \
                   --out "${OUT_PREFIX}_chr${CHR}" \
                   2>/dev/null || echo "     chr${CHR}: FALLÓ ambos intentos"
    }
done

# Consolidar scores
echo "     Consolidando scores..."
SCOREEOF
chmod +x "$SCORE_SCRIPT"

# Correr scoring para T2DGGI
bash "$SCORE_SCRIPT" "$PLINK_DIR" "$W_T2DGGI" "${SCORE_DIR}/sbrc_t2dggi" "T2DGGI" "$SIF"

# Correr scoring para MCPS
bash "$SCORE_SCRIPT" "$PLINK_DIR" "$W_MCPS" "${SCORE_DIR}/sbrc_mcps" "MCPS" "$SIF"

# ================================================================
# FASE 2: CONSOLIDAR + VALIDACIÓN CLÍNICA EN R
# ================================================================
echo ""
echo "============================================================"
echo "  FASE 2: Consolidación y Validación Clínica Comparativa"
echo "============================================================"

R_SCRIPT="$TMPDIR/validate_sbrc_comparison.R"
cat << 'EOF' > "$R_SCRIPT"
suppressMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(pROC)
  library(scales)
  library(patchwork)
})

score_dir <- "/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes/sbayesrc/scoring"
base_dir  <- "/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes"
fig_dir   <- "/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes/sbayesrc/figures_validation"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

col_t2d <- "#4DBBD5"   # Azul — T2DGGI
col_mcp <- "#E64B35"   # Rojo — MCPS

# ================================================================
# 1. CONSOLIDAR SCORES POR CROMOSOMA
# ================================================================
cat(">> Fase 1: Consolidando scores plink2 por cromosoma...\n")

consolidate_scores <- function(prefix, label) {
  score_files <- list.files(score_dir, pattern = paste0("^", prefix, "_chr[0-9]+\\.sscore$"),
                            full.names = TRUE)
  if (length(score_files) == 0) {
    # Intentar .sscore.zst
    score_files <- list.files(score_dir, pattern = paste0("^", prefix, "_chr[0-9]+\\.sscore"),
                              full.names = TRUE)
  }

  cat(sprintf("  %s: encontrados %d archivos de score\n", label, length(score_files)))
  if (length(score_files) == 0) stop(sprintf("No se encontraron scores para %s", label))

  # Leer el primer archivo para obtener columnas
  df_first <- fread(score_files[1])
  cat(sprintf("  Columnas: %s\n", paste(names(df_first), collapse = ", ")))

  # Detectar columnas de ID y score
  id_col <- intersect(c("#IID", "IID", "#FID"), names(df_first))[1]
  # SCORE1_SUM o SCORE1_AVG — preferir SUM si existe
  score_col <- intersect(c("SCORE1_SUM", "SCORE1_AVG", "SCORE_SUM", "SCORE_AVG"), names(df_first))[1]
  if (is.na(score_col)) {
    # Buscar cualquier columna con SCORE
    score_col <- grep("SCORE", names(df_first), value = TRUE)[1]
  }
  cat(sprintf("  ID col: %s, Score col: %s\n", id_col, score_col))

  # Consolidar sumando todos los cromosomas
  all_scores <- NULL
  for (f in score_files) {
    df_chr <- fread(f, select = c(id_col, score_col))
    setnames(df_chr, c(id_col, score_col), c("IID", "PRS_CHR"))

    if (is.null(all_scores)) {
      all_scores <- df_chr
    } else {
      all_scores <- merge(all_scores, df_chr, by = "IID", all = TRUE, suffixes = c("", "_new"))
      all_scores[, PRS_CHR := fifelse(is.na(PRS_CHR), 0, PRS_CHR) +
                               fifelse(is.na(PRS_CHR_new), 0, PRS_CHR_new)]
      all_scores[, PRS_CHR_new := NULL]
    }
  }
  setnames(all_scores, "PRS_CHR", paste0("PRS_", label))
  cat(sprintf("  %s: %s muestras consolidadas\n", label, format(nrow(all_scores), big.mark = ",")))
  return(all_scores)
}

# Intentar consolidar; si no hay archivos, buscar si ya existe un TSV previo
tryCatch({
  df_t2d <- consolidate_scores("sbrc_t2dggi", "T2DGGI")
  df_mcp <- consolidate_scores("sbrc_mcps", "MCPS")
  df_prs <- merge(df_t2d, df_mcp, by = "IID", all = TRUE)
}, error = function(e) {
  cat(sprintf("\n⚠ Error consolidando scores: %s\n", e$message))
  cat("  Intentando cargar PRS existente como fallback...\n")

  # Fallback: si ya existe un TSV previo con PRS, usarlo
  fallback <- file.path(base_dir, "prs_diabetes_mcps_gwas_v4_final.tsv")
  if (file.exists(fallback)) {
    df_prev <- fread(fallback)
    cat(sprintf("  Cargado fallback: %s (%s filas)\n", fallback, nrow(df_prev)))
    # Crear PRS_T2DGGI como el PRS existente, PRS_MCPS como NA
    df_prs <<- data.table(IID = df_prev$IID, PRS_T2DGGI = df_prev$PRS_TOTAL, PRS_MCPS = NA_real_)
    warning("Solo se pudo cargar un PRS. Las comparaciones serán limitadas.")
  } else {
    stop("No hay scores ni fallback disponible.")
  }
})

# ================================================================
# 2. LINK IID → PATID Y MERGE CON CLÍNICA
# ================================================================
cat("\n>> Fase 2: Merge con datos clínicos...\n")

link_file <- file.path(base_dir, "RGN_LINK_IID.csv")
clin_file <- file.path(base_dir, "MCPS BASELINE.csv")

df_link <- fread(link_file)
df_clin <- fread(clin_file)

# Cruzar IID -> PATID
df_prs <- merge(df_prs, df_link[, .(IID, PATID)], by = "IID", all.x = TRUE)

# Merge con clínica
cols_clin <- c("PATID", "AGE", "MALE", "WEIGHT", "HEIGHT",
               "BASE_DIABETES", "BASE_HBA1C",
               "DRUG_D1", "DRUG_D2", "DRUG_D3", "DRUG_D4")
cols_avail <- intersect(cols_clin, names(df_clin))
df <- merge(df_prs, df_clin[, ..cols_avail], by = "PATID", all.x = FALSE)

cat(sprintf("  Pacientes cruzados (PRS + clínica): %s\n", format(nrow(df), big.mark = ",")))

# ================================================================
# 3. FENOTIPO (idéntico a validate_prs_clinical_R.sbatch)
# ================================================================
cat(">> Fase 3: Construyendo Fenotipo para Diabetes T2...\n")

for (med in c("DRUG_D1", "DRUG_D2", "DRUG_D3", "DRUG_D4")) {
  if (med %in% names(df)) {
    df[[med]][is.na(df[[med]])] <- 0
  } else {
    df[[med]] <- 0
  }
}

df <- df %>%
  mutate(
    HEIGHT_M = HEIGHT / 100,
    BMI = WEIGHT / (HEIGHT_M^2),
    dx_previo = (BASE_DIABETES == 1) | (DRUG_D1 == 1) | (DRUG_D2 == 1) | (DRUG_D3 == 1) | (DRUG_D4 == 1),
    hba1c_alta = (BASE_DIABETES == 0) & (!is.na(BASE_HBA1C)) & (BASE_HBA1C >= 6.5),
    control_puro = (BASE_DIABETES == 0) & (!is.na(BASE_HBA1C)) & (BASE_HBA1C < 6.5) & !dx_previo,
    T2D_STATUS = case_when(
      dx_previo | hba1c_alta ~ 1,
      control_puro ~ 0,
      TRUE ~ NA_real_
    )
  ) %>%
  filter(!is.na(AGE), !is.na(MALE), !is.na(BMI), !is.na(T2D_STATUS))

# Estandarizar ambos PRS (Z-score)
if (!all(is.na(df$PRS_T2DGGI))) {
  df$PRS_Z_T2DGGI <- as.numeric(scale(df$PRS_T2DGGI))
}
if (!all(is.na(df$PRS_MCPS))) {
  df$PRS_Z_MCPS <- as.numeric(scale(df$PRS_MCPS))
}

# Filtrar que al menos uno no sea NA
df <- df %>% filter(!is.na(PRS_Z_T2DGGI) | !is.na(PRS_Z_MCPS))

n_total    <- nrow(df)
n_cases    <- sum(df$T2D_STATUS == 1)
n_controls <- sum(df$T2D_STATUS == 0)

cat(sprintf("  Pacientes con datos completos: %s\n", format(n_total, big.mark = ",")))
cat(sprintf("  Casos T2D: %s | Controles: %s\n",
    format(n_cases, big.mark = ","), format(n_controls, big.mark = ",")))

# Subset con AMBOS PRS disponibles (para comparaciones directas)
df_both <- df %>% filter(!is.na(PRS_Z_T2DGGI) & !is.na(PRS_Z_MCPS))
cat(sprintf("  Pacientes con AMBOS PRS: %s\n", format(nrow(df_both), big.mark = ",")))

# ================================================================
# 4. REGRESIÓN LOGÍSTICA — TRES MODELOS
# ================================================================
cat("\n>> Fase 4: Regresión Logística...\n")

# Función Pseudo-R² de Nagelkerke
calc_nagelkerke <- function(model, null_model) {
  n <- nobs(model)
  ll_mod <- as.numeric(logLik(model))
  ll_null <- as.numeric(logLik(null_model))
  r2_cs <- 1 - exp(-2/n * (ll_mod - ll_null))
  r2_max <- 1 - exp(2/n * ll_null)
  return((r2_cs / r2_max) * 100)
}

mod_null <- glm(T2D_STATUS ~ 1, data = df_both, family = binomial)
mod_base <- glm(T2D_STATUS ~ AGE + MALE + BMI, data = df_both, family = binomial)
mod_t2d  <- glm(T2D_STATUS ~ AGE + MALE + BMI + PRS_Z_T2DGGI, data = df_both, family = binomial)
mod_mcp  <- glm(T2D_STATUS ~ AGE + MALE + BMI + PRS_Z_MCPS, data = df_both, family = binomial)
mod_full <- glm(T2D_STATUS ~ AGE + MALE + BMI + PRS_Z_T2DGGI + PRS_Z_MCPS, data = df_both, family = binomial)

r2_base <- calc_nagelkerke(mod_base, mod_null)
r2_t2d  <- calc_nagelkerke(mod_t2d, mod_null)
r2_mcp  <- calc_nagelkerke(mod_mcp, mod_null)
r2_full <- calc_nagelkerke(mod_full, mod_null)

cat("\n[!] VARIANZA EXPLICADA (R² de Nagelkerke):\n")
cat(sprintf("    Clínica Base:                  %.2f%%\n", r2_base))
cat(sprintf("    Clínica + PRS T2DGGI:          %.2f%% (+%.2f%%)\n", r2_t2d, r2_t2d - r2_base))
cat(sprintf("    Clínica + PRS MCPS:            %.2f%% (+%.2f%%)\n", r2_mcp, r2_mcp - r2_base))
cat(sprintf("    Clínica + Ambos PRS:           %.2f%% (+%.2f%%)\n", r2_full, r2_full - r2_base))

# OR e IC del PRS en cada modelo
for (mod_name in c("T2DGGI", "MCPS")) {
  mod <- if (mod_name == "T2DGGI") mod_t2d else mod_mcp
  prs_var <- paste0("PRS_Z_", mod_name)
  coefs <- summary(mod)$coefficients
  if (prs_var %in% rownames(coefs)) {
    or_val <- exp(coefs[prs_var, "Estimate"])
    ci_lo  <- exp(coefs[prs_var, "Estimate"] - 1.96 * coefs[prs_var, "Std. Error"])
    ci_hi  <- exp(coefs[prs_var, "Estimate"] + 1.96 * coefs[prs_var, "Std. Error"])
    p_val  <- coefs[prs_var, "Pr(>|z|)"]
    cat(sprintf("\n[!] OR del PRS %s: %.3f (95%% CI: %.3f - %.3f), p = %.2e\n",
                mod_name, or_val, ci_lo, ci_hi, p_val))
  }
}

cat("\n--- Modelo Completo (Ambos PRS) ---\n")
print(summary(mod_full))

# ================================================================
# GRÁFICA 1: CURVAS ROC COMPARATIVAS
# ================================================================
cat("\n>> Generando gráficas de validación...\n")
cat("   Gráfica 1: Curvas ROC...\n")

roc_base <- roc(df_both$T2D_STATUS, fitted(mod_base), quiet = TRUE)
roc_t2d  <- roc(df_both$T2D_STATUS, fitted(mod_t2d),  quiet = TRUE)
roc_mcp  <- roc(df_both$T2D_STATUS, fitted(mod_mcp),  quiet = TRUE)
roc_full <- roc(df_both$T2D_STATUS, fitted(mod_full),  quiet = TRUE)

auc_base <- sprintf("%.3f", auc(roc_base))
auc_t2d  <- sprintf("%.3f", auc(roc_t2d))
auc_mcp  <- sprintf("%.3f", auc(roc_mcp))
auc_full <- sprintf("%.3f", auc(roc_full))

cat(sprintf("   AUC Base:     %s\n", auc_base))
cat(sprintf("   AUC +T2DGGI:  %s\n", auc_t2d))
cat(sprintf("   AUC +MCPS:    %s\n", auc_mcp))
cat(sprintf("   AUC +Ambos:   %s\n", auc_full))

# DeLong test
delong_test <- tryCatch({
  roc.test(roc_t2d, roc_mcp, method = "delong")
}, error = function(e) NULL)

delong_p <- if (!is.null(delong_test)) sprintf("p = %.4f", delong_test$p.value) else "N/A"
cat(sprintf("   DeLong test T2DGGI vs MCPS: %s\n", delong_p))

df_roc <- rbind(
  data.frame(Spec = roc_base$specificities, Sens = roc_base$sensitivities,
             Model = sprintf("Clínica Base (AUC = %s)", auc_base)),
  data.frame(Spec = roc_t2d$specificities, Sens = roc_t2d$sensitivities,
             Model = sprintf("+ PRS T2DGGI (AUC = %s)", auc_t2d)),
  data.frame(Spec = roc_mcp$specificities, Sens = roc_mcp$sensitivities,
             Model = sprintf("+ PRS MCPS (AUC = %s)", auc_mcp)),
  data.frame(Spec = roc_full$specificities, Sens = roc_full$sensitivities,
             Model = sprintf("+ Ambos PRS (AUC = %s)", auc_full))
)
df_roc$Model <- factor(df_roc$Model, levels = unique(df_roc$Model))

p_roc <- ggplot(df_roc, aes(x = 1 - Spec, y = Sens, color = Model)) +
  geom_line(linewidth = 1.1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray60") +

  scale_color_manual(values = c("gray50", col_t2d, col_mcp, "#7B2D8E")) +

  labs(
    title = "Curvas ROC: Predicción de Diabetes Tipo 2",
    subtitle = sprintf("Cohorte MCPS (n = %s) | Fenotipo de Oro: HbA1c ≥ 6.5%% + Dx + Meds | DeLong: %s",
      format(nrow(df_both), big.mark = ","), delong_p),
    x = "Tasa de Falsos Positivos (1 - Especificidad)",
    y = "Tasa de Verdaderos Positivos (Sensibilidad)"
  ) +

  theme_minimal(base_size = 14) +
  theme(
    plot.title    = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(color = "gray40", size = 10),
    plot.margin   = margin(15, 20, 15, 15),
    legend.position = c(0.72, 0.25),
    legend.background = element_rect(fill = "white", color = "gray80", linewidth = 0.3),
    legend.title = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(fig_dir, "01_roc_comparison.png"), p_roc, width = 9, height = 8, dpi = 300)
ggsave(file.path(fig_dir, "01_roc_comparison.pdf"), p_roc, width = 9, height = 8)

# ================================================================
# GRÁFICA 2: DECILES DE RIESGO (lado a lado)
# ================================================================
cat("   Gráfica 2: Deciles de riesgo...\n")

df_both <- df_both %>%
  mutate(
    Decile_T2DGGI = ntile(PRS_Z_T2DGGI, 10),
    Decile_MCPS   = ntile(PRS_Z_MCPS, 10)
  )

dec_t2d <- df_both %>%
  group_by(Decile = Decile_T2DGGI) %>%
  summarise(Prevalence = mean(T2D_STATUS) * 100, n = n(), .groups = "drop") %>%
  mutate(GWAS = "T2DGGI")

dec_mcp <- df_both %>%
  group_by(Decile = Decile_MCPS) %>%
  summarise(Prevalence = mean(T2D_STATUS) * 100, n = n(), .groups = "drop") %>%
  mutate(GWAS = "MCPS")

df_dec <- bind_rows(dec_t2d, dec_mcp)
df_dec$GWAS <- factor(df_dec$GWAS, levels = c("T2DGGI", "MCPS"))

# Fold-enrichment (decil 10 / decil 1)
fold_t2d <- dec_t2d$Prevalence[dec_t2d$Decile == 10] / dec_t2d$Prevalence[dec_t2d$Decile == 1]
fold_mcp <- dec_mcp$Prevalence[dec_mcp$Decile == 10] / dec_mcp$Prevalence[dec_mcp$Decile == 1]
cat(sprintf("   Fold-enrichment D10/D1: T2DGGI=%.2fx, MCPS=%.2fx\n", fold_t2d, fold_mcp))

p_dec <- ggplot(df_dec, aes(x = factor(Decile), y = Prevalence, fill = GWAS)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8),
           width = 0.7, color = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.1f%%", Prevalence)),
            position = position_dodge(width = 0.8), vjust = -0.5, size = 3, fontface = "bold") +

  scale_fill_manual(values = c("T2DGGI" = col_t2d, "MCPS" = col_mcp), name = "PRS (SBayesRC)") +

  labs(
    title = "Prevalencia de Diabetes T2 por Decil de Riesgo Genético",
    subtitle = sprintf("Decil 1 = menor riesgo, Decil 10 = mayor | Fold D10/D1: T2DGGI=%.1fx, MCPS=%.1fx",
      fold_t2d, fold_mcp),
    x = "Decil de Riesgo Poligénico",
    y = "Prevalencia de Diabetes (%)",
    caption = "Fenotipo de Oro: HbA1c ≥ 6.5% | Diagnóstico previo | Medicación antidiabética"
  ) +

  theme_minimal(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(color = "gray40", size = 11),
    plot.caption  = element_text(color = "gray50", size = 9),
    plot.margin   = margin(15, 20, 15, 15),
    legend.position = c(0.15, 0.88),
    legend.background = element_rect(fill = "white", color = "gray80", linewidth = 0.3),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12)))

ggsave(file.path(fig_dir, "02_deciles_comparison.png"), p_dec, width = 12, height = 7, dpi = 300)
ggsave(file.path(fig_dir, "02_deciles_comparison.pdf"), p_dec, width = 12, height = 7)

# ================================================================
# GRÁFICA 3: DENSIDAD Casos vs Controles (facetado por PRS)
# ================================================================
cat("   Gráfica 3: Densidades Casos vs Controles...\n")

df_both$Status_Label <- ifelse(df_both$T2D_STATUS == 1, "Casos (Diabetes)", "Controles (Sanos)")

df_dens <- rbind(
  data.frame(PRS_Z = df_both$PRS_Z_T2DGGI, Status = df_both$Status_Label, PRS_Type = "PRS T2DGGI (Consorcio Global)"),
  data.frame(PRS_Z = df_both$PRS_Z_MCPS,   Status = df_both$Status_Label, PRS_Type = "PRS MCPS (GWAS Propio)")
)

p_dens <- ggplot(df_dens, aes(x = PRS_Z, fill = Status, color = Status)) +
  geom_density(alpha = 0.35, linewidth = 0.8) +
  facet_wrap(~PRS_Type, ncol = 1, scales = "free_y") +

  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +

  scale_fill_manual(values = c("Casos (Diabetes)" = "#E64B35", "Controles (Sanos)" = "#4DBBD5")) +
  scale_color_manual(values = c("Casos (Diabetes)" = "#E64B35", "Controles (Sanos)" = "#4DBBD5")) +

  labs(
    title = "Distribución del Riesgo Genético: Casos vs Controles",
    subtitle = "Comparación de discriminación entre PRS con GWAS global vs GWAS MCPS",
    x = "PRS Estandarizado (Z-Score)",
    y = "Densidad"
  ) +

  theme_minimal(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(color = "gray40", size = 11),
    plot.margin   = margin(15, 20, 15, 15),
    legend.position = "top",
    legend.title = element_blank(),
    strip.text = element_text(face = "bold", size = 12),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(fig_dir, "03_density_cases_controls.png"), p_dens, width = 10, height = 9, dpi = 300)
ggsave(file.path(fig_dir, "03_density_cases_controls.pdf"), p_dens, width = 10, height = 9)

# ================================================================
# GRÁFICA 4: CALIBRACIÓN COMPARATIVA
# ================================================================
cat("   Gráfica 4: Calibración...\n")

calib_fn <- function(model, label) {
  pred <- predict(model, type = "response")
  dec  <- ntile(pred, 10)
  data.frame(
    pred_decile = 1:10,
    mean_pred = tapply(pred, dec, mean),
    obs_prob  = tapply(df_both$T2D_STATUS, dec, mean),
    se_obs    = tapply(df_both$T2D_STATUS, dec, function(x) sqrt(mean(x)*(1-mean(x))/length(x))),
    Model = label
  )
}

df_calib <- rbind(
  calib_fn(mod_t2d, "Clínica + PRS T2DGGI"),
  calib_fn(mod_mcp, "Clínica + PRS MCPS")
)

p_calib <- ggplot(df_calib, aes(x = mean_pred, y = obs_prob, color = Model)) +
  geom_point(size = 3.5) +
  geom_errorbar(aes(ymin = obs_prob - 1.96*se_obs, ymax = obs_prob + 1.96*se_obs),
                width = 0.008, linewidth = 0.7) +
  geom_line(linewidth = 0.6, alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray30", linewidth = 0.8) +

  scale_color_manual(values = c("Clínica + PRS T2DGGI" = col_t2d,
                                 "Clínica + PRS MCPS"   = col_mcp)) +

  labs(
    title = "Calibración de Modelos Predictivos",
    subtitle = "Probabilidad predicha vs observada por decil de riesgo",
    x = "Probabilidad Predicha de Diabetes T2",
    y = "Probabilidad Observada (Casos Reales)",
    caption = "Línea punteada: calibración perfecta | Barras: IC 95%"
  ) +

  theme_minimal(base_size = 14) +
  theme(
    plot.title    = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(color = "gray40", size = 11),
    plot.caption  = element_text(color = "gray50", size = 9),
    plot.margin   = margin(15, 20, 15, 15),
    legend.position = c(0.25, 0.85),
    legend.background = element_rect(fill = "white", color = "gray80", linewidth = 0.3),
    legend.title = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(fig_dir, "04_calibration_comparison.png"), p_calib, width = 8, height = 7, dpi = 300)
ggsave(file.path(fig_dir, "04_calibration_comparison.pdf"), p_calib, width = 8, height = 7)

# ================================================================
# GRÁFICA 5: ESTRATIFICACIÓN — OR por decil
# ================================================================
cat("   Gráfica 5: OR por decil (referencia = Decil 5)...\n")

calc_or_by_decile <- function(data, prs_col, label) {
  data$Dec <- ntile(data[[prs_col]], 10)
  data$Dec_F <- relevel(as.factor(data$Dec), ref = "5")
  mod <- glm(T2D_STATUS ~ AGE + MALE + BMI + Dec_F, data = data, family = binomial)
  coefs <- summary(mod)$coefficients
  dec_rows <- grep("Dec_F", rownames(coefs))

  or_df <- data.frame(
    Decile = as.integer(gsub("Dec_F", "", rownames(coefs)[dec_rows])),
    OR = exp(coefs[dec_rows, "Estimate"]),
    CI_lo = exp(coefs[dec_rows, "Estimate"] - 1.96 * coefs[dec_rows, "Std. Error"]),
    CI_hi = exp(coefs[dec_rows, "Estimate"] + 1.96 * coefs[dec_rows, "Std. Error"]),
    GWAS = label
  )
  # Agregar decil 5 (referencia: OR = 1)
  or_df <- rbind(or_df, data.frame(Decile = 5, OR = 1, CI_lo = 1, CI_hi = 1, GWAS = label))
  or_df <- or_df[order(or_df$Decile), ]
  return(or_df)
}

or_t2d <- calc_or_by_decile(df_both, "PRS_Z_T2DGGI", "T2DGGI")
or_mcp <- calc_or_by_decile(df_both, "PRS_Z_MCPS",   "MCPS")
df_or  <- rbind(or_t2d, or_mcp)
df_or$GWAS <- factor(df_or$GWAS, levels = c("T2DGGI", "MCPS"))

cat(sprintf("   OR Decil 10 vs 5: T2DGGI=%.2f, MCPS=%.2f\n",
    or_t2d$OR[or_t2d$Decile == 10], or_mcp$OR[or_mcp$Decile == 10]))

p_or <- ggplot(df_or, aes(x = factor(Decile), y = OR, color = GWAS, group = GWAS)) +
  geom_pointrange(aes(ymin = CI_lo, ymax = CI_hi),
                  position = position_dodge(width = 0.4),
                  size = 0.7, linewidth = 0.7) +
  geom_line(position = position_dodge(width = 0.4), linewidth = 0.6, alpha = 0.6) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40", linewidth = 0.5) +

  scale_color_manual(values = c("T2DGGI" = col_t2d, "MCPS" = col_mcp), name = "PRS (SBayesRC)") +

  labs(
    title = "Odds Ratio de Diabetes T2 por Decil Genético",
    subtitle = "Referencia: Decil 5 (riesgo promedio) | Ajustado por Edad, Sexo, IMC",
    x = "Decil de Riesgo Poligénico",
    y = "Odds Ratio (IC 95%)",
    caption = "Fenotipo de Oro: HbA1c ≥ 6.5% + Diagnóstico + Medicación"
  ) +

  theme_minimal(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(color = "gray40", size = 11),
    plot.caption  = element_text(color = "gray50", size = 9),
    plot.margin   = margin(15, 20, 15, 15),
    legend.position = c(0.15, 0.88),
    legend.background = element_rect(fill = "white", color = "gray80", linewidth = 0.3),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(fig_dir, "05_or_by_decile.png"), p_or, width = 11, height = 7, dpi = 300)
ggsave(file.path(fig_dir, "05_or_by_decile.pdf"), p_or, width = 11, height = 7)

# ================================================================
# GRÁFICA 6: SCATTER PRS_T2DGGI vs PRS_MCPS coloreado por T2D
# ================================================================
cat("   Gráfica 6: Scatter PRS T2DGGI vs MCPS...\n")

# Subsample para visualización
set.seed(42)
n_plot <- min(nrow(df_both), 50000)
df_sc <- df_both[sample(nrow(df_both), n_plot), ]

r_prs <- cor(df_both$PRS_Z_T2DGGI, df_both$PRS_Z_MCPS, use = "complete.obs")
cat(sprintf("   Correlación PRS T2DGGI vs MCPS: r = %.4f\n", r_prs))

p_sc <- ggplot(df_sc, aes(x = PRS_Z_T2DGGI, y = PRS_Z_MCPS, color = factor(T2D_STATUS))) +
  geom_point(alpha = 0.15, size = 0.5) +
  geom_smooth(method = "lm", color = "black", linewidth = 0.8, se = FALSE, linetype = "dashed") +

  scale_color_manual(
    values = c("0" = "#4DBBD5", "1" = "#E64B35"),
    labels = c("Control", "Caso T2D"),
    name = "Estatus"
  ) +

  labs(
    title = "Concordancia entre PRS: T2DGGI vs MCPS",
    subtitle = sprintf("r = %.4f | n = %s | Coloreado por estatus T2D",
      r_prs, format(nrow(df_both), big.mark = ",")),
    x = "PRS T2DGGI (Z-Score)",
    y = "PRS MCPS (Z-Score)"
  ) +

  theme_minimal(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(color = "gray40", size = 11),
    plot.margin   = margin(15, 20, 15, 15),
    legend.position = c(0.15, 0.88),
    legend.background = element_rect(fill = "white", color = "gray80", linewidth = 0.3),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(fig_dir, "06_prs_concordance_scatter.png"), p_sc, width = 9, height = 8, dpi = 300)
ggsave(file.path(fig_dir, "06_prs_concordance_scatter.pdf"), p_sc, width = 9, height = 8)

# ================================================================
# GRÁFICA 7: BARPLOT RESUMEN MÉTRICAS (R², AUC, OR)
# ================================================================
cat("   Gráfica 7: Resumen de métricas...\n")

# Panel A: R² incremental
df_r2 <- data.frame(
  Model = factor(c("Base", "+ T2DGGI", "+ MCPS", "+ Ambos"),
                 levels = c("Base", "+ T2DGGI", "+ MCPS", "+ Ambos")),
  R2 = c(r2_base, r2_t2d, r2_mcp, r2_full)
)

p7a <- ggplot(df_r2, aes(x = Model, y = R2, fill = Model)) +
  geom_bar(stat = "identity", width = 0.6, color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f%%", R2)), vjust = -0.5, size = 4, fontface = "bold") +
  scale_fill_manual(values = c("gray60", col_t2d, col_mcp, "#7B2D8E")) +
  labs(title = "R² de Nagelkerke", x = "", y = "R² (%)") +
  ylim(0, max(df_r2$R2) * 1.2) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

# Panel B: AUC
df_auc <- data.frame(
  Model = factor(c("Base", "+ T2DGGI", "+ MCPS", "+ Ambos"),
                 levels = c("Base", "+ T2DGGI", "+ MCPS", "+ Ambos")),
  AUC = c(as.numeric(auc_base), as.numeric(auc_t2d), as.numeric(auc_mcp), as.numeric(auc_full))
)

p7b <- ggplot(df_auc, aes(x = Model, y = AUC, fill = Model)) +
  geom_bar(stat = "identity", width = 0.6, color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.3f", AUC)), vjust = -0.5, size = 4, fontface = "bold") +
  scale_fill_manual(values = c("gray60", col_t2d, col_mcp, "#7B2D8E")) +
  labs(title = "AUC (Curva ROC)", x = "", y = "AUC") +
  coord_cartesian(ylim = c(0.5, max(df_auc$AUC) * 1.05)) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

p7 <- p7a + p7b +
  plot_annotation(
    title = "Resumen de Métricas Predictivas — SBayesRC T2DGGI vs MCPS",
    subtitle = sprintf("Cohorte MCPS: %s casos, %s controles | Fenotipo de Oro",
      format(n_cases, big.mark = ","), format(n_controls, big.mark = ",")),
    theme = theme(
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(color = "gray40", size = 12)
    )
  )

ggsave(file.path(fig_dir, "07_metrics_summary.png"), p7, width = 13, height = 6, dpi = 300)
ggsave(file.path(fig_dir, "07_metrics_summary.pdf"), p7, width = 13, height = 6)

# ================================================================
# GRÁFICA 8: NET RECLASSIFICATION (NRI simplificado)
# ================================================================
cat("   Gráfica 8: Reclasificación neta...\n")

# Clasificar en riesgo alto (top 20%) con cada PRS
df_both$high_risk_t2d <- as.integer(df_both$PRS_Z_T2DGGI >= quantile(df_both$PRS_Z_T2DGGI, 0.8))
df_both$high_risk_mcp <- as.integer(df_both$PRS_Z_MCPS   >= quantile(df_both$PRS_Z_MCPS, 0.8))

# Tabla de reclasificación
reclass <- table(T2DGGI = df_both$high_risk_t2d, MCPS = df_both$high_risk_mcp, T2D = df_both$T2D_STATUS)

df_reclass <- df_both %>%
  mutate(
    Reclassified = case_when(
      high_risk_t2d == 1 & high_risk_mcp == 1 ~ "Ambos: Alto riesgo",
      high_risk_t2d == 0 & high_risk_mcp == 0 ~ "Ambos: Bajo riesgo",
      high_risk_t2d == 1 & high_risk_mcp == 0 ~ "Solo T2DGGI alto",
      high_risk_t2d == 0 & high_risk_mcp == 1 ~ "Solo MCPS alto"
    ),
    Status_Label = ifelse(T2D_STATUS == 1, "Caso T2D", "Control")
  )

df_reclass_summ <- df_reclass %>%
  group_by(Reclassified, Status_Label) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Status_Label) %>%
  mutate(pct = n / sum(n) * 100)

df_reclass_summ$Reclassified <- factor(df_reclass_summ$Reclassified,
  levels = c("Ambos: Alto riesgo", "Solo T2DGGI alto", "Solo MCPS alto", "Ambos: Bajo riesgo"))

p_nri <- ggplot(df_reclass_summ, aes(x = Reclassified, y = pct, fill = Status_Label)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.75),
           width = 0.65, color = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.1f%%\n(n=%s)", pct, format(n, big.mark = ","))),
            position = position_dodge(width = 0.75), vjust = -0.3, size = 3.2, fontface = "bold") +

  scale_fill_manual(values = c("Caso T2D" = "#E64B35", "Control" = "#4DBBD5"), name = "") +

  labs(
    title = "Reclasificación de Riesgo: T2DGGI vs MCPS",
    subtitle = "Umbral: Top 20% del PRS = Alto riesgo | Comparación de clasificación",
    x = "",
    y = "Porcentaje dentro del grupo",
    caption = "Individuos clasificados como alto riesgo por un PRS pero no por el otro"
  ) +

  theme_minimal(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(color = "gray40", size = 10),
    plot.caption  = element_text(color = "gray50", size = 9),
    plot.margin   = margin(15, 20, 15, 15),
    legend.position = "top",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 10)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18)))

ggsave(file.path(fig_dir, "08_reclassification.png"), p_nri, width = 12, height = 7, dpi = 300)
ggsave(file.path(fig_dir, "08_reclassification.pdf"), p_nri, width = 12, height = 7)

# ================================================================
# GUARDAR TABLA RESUMEN
# ================================================================
cat("\n>> Exportando tabla resumen...\n")

summary_table <- data.frame(
  Métrica = c("N total", "Casos T2D", "Controles",
              "AUC Base", "AUC + T2DGGI", "AUC + MCPS", "AUC + Ambos",
              "R² Base", "R² + T2DGGI", "R² + MCPS", "R² + Ambos",
              "R² incremental T2DGGI", "R² incremental MCPS",
              "Fold D10/D1 T2DGGI", "Fold D10/D1 MCPS",
              "Correlación PRS T2DGGI vs MCPS",
              "DeLong p-value"),
  Valor = c(
    format(n_total, big.mark = ","), format(n_cases, big.mark = ","), format(n_controls, big.mark = ","),
    auc_base, auc_t2d, auc_mcp, auc_full,
    sprintf("%.2f%%", r2_base), sprintf("%.2f%%", r2_t2d), sprintf("%.2f%%", r2_mcp), sprintf("%.2f%%", r2_full),
    sprintf("+%.2f%%", r2_t2d - r2_base), sprintf("+%.2f%%", r2_mcp - r2_base),
    sprintf("%.2fx", fold_t2d), sprintf("%.2fx", fold_mcp),
    sprintf("%.4f", r_prs),
    delong_p
  )
)

write.csv(summary_table, file.path(fig_dir, "validation_summary_table.csv"), row.names = FALSE)

# ================================================================
# RESUMEN FINAL
# ================================================================
cat("\n============================================================\n")
cat("  RESUMEN — VALIDACIÓN CLÍNICA COMPARATIVA\n")
cat("============================================================\n")
cat(sprintf("  Cohorte: MCPS (n = %s)\n", format(n_total, big.mark = ",")))
cat(sprintf("  Casos T2D: %s | Controles: %s\n", format(n_cases, big.mark = ","), format(n_controls, big.mark = ",")))
cat(sprintf("\n  AUC:  Base=%-7s  +T2DGGI=%-7s  +MCPS=%-7s  +Ambos=%s\n", auc_base, auc_t2d, auc_mcp, auc_full))
cat(sprintf("  R²:   Base=%-7.2f%%  +T2DGGI=%-7.2f%%  +MCPS=%-7.2f%%  +Ambos=%.2f%%\n", r2_base, r2_t2d, r2_mcp, r2_full))
cat(sprintf("  Fold: T2DGGI D10/D1=%.2fx | MCPS D10/D1=%.2fx\n", fold_t2d, fold_mcp))
cat(sprintf("  Correlación PRS: r = %.4f\n", r_prs))
cat(sprintf("  DeLong test: %s\n", delong_p))
cat(sprintf("\n  8 gráficas + tabla resumen en: %s/\n", fig_dir))
cat("============================================================\n")
EOF

apptainer exec \
    --bind /mnt/cephfs:/mnt/cephfs \
    --bind "$TMPDIR:/tmp,$TMPDIR:/var/tmp" \
    "$SIF" \
    /opt/miniforge3/envs/r-bio/bin/Rscript "$R_SCRIPT"

rm -rf "$TMPDIR"
echo ""
echo "==== VALIDACIÓN COMPARATIVA COMPLETADA ===="
