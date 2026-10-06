#!/bin/bash
#SBATCH --job-name=val_prs_R
#SBATCH --account=researchers
#SBATCH --qos=vip
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=/mnt/cephfs/orgs/home/angel.pacheco/logs/val_prs_R_%j.out
#SBATCH --error=/mnt/cephfs/orgs/home/angel.pacheco/logs/val_prs_R_%j.err

set -euo pipefail

SIF="/mnt/cephfs/biocontainers/images/jupyter-biotools-1.2.sif"
TMPDIR="/mnt/cephfs/scratch/$USER/val_tmp/$SLURM_JOB_ID"
mkdir -p "$TMPDIR"
chmod 700 "$TMPDIR"
export TMPDIR TMP="$TMPDIR" TEMP="$TMPDIR"

echo "============================================================"
echo "VALIDACIÓN CLÍNICA DEL PRS DIABETES T2 (MCPS) - VERSIÓN R"
echo "Usando Fenotipo de Oro: HbA1c >= 6.5% + Diagnóstico + Meds"
echo "============================================================"

R_SCRIPT="$TMPDIR/validate_prs.R"
cat << 'EOF' > "$R_SCRIPT"
# Cargar librerías necesarias
library(dplyr)
library(ggplot2)
library(pROC)
library(scales)

out_dir <- "/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes"

# ================================================================
# 1. CARGA DE DATOS Y MERGE
# ================================================================
cat(">> Cargando datos...\n")
df_prs <- read.delim(file.path(out_dir, "prs_diabetes_mcps_gwas_v4_final.tsv"), stringsAsFactors = FALSE)
df_clin <- read.csv(file.path(out_dir, "MCPS BASELINE.csv"), stringsAsFactors = FALSE)

# Filtrar columnas clínicas relevantes
cols_clin <- c("PATID", "AGE", "MALE", "WEIGHT", "HEIGHT", 
               "BASE_DIABETES", "BASE_HBA1C", 
               "DRUG_D1", "DRUG_D2", "DRUG_D3", "DRUG_D4")
df_clin <- df_clin[, intersect(names(df_clin), cols_clin)]

# Cruce
df <- merge(df_prs, df_clin, by = "PATID")
cat(sprintf("  Pacientes cruzados inicialmente: %s\n", format(nrow(df), big.mark=",")))

# ================================================================
# 2. CONSTRUCCIÓN DEL FENOTIPO Y LIMPIEZA
# ================================================================
cat(">> Construyendo Fenotipo de Oro para Diabetes...\n")

# Rellenar NAs en medicamentos con 0
for (med in c("DRUG_D1", "DRUG_D2", "DRUG_D3", "DRUG_D4")) {
  if (med %in% names(df)) {
    df[[med]][is.na(df[[med]])] <- 0
  } else {
    df[[med]] <- 0
  }
}

# Crear el estatus de T2D y calcular IMC
df <- df %>%
  mutate(
    # IMC
    HEIGHT_M = HEIGHT / 100,
    BMI = WEIGHT / (HEIGHT_M^2),
    
    # Condiciones Fenotipo
    dx_previo = (BASE_DIABETES == 1) | (DRUG_D1 == 1) | (DRUG_D2 == 1) | (DRUG_D3 == 1) | (DRUG_D4 == 1),
    hba1c_alta = (BASE_DIABETES == 0) & (!is.na(BASE_HBA1C)) & (BASE_HBA1C >= 6.5),
    control_puro = (BASE_DIABETES == 0) & (!is.na(BASE_HBA1C)) & (BASE_HBA1C < 6.5) & !dx_previo,
    
    T2D_STATUS = case_when(
      dx_previo | hba1c_alta ~ 1,
      control_puro ~ 0,
      TRUE ~ NA_real_
    )
  ) %>%
  filter(!is.na(AGE), !is.na(MALE), !is.na(BMI), !is.na(T2D_STATUS), !is.na(PRS_Z))

cat(sprintf("  Pacientes finales tras limpieza: %s\n", format(nrow(df), big.mark=",")))
cat(sprintf("  Casos confirmados: %s\n", format(sum(df$T2D_STATUS == 1), big.mark=",")))
cat(sprintf("  Controles verdaderos: %s\n", format(sum(df$T2D_STATUS == 0), big.mark=",")))

# ================================================================
# 3. REGRESIÓN LOGÍSTICA
# ================================================================
cat("\n>> Evaluando Modelos Predictivos (Regresión Logística)...\n")

# Modelo Base (Clínica)
mod_base <- glm(T2D_STATUS ~ AGE + MALE + BMI, data = df, family = binomial(link = "logit"))

# Modelo Completo (Clínica + PRS)
mod_prs <- glm(T2D_STATUS ~ AGE + MALE + BMI + PRS_Z, data = df, family = binomial(link = "logit"))

# Extraer OR e Intervalos de Confianza para PRS_Z
prs_coef <- summary(mod_prs)$coefficients["PRS_Z", ]
or_prs <- exp(prs_coef["Estimate"])
ci_lower <- exp(prs_coef["Estimate"] - 1.96 * prs_coef["Std. Error"])
ci_upper <- exp(prs_coef["Estimate"] + 1.96 * prs_coef["Std. Error"])
p_val <- prs_coef["Pr(>|z|)"]

cat("\n--- Resultados del Modelo Completo ---\n")
print(summary(mod_prs))
cat(sprintf("\n[!] ODDS RATIO DEL PRS: %.3f (95%% CI: %.3f - %.3f), p-value: %e\n", 
            or_prs, ci_lower, ci_upper, p_val))

# ================================================================
# 4. GRÁFICA 1: CURVA ROC
# ================================================================
cat("\n>> Generando gráficas de validación...\n")

roc_base <- roc(df$T2D_STATUS, fitted(mod_base), quiet = TRUE)
roc_prs <- roc(df$T2D_STATUS, fitted(mod_prs), quiet = TRUE)

auc_base <- format(round(auc(roc_base), 3), nsmall=3)
auc_prs <- format(round(auc(roc_prs), 3), nsmall=3)

# Crear data.frame para ggplot
df_roc <- data.frame(
  Specificity = c(roc_base$specificities, roc_prs$specificities),
  Sensitivity = c(roc_base$sensitivities, roc_prs$sensitivities),
  Model = factor(
    c(rep(sprintf("Clínica Base (AUC = %s)", auc_base), length(roc_base$sensitivities)),
      rep(sprintf("+ PRS Genético (AUC = %s)", auc_prs), length(roc_prs$sensitivities))),
    levels = c(sprintf("Clínica Base (AUC = %s)", auc_base), 
               sprintf("+ PRS Genético (AUC = %s)", auc_prs))
  )
)

p_roc <- ggplot(df_roc, aes(x = 1 - Specificity, y = Sensitivity, color = Model)) +
  geom_line(linewidth = 1.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  scale_color_manual(values = c("gray50", "#E64B35")) +
  labs(
    title = "Curva ROC: Predicción de Diabetes Tipo 2",
    subtitle = "Cohorte MCPS",
    x = "Tasa de Falsos Positivos (1 - Especificidad)",
    y = "Tasa de Verdaderos Positivos (Sensibilidad)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = c(0.75, 0.2),
    legend.background = element_rect(fill = "white", color = "gray80"),
    legend.title = element_blank()
  )

ggsave(file.path(out_dir, "prs_validation_roc.png"), p_roc, width = 8, height = 7, dpi = 300)

# ================================================================
# 5. GRÁFICA 2: DECILES DE RIESGO
# ================================================================
df <- df %>%
  mutate(PRS_Decile = ntile(PRS_TOTAL, 10))

df_deciles <- df %>%
  group_by(PRS_Decile) %>%
  summarise(Prevalence = mean(T2D_STATUS) * 100, n = n())

p_dec <- ggplot(df_deciles, aes(x = factor(PRS_Decile), y = Prevalence, fill = Prevalence)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
  scale_fill_gradient(low = "#FFCDD2", high = "#B71C1C") +
  geom_text(aes(label = sprintf("%.1f%%", Prevalence)), vjust = -0.5, size = 4, fontface = "bold") +
  labs(
    title = "Prevalencia de Diabetes T2 por Decil Genético",
    x = "Decil de Riesgo Poligénico (1 = Menor, 10 = Mayor)",
    y = "Prevalencia de Diabetes (%)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "none",
    panel.grid.major.x = element_blank()
  ) +
  ylim(0, max(df_deciles$Prevalence) * 1.15)

ggsave(file.path(out_dir, "prs_validation_deciles.png"), p_dec, width = 10, height = 6, dpi = 300)

# ================================================================
# 6. GRÁFICA 3: DENSIDAD (Casos vs Controles)
# ================================================================
df$Status_Label <- ifelse(df$T2D_STATUS == 1, "Casos (Diabetes)", "Controles (Sanos)")

p_den <- ggplot(df, aes(x = PRS_Z, fill = Status_Label, color = Status_Label)) +
  geom_density(alpha = 0.4, linewidth = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  scale_fill_manual(values = c("Casos (Diabetes)" = "#E64B35", "Controles (Sanos)" = "#4DBBD5")) +
  scale_color_manual(values = c("Casos (Diabetes)" = "#E64B35", "Controles (Sanos)" = "#4DBBD5")) +
  labs(
    title = "Distribución del Riesgo Genético: Casos vs Controles",
    x = "PRS Estandarizado (Z-Score)",
    y = "Densidad"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "top",
    legend.title = element_blank()
  )

ggsave(file.path(out_dir, "prs_validation_density.png"), p_den, width = 10, height = 6, dpi = 300)

cat(sprintf("  Gráficas exportadas en: %s/\n", out_dir))
cat("==== PIPELINE DE VALIDACIÓN COMPLETADO ====\n")
EOF

apptainer exec \
    --bind /mnt/cephfs:/mnt/cephfs \
    --bind "$TMPDIR:/tmp,$TMPDIR:/var/tmp" \
    "$SIF" \
    /opt/miniforge3/envs/r-bio/bin/Rscript "$R_SCRIPT"

rm -rf "$TMPDIR"
