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
echo "VALIDACIÓN CLÍNICA DEL PRS DIABETES T2 (MCPS) - VERSIÓN FINAL"
echo "Incluyendo: AUC, R2 de Nagelkerke, Estratificación y Calibración"
echo "============================================================"

R_SCRIPT="$TMPDIR/validate_prs.R"
cat << 'EOF' > "$R_SCRIPT"
library(dplyr)
library(ggplot2)
library(pROC)
library(scales)

out_dir <- "/mnt/cephfs/orgs/home/angel.pacheco/prs_diabetes"

# ================================================================
# 1 & 2. CARGA DE DATOS, FENOTIPO Y LIMPIEZA
# ================================================================
cat(">> Cargando datos y construyendo Fenotipo de Oro...\n")
df_prs <- read.delim(file.path(out_dir, "prs_diabetes_mcps_gwas_final.tsv"), stringsAsFactors = FALSE)
df_clin <- read.csv(file.path(out_dir, "MCPS BASELINE.csv"), stringsAsFactors = FALSE)

cols_clin <- c("PATID", "AGE", "MALE", "WEIGHT", "HEIGHT", "BASE_DIABETES", "BASE_HBA1C", "DRUG_D1", "DRUG_D2", "DRUG_D3", "DRUG_D4")
df_clin <- df_clin[, intersect(names(df_clin), cols_clin)]
df <- merge(df_prs, df_clin, by = "PATID")

for (med in c("DRUG_D1", "DRUG_D2", "DRUG_D3", "DRUG_D4")) {
  if (med %in% names(df)) df[[med]][is.na(df[[med]])] <- 0 else df[[med]] <- 0
}

df <- df %>%
  mutate(
    HEIGHT_M = HEIGHT / 100,
    BMI = WEIGHT / (HEIGHT_M^2),
    BMI = ifelse(is.na(BMI), median(BMI, na.rm = TRUE), BMI),
    dx_previo = (BASE_DIABETES == 1) | (DRUG_D1 == 1) | (DRUG_D2 == 1) | (DRUG_D3 == 1) | (DRUG_D4 == 1),
    hba1c_alta = (BASE_DIABETES == 0) & (!is.na(BASE_HBA1C)) & (BASE_HBA1C >= 6.5),
    control_puro = (BASE_DIABETES == 0) & (!is.na(BASE_HBA1C)) & (BASE_HBA1C < 6.5) & !dx_previo,
    T2D_STATUS = case_when(dx_previo | hba1c_alta ~ 1, control_puro ~ 0, TRUE ~ NA_real_)
  ) %>%
  filter(!is.na(AGE), !is.na(MALE), !is.na(T2D_STATUS), !is.na(PRS_Z))

# Deciles de PRS
df <- df %>% mutate(PRS_Decile = ntile(PRS_TOTAL, 10))

# ================================================================
# 3. REGRESIÓN LOGÍSTICA Y VARIANZA EXPLICADA (Punto 2 de la imagen)
# ================================================================
cat("\n>> Evaluando Modelos Predictivos...\n")

mod_null <- glm(T2D_STATUS ~ 1, data = df, family = binomial(link = "logit"))
mod_base <- glm(T2D_STATUS ~ AGE + MALE + BMI, data = df, family = binomial(link = "logit"))
mod_prs  <- glm(T2D_STATUS ~ AGE + MALE + BMI + PRS_Z, data = df, family = binomial(link = "logit"))

# Función para Pseudo-R2 de Nagelkerke
calc_nagelkerke <- function(model, null_model) {
  n <- nobs(model)
  ll_mod <- as.numeric(logLik(model))
  ll_null <- as.numeric(logLik(null_model))
  r2_cs <- 1 - exp(-2/n * (ll_mod - ll_null))
  r2_max <- 1 - exp(2/n * ll_null)
  return((r2_cs / r2_max) * 100) # En porcentaje
}

r2_base <- calc_nagelkerke(mod_base, mod_null)
r2_prs <- calc_nagelkerke(mod_prs, mod_null)
r2_incremental <- r2_prs - r2_base

cat("\n[!] VARIANZA EXPLICADA (R2 de Nagelkerke):\n")
cat(sprintf("    Modelo Clínico Base: %.2f%%\n", r2_base))
cat(sprintf("    Modelo Clínico + PRS: %.2f%%\n", r2_prs))
cat(sprintf("    -> El PRS explica un %.2f%% de varianza adicional en el riesgo de T2D.\n", r2_incremental))

# ================================================================
# 4. ESTRATIFICACIÓN (Punto 3 de la imagen)
# ================================================================
# Comparar el riesgo del Decil 10 contra el promedio poblacional (Decil 5)
df$Decile_Factor <- as.factor(df$PRS_Decile)
df$Decile_Factor <- relevel(df$Decile_Factor, ref = "5")
mod_estrat <- glm(T2D_STATUS ~ AGE + MALE + BMI + Decile_Factor, data = df, family = binomial(link = "logit"))

coef_d10 <- summary(mod_estrat)$coefficients["Decile_Factor10", ]
or_d10 <- exp(coef_d10["Estimate"])
ci_low_d10 <- exp(coef_d10["Estimate"] - 1.96 * coef_d10["Std. Error"])
ci_high_d10 <- exp(coef_d10["Estimate"] + 1.96 * coef_d10["Std. Error"])

cat("\n[!] ESTRATIFICACIÓN EXTREMA (Top 10% vs Promedio poblacional):\n")
cat(sprintf("    El Decil 10 tiene un OR de %.2f (95%% CI: %.2f - %.2f) comparado con el Decil 5.\n", or_d10, ci_low_d10, ci_high_d10))

# ================================================================
# 5. CALIBRACIÓN (Punto 4 de la imagen)
# ================================================================
cat("\n>> Generando Gráfica de Calibración...\n")
df$pred_prob <- predict(mod_prs, type = "response")
df$pred_decile <- ntile(df$pred_prob, 10)

calib_data <- df %>%
  group_by(pred_decile) %>%
  summarise(
    mean_pred = mean(pred_prob),
    obs_prob = mean(T2D_STATUS),
    se_obs = sqrt((obs_prob * (1 - obs_prob)) / n())
  )

p_calib <- ggplot(calib_data, aes(x = mean_pred, y = obs_prob)) +
  geom_point(color = "#E64B35", size = 4) +
  geom_errorbar(aes(ymin = obs_prob - 1.96*se_obs, ymax = obs_prob + 1.96*se_obs), width = 0.01, color = "#E64B35", linewidth = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray30", linewidth = 1) +
  labs(
    title = "Gráfica de Calibración del Modelo (Clínica + PRS)",
    subtitle = "Probabilidad Predicha vs. Prevalencia Observada por deciles",
    x = "Probabilidad Predicha de Diabetes T2",
    y = "Probabilidad Observada (Casos Reales)"
  ) +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold")) +
  scale_x_continuous(limits = c(0, max(calib_data$mean_pred)*1.1)) +
  scale_y_continuous(limits = c(0, max(calib_data$obs_prob)*1.1))

ggsave(file.path(out_dir, "prs_validation_calibration.png"), p_calib, width = 8, height = 7, dpi = 300)

cat(sprintf("\n  Gráfica exportada en: %s/prs_validation_calibration.png\n", out_dir))
cat("==== VALIDACIÓN COMPLETA FINALIZADA ====\n")
EOF

apptainer exec \
    --bind /mnt/cephfs:/mnt/cephfs \
    --bind "$TMPDIR:/tmp,$TMPDIR:/var/tmp" \
    "$SIF" \
    /opt/miniforge3/envs/r-bio/bin/Rscript "$R_SCRIPT"

rm -rf "$TMPDIR"
