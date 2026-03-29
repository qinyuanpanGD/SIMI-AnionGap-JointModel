rm(list = ls())
library(data.table)
library(dplyr)
library(tableone)

out_path <- "./"
mimic_path <- "./"

# 1. 读取 MIMIC-IV 队列
mimic_df <- fread(paste0(mimic_path, "cohort_final_baseline.csv"))
mimic_sub <- mimic_df[, .(
  admission_age, gender, SOFA_score, SAPS_II, charlson_comorbidity_index,
  Heart_rate, MAP, SpO2, WBC, Creatinine, Lactate, baseline_AG, fluid_input_24h,
  Vasopressor_use, Mechanical_ventilation_use, RRT_use, Sedative_use,
  mortality_28d = death_28_day
)]
mimic_sub[, Cohort := "MIMIC-IV (Primary)"]

# 2. 读取 Zigong 队列
zigong_df <- fread(paste0(out_path, "zigong_baseline_features.csv"))
zigong_sub <- zigong_df[, .(
  admission_age, gender, SOFA_score, SAPS_II, charlson_comorbidity_index,
  Heart_rate, MAP, SpO2, WBC, Creatinine, Lactate = lactate, baseline_AG = aniongap, fluid_input_24h,
  Vasopressor_use, Mechanical_ventilation_use, RRT_use, Sedative_use,
  mortality_28d
)]
zigong_sub[, Cohort := "Zigong (Validation)"]
zigong_sub[gender == "Male", gender := "M"]
zigong_sub[gender == "Female", gender := "F"]

# 3. 合并数据集
combined_df <- rbind(mimic_sub, zigong_sub, fill = TRUE)
combined_df$Cohort <- factor(combined_df$Cohort, levels = c("MIMIC-IV (Primary)", "Zigong (Validation)"))

binary_vars <- c("Vasopressor_use", "Mechanical_ventilation_use", "RRT_use", "Sedative_use", "mortality_28d")
for (v in binary_vars) {
  combined_df[[v]] <- factor(combined_df[[v]], levels = c(0, 1), labels = c("No", "Yes"))
}

# 4. 生成基线对比表
vars <- c("admission_age", "gender", "SOFA_score", "SAPS_II", "charlson_comorbidity_index",
          "Heart_rate", "MAP", "SpO2", "WBC", "Creatinine", "Lactate", "baseline_AG", "fluid_input_24h",
          "Vasopressor_use", "Mechanical_ventilation_use", "RRT_use", "Sedative_use", "mortality_28d")

catVars <- c("gender", "Vasopressor_use", "Mechanical_ventilation_use", "RRT_use", "Sedative_use", "mortality_28d")
nonNormalVars <- c("admission_age", "SOFA_score", "SAPS_II", "charlson_comorbidity_index", 
                   "Heart_rate", "MAP", "SpO2", "WBC", "Creatinine", "Lactate", "baseline_AG", "fluid_input_24h")

tab1 <- CreateTableOne(vars = vars, strata = "Cohort", data = combined_df, factorVars = catVars, test = TRUE)

tab1_mat <- print(tab1, nonnormal = nonNormalVars, quote = FALSE, noSpaces = TRUE, printToggle = FALSE, showAllLevels = TRUE)

tab1_df <- as.data.frame(tab1_mat)
tab1_df$Variable <- rownames(tab1_mat)
if("level" %in% colnames(tab1_df)) {
  for(i in 2:nrow(tab1_df)) {
    if(tab1_df$Variable[i] == "" || is.na(tab1_df$Variable[i])) {
      tab1_df$Variable[i] <- tab1_df$Variable[i-1]
    }
  }
  tab1_df$Variable <- ifelse(tab1_df$level == "" | is.na(tab1_df$level), 
                             tab1_df$Variable, 
                             paste0(tab1_df$Variable, " = ", tab1_df$level))
}

tab1_df <- tab1_df[, c("Variable", "MIMIC-IV (Primary)", "Zigong (Validation)", "p", "test")]

write.csv(tab1_df, paste0(out_path, "Table_S7_Baseline_Comparison.csv"), row.names = FALSE)
