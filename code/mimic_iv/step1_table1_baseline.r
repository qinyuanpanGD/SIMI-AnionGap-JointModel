library(tidyverse)
library(tableone)

# 1. 读取数据 (使用原始未插补数据)
cohort <- read_csv("cohort_final_baseline.csv")

# 2. 数据清洗与变量转换

cat_vars <- c("gender", "hospital_expire_flag", 
              "Diabetes", "Renal", "Liver", "CAD", "Stroke", "Malignancy", "CHF",
              "RRT_use", "Mechanical_ventilation_use", "Vasopressor_use", "Sedative_use")

cont_vars <- c("admission_age", "SOFA_score", "baseline_AG", "los_icu", "baseline_Troponin",
               "WBC", "Hemoglobin", "Platelet",
               "PT", "APTT", "BUN", "Creatinine",
               "Glucose", "Lactate", "pH", "pCO2",
               "Heart_rate", "MAP", "Temperature", "Respiratory_rate", "SpO2",
               "SAPS_II", "APS_III", "charlson_comorbidity_index", "fluid_input_24h")

cohort <- cohort %>%
  mutate(
    hospital_expire_flag = factor(hospital_expire_flag, levels = c(0, 1), labels = c("Survivors", "Non-Survivors")),
    gender = factor(gender),
    across(all_of(cat_vars[-c(1,2)]), factor) # 其他二分类变量
  )

# 3. 构建 Table 1

all_vars <- c(cont_vars, cat_vars[-2]) # 排除分组变量本身

non_normal_vars <- c("los_icu", "baseline_Troponin", "SOFA_score", "baseline_AG", 
                     "WBC", "Platelet", "BUN", "Creatinine", "Lactate", "pCO2",
                     "SAPS_II", "APS_III", "charlson_comorbidity_index", "fluid_input_24h")

tab1 <- CreateTableOne(vars = all_vars, strata = "hospital_expire_flag", data = cohort, addOverall = TRUE)

# 4. 打印结果 (包含 P 值)

tab1_mat <- print(tab1, nonnormal = non_normal_vars, quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
write.csv(tab1_mat, "Table_1_Baseline.csv")

# 5. 简单展示部分结果
