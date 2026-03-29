library(tidyverse)
library(missForest)
library(corrplot)

# 1. 读取数据
cohort <- read_csv("cohort_final_baseline.csv")

# 2. 缺失值处理

impute_vars <- c("baseline_AG", "fluid_input_24h", "charlson_comorbidity_index", 
                 "admission_age", "gender", "WBC", "Vasopressor_use", 
                 "MAP", "SOFA_score", "Mechanical_ventilation_use", 
                 "RRT_use", "Creatinine", 
                 "Sedative_use",
                 "SAPS_II", "APS_III", "BUN")

missing_cols <- setdiff(impute_vars, colnames(cohort))
if(length(missing_cols) > 0) {
  warning("以下变量在数据中未找到，将跳过: ", paste(missing_cols, collapse = ", "))
  impute_vars <- intersect(impute_vars, colnames(cohort))
}

# 2.1 计算缺失率并保存 (已在 step1_missing.r 中完成，此处跳过以防覆盖)

# 2.2 准备插补数据
data_to_impute <- cohort %>%
  select(all_of(impute_vars)) %>%
  as.data.frame() # missForest 需要 data.frame

data_to_impute <- data_to_impute %>%
  mutate(across(where(is.character), as.factor)) %>%
  mutate(across(where(is.numeric), as.numeric))

set.seed(123) # 确保可重复性
rf_imputed <- missForest(data_to_impute, verbose = TRUE)

data_imputed_df <- rf_imputed$ximp

cohort_imputed <- cohort
for(v in impute_vars) {
  cohort_imputed[[v]] <- data_imputed_df[[v]]
}

write_csv(cohort_imputed, "cohort_imputed.csv")

# 3. 共线性分析 (Correlation Matrix)

cor_vars <- impute_vars[sapply(cohort_imputed[impute_vars], is.numeric)]
cor_data <- cohort_imputed %>% select(all_of(cor_vars))

M <- cor(cor_data, use = "pairwise.complete.obs")

pdf("figure_correlation_matrix.pdf", width = 10, height = 10)
corrplot(M, method = "color", type = "upper", order = "hclust", 
         addCoef.col = "black", # 添加相关系数
         tl.col = "black", tl.srt = 45, # 标签颜色和旋转
         number.cex = 0.7, # 数字大小
         diag = FALSE) # 不显示对角线
dev.off()

