library(tidyverse)
library(survival)
library(survminer)
library(Hmisc)
library(missForest)

# 1. 读取数据 (优先读取插补后的数据)
if(file.exists("cohort_imputed.csv")) {
  cohort <- read_csv("cohort_imputed.csv")
} else {
  warning("cohort_imputed.csv not found, using raw baseline data (may contain missing values).")
  cohort <- read_csv("cohort_final_baseline.csv")
}

# 2. 数据准备
data_cox <- cohort %>%
  filter(time_28_day > 0) %>%
  mutate(
    status = death_28_day, # 0=Alive, 1=Dead (28-day)
    time = time_28_day,    # Censored at 28 days
    fluid_input_24h_L = fluid_input_24h / 1000 # Convert mL to L
  )

exposure_var <- "baseline_AG" # 修正变量名
covariates <- c(
  "admission_age", 
  "gender", 
  "charlson_comorbidity_index", 
  "SOFA_score", 
  "fluid_input_24h_L", 
  "WBC", 
  "Vasopressor_use", 
  "MAP", 
  "Mechanical_ventilation_use", 
  "RRT_use", 
  "Creatinine", 
  "Sedative_use",
  "SAPS_II" # 新增 SAPS II
)

all_vars <- c(exposure_var, covariates)

# 3. 缺失值插补 (使用 missForest, 保持与 Step 2 一致)
set.seed(123)
data_to_impute <- data_cox %>%
  select(all_of(all_vars)) %>%
  as.data.frame() %>%
  mutate(across(where(is.character), as.factor))

if(any(sapply(data_to_impute, function(x) sum(is.na(x)) > 0))) {
  rf_imputed <- missForest(data_to_impute, verbose = FALSE)
  data_imputed_df <- rf_imputed$ximp

  for(v in all_vars) {
    data_cox[[v]] <- data_imputed_df[[v]]
  }
} else {
}

# 4. 单因素 Cox 回归 (Univariate)
uni_res <- list()

for(var in all_vars) {
  f <- as.formula(paste("Surv(time, status) ~", var))
  fit <- coxph(f, data = data_cox)
  s <- summary(fit)

  hr <- s$conf.int[1, "exp(coef)"]
  lower <- s$conf.int[1, "lower .95"]
  upper <- s$conf.int[1, "upper .95"]
  p <- s$coefficients[1, "Pr(>|z|)"]

  uni_res[[var]] <- data.frame(
    Variable = var,
    Uni_HR = sprintf("%.3f", hr),
    Uni_CI = sprintf("(%.3f-%.3f)", lower, upper),
    Uni_P = ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)),
    stringsAsFactors = FALSE
  )
}
uni_df <- do.call(rbind, uni_res)

# 5. 多因素 Cox 回归 (Multivariate)
multi_formula <- as.formula(paste("Surv(time, status) ~", paste(all_vars, collapse = " + ")))
fit_multi <- coxph(multi_formula, data = data_cox)
s_multi <- summary(fit_multi)

multi_res <- data.frame(
  Variable = rownames(s_multi$conf.int),
  Multi_HR = sprintf("%.3f", s_multi$conf.int[, "exp(coef)"]),
  Multi_CI = sprintf("(%.3f-%.3f)", s_multi$conf.int[, "lower .95"], s_multi$conf.int[, "upper .95"]),
  Multi_P = ifelse(s_multi$coefficients[, "Pr(>|z|)"] < 0.001, "<0.001", sprintf("%.3f", s_multi$coefficients[, "Pr(>|z|)"])),
  stringsAsFactors = FALSE
)

multi_res$Variable <- gsub("genderM", "gender", multi_res$Variable) 

# 6. 合并表格 (Table S3)
table_s3 <- left_join(uni_df, multi_res, by = "Variable")

final_table <- table_s3 %>%
  select(Variable, 
         `Univariate HR (95% CI)` = Uni_HR, `Univariate CI` = Uni_CI, `Univariate P` = Uni_P,
         `Multivariate HR (95% CI)` = Multi_HR, `Multivariate CI` = Multi_CI, `Multivariate P` = Multi_P) %>%
  mutate(
    `Univariate Analysis` = paste0(`Univariate HR (95% CI)`, " ", `Univariate CI`),
    `Multivariate Analysis` = paste0(`Multivariate HR (95% CI)`, " ", `Multivariate CI`)
  ) %>%
  select(Variable, `Univariate Analysis`, `Univariate P`, `Multivariate Analysis`, `Multivariate P`)

write_csv(final_table, "Table_S3_Cox_Regression.csv")
