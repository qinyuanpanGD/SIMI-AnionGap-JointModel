library(tidyverse)
library(survival)
library(survminer)
library(Hmisc)
library(missForest)

# 1. 读取数据 (使用插补后数据优先)
if(file.exists("cohort_imputed.csv")) {
  cohort <- read_csv("cohort_imputed.csv")
} else {
  cohort <- read_csv("cohort_final_baseline.csv")
}

# 2. 数据准备 (28-day Mortality)
data_cox <- cohort %>%
  filter(time_28_day > 0) %>%
  mutate(
    status = death_28_day, 
    time = time_28_day,
    fluid_input_24h_L = fluid_input_24h / 1000 # Convert mL to L
  ) %>%
  filter(!is.na(Lactate)) 

exposure_var <- "Lactate"
covariates <- c(
  "admission_age", "gender", "charlson_comorbidity_index", 
  "SOFA_score", "fluid_input_24h_L", "WBC", 
  "Vasopressor_use", "MAP", "Mechanical_ventilation_use", 
  "RRT_use", "Creatinine", "Sedative_use", "SAPS_II"
)

data_cox <- data_cox %>%
  mutate(
    gender = as.factor(gender),
    Vasopressor_use = as.factor(Vasopressor_use),
    Mechanical_ventilation_use = as.factor(Mechanical_ventilation_use),
    RRT_use = as.factor(RRT_use),
    Sedative_use = as.factor(Sedative_use)
  )

all_vars <- c(exposure_var, covariates)

# 3. 缺失值插补 (如果未插补)
if(any(is.na(data_cox[, all_vars]))) {
  set.seed(123)
  data_to_impute <- data_cox %>% select(all_of(all_vars)) %>% as.data.frame() %>% mutate(across(where(is.character), as.factor))
  rf_imputed <- missForest(data_to_impute, verbose = FALSE)
  data_imputed_df <- rf_imputed$ximp
  for(v in all_vars) data_cox[[v]] <- data_imputed_df[[v]]
}

# 4. 单因素 Cox
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

# 5. 多因素 Cox
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
multi_res$Variable <- gsub("Vasopressor_use1", "Vasopressor_use", multi_res$Variable) 
multi_res$Variable <- gsub("Mechanical_ventilation_use1", "Mechanical_ventilation_use", multi_res$Variable) 
multi_res$Variable <- gsub("RRT_use1", "RRT_use", multi_res$Variable) 
multi_res$Variable <- gsub("Sedative_use1", "Sedative_use", multi_res$Variable)

# 6. 合并表格
table_s3 <- left_join(uni_df, multi_res, by = "Variable")
final_table <- table_s3 %>%
  select(Variable, 
         `Univariate HR (95% CI)` = Uni_HR, `Univariate CI` = Uni_CI, `Univariate P` = Uni_P,
         `Multivariate HR (95% CI)` = Multi_HR, `Multivariate CI` = Multi_CI, `Multivariate P` = Multi_P)

write_csv(final_table, "script_lac/Table_S3_Cox_Regression_Lac.csv")
