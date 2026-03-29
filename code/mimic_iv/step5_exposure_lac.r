source("renv/activate.R")

library(tidyverse)
library(JMbayes2)
library(survival)
library(nlme)
library(splines)
library(GLMMadaptive)

# 1. 读取数据
if(file.exists("../cohort_imputed.csv")) {
  cohort <- read_csv("cohort_imputed.csv")
} else {
  cohort <- read_csv("cohort_final_baseline.csv")
}
long_data <- read_csv("long_data_daily_lac.csv")

# 2. 数据准备 (28-day)
id_data <- cohort %>%
  filter(time_28_day > 0) %>%
  mutate(
    status = death_28_day,
    time = time_28_day,
    fluid_input_24h = fluid_input_24h / 1000 # Convert mL to L
  )

base_vars <- c("admission_age", "gender", "charlson_comorbidity_index", 
               "SOFA_score", "fluid_input_24h", "WBC", 
               "Vasopressor_use", "MAP", "Mechanical_ventilation_use", 
               "RRT_use", "Creatinine", "Sedative_use", "SAPS_II")

if(any(is.na(id_data[, base_vars]))) {
  set.seed(123)
  data_to_impute <- id_data %>% select(all_of(c("stay_id", "time", "status", base_vars))) %>% as.data.frame() %>% mutate(across(where(is.character), as.factor))
  rf_imputed <- missForest(data_to_impute %>% select(-stay_id, -time, -status), verbose = FALSE)
  imputed_df <- rf_imputed$ximp
  for(v in base_vars) id_data[[v]] <- imputed_df[[v]]
}

long_df <- long_data %>%
  inner_join(id_data %>% select(stay_id), by = "stay_id") %>%
  mutate(
    visit_time = as.numeric(visit_day),
    fluid_input_daily = fluid_input_daily / 1000 # Convert mL to L
  ) %>%
  filter(!is.na(visit_time)) %>%
  filter(stay_id %in% id_data$stay_id) %>%
  left_join(id_data %>% select(stay_id, time), by = "stay_id") %>%
  filter(visit_time <= time & visit_time <= 28)

cox_fit <- coxph(Surv(time, status) ~ admission_age + gender + charlson_comorbidity_index + 
                   SOFA_score + fluid_input_24h + WBC + Vasopressor_use + 
                   MAP + Mechanical_ventilation_use + RRT_use + 
                   Creatinine + Sedative_use + SAPS_II,
                 data = id_data)

# 3. 阈值设定
if(file.exists("lac_threshold_value.txt")) {
  lac_threshold_str <- readLines("lac_threshold_value.txt", warn=FALSE)
  lac_threshold <- as.numeric(lac_threshold_str[1])
} else {
  stop("Error: 'lac_threshold_value.txt' not found. Please run step2_rcs_lac.r first to generate the threshold.")
}

# 4. 数据准备：一次性创建所有需要的列 (high_lac, cum_lac)

# 4.1 High Lactate
long_df$high_lac <- ifelse(long_df$lactate > lac_threshold, 1, 0)

# 4.2 Cumulative Lactate
long_df <- long_df %>%
  arrange(stay_id, visit_time) %>%
  group_by(stay_id) %>%
  mutate(cum_lac = cumsum(ifelse(lactate > lac_threshold, lactate - lac_threshold, 0))) %>%
  ungroup()

# 5. 统一拟合子模型 (使用完全相同的 long_df 对象)

lme_fluid <- lme(fluid_input_daily ~ visit_time,
                 random = ~ visit_time | stay_id,
                 data = long_df,
                 na.action = na.exclude,
                 control = lmeControl(opt = "nlminb", msMaxIter = 200))

lme_bin <- mixed_model(high_lac ~ visit_time, random = ~ visit_time | stay_id, 
                       data = long_df, family = binomial())

ctrl <- lmeControl(opt = "nlminb", msMaxIter = 200)
lme_cum <- lme(cum_lac ~ ns(visit_time, 3), 
               random = ~ ns(visit_time, 3) | stay_id, 
               data = long_df,
               na.action = na.exclude,
               control = ctrl)

# 6. Threshold Exposure Analysis
if(file.exists("mcmc_params.txt")) {
  params <- readLines("mcmc_params.txt")
  iter_line <- grep("n_iter", params, value=TRUE)
  burnin_line <- grep("n_burnin", params, value=TRUE)

  if(length(iter_line) > 0 && length(burnin_line) > 0) {
    n_iter_val <- as.integer(sub(".*=", "", iter_line))
    n_burnin_val <- as.integer(sub(".*=", "", burnin_line))
  } else {
    n_iter_val <- 15000L
    n_burnin_val <- 3000L
  }
} else {
  n_iter_val <- 15000L
  n_burnin_val <- 3000L
}

jm_fit_bin <- jm(cox_fit, list(lme_bin, lme_fluid), time_var = "visit_time",
                 n_iter = n_iter_val, n_burnin = n_burnin_val)

# 7. Cumulative Exposure Analysis

jm_fit_cum <- jm(cox_fit, list(lme_cum, lme_fluid), time_var = "visit_time",
                 n_iter = n_iter_val, n_burnin = n_burnin_val) # 自动继承

# 8. Output formatting (Mimic the image)

format_table <- function(jm_model, tv_var_name, tv_display_name) {
  surv_summary <- summary(jm_model)$Survival

  tv_idx <- grep("value\\(", rownames(surv_summary))

  base_df <- data.frame(
    Category = "Baseline characteristics",
    Variable = rownames(surv_summary)[-tv_idx],
    `HR (95% CI)` = sprintf("%.3f (%.3f-%.3f)", 
                            exp(surv_summary[-tv_idx, "Mean"]),
                            exp(surv_summary[-tv_idx, "2.5%"]),
                            exp(surv_summary[-tv_idx, "97.5%"])),
    `P value` = ifelse(surv_summary[-tv_idx, "P"] < 0.001, "<0.001", sprintf("%.3f", surv_summary[-tv_idx, "P"])),
    check.names = FALSE
  )

  tv_df <- data.frame(
    Category = "Time-varying variables",
    Variable = rownames(surv_summary)[tv_idx],
    `HR (95% CI)` = sprintf("%.3f (%.3f-%.3f)", 
                            exp(surv_summary[tv_idx, "Mean"]),
                            exp(surv_summary[tv_idx, "2.5%"]),
                            exp(surv_summary[tv_idx, "97.5%"])),
    `P value` = ifelse(surv_summary[tv_idx, "P"] < 0.001, "<0.001", sprintf("%.3f", surv_summary[tv_idx, "P"])),
    check.names = FALSE
  )

  tv_df$Variable <- ifelse(grepl(tv_var_name, tv_df$Variable), tv_display_name, tv_df$Variable)
  tv_df$Variable <- ifelse(grepl("fluid_input_daily", tv_df$Variable), "Fluid input daily", tv_df$Variable)

  rbind(base_df, tv_df)
}

table_s4 <- format_table(jm_fit_bin, "high_lac", paste0("Any exposure to high Lactate (> ", lac_threshold, " mmol/L)"))
write_csv(table_s4, "Table_Threshold_Exposure_Lac.csv")

table_s5 <- format_table(jm_fit_cum, "cum_lac", paste0("Cumulative Lactate burden (excess > ", lac_threshold, " mmol/L)"))
write_csv(table_s5, "Table_Cumulative_Exposure_Lac.csv")

