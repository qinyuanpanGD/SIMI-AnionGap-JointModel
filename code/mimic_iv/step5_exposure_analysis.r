library(tidyverse)
library(JMbayes2)
library(survival)
library(nlme)
library(splines)
library(missForest)
library(GLMMadaptive) # For mixed_model

# 1. 读取数据
if(file.exists("cohort_imputed.csv")) {
  cohort <- read_csv("cohort_imputed.csv")
} else {
  cohort <- read_csv("cohort_final_baseline.csv")
}
long_data <- read_csv("long_data_daily.csv")

# 2. 数据准备 - 生存数据 (id_data)
id_data <- cohort %>%
  filter(time_28_day > 0) %>%
  mutate(
    status = death_28_day, # 0=Alive, 1=Dead (28-day)
    time = time_28_day,    # Censored at 28 days
    fluid_input_24h = fluid_input_24h / 1000 # Convert mL to L
  )

base_vars <- c("admission_age", "gender", "charlson_comorbidity_index", 
               "SOFA_score", "fluid_input_24h", "WBC", 
               "Vasopressor_use", "MAP", "Mechanical_ventilation_use", 
               "RRT_use", "Creatinine", "Sedative_use",
               "SAPS_II")

set.seed(123)
data_to_impute <- id_data %>%
  select(all_of(c("stay_id", "time", "status", base_vars))) %>%
  as.data.frame() %>%
  mutate(across(where(is.character), as.factor))

if(any(sapply(data_to_impute, function(x) sum(is.na(x)) > 0))) {
  rf_imputed <- missForest(data_to_impute %>% select(-stay_id, -time, -status), verbose = FALSE)
  imputed_df <- rf_imputed$ximp
  for(v in base_vars) {
    id_data[[v]] <- imputed_df[[v]]
  }
}

# 3. 数据准备 - 纵向数据 (long_data)
long_df <- long_data %>%
  inner_join(id_data %>% select(stay_id), by = "stay_id") %>%
  mutate(
    visit_time = as.numeric(visit_day),
    fluid_input_daily = fluid_input_daily / 1000 # Convert mL to L
  ) %>%
  filter(!is.na(visit_time)) %>%
  filter(stay_id %in% id_data$stay_id) %>%
  left_join(id_data %>% select(stay_id, time), by = "stay_id") %>%
  filter(visit_time <= time) %>%
  filter(visit_time <= 28)

# 4. 生存子模型 (Cox Model)
cox_fit <- coxph(Surv(time, status) ~ admission_age + gender + charlson_comorbidity_index + 
                   SOFA_score + fluid_input_24h + WBC + Vasopressor_use + 
                   MAP + Mechanical_ventilation_use + RRT_use + 
                   Creatinine + Sedative_use + SAPS_II, 
                 data = id_data)

# 5. Step 6: 阈值暴露分析 (Table 3)
if(file.exists("ag_threshold_value.txt")) {
  ag_threshold_str <- readLines("ag_threshold_value.txt", warn=FALSE)
  ag_threshold <- as.numeric(ag_threshold_str[1])
} else {
  stop("Error: 'ag_threshold_value.txt' not found. Please run step2_rcs.r first.")
}

# 1. 一次性准备所有衍生列
long_df$high_ag <- ifelse(long_df$aniongap > ag_threshold, 1, 0)
long_df <- long_df %>%
  arrange(stay_id, visit_time) %>%
  group_by(stay_id) %>%
  mutate(cum_ag = cumsum(ifelse(aniongap > ag_threshold, aniongap - ag_threshold, 0))) %>%
  ungroup()

lme_fluid <- lme(fluid_input_daily ~ visit_time,
                 random = ~ visit_time | stay_id,
                 data = long_df,
                 na.action = na.exclude,
                 control = lmeControl(opt = "nlminb", msMaxIter = 200))

lme_bin <- mixed_model(high_ag ~ visit_time, random = ~ visit_time | stay_id,
                       data = long_df, family = binomial())

ctrl <- lmeControl(opt = "nlminb", msMaxIter = 200)
lme_cum <- lme(cum_ag ~ ns(visit_time, 3),
               random = ~ ns(visit_time, 3) | stay_id,
               data = long_df,
               na.action = na.exclude,
               control = ctrl)

if(file.exists("mcmc_params.txt")) {
  params <- readLines("mcmc_params.txt")
  iter_line <- grep("n_iter", params, value = TRUE)
  burnin_line <- grep("n_burnin", params, value = TRUE)
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

fit_jm_with_convergence <- function(surv_model, mixed_models, time_var, n_iter, n_burnin, target_rhat = 1.05, max_retries = 5) {
  current_iter <- n_iter
  current_burnin <- n_burnin

  fit <- jm(surv_model, mixed_models, time_var = time_var,
            n_iter = current_iter, n_burnin = current_burnin)

  for(i in 1:max_retries) {
    sum_obj <- summary(fit)
    rhats_surv <- sum_obj$Survival[, "Rhat"]
    rhats_long <- unlist(lapply(sum_obj$Outcome, function(x) x[, "Rhat"]))
    max_rhat <- max(c(rhats_surv, rhats_long), na.rm = TRUE)

    if(max_rhat <= target_rhat) {
      break
    } else {
      if(i == max_retries) {
      } else {
        current_iter <- current_iter + 5000L
        current_burnin <- current_burnin + 1000L
        fit <- jm(surv_model, mixed_models, time_var = time_var,
                  n_iter = current_iter, n_burnin = current_burnin)
      }
    }
  }
  return(fit)
}

jm_fit_bin <- fit_jm_with_convergence(cox_fit, list(lme_bin, lme_fluid), "visit_time",
                                      n_iter_val, n_burnin_val, target_rhat = 1.05)

surv_summary_bin <- summary(jm_fit_bin)$Survival
row_idx <- grep("^value\\(high_ag\\)$", rownames(surv_summary_bin))
hr_high_ag <- as.numeric(exp(surv_summary_bin[row_idx, "Mean"]))
ci_high_ag <- as.numeric(exp(surv_summary_bin[row_idx, c("2.5%", "97.5%")]))
p_high_ag <- as.numeric(surv_summary_bin[row_idx, "P"])

            ag_threshold, hr_high_ag, ci_high_ag[1], ci_high_ag[2], p_high_ag))

# 6. Step 7: 累积暴露分析 (Table S6)

# 1. 联合建模
jm_fit_cum <- fit_jm_with_convergence(cox_fit, list(lme_cum, lme_fluid), "visit_time",
                                      n_iter_val, n_burnin_val, target_rhat = 1.05)

surv_summary_cum <- summary(jm_fit_cum)$Survival
row_idx_cum <- grep("^value\\(cum_ag\\)$", rownames(surv_summary_cum))
hr_cum <- as.numeric(exp(surv_summary_cum[row_idx_cum, "Mean"]))
ci_cum <- as.numeric(exp(surv_summary_cum[row_idx_cum, c("2.5%", "97.5%")]))
p_cum <- as.numeric(surv_summary_cum[row_idx_cum, "P"])

            hr_cum, ci_cum[1], ci_cum[2], p_cum))

# 7. 格式化并保存结果为两个独立的 CSV 文件 (模仿乳酸格式)
format_table <- function(jm_fit_obj, tv_var_name, tv_display_name) {
  surv_summary <- summary(jm_fit_obj)$Survival

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

table_threshold <- format_table(jm_fit_bin, "high_ag", paste0("Any exposure to high AG (> ", ag_threshold, ")"))
write_csv(table_threshold, "Table_Threshold_Exposure_AG.csv")

table_cumulative <- format_table(jm_fit_cum, "cum_ag", paste0("Cumulative AG burden (excess > ", ag_threshold, ")"))
write_csv(table_cumulative, "Table_Cumulative_Exposure_AG.csv")
