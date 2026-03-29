library(tidyverse)
library(JMbayes2)
library(survival)
library(nlme)
library(splines)

# 1. 读取数据
cohort <- read_csv("cohort_final_baseline.csv")
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
               "SAPS_II", "baseline_AG") # 必须包含 baseline_AG

n_total <- nrow(id_data)

id_data <- id_data %>%
  drop_na(all_of(base_vars))

n_complete_baseline <- nrow(id_data)

# 3. 数据准备 - 纵向数据 (long_data)
# 1. 排除完全没有纵向测量的患者。
# 2. 联合模型 (Joint Model) 天然可以处理非平衡数据（即不需要每个人都有28天的完整数据）。

long_df_early <- long_data %>%
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

valid_stay_ids <- unique(long_df_early$stay_id)
n_before_long_filter <- nrow(id_data)
id_data <- id_data %>% filter(stay_id %in% valid_stay_ids)
n_final <- nrow(id_data)

# 4. 生存子模型 (Cox Model)
cox_fit <- coxph(Surv(time, status) ~ admission_age + gender + charlson_comorbidity_index + 
                   SOFA_score + fluid_input_24h + WBC + Vasopressor_use + 
                   MAP + Mechanical_ventilation_use + RRT_use + 
                   Creatinine + Sedative_use + SAPS_II, 
                 data = id_data)

# 5. 纵向子模型 (LME Submodels)
ctrl <- lmeControl(opt = "nlminb", msMaxIter = 200)

lme_ag <- lme(aniongap ~ ns(visit_time, 3),
              random = ~ ns(visit_time, 3) | stay_id,
              data = long_df_early,
              na.action = na.exclude,
              control = ctrl)

lme_fluid <- lme(fluid_input_daily ~ visit_time,
                 random = ~ visit_time | stay_id,
                 data = long_df_early,
                 na.action = na.exclude,
                 control = ctrl)

# 6. 联合模型 (Joint Model) - JMbayes2

fit_jm_with_convergence <- function(surv_model, mixed_models, time_var, n_iter, n_burnin, target_rhat = 1.05, max_retries = 5) {
  current_iter <- n_iter
  current_burnin <- n_burnin

  fit <- jm(surv_model, mixed_models, time_var = time_var,
            n_iter = current_iter, n_burnin = current_burnin)

  for(i in 1:max_retries) {
    sum_obj <- summary(fit)
    rhats_surv <- sum_obj$Survival[, "Rhat"]
    rhats_long <- unlist(lapply(sum_obj$Outcome, function(x) x[, "Rhat"]))

    all_rhats <- c(rhats_surv, rhats_long)
    max_rhat <- max(all_rhats, na.rm = TRUE)

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

jm_fit_early <- fit_jm_with_convergence(cox_fit, list(lme_ag, lme_fluid), 
                                      time_var = "visit_time",
                                      n_iter = 5000L, n_burnin = 1000L,
                                      target_rhat = 1.05, max_retries = 5) 

# 7. 提取并保存结果
surv_summary <- summary(jm_fit_early)$Survival
results_sensitivity <- data.frame(
  Variable = rownames(surv_summary),
  HR = exp(surv_summary[, "Mean"]),
  Lower_CI = exp(surv_summary[, "2.5%"]),
  Upper_CI = exp(surv_summary[, "97.5%"]),
  P_Value = surv_summary[, "P"]
)

results_formatted <- results_sensitivity %>%
  mutate(
    `HR (95% CI)` = sprintf("%.3f (%.3f - %.3f)", HR, Lower_CI, Upper_CI),
    `P Value` = ifelse(P_Value < 0.001, "<0.001", sprintf("%.3f", P_Value))
  ) %>%
  select(Variable, `HR (95% CI)`, `P Value`)

write_csv(results_formatted, "Table_S6_Sensitivity_CompleteCase_Day28.csv")
