library(tidyverse)
library(JMbayes2)
library(survival)
library(nlme)
library(splines)

cohort_base <- read_csv("zigong_baseline_features.csv")
simi_cohort <- read_csv("zigong_simi_cohort.csv") %>% select(INP_NO, death_time_icu, follow_time_icu)
cohort <- cohort_base %>% left_join(simi_cohort, by = "INP_NO")
long_data <- read_csv("zigong_longitudinal_features.csv")

# 1. 准备生存数据 (id_data)
id_data <- cohort %>%
  mutate(
    time = ifelse(!is.na(death_time_icu), death_time_icu / 24, follow_time_icu / 24),
    status = ifelse(!is.na(death_time_icu) & death_time_icu <= 28*24, 1, 0),
    time = ifelse(time > 28, 28, time),
    fluid_input_24h = fluid_input_24h / 1000,
    stay_id = INP_NO
  ) %>%
  filter(time > 0) %>%
  drop_na(admission_age, gender, SOFA_score, charlson_comorbidity_index, fluid_input_24h, WBC, 
          Vasopressor_use, MAP, Mechanical_ventilation_use, 
          RRT_use, Creatinine, Sedative_use)

# 2. 准备纵向数据 (long_data)
long_df <- long_data %>%
  rename(stay_id = INP_NO, visit_time = day) %>%
  inner_join(id_data %>% select(stay_id, admission_age, gender, time), by = "stay_id") %>%
  mutate(fluid_input_daily = fluid_input_L) %>%
  filter(!is.na(visit_time)) %>%
  filter(visit_time <= time) %>%
  filter(visit_time <= 28)

# 3. 拟合子模型
cox_fit <- coxph(Surv(time, status) ~ admission_age + gender + 
                   SOFA_score + charlson_comorbidity_index + fluid_input_24h + WBC + Vasopressor_use + 
                   MAP + Mechanical_ventilation_use + RRT_use + 
                   Creatinine + Sedative_use,
                 data = id_data)

ctrl <- lmeControl(opt = "nlminb", msMaxIter = 200)

lme_lac <- lme(lactate ~ ns(visit_time, 3),
              random = ~ visit_time | stay_id,
              data = long_df,
              na.action = na.exclude,
              control = lmeControl(opt = "optim", msMaxIter = 500))

lme_fluid <- lme(fluid_input_daily ~ visit_time,
                 random = ~ visit_time | stay_id,
                 data = long_df,
                 na.action = na.exclude,
                 control = ctrl)

# 4. 联合模型
current_n_iter <- 5000L
current_n_burnin <- 1000L

jm_fit <- jm(cox_fit, list(lme_lac, lme_fluid),
             time_var = "visit_time",
             n_iter = current_n_iter, n_burnin = current_n_burnin)

# 5. 输出结果 Table 2
surv_summary <- summary(jm_fit)$Survival
table2 <- data.frame(
  Variable = rownames(surv_summary),
  HR = exp(surv_summary[, "Mean"]),
  Lower_CI = exp(surv_summary[, "2.5%"]),
  Upper_CI = exp(surv_summary[, "97.5%"]),
  P_Value = surv_summary[, "P"]
)
table2_formatted <- table2 %>%
  mutate(
    `HR (95% CI)` = sprintf("%.3f (%.3f-%.3f)", HR, Lower_CI, Upper_CI),
    `P Value` = ifelse(P_Value < 0.001, "<0.001", sprintf("%.3f", P_Value))
  ) %>% select(Variable, `HR (95% CI)`, `P Value`)

write_csv(table2_formatted, "Table_2_JMbayes_Lac_Zigong.csv")

# 6. 时变模型与绘图
knots_28 <- c(7, 14, 21)
form_splines <- list(
  "lactate" = ~ value(lactate) * ns(visit_time, k = knots_28, B = c(0, 28))
)
jm_fit_tv <- update(jm_fit, functional_forms = form_splines,
                    n_iter = current_n_iter, n_burnin = current_n_burnin)

time_seq <- seq(0, 28, length.out = 100)
W_time <- ns(time_seq, k = knots_28, B = c(0, 28))
mcmc_alphas <- jm_fit_tv$mcmc$alphas
alphas_all <- do.call(rbind, mcmc_alphas)

alpha_main_col <- grep("^value\\(lactate\\)$", colnames(alphas_all), value = TRUE)
alpha_inter_cols <- grep("value\\(lactate\\):ns", colnames(alphas_all), value = TRUE)

alphas_main <- if(length(alpha_main_col) > 0) alphas_all[, alpha_main_col] else rep(0, nrow(alphas_all))

if(length(alpha_inter_cols) > 0) {
  alphas_inter <- alphas_all[, alpha_inter_cols]
  if(ncol(alphas_inter) < ncol(W_time)) W_time <- W_time[, 1:ncol(alphas_inter)]
  inter_term <- W_time %*% t(alphas_inter)
} else {
  inter_term <- matrix(0, nrow = 100, ncol = nrow(alphas_all))
}

pred_log_hr <- sweep(inter_term, 2, alphas_main, "+")
hr_summary <- apply(pred_log_hr, 1, function(x) c(Mean = mean(x), Lower = quantile(x, 0.025), Upper = quantile(x, 0.975)))

plot_df <- as.data.frame(t(hr_summary))
plot_df$Time = as.numeric(time_seq)
plot_df$HR = as.numeric(exp(plot_df$Mean))
plot_df$Lower_CI = as.numeric(exp(plot_df$Lower))
plot_df$Upper_CI = as.numeric(exp(plot_df$Upper))

p_tv <- ggplot(plot_df, aes(x = Time, y = HR)) +
  geom_line(color = "#E74C3C", linewidth = 1) +
  geom_ribbon(aes(ymin = Lower_CI, ymax = Upper_CI), fill = "#E74C3C", alpha = 0.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
  labs(title = "Zigong Cohort: Time-Varying HR of Lactate", x = "Days since ICU Admission", y = "Hazard Ratio (95% CI)") +
  theme_classic(base_size = 14)

ggsave("Figure_2_TimeVarying_HR_Lac_Zigong.pdf", p_tv, width = 8, height = 6)
