library(tidyverse)
library(survival)
library(rms)
library(ggplot2)

# 1. 读取数据 (已插补)
if(file.exists("cohort_imputed.csv")) {
  data_rcs <- read_csv("cohort_imputed.csv")
} else {
  warning("cohort_imputed.csv not found, using raw data (may have missing values). Please run step1_missing_collinearity.r first.")
  data_rcs <- read_csv("cohort_final_baseline.csv")
}

# 2. 数据准备
data_rcs <- data_rcs %>%
  filter(time_28_day > 0) %>%
  mutate(
    status = death_28_day, # 0=Alive, 1=Dead (28-day)
    time = time_28_day     # Censored at 28 days
  ) %>%
  as.data.frame()

# 3. Set RMS environment
dd <- datadist(data_rcs)
options(datadist = 'dd')

# 4. Remove highly correlated variables (manual decision)
vars_to_remove <- c("APS_III", "BUN")
data_rcs <- data_rcs %>% select(-any_of(vars_to_remove))

# 5. 建立 Cox 模型 + RCS (Anion Gap)

fit_rcs <- cph(Surv(time, status) ~ rcs(baseline_AG, 4) + 
                 admission_age + gender + charlson_comorbidity_index + 
                 SOFA_score + fluid_input_24h + 
                 WBC + Vasopressor_use + MAP + 
                 Mechanical_ventilation_use + RRT_use + 
                 Creatinine + 
                 Sedative_use + SAPS_II, # 保留 SAPS II
               data = data_rcs, x=TRUE, y=TRUE)

# 6. 预测并画图
ag_limits <- quantile(data_rcs$baseline_AG, probs = c(0.01, 0.99), na.rm = TRUE)
ag_range <- seq(ag_limits[1], ag_limits[2], length.out=500)

pred_data <- Predict(fit_rcs, baseline_AG = ag_range, fun = exp)

min_hr_idx <- which.min(pred_data$yhat)
ag_nadir <- ag_range[min_hr_idx]

log_hazard <- pred_data$yhat
log_hr_relative <- log_hazard - log_hazard[min_hr_idx]
log_hr_lower <- pred_data$lower - log_hazard[min_hr_idx]
log_hr_upper <- pred_data$upper - log_hazard[min_hr_idx]

res_df <- data.frame(
  AG = ag_range,
  HR = exp(log_hr_relative),
  Lower = exp(log_hr_lower),
  Upper = exp(log_hr_upper)
)

sig_start_idx <- which(res_df$AG > ag_nadir & res_df$Lower > 1)[1]

if(!is.na(sig_start_idx)) {
  ag_sig_threshold <- res_df$AG[sig_start_idx]

  ag_threshold_auto <- round(ag_sig_threshold * 2) / 2
} else {
  ag_threshold_auto <- round(ag_nadir * 2) / 2
}

writeLines(as.character(ag_threshold_auto), "ag_threshold_value.txt")

plot_df <- res_df

# 7. 绘图 (ggplot2)
p <- ggplot(plot_df, aes(x = AG, y = HR)) +
  geom_line(color = "#E74C3C", linewidth = 1) + # 红色实线
  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.1, fill = "#E74C3C") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = ag_threshold_auto, linetype = "dotted", color = "#2980B9", linewidth=0.8) +
  annotate("text", x = ag_threshold_auto, y = 1.1, 
           label = paste0("Threshold: ", ag_threshold_auto), 
           color = "#2980B9", hjust = -0.1, vjust = 0) +
  labs(
    title = "Restricted Cubic Spline of Anion Gap vs. 28-Day Mortality",
    subtitle = paste("Reference: Nadir at AG =", round(ag_nadir, 1), "| Significant Risk >", ag_threshold_auto),
    x = "Baseline Anion Gap (mEq/L)",
    y = "Hazard Ratio (95% CI)"
  ) +
  theme_classic(base_size = 14) +
  coord_cartesian(ylim = c(0.1, 10))

ggsave("figure_rcs_ag.pdf", p, width = 8, height = 6)
