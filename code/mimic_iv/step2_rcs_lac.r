library(tidyverse)
library(survival)
library(rms)
library(ggplot2)

# 1. 读取数据 (使用插补后的基线文件，如果存在)

if(file.exists("../cohort_imputed.csv")) {
  data_rcs <- read_csv("cohort_imputed.csv")
} else if(file.exists("../cohort_final_baseline.csv")) {
  data_rcs <- read_csv("cohort_final_baseline.csv")
} else {
  stop("Data file not found.")
}

# 2. 数据准备 (28-day Mortality)

  if(file.exists("long_data_daily_lac.csv")) {
    long_lac <- read_csv("long_data_daily_lac.csv")

    day0_lac <- long_lac %>%
      filter(visit_day == 0) %>%
      select(stay_id, lactate_day0 = lactate) %>%
      filter(!is.na(lactate_day0))

    data_rcs <- data_rcs %>%
      left_join(day0_lac, by = "stay_id") %>%
      mutate(
        Lactate = ifelse(!is.na(lactate_day0), lactate_day0, Lactate)
      ) %>%
      select(-lactate_day0)

  } else {
  warning("long_data_daily_lac.csv not found. Using original baseline Lactate.")
}

data_rcs <- data_rcs %>%
  filter(time_28_day > 0) %>%
  mutate(
    status = death_28_day,
    time = time_28_day
  ) %>%
  filter(!is.na(Lactate) & Lactate > 0 & Lactate < 30)

# 3. 设置 RMS 环境
dd <- datadist(data_rcs)
options(datadist = 'dd')

# 4. 建立 Cox 模型 + RCS (Lactate)

fit_rcs <- cph(Surv(time, status) ~ rcs(Lactate, 4) + 
                 admission_age + gender + charlson_comorbidity_index + 
                 SOFA_score + rcs(fluid_input_24h, 3) + 
                 WBC + Vasopressor_use + MAP + 
                 Mechanical_ventilation_use + RRT_use + 
                 Creatinine + 
                 Sedative_use + SAPS_II,
               data = data_rcs, x=TRUE, y=TRUE)

# 5. 预测并画图
lac_limits <- quantile(data_rcs$Lactate, probs = c(0.01, 0.99), na.rm = TRUE)
lac_max_plot <- min(lac_limits[2], 15) 
lac_range <- seq(lac_limits[1], lac_max_plot, length.out=500)

pred_data <- Predict(fit_rcs, Lactate = lac_range, fun = exp)

min_hr_idx <- which.min(pred_data$yhat)
lac_nadir <- lac_range[min_hr_idx]

log_hazard <- pred_data$yhat
log_hr_relative <- log_hazard - log_hazard[min_hr_idx]
log_hr_lower <- pred_data$lower - log_hazard[min_hr_idx]
log_hr_upper <- pred_data$upper - log_hazard[min_hr_idx]

res_df <- data.frame(
  Lactate = lac_range,
  HR = exp(log_hr_relative),
  Lower = exp(log_hr_lower),
  Upper = exp(log_hr_upper)
)

sig_start_idx <- which(res_df$Lactate > lac_nadir & res_df$Lower > 1)[1]

if(!is.na(sig_start_idx)) {
  lac_sig_threshold <- res_df$Lactate[sig_start_idx]
  lac_threshold_auto <- round(lac_sig_threshold, 0)
} else {
  lac_threshold_auto <- round(lac_nadir, 0)
}

writeLines(as.character(lac_threshold_auto), "lac_threshold_value.txt")

# 6. 绘图
p <- ggplot(res_df, aes(x = Lactate, y = HR)) +
  geom_line(color = "#E74C3C", linewidth = 1) + 
  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.1, fill = "#E74C3C") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = lac_threshold_auto, linetype = "dotted", color = "#2980B9", linewidth=0.8) +
  annotate("text", x = lac_threshold_auto, y = max(res_df$Upper, na.rm=T) * 0.8, 
           label = paste0("Threshold: ", lac_threshold_auto), 
           color = "#2980B9", hjust = -0.1, vjust = 0) +
  labs(
    title = "Restricted Cubic Spline of Lactate vs. 28-Day Mortality",
    subtitle = paste("Reference: Nadir at", round(lac_nadir, 1)),
    x = "Baseline Lactate (mmol/L)",
    y = "Hazard Ratio (95% CI)"
  ) +
  theme_classic(base_size = 14) +
  coord_cartesian(ylim = c(0.1, 10))

ggsave("figure_rcs_lac.pdf", p, width = 8, height = 6)
