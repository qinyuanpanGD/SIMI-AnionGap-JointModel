if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, data.table, WeightIt, survey, broom, survival)

# 1. 读取插补后的数据 ------------------------------------------------------------
if(!file.exists("cohort_imputed_causal.csv")) {
  stop("Error: 'cohort_imputed_causal.csv' not found. Please run Step 7 first.")
}
data_imputed <- read_csv("cohort_imputed_causal.csv", show_col_types = FALSE)

# 2. 数据准备 ------------------------------------------------------------------
base_vars <- c("admission_age", "gender", "charlson_comorbidity_index", 
               "SOFA_score", "fluid_input_24h", "WBC", 
               "Vasopressor_use", "MAP", "Mechanical_ventilation_use", 
               "RRT_use", "Creatinine", "Sedative_use",
               "SAPS_II")

if(!file.exists("ag_threshold_value.txt")) {
  stop("Error: 'ag_threshold_value.txt' not found. 请先运行 RCS 脚本生成阈值。")
}
ag_threshold <- as.numeric(readLines("ag_threshold_value.txt", warn=FALSE)[1])

if("baseline_AG" %in% names(data_imputed)) {
   data_imputed <- data_imputed %>%
     mutate(exposure = ifelse(baseline_AG >= ag_threshold, 1, 0))
}

data_imputed <- data_imputed %>%
  mutate(
    gender = as.factor(gender),
    Vasopressor_use = as.factor(Vasopressor_use),
    Mechanical_ventilation_use = as.factor(Mechanical_ventilation_use),
    RRT_use = as.factor(RRT_use),
    Sedative_use = as.factor(Sedative_use),
    outcome = as.factor(outcome),
    exposure = as.factor(exposure) # 0为对照，1为处理
  )

# 3. 计算 Overlap Weights (OW) -------------------------------------------------

covariates <- base_vars
data_for_weight <- data_imputed %>%
  select(exposure, outcome, all_of(covariates))

new_names <- paste0("V", 1:length(covariates))
names(data_for_weight)[3:ncol(data_for_weight)] <- new_names

var_map <- data.frame(Original = covariates, New = new_names)

ps_formula <- as.formula(paste("exposure ~", paste(new_names, collapse = " + ")))

W.out <- weightit(ps_formula, data = data_for_weight, 
                 estimand = "ATO", 
                 method = "gbm", 
                 stop.method = "es.mean")

summary(W.out)

# 4. 双重稳健估计 (Outcome Model) ----------------------------------------------

design_dr <- svydesign(ids = ~1, weights = ~W.out$weights, data = data_imputed)

covariates_str <- paste(covariates, collapse = " + ")
dr_formula_str <- paste("outcome ~ exposure +", covariates_str)
dr_formula <- as.formula(dr_formula_str)

fit_dr <- svyglm(dr_formula, design = design_dr, family = quasibinomial())

# 5. 提取并保存结果 ------------------------------------------------------------------
res_dr <- tidy(fit_dr, conf.int = TRUE, exponentiate = TRUE) %>%
  filter(term == "exposure1") %>% # 提取暴露因子为1 (即高AG组) 的效应
  select(estimate, conf.low, conf.high, p.value) %>%
  rename(OR = estimate, CI_Lower = conf.low, CI_Upper = conf.high, P_Value = p.value) %>%
  mutate(Method = "Doubly Robust (OW + Adjustment)")

if(nrow(res_dr) == 0) {
  res_dr <- tidy(fit_dr, conf.int = TRUE, exponentiate = TRUE) %>%
    filter(term == "exposure") %>% 
    select(estimate, conf.low, conf.high, p.value) %>%
    rename(OR = estimate, CI_Lower = conf.low, CI_Upper = conf.high, P_Value = p.value) %>%
    mutate(Method = "Doubly Robust (OW + Adjustment)")
}

write_csv(res_dr, "Table_8_DR_Analysis_Result.csv")

# 6. 生成森林图展示 DR 结果 ------------------------------------------------------

plot_data <- res_dr %>%
  mutate(
    Label = "High vs Low Anion Gap",
    Significance = case_when(
      P_Value < 0.001 ~ "***",
      P_Value < 0.01  ~ "**",
      P_Value < 0.05  ~ "*",
      TRUE            ~ ""
    ),
    OR_Text = sprintf("%.2f (%.2f - %.2f)%s", OR, CI_Lower, CI_Upper, Significance)
  )

p_forest <- ggplot(plot_data, aes(y = Label, x = OR)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50", linewidth = 1) +
  geom_errorbarh(aes(xmin = CI_Lower, xmax = CI_Upper), height = 0.2, color = "#0072B2", linewidth = 1.2) +
  geom_point(color = "#D55E00", size = 4) +
  geom_text(aes(label = OR_Text), vjust = -1.5, size = 5) +
  scale_x_continuous(limits = c(min(plot_data$CI_Lower) * 0.8, max(plot_data$CI_Upper) * 1.2)) +
  labs(
    title = "Forest Plot of 28-day Mortality",
    subtitle = "Doubly Robust Estimation (Overlap Weighting + Adjusted Regression)",
    x = "Odds Ratio (95% CI)",
    y = ""
  ) +
  theme_classic() + # 使用经典的学术主题
  theme(
    text = element_text(size = 14, family = "sans"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, color = "gray30"),
    axis.text.y = element_text(face = "bold", size = 12),
    axis.line = element_line(linewidth = 0.8),
    axis.ticks = element_line(linewidth = 0.8)
  )

pdf("Figure_8_DR_Forest_Plot.pdf", width = 8, height = 4)
dev.off()

