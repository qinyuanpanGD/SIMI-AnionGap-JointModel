if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, data.table, WeightIt, cobalt, survey, broom, survival, missForest)

# 1. 读取数据 ------------------------------------------------------------
cohort <- read_csv("cohort_final_baseline.csv", show_col_types = FALSE)

# 2. 数据准备 ------------------------------------------------------------------
base_vars <- c("admission_age", "gender", "charlson_comorbidity_index", 
               "SOFA_score", "fluid_input_24h", "WBC", 
               "Vasopressor_use", "MAP", "Mechanical_ventilation_use", 
               "RRT_use", "Creatinine", "Sedative_use",
               "SAPS_II")

clean_var_names <- c(
  "admission_age" = "Age (years)",
  "gender" = "Gender",
  "gender_M" = "Male Gender",
  "gender_F" = "Female Gender",
  "charlson_comorbidity_index" = "Charlson Comorbidity Index",
  "SOFA_score" = "SOFA Score",
  "fluid_input_24h" = "24h Fluid Input (L)",
  "WBC" = "White Blood Cell",
  "Vasopressor_use" = "Vasopressor Use",
  "MAP" = "Mean Arterial Pressure",
  "Mechanical_ventilation_use" = "Mechanical Ventilation",
  "RRT_use" = "Renal Replacement Therapy",
  "Creatinine" = "Creatinine",
  "Sedative_use" = "Sedative Use",
  "SAPS_II" = "SAPS II Score",
  "prop.score" = "Propensity Score"
)

if(!file.exists("ag_threshold_value.txt")) {
  stop("Error: 'ag_threshold_value.txt' not found. 请先运行 RCS 脚本生成阈值。")
}
ag_threshold <- as.numeric(readLines("ag_threshold_value.txt", warn=FALSE)[1])

if(!"baseline_AG" %in% names(cohort)) {
  stop("Error: 'baseline_AG' not found in cohort data.")
}

data_analysis <- cohort %>%
  filter(!is.na(baseline_AG)) %>%
  filter(!is.na(death_28_day)) %>%
  mutate(fluid_input_24h = fluid_input_24h / 1000) %>% 
  select(stay_id, baseline_AG, death_28_day, all_of(base_vars)) %>%
  mutate(
    exposure = ifelse(baseline_AG >= ag_threshold, 1, 0),
    outcome = death_28_day,
    gender = as.factor(gender),
    Vasopressor_use = as.factor(Vasopressor_use),
    Mechanical_ventilation_use = as.factor(Mechanical_ventilation_use),
    RRT_use = as.factor(RRT_use),
    Sedative_use = as.factor(Sedative_use)
  )

# 3. 缺失值插补 (Random Forest) ------------------------------------------------
set.seed(123) # 设置随机种子以保证结果可重复

data_to_impute <- data_analysis %>% select(-stay_id) %>% as.data.frame()

if(any(sapply(data_to_impute, function(x) sum(is.na(x)) > 0))) {
  rf_imputed <- missForest(data_to_impute, verbose = FALSE)
  data_imputed <- rf_imputed$ximp
  data_imputed$stay_id <- data_analysis$stay_id
} else {
  data_imputed <- data_analysis
}

write_csv(data_imputed, "cohort_imputed_causal.csv")

# 4. 计算 Overlap Weights (OW) -------------------------------------------------

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

# 5. 平衡性检验 (可视化 & 表格) -------------------------------------------------------

var_map_list <- setNames(covariates, new_names)

clean_var_names_list <- sapply(var_map_list, function(x) {
  if (x %in% names(clean_var_names)) return(clean_var_names[[x]])
  return(x)
})
names(clean_var_names_list) <- names(var_map_list)

new_var_names_df <- data.frame(
  old = c(new_names, paste0(new_names[2], "_M")), # 针对gender_M
  new = c(clean_var_names_list, clean_var_names["gender_M"])
)

p_love <- love.plot(W.out, 
                    thresholds = c(m = .1), # 添加 0.1 阈值线
                    binary = "std",
                    abs = TRUE, # 取绝对值
                    var.names = new_var_names_df, # 使用美化后的变量名
                    line = TRUE,
                    stars = "raw",
                    colors = c("#0072B2", "#D55E00"), # 采用更学术配色的红蓝
                    shapes = c("circle", "triangle"),
                    sample.names = c("Unadjusted", "OW Weighted"),
                    title = "Covariate Balance (Love Plot) - Overlap Weighting") +
  theme_bw() + # 使用白底黑框主题，更适合学术论文
  theme(
    legend.position = "top",
    text = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

pdf("Figure_7_OW_Balance_Plot.pdf", width = 10, height = 8)
dev.off()

p_bal <- bal.plot(W.out, var.name = "prop.score", which = "both") +
  labs(title = "Propensity Score Distribution (Before vs After OW)", 
       x = "Propensity Score") +
  scale_fill_manual(values = c("#0072B2", "#D55E00"), name = "Group") +
  theme_bw() +
  theme(
    legend.position = "top",
    text = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

pdf("Figure_7_OW_PS_Distribution.pdf", width = 8, height = 6)
dev.off()

bal_tab <- bal.tab(W.out, un = TRUE, 
                   stats = c("m"), 
                   disp = c("means", "sds"),
                   binary = "std")

bal_df <- bal_tab$Balance

bal_df$Variable <- rownames(bal_df)

for(i in 1:nrow(var_map)) {
  bal_df$Variable[bal_df$Variable == var_map$New[i]] <- var_map$Original[i]
  bal_df$Variable <- gsub(paste0("^", var_map$New[i], "_"), paste0(var_map$Original[i], "_"), bal_df$Variable)
}

bal_df$Variable <- sapply(bal_df$Variable, function(x) {
  if(x %in% names(clean_var_names)) return(clean_var_names[[x]])
  return(x)
})

weighted_table <- bal_df %>%
  select(Variable, 
         Mean_Control_Un = M.0.Un, SD_Control_Un = SD.0.Un,
         Mean_Treated_Un = M.1.Un, SD_Treated_Un = SD.1.Un,
         SMD_Un = Diff.Un,
         Mean_Control_W = M.0.Adj, SD_Control_W = SD.0.Adj,
         Mean_Treated_W = M.1.Adj, SD_Treated_W = SD.1.Adj,
         SMD_W = Diff.Adj) %>%
  mutate(
    `Control (Unadj)` = sprintf("%.2f (%.2f)", Mean_Control_Un, SD_Control_Un),
    `Treated (Unadj)` = sprintf("%.2f (%.2f)", Mean_Treated_Un, SD_Treated_Un),
    `SMD (Unadj)` = sprintf("%.3f", abs(SMD_Un)),

    `Control (Weighted)` = sprintf("%.2f (%.2f)", Mean_Control_W, SD_Control_W),
    `Treated (Weighted)` = sprintf("%.2f (%.2f)", Mean_Treated_W, SD_Treated_W),
    `SMD (Weighted)` = sprintf("%.3f", abs(SMD_W))
  ) %>%
  select(Variable, `Control (Unadj)`, `Treated (Unadj)`, `SMD (Unadj)`,
         `Control (Weighted)`, `Treated (Weighted)`, `SMD (Weighted)`)

write_csv(weighted_table, "Table_7_Weighted_Baseline.csv")

# 6. 结局分析 (加权回归) -------------------------------------------------------

design_ow <- svydesign(ids = ~1, weights = ~W.out$weights, data = data_imputed)

fit_ow <- svyglm(outcome ~ exposure, design = design_ow, family = quasibinomial())

res_ow <- tidy(fit_ow, conf.int = TRUE, exponentiate = TRUE) %>%
  filter(term == "exposure1") %>% # 注意提取因子变量为1时的结果
  select(estimate, conf.low, conf.high, p.value) %>%
  rename(OR = estimate, CI_Lower = conf.low, CI_Upper = conf.high, P_Value = p.value) %>%
  mutate(Method = "Overlap Weighting (OW)")

if(nrow(res_ow) == 0) {
  res_ow <- tidy(fit_ow, conf.int = TRUE, exponentiate = TRUE) %>%
    filter(term == "exposure") %>% 
    select(estimate, conf.low, conf.high, p.value) %>%
    rename(OR = estimate, CI_Lower = conf.low, CI_Upper = conf.high, P_Value = p.value) %>%
    mutate(Method = "Overlap Weighting (OW)")
}

write_csv(res_ow, "Table_7_OW_Analysis_Result.csv")
