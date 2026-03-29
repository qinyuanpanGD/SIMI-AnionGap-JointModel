library(tidyverse)
library(JMbayes2)
library(survival)
library(nlme)
library(splines)
library(missForest)

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
               "SAPS_II") # 新增 SAPS II

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
  inner_join(id_data %>% select(stay_id, admission_age, gender), by = "stay_id") %>%
  mutate(
    visit_time = as.numeric(visit_day), # Use visit_day
    fluid_input_daily = fluid_input_daily / 1000 # Convert mL to L
  ) %>%
  filter(!is.na(visit_time)) %>%
  filter(stay_id %in% id_data$stay_id) %>%
  left_join(id_data %>% select(stay_id, time), by = "stay_id") %>%
  filter(visit_time <= time) %>%
  filter(visit_time <= 28)

# 4. 构建模型

# 4.1 生存子模型 (Cox Model)
cox_fit <- coxph(Surv(time, status) ~ admission_age + gender + charlson_comorbidity_index + 
                   SOFA_score + fluid_input_24h + WBC + Vasopressor_use + 
                   MAP + Mechanical_ventilation_use + RRT_use + 
                   Creatinine + Sedative_use + SAPS_II, # 加入 SAPS II
                 data = id_data)

# 4.2 纵向子模型 (LME Submodels)

ctrl <- lmeControl(opt = "nlminb", msMaxIter = 200)

lme_ag <- lme(aniongap ~ ns(visit_time, 3),
              random = ~ ns(visit_time, 3) | stay_id,
              data = long_df,
              na.action = na.exclude,
              control = ctrl)

lme_fluid <- lme(fluid_input_daily ~ visit_time,
                 random = ~ visit_time | stay_id,
                 data = long_df,
                 na.action = na.exclude,
                 control = ctrl)

# 5. 联合模型 (Joint Model) - JMbayes2

current_n_iter <- 5000L
current_n_burnin <- 1000L
max_retries <- 3
target_rhat <- 1.05 # 收紧至 1.05

jm_fit <- jm(cox_fit, list(lme_ag, lme_fluid),
             time_var = "visit_time",
             n_iter = current_n_iter, n_burnin = current_n_burnin)

for(i in 1:max_retries) {
  sum_obj <- summary(jm_fit)
  rhats_surv <- sum_obj$Survival[, "Rhat"]
  rhats_long <- unlist(lapply(sum_obj$Outcome, function(x) x[, "Rhat"]))
  all_rhats <- c(rhats_surv, rhats_long)
  max_rhat <- max(all_rhats, na.rm = TRUE)

  if(max_rhat <= target_rhat) {
    break
  } else {
    if(i == max_retries) {
    } else {
      current_n_iter <- current_n_iter + 5000L
      current_n_burnin <- current_n_burnin + 1000L
      jm_fit <- jm(cox_fit, list(lme_ag, lme_fluid),
                   time_var = "visit_time",
                   n_iter = current_n_iter, n_burnin = current_n_burnin)
    }
  }
}

# 6. 输出结果 (Table 2 - 28 Days)

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
  ) %>%
  select(Variable, `HR (95% CI)`, `P Value`)

write_csv(table2_formatted, "Table_2_JMbayes_Main_28Days.csv")

# 7. 时变系数模型 (Figure 2 - Time-varying Effect)

knots_28 <- c(7, 14, 21)

form_splines <- list(
  "aniongap" = ~ value(aniongap) * ns(visit_time, k = knots_28, B = c(0, 28))
)

jm_fit_tv <- update(jm_fit, functional_forms = form_splines,
                    n_iter = current_n_iter, n_burnin = current_n_burnin)

# 8. 绘制时变 HR

# 1. 准备预测时间点
time_seq <- seq(0, 28, length.out = 100)

# 2. 构建样条基函数矩阵
W_time <- ns(time_seq, k = knots_28, B = c(0, 28))

# 3. 提取 MCMC 样本 (Alphas)
mcmc_alphas <- jm_fit_tv$mcmc$alphas
chains <- length(mcmc_alphas)

alpha_main_col <- grep("^value\\(aniongap\\)$", colnames(mcmc_alphas[[1]]), value = TRUE)
alpha_inter_cols <- grep("value\\(aniongap\\):ns", colnames(mcmc_alphas[[1]]), value = TRUE)

# 4. 计算每个时间点的 HR
alphas_all <- do.call(rbind, mcmc_alphas)

if(length(alpha_main_col) > 0) {
  alphas_main <- alphas_all[, alpha_main_col]
} else {
  alphas_main <- rep(0, nrow(alphas_all))
}

if(length(alpha_inter_cols) > 0) {
  alphas_inter <- alphas_all[, alpha_inter_cols]
  if(ncol(alphas_inter) != ncol(W_time)) {
    if(ncol(alphas_inter) < ncol(W_time)) W_time <- W_time[, 1:ncol(alphas_inter)]
  }
  inter_term <- W_time %*% t(alphas_inter)
} else {
  inter_term <- matrix(0, nrow = 100, ncol = nrow(alphas_all))
}

pred_log_hr <- sweep(inter_term, 2, alphas_main, "+")

# 5. 计算汇总统计量
hr_summary <- apply(pred_log_hr, 1, function(x) {
  c(Mean = mean(x), 
    Lower = quantile(x, 0.025), 
    Upper = quantile(x, 0.975))
})

# 6. 绘图
plot_df <- as.data.frame(t(hr_summary))
plot_df$Time = as.numeric(time_seq)
plot_df$HR = as.numeric(exp(plot_df$Mean))
plot_df$Lower_CI = as.numeric(exp(plot_df$Lower))
plot_df$Upper_CI = as.numeric(exp(plot_df$Upper))

p_tv <- ggplot(plot_df, aes(x = Time, y = HR)) +
  geom_line(color = "#E74C3C", linewidth = 1) +
  geom_ribbon(aes(ymin = Lower_CI, ymax = Upper_CI), fill = "#E74C3C", alpha = 0.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
  labs(title = "Time-Varying Hazard Ratio of Anion Gap (28 Days)",
       x = "Days since ICU Admission",
       y = "Hazard Ratio (95% CI)") +
  theme_classic(base_size = 14)

ggsave("Figure_2_TimeVarying_HR_AG.pdf", p_tv, width = 8, height = 6)

# 9. MCMC Diagnostics
pdf("Figure_S4_Diagnostics_AG.pdf", width = 12, height = 6)
alpha_name <- "value(aniongap)"
chain_values <- lapply(jm_fit$mcmc$alphas, function(x) as.numeric(x[, alpha_name]))
chain_colors <- c("#FF4D4D", "#4A90E2", "#7ED957")
fill_colors <- grDevices::adjustcolor(chain_colors, alpha.f = 0.35)
n_iter <- length(chain_values[[1]])

par(mfrow = c(1, 2), mar = c(5, 4, 4, 1) + 0.1)

plot(1:n_iter, chain_values[[1]], type = "l", col = chain_colors[1], lwd = 1,
     xlab = "Iteration", ylab = "Value", main = paste0("Traceplot of ", alpha_name))
if (length(chain_values) >= 2) lines(1:n_iter, chain_values[[2]], col = chain_colors[2], lwd = 1)
if (length(chain_values) >= 3) lines(1:n_iter, chain_values[[3]], col = chain_colors[3], lwd = 1)
legend("top",
       legend = c("Chain", "Chain 1", "Chain 2", "Chain 3"),
       col = c(NA, chain_colors),
       lty = c(NA, 1, 1, 1),
       lwd = c(NA, 1, 1, 1),
       bty = "n",
       horiz = TRUE)

density_list <- lapply(chain_values, density)
x_lim <- range(unlist(lapply(density_list, `[[`, "x")))
y_lim <- c(0, max(unlist(lapply(density_list, `[[`, "y"))) * 1.05)
plot(density_list[[1]], type = "n", xlim = x_lim, ylim = y_lim,
     xlab = "Value", ylab = "Density", main = paste0("Density plot of ", alpha_name))
for (i in seq_along(density_list)) {
  polygon(c(density_list[[i]]$x, rev(density_list[[i]]$x)),
          c(density_list[[i]]$y, rep(0, length(density_list[[i]]$y))),
          col = fill_colors[i], border = NA)
}
for (i in seq_along(density_list)) {
  lines(density_list[[i]], col = chain_colors[i], lwd = 1)
}
legend("top",
       legend = c("Chain", "Chain 1", "Chain 2", "Chain 3"),
       pch = c(NA, 15, 15, 15),
       col = c(NA, chain_colors),
       pt.cex = c(NA, 1.2, 1.2, 1.2),
       bty = "n",
       horiz = TRUE)
dev.off()
