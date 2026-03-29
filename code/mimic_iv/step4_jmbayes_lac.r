source("renv/activate.R")

required_packages <- c("dplyr", "readr", "ggplot2", "JMbayes2", "survival", "nlme", "splines", "missForest")

Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1"
)

missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if(length(missing_packages) > 0) {
  stop(sprintf(
    "Missing packages: %s. Please install them once in advance (e.g., renv::install()).",
    paste(missing_packages, collapse = ", ")
  ))
}

library(dplyr)
library(readr)
library(ggplot2)
library(JMbayes2)
library(survival)
library(nlme)
library(splines)
library(missForest)

# 1. 读取数据
if(file.exists("../cohort_imputed.csv")) {
  cohort <- read_csv("cohort_imputed.csv")
} else {
  cohort <- read_csv("cohort_final_baseline.csv")
}

if(file.exists("long_data_daily_lac.csv")) {
  long_data <- read_csv("long_data_daily_lac.csv")
} else {
  stop("long_data_daily_lac.csv not found.")
}

# 2. 数据准备 - 生存数据 (28-day)
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

# 3. 数据准备 - 纵向数据 (28 days)
long_df <- long_data %>%
  inner_join(id_data %>% select(stay_id, admission_age, gender), by = "stay_id") %>%
  mutate(visit_time = as.numeric(visit_day),
         fluid_input_daily = fluid_input_daily / 1000) %>% # Convert mL to L
  filter(!is.na(visit_time)) %>%
  filter(stay_id %in% id_data$stay_id) %>%
  left_join(id_data %>% select(stay_id, time), by = "stay_id") %>%
  filter(visit_time <= time) %>%
  filter(visit_time <= 28) # 改为 28 天

# 4. 构建模型

# 4.1 Cox Submodel
cox_fit <- coxph(Surv(time, status) ~ admission_age + gender + charlson_comorbidity_index + 
                   SOFA_score + fluid_input_24h + WBC + Vasopressor_use + 
                   MAP + Mechanical_ventilation_use + RRT_use + 
                   Creatinine + Sedative_use + SAPS_II,
                 data = id_data, model = TRUE)

# 4.2 LME Submodels
ctrl <- lmeControl(opt = "optim", msMaxIter = 500) # 更换为更稳健的 optim

lme_fluid <- lme(fluid_input_daily ~ visit_time,
                 random = ~ visit_time | stay_id,
                 data = long_df,
                 na.action = na.exclude,
                 control = ctrl)

# 5. Joint Model

current_n_iter <- 5000L
current_n_burnin <- 1000L
max_retries <- 5
target_rhat <- 1.05 # 目标 Rhat 阈值

fit_jm_with_retry <- function(long_data, cox_model, init_iter, init_burnin, max_r, t_rhat) {

  try_lme <- function(k_fixed, k_random) {

    if (k_fixed == 1) {
      form_fixed <- as.formula("lactate ~ visit_time")
    } else {
      form_fixed <- as.formula(paste0("lactate ~ ns(visit_time, ", k_fixed, ")"))
    }

    form_random <- as.formula(paste0("~ ns(visit_time, ", k_random, ") | stay_id"))
    if (k_random == 0) form_random <- as.formula("~ 1 | stay_id")
    if (k_random == 1) form_random <- as.formula("~ visit_time | stay_id")

    tryCatch({
      lme(form_fixed, random = form_random, data = long_data, na.action = na.exclude, control = ctrl)
    }, error = function(e) {
      return(NULL)
    })
  }

  run_jm_loop <- function(lme_mod) {
    c_iter <- init_iter
    c_burnin <- init_burnin

    jm_fit <- tryCatch({
      jm(cox_model, list(lme_mod, lme_fluid), time_var = "visit_time", n_iter = c_iter, n_burnin = c_burnin, cores = 3)
    }, error = function(e) {
      return(NULL)
    })

    if(is.null(jm_fit)) return(list(fit = NULL, success = FALSE, iter = c_iter, burnin = c_burnin, lme_mod = lme_mod))

    for(i in 1:max_r) {
      sum_obj <- summary(jm_fit)
      rhats_surv <- sum_obj$Survival[, "Rhat"]
      rhats_long <- unlist(lapply(sum_obj$Outcome, function(x) x[, "Rhat"]))
      max_rhat <- max(c(rhats_surv, rhats_long), na.rm = TRUE)

      if(max_rhat <= t_rhat) {
        return(list(fit = jm_fit, success = TRUE, iter = c_iter, burnin = c_burnin, lme_mod = lme_mod))
      }

      if(i < max_r) {
        c_iter <- c_iter + 5000L
        c_burnin <- c_burnin + 1000L
        jm_fit <- tryCatch({
          jm(cox_model, list(lme_mod, lme_fluid), time_var = "visit_time", n_iter = c_iter, n_burnin = c_burnin, cores = 3)
        }, error = function(e) {
          return(jm_fit) # 返回上一次的拟合结果
        })
      }
    }
    return(list(fit = jm_fit, success = FALSE, iter = c_iter, burnin = c_burnin, lme_mod = lme_mod))
  }

  strategies <- list(
    c(3, 3), c(3, 2), c(3, 1), c(2, 2), c(2, 1), c(2, 0), c(1, 1), c(1, 0)
  )

  for(strat in strategies) {
    lme_lac_tmp <- try_lme(strat[1], strat[2])
    if(is.null(lme_lac_tmp)) next # 如果 lme 直接报错，尝试下一个降级策略

    jm_res <- run_jm_loop(lme_lac_tmp)

    if(jm_res$success) {
      return(jm_res)
    } else {
    }
  }

  return(jm_res) # 如果全都不行，返回最后一个模型
}

res_list <- fit_jm_with_retry(long_df, cox_fit, current_n_iter, current_n_burnin, max_retries, target_rhat)
jm_fit <- res_list$fit
lme_lac_final <- res_list$lme_mod
current_n_iter <- res_list$iter
current_n_burnin <- res_list$burnin

# 6. 输出结果

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

write_csv(table2_formatted, "Table_2_JMbayes_Lac_28Days.csv")

# 7. Time-varying Effect (Figure)

knots_28 <- c(7, 14, 21)
form_splines <- list(
  "lactate" = ~ value(lactate) * ns(visit_time, k = knots_28, B = c(0, 28))
)

jm_fit_tv <- update(jm_fit, functional_forms = form_splines,
                    n_iter = current_n_iter, n_burnin = current_n_burnin, cores = 3,
                    Surv_object = cox_fit, Mixed_objects = list(lme_lac_final, lme_fluid))

time_seq <- seq(0, 28, length.out = 100)
W_time <- ns(time_seq, k = knots_28, B = c(0, 28))

mcmc_alphas <- jm_fit_tv$mcmc$alphas
alpha_main_col <- grep("^value\\(lactate\\)$", colnames(mcmc_alphas[[1]]), value = TRUE)
alpha_inter_cols <- grep("value\\(lactate\\):ns", colnames(mcmc_alphas[[1]]), value = TRUE)

alphas_all <- do.call(rbind, mcmc_alphas)

if(length(alpha_main_col) > 0) alphas_main <- alphas_all[, alpha_main_col] else alphas_main <- rep(0, nrow(alphas_all))

if(length(alpha_inter_cols) > 0) {
  alphas_inter <- alphas_all[, alpha_inter_cols]
  if(ncol(alphas_inter) < ncol(W_time)) W_time <- W_time[, 1:ncol(alphas_inter)]
  inter_term <- W_time %*% t(alphas_inter)
} else {
  inter_term <- matrix(0, nrow = 100, ncol = nrow(alphas_all))
}

pred_log_hr <- sweep(inter_term, 2, alphas_main, "+")
hr_summary <- apply(pred_log_hr, 1, function(x) {
  c(Mean = mean(x), Lower = quantile(x, 0.025), Upper = quantile(x, 0.975))
})

plot_df <- as.data.frame(t(hr_summary))
plot_df$Time = as.numeric(time_seq)
plot_df$HR = as.numeric(exp(plot_df$Mean))
plot_df$Lower_CI = as.numeric(exp(plot_df$Lower))
plot_df$Upper_CI = as.numeric(exp(plot_df$Upper))

p_tv <- ggplot(plot_df, aes(x = Time, y = HR)) +
  geom_line(color = "#E74C3C", linewidth = 1) +
  geom_ribbon(aes(ymin = Lower_CI, ymax = Upper_CI), fill = "#E74C3C", alpha = 0.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
  labs(title = "Time-Varying Hazard Ratio of Lactate (28 Days)",
       x = "Days since ICU Admission",
       y = "Hazard Ratio (95% CI)") +
  theme_classic(base_size = 14)

ggsave("Figure_2_TimeVarying_HR_Lac_28Days.pdf", p_tv, width = 8, height = 6)

# 8. MCMC Diagnostics
pdf("Figure_S4_Diagnostics_Lac.pdf", width = 12, height = 6)
alpha_name <- "value(lactate)"
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
