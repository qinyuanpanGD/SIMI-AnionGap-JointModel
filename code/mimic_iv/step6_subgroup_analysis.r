library(tidyverse)
library(JMbayes2)
library(survival)
library(nlme)
library(splines)
library(missForest)

# 1. 读取数据
cohort <- read_csv("cohort_final_baseline.csv")
long_data <- read_csv("long_data_daily.csv")

# 2. 数据准备与插补
id_data <- cohort %>%
  filter(!is.na(los_icu) & los_icu > 0) %>%
  mutate(
    status = hospital_expire_flag, 
    time = los_icu,
    fluid_input_24h = fluid_input_24h / 1000 # Convert mL to L
  )

base_vars <- c("admission_age", "gender", "charlson_comorbidity_index", 
               "SOFA_score", "fluid_input_24h", "WBC", 
               "Vasopressor_use", "MAP", "Mechanical_ventilation_use", 
               "RRT_use", "Creatinine", "Sedative_use", "SAPS_II")

set.seed(123)
data_to_impute <- id_data %>%
  select(all_of(c("stay_id", "time", "status", base_vars))) %>%
  as.data.frame() %>%
  mutate(across(where(is.character), as.factor))

if(any(sapply(data_to_impute, function(x) sum(is.na(x)) > 0))) {
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
  filter(!is.na(visit_time) & stay_id %in% id_data$stay_id) %>%
  left_join(id_data %>% select(stay_id, time), by = "stay_id") %>%
  filter(visit_time <= time & visit_time <= 7)

# 3. 定义亚组 (根据参考图片)

age_median <- round(median(id_data$admission_age, na.rm = TRUE))
id_data$age_group <- ifelse(id_data$admission_age <= age_median, paste0("<=", age_median), paste0(">", age_median))

sofa_median <- median(id_data$SOFA_score, na.rm = TRUE)
id_data$sofa_group <- ifelse(id_data$SOFA_score <= sofa_median, paste0("<=", sofa_median), paste0(">", sofa_median))

sapsii_median <- round(median(id_data$SAPS_II, na.rm = TRUE))
id_data$sapsii_group <- ifelse(id_data$SAPS_II <= sapsii_median, paste0("<=", sapsii_median), paste0(">", sapsii_median))

id_data$vaso_group <- ifelse(id_data$Vasopressor_use == 1, "Yes", "No")

id_data$mv_group <- ifelse(id_data$Mechanical_ventilation_use == 1, "Yes", "No")

id_data$rrt_group <- ifelse(id_data$RRT_use == 1, "Yes", "No")

# 4. 定义亚组分析函数（单模型显式交互项）
mcmc_params <- tryCatch({
  lines <- readLines("mcmc_params.txt")
  list(
    n_iter = as.integer(sub("n_iter=", "", lines[grep("n_iter=", lines)])),
    n_burnin = as.integer(sub("n_burnin=", "", lines[grep("n_burnin=", lines)]))
  )
}, error = function(e) {
  list(n_iter = 3000L, n_burnin = 1000L)
})

summarise_draws <- function(draws) {
  m <- mean(draws, na.rm = TRUE)
  s <- sd(draws, na.rm = TRUE)
  l <- unname(quantile(draws, 0.025, na.rm = TRUE))
  u <- unname(quantile(draws, 0.975, na.rm = TRUE))
  p <- ifelse(is.na(s) || s == 0, NA_real_, 2 * pnorm(-abs(m / s)))
  c(Mean = m, SD = s, Lower = l, Upper = u, P = p)
}

run_subgroup_jm_interaction <- function(subgroup_name, ref_value, alt_value, label_ref, label_alt, remove_var) {

  data_model <- id_data %>%
    filter(!is.na(.data[[subgroup_name]]), .data[[subgroup_name]] %in% c(ref_value, alt_value))
  if(nrow(data_model) == 0) {
    return(data.frame(
      GroupVar = subgroup_name,
      Subgroup = c(label_ref, label_alt),
      N = c(0, 0),
      Events = c(0, 0),
      HR = NA, Lower = NA, Upper = NA, P_Value = NA, P_Interaction = NA_real_
    ))
  }

  data_model <- data_model %>%
    mutate(
      subgroup_bin = ifelse(.data[[subgroup_name]] == alt_value, 1, 0)
    )

  n_ref <- sum(data_model[[subgroup_name]] == ref_value, na.rm = TRUE)
  n_alt <- sum(data_model[[subgroup_name]] == alt_value, na.rm = TRUE)
  e_ref <- sum(data_model$status[data_model[[subgroup_name]] == ref_value], na.rm = TRUE)
  e_alt <- sum(data_model$status[data_model[[subgroup_name]] == alt_value], na.rm = TRUE)

  if(e_ref < 10 || e_alt < 10) {
    return(data.frame(
      GroupVar = subgroup_name,
      Subgroup = c(label_ref, label_alt),
      N = c(n_ref, n_alt),
      Events = c(e_ref, e_alt),
      HR = NA, Lower = NA, Upper = NA, P_Value = NA, P_Interaction = NA_real_
    ))
  }

  vars_to_adj <- setdiff(base_vars, remove_var)
  f <- as.formula(paste("Surv(time, status) ~ subgroup_bin +", paste(vars_to_adj, collapse = " + ")))

  tryCatch({
    cox_fit <- coxph(f, data = data_model)

    long_model <- long_df %>%
      inner_join(data_model %>% select(stay_id, subgroup_bin), by = "stay_id")

    ctrl <- lmeControl(opt = "optim", msMaxIter = 500)

    lme_ag <- tryCatch({
      lme(aniongap ~ ns(visit_time, 3),
          random = ~ ns(visit_time, 3) | stay_id,
          data = long_model,
          na.action = na.exclude,
          control = ctrl)
    }, error = function(e) {
      tryCatch({
        lme(aniongap ~ ns(visit_time, 2),
            random = ~ ns(visit_time, 2) | stay_id,
            data = long_model,
            na.action = na.exclude,
            control = ctrl)
      }, error = function(e2) {
        lme(aniongap ~ visit_time,
            random = ~ visit_time | stay_id,
            data = long_model,
            na.action = na.exclude,
            control = ctrl)
      })
    })

    lme_fluid <- tryCatch({
      lme(fluid_input_daily ~ ns(visit_time, 3),
          random = ~ ns(visit_time, 3) | stay_id,
          data = long_model,
          na.action = na.exclude,
          control = ctrl)
    }, error = function(e) {
      tryCatch({
        lme(fluid_input_daily ~ ns(visit_time, 2),
            random = ~ ns(visit_time, 2) | stay_id,
            data = long_model,
            na.action = na.exclude,
            control = ctrl)
      }, error = function(e2) {
        lme(fluid_input_daily ~ visit_time,
            random = ~ visit_time | stay_id,
            data = long_model,
            na.action = na.exclude,
            control = ctrl)
      })
    })

    jm_fit <- jm(
      cox_fit,
      list(lme_ag, lme_fluid),
      time_var = "visit_time",
      n_iter = mcmc_params$n_iter, n_burnin = mcmc_params$n_burnin, cores = 4
    )

    jm_fit_inter <- update(
      jm_fit,
      functional_forms = list("aniongap" = ~ value(aniongap) * subgroup_bin),
      n_iter = mcmc_params$n_iter, n_burnin = mcmc_params$n_burnin, cores = 4
    )

    alphas_all <- do.call(rbind, jm_fit_inter$mcmc$alphas)
    alpha_names <- colnames(alphas_all)

    main_col <- grep("^value\\(aniongap\\)$", alpha_names, value = TRUE)
    inter_col <- grep("value\\(aniongap\\).*subgroup_bin|subgroup_bin.*value\\(aniongap\\)", alpha_names, value = TRUE)

    if(length(main_col) == 0 || length(inter_col) == 0) {
      stop("Cannot extract interaction coefficients from JM alphas")
    }

    loghr_ref <- alphas_all[, main_col[1]]
    loghr_alt <- alphas_all[, main_col[1]] + alphas_all[, inter_col[1]]

    ref_stats <- summarise_draws(loghr_ref)
    alt_stats <- summarise_draws(loghr_alt)

    surv_summ <- summary(jm_fit_inter)$Survival
    inter_row <- grep("value\\(aniongap\\).*subgroup_bin|subgroup_bin.*value\\(aniongap\\)", rownames(surv_summ))
    p_int <- if(length(inter_row) == 0) NA_real_ else surv_summ[inter_row[1], "P"]

    data.frame(
      GroupVar = subgroup_name,
      Subgroup = c(label_ref, label_alt),
      N = c(n_ref, n_alt),
      Events = c(e_ref, e_alt),
      HR = exp(c(ref_stats["Mean"], alt_stats["Mean"])),
      Lower = exp(c(ref_stats["Lower"], alt_stats["Lower"])),
      Upper = exp(c(ref_stats["Upper"], alt_stats["Upper"])),
      P_Value = c(ref_stats["P"], alt_stats["P"]),
      P_Interaction = c(p_int, p_int)
    )
  }, error = function(e) {
    data.frame(
      GroupVar = subgroup_name,
      Subgroup = c(label_ref, label_alt),
      N = c(n_ref, n_alt),
      Events = c(e_ref, e_alt),
      HR = NA, Lower = NA, Upper = NA, P_Value = NA, P_Interaction = NA_real_
    )
  })
}

# 5. 执行所有亚组分析
res_age <- run_subgroup_jm_interaction("age_group", paste0("<=", age_median), paste0(">", age_median), paste0("Age <= ", age_median), paste0("Age > ", age_median), "admission_age")
res_gender <- run_subgroup_jm_interaction("gender", "M", "F", "Male", "Female", "gender")
res_sofa <- run_subgroup_jm_interaction("sofa_group", paste0("<=", sofa_median), paste0(">", sofa_median), paste0("SOFA <= ", sofa_median), paste0("SOFA > ", sofa_median), "SOFA_score")
res_vaso <- run_subgroup_jm_interaction("vaso_group", "No", "Yes", "Vasopressor: No", "Vasopressor: Yes", "Vasopressor_use")
res_rrt <- run_subgroup_jm_interaction("rrt_group", "No", "Yes", "RRT: No", "RRT: Yes", "RRT_use")
res_mv <- run_subgroup_jm_interaction("mv_group", "No", "Yes", "MV: No", "MV: Yes", "Mechanical_ventilation_use")
res_sapsii <- run_subgroup_jm_interaction("sapsii_group", paste0("<=", sapsii_median), paste0(">", sapsii_median), paste0("SAPS II <= ", sapsii_median), paste0("SAPS II > ", sapsii_median), "SAPS_II")

final_res <- bind_rows(
  res_age %>% mutate(ord = match(Subgroup, c(paste0("Age <= ", age_median), paste0("Age > ", age_median)))),
  res_gender %>% mutate(ord = match(Subgroup, c("Male", "Female"))),
  res_sofa %>% mutate(ord = match(Subgroup, c(paste0("SOFA <= ", sofa_median), paste0("SOFA > ", sofa_median)))),
  res_vaso %>% mutate(ord = match(Subgroup, c("Vasopressor: Yes", "Vasopressor: No"))),
  res_rrt %>% mutate(ord = match(Subgroup, c("RRT: Yes", "RRT: No"))),
  res_mv %>% mutate(ord = match(Subgroup, c("MV: Yes", "MV: No"))),
  res_sapsii %>% mutate(ord = match(Subgroup, c(paste0("SAPS II <= ", sapsii_median), paste0("SAPS II > ", sapsii_median))))
) %>%
  arrange(factor(GroupVar, levels = c("age_group", "gender", "sofa_group", "vaso_group", "rrt_group", "mv_group", "sapsii_group")), ord) %>%
  select(-ord)

final_res_formatted <- final_res %>%
  mutate(
    `HR (95% CI)` = sprintf("%.3f (%.3f - %.3f)", HR, Lower, Upper),
    `P Value` = ifelse(P_Value < 0.001, "<0.001", sprintf("%.3f", P_Value)),
    `P for interaction` = ifelse(
      is.na(P_Interaction),
      "",
      ifelse(P_Interaction < 0.001, "<0.001", sprintf("%.3f", P_Interaction))
    )
  ) %>%
  select(Subgroup, N, Events, `HR (95% CI)`, `P Value`, `P for interaction`)

write_csv(final_res_formatted, "Table_S_Subgroup_Analysis.csv")
