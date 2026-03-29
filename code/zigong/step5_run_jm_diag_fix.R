library(JMbayes2)
library(tidyverse)
library(survival)
library(nlme)
library(splines)

cohort_base <- read_csv("zigong_baseline_features.csv", show_col_types = FALSE)
simi_cohort <- read_csv("zigong_simi_cohort.csv", show_col_types = FALSE) %>% select(INP_NO, death_time_icu, follow_time_icu)
cohort <- cohort_base %>% left_join(simi_cohort, by = "INP_NO")
long_data <- read_csv("zigong_longitudinal_features.csv", show_col_types = FALSE)

id_data <- cohort %>%
  mutate(
    time = ifelse(!is.na(death_time_icu), death_time_icu / 24, follow_time_icu / 24),
    status = ifelse(!is.na(death_time_icu) & death_time_icu <= 28*24, 1, 0),
    time = ifelse(time > 28, 28, time),
    fluid_input_24h = fluid_input_24h / 1000,
    stay_id = INP_NO
  ) %>%
  filter(time > 0) %>%
  drop_na(admission_age, gender, SOFA_score, fluid_input_24h, WBC, 
          Vasopressor_use, MAP, Mechanical_ventilation_use, 
          RRT_use, Creatinine, Sedative_use)

long_df <- long_data %>%
  rename(stay_id = INP_NO, visit_time = day) %>%
  inner_join(id_data %>% select(stay_id, admission_age, gender, time), by = "stay_id") %>%
  mutate(fluid_input_daily = fluid_input_L) %>%
  filter(!is.na(visit_time)) %>%
  filter(visit_time <= time) %>%
  filter(visit_time <= 28)

cox_fit <- coxph(Surv(time, status) ~ admission_age + gender + 
                   SOFA_score + fluid_input_24h + WBC + Vasopressor_use + 
                   MAP + Mechanical_ventilation_use + RRT_use + 
                   Creatinine + Sedative_use,
                 data = id_data)

ctrl <- lmeControl(opt = "optim", msMaxIter = 1000)

lme_ag <- lme(aniongap ~ ns(visit_time, 3),
              random = ~ visit_time | stay_id,
              data = long_df,
              na.action = na.exclude,
              control = ctrl)
lme_fluid <- lme(fluid_input_daily ~ visit_time,
                 random = ~ visit_time | stay_id,
                 data = long_df,
                 na.action = na.exclude,
                 control = ctrl)

jmf_ag <- jm(cox_fit, list(lme_ag, lme_fluid), time_var = "visit_time", n_iter = 1000L, n_burnin = 200L, save_random_effects = TRUE)

pdf("Figure_S11_Diagnostics_AG_Zigong.pdf", width = 10, height = 8)
par(mfrow = c(2, 2))
JMbayes2::traceplot(jmf_ag, parm = "alphas")
JMbayes2::ggdensityplot(jmf_ag, parm = "alphas")
dev.off()

lme_lac <- lme(lactate ~ ns(visit_time, 3),
              random = ~ ns(visit_time, 3) | stay_id,
              data = long_df,
              na.action = na.exclude,
              control = ctrl)

jmf_lac <- jm(cox_fit, list(lme_lac, lme_fluid), time_var = "visit_time", n_iter = 1000L, n_burnin = 200L, save_random_effects = TRUE)

pdf("Figure_S12_Diagnostics_Lac_Zigong.pdf", width = 10, height = 8)
par(mfrow = c(2, 2))
JMbayes2::traceplot(jmf_lac, parm = "alphas")
JMbayes2::ggdensityplot(jmf_lac, parm = "alphas")
dev.off()

