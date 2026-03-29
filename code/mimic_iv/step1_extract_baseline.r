library(tidyverse)
library(DBI)
library(odbc)
library(dbplyr)
library(lubridate)

# 1. 数据库连接
db = "mimic4_v31"
con <- DBI::dbConnect(odbc::odbc(),
                      dsn = db,
                      database = db,
                      uid = "your_username",
                      pwd = "your_password",
                      server   = "localhost",
                      port     = 5432)

# 2. 读取队列文件
if(file.exists("cohort_test1.csv")) {
  cohort <- read_csv("cohort_test1.csv")
} else {
  stop("找不到 cohort_test1.csv，请先运行 step0_cohort_extraction.r")
}

# 3. 准备时间窗口 (D0 定义)
icu_times <- tbl(con, I("mimiciv_icu.icustays")) %>%
  filter(stay_id %in% !!cohort$stay_id) %>%
  select(stay_id, intime) %>%
  collect()

cohort_times <- cohort %>%
  inner_join(icu_times, by = "stay_id") %>%
  mutate(
    intime = as.POSIXct(intime),
    intime_hour = hour(intime),
    d0_end = if_else(intime_hour < 7,
                     as.POSIXct(paste0(as.Date(intime), " 07:00:00"), tz = "UTC"),
                     as.POSIXct(paste0(as.Date(intime) + 1, " 07:00:00"), tz = "UTC")),
    baseline_start = intime - hours(24)
  ) %>%
  select(subject_id, stay_id, hadm_id, intime, d0_end, baseline_start)

# 4. 提取基线 AG (Strict Definition)

item_na <- c(50983, 50821)
item_k <- c(50971, 50822)
item_cl <- c(50902, 50806)
item_hco3 <- c(50882, 50803)

ag_items <- c(item_na, item_k, item_cl, item_hco3)

lab_raw <- tbl(con, I("mimiciv_hosp.labevents")) %>%
  filter(hadm_id %in% !!cohort$hadm_id) %>%
  filter(itemid %in% ag_items) %>%
  select(hadm_id, charttime, itemid, valuenum) %>%
  collect()

ag_calc <- lab_raw %>%
  mutate(
    label = case_when(
      itemid %in% item_na ~ "Na",
      itemid %in% item_k ~ "K",
      itemid %in% item_cl ~ "Cl",
      itemid %in% item_hco3 ~ "HCO3"
    )
  ) %>%
  select(-itemid) %>%
  filter(!is.na(valuenum)) %>%
  pivot_wider(names_from = label, values_from = valuenum, values_fn = mean) %>% # mean in case of duplicates (rare corrections)

  # 1. 源头质控：严格按照用户指定范围过滤 (Outliers)
  filter(
    (Na >= 110 & Na <= 180) &      # Na: 110-180
    (K >= 1 & K <= 10) &           # K: 1-10
    (Cl >= 80 & Cl <= 160) &       # Cl: 80-160
    (HCO3 >= 5 & HCO3 <= 50)       # HCO3: 5-50
  ) %>%

  filter(!is.na(Na) & !is.na(K) & !is.na(Cl) & !is.na(HCO3)) %>%

  # 2. 计算 AG (含钾公式，保留钾离子)
  mutate(AG_strict = Na + K - Cl - HCO3) %>%

  # 3. 再次对计算结果进行合理性质控 (AG < 0 视为异常/缺失)
  filter(AG_strict >= 0 & AG_strict < 100) %>%

  select(hadm_id, charttime, AG_strict)

baseline_ag <- ag_calc %>%
  inner_join(cohort_times, by = "hadm_id") %>%
  filter(charttime >= baseline_start & charttime <= d0_end) %>%
  group_by(stay_id) %>%
  summarise(baseline_AG = max(AG_strict, na.rm = TRUE), .groups = "drop") %>%
  mutate(baseline_AG = ifelse(is.infinite(baseline_AG), NA, baseline_AG))

# 5. 提取基线 Troponin (Max in [intime - 24h, intime + 24h])
item_trop <- c(51003) # Troponin T

trop_raw <- tbl(con, I("mimiciv_hosp.labevents")) %>%
  filter(hadm_id %in% !!cohort$hadm_id) %>%
  filter(itemid %in% item_trop) %>%
  select(hadm_id, charttime, valuenum) %>%
  collect()

baseline_trop <- trop_raw %>%
  inner_join(cohort_times, by = "hadm_id") %>%
  filter(charttime >= intime - hours(24) & charttime <= intime + hours(24)) %>%
  group_by(stay_id) %>%
  summarise(baseline_Troponin = max(valuenum, na.rm = TRUE), .groups = "drop")

# 6. 提取基线 Fluid Input (Total in D0)
input_raw <- tbl(con, I("mimiciv_icu.inputevents")) %>%
  filter(stay_id %in% !!cohort$stay_id) %>%
  filter(amountuom %in% c("ml", "mL", "L")) %>%
  select(stay_id, starttime, amount, amountuom) %>%
  collect() %>%
  mutate(amount_ml = ifelse(amountuom == "L", amount * 1000, amount))

baseline_fluid <- input_raw %>%
  inner_join(cohort_times, by = "stay_id") %>%
  filter(starttime >= intime & starttime <= d0_end) %>%
  group_by(stay_id) %>%
  summarise(baseline_Fluid = sum(amount_ml, na.rm = TRUE), .groups = "drop")

# 7. 其他常规基线变量 (使用 First Day Derived Tables 近似，因为差异通常不大且提取成本低)

lab_derived <- tbl(con, I("mimiciv_derived.first_day_lab")) %>%
  filter(stay_id %in% !!cohort$stay_id) %>%
  select(stay_id, 
         wbc_max, hemoglobin_min, platelets_min, 
         creatinine_max, bun_max, glucose_max, pt_max, ptt_max) %>%
  collect()

bg_derived <- tbl(con, I("mimiciv_derived.first_day_bg")) %>%
  filter(stay_id %in% !!cohort$stay_id) %>%
  select(stay_id, lactate_max, ph_min, pco2_max) %>%
  collect()

vitals_derived <- tbl(con, I("mimiciv_derived.first_day_vitalsign")) %>%
  filter(stay_id %in% !!cohort$stay_id) %>%
  select(stay_id, mbp_mean, heart_rate_mean, spo2_mean, temperature_mean, resp_rate_mean) %>%
  collect()

sapsii <- tbl(con, I("mimiciv_derived.sapsii")) %>%
  filter(stay_id %in% !!cohort$stay_id) %>%
  select(stay_id, sapsii) %>%
  collect()

apsiii <- tbl(con, I("mimiciv_derived.apsiii")) %>%
  filter(stay_id %in% !!cohort$stay_id) %>%
  select(stay_id, apsiii) %>%
  collect()

charlson <- tbl(con, I("mimiciv_derived.charlson")) %>%
  filter(hadm_id %in% !!cohort$hadm_id) %>%
  select(hadm_id, myocardial_infarct, congestive_heart_failure, 
         diabetes_with_cc, diabetes_without_cc, 
         renal_disease, mild_liver_disease, severe_liver_disease,
         cerebrovascular_disease, malignant_cancer, metastatic_solid_tumor,
         charlson_comorbidity_index) %>%
  collect() %>%
  mutate(
    icd_cad = myocardial_infarct,
    icd_chf = congestive_heart_failure,
    icd_diabetes = ifelse(diabetes_with_cc==1 | diabetes_without_cc==1, 1, 0),
    icd_renal = renal_disease,
    icd_liver = ifelse(mild_liver_disease==1 | severe_liver_disease==1, 1, 0),
    icd_stroke = cerebrovascular_disease,
    icd_malignancy = ifelse(malignant_cancer==1 | metastatic_solid_tumor==1, 1, 0)
  ) %>%
  select(hadm_id, icd_cad, icd_chf, icd_diabetes, icd_renal, icd_liver, icd_stroke, icd_malignancy, charlson_comorbidity_index)

vent <- tbl(con, I("mimiciv_derived.ventilation")) %>%
  filter(stay_id %in% !!cohort$stay_id) %>%
  collect() %>%
  inner_join(cohort_times, by = "stay_id") %>%
  filter(starttime <= intime + hours(24)) %>%
  distinct(stay_id) %>%
  mutate(Mechanical_ventilation_use = 1)

vaso <- tbl(con, I("mimiciv_derived.vasoactive_agent")) %>%
  filter(stay_id %in% !!cohort$stay_id) %>%
  collect() %>%
  inner_join(cohort_times, by = "stay_id") %>%
  filter(starttime <= intime + hours(24)) %>%
  distinct(stay_id) %>%
  mutate(Vasopressor_use = 1)

rrt <- tbl(con, I("mimiciv_derived.first_day_rrt")) %>%
  filter(stay_id %in% !!cohort$stay_id) %>%
  select(stay_id, dialysis_present) %>%
  collect() %>%
  mutate(RRT_use = dialysis_present) %>%
  select(stay_id, RRT_use)

sedative <- tbl(con, I("mimiciv_icu.inputevents")) %>%
  filter(stay_id %in% !!cohort$stay_id) %>%
  filter(itemid %in% c(221668, 222168, 221744)) %>%
  select(stay_id, starttime) %>%
  collect() %>%
  inner_join(cohort_times, by = "stay_id") %>%
  filter(starttime <= intime + hours(24)) %>%
  distinct(stay_id) %>%
  mutate(Sedative_use = 1)

# 8. 合并最终基线表
final_df <- cohort %>%
  select(stay_id, subject_id, hadm_id, admission_age, gender, sofa_score, los_icu, hospital_expire_flag, dod) %>%
  inner_join(icu_times, by = "stay_id") %>% # Add intime here explicitly
  left_join(baseline_ag, by = "stay_id") %>%
  left_join(baseline_trop, by = "stay_id") %>%
  left_join(baseline_fluid, by = "stay_id") %>%
  left_join(lab_derived, by = "stay_id") %>%
  left_join(bg_derived, by = "stay_id") %>%
  left_join(vitals_derived, by = "stay_id") %>%
  left_join(sapsii, by = "stay_id") %>%
  left_join(apsiii, by = "stay_id") %>%
  left_join(charlson, by = "hadm_id") %>%
  left_join(vent, by = "stay_id") %>%
  left_join(vaso, by = "stay_id") %>%
  left_join(rrt, by = "stay_id") %>%
  left_join(sedative, by = "stay_id") %>%
  mutate(
    Mechanical_ventilation_use = replace_na(Mechanical_ventilation_use, 0),
    Vasopressor_use = replace_na(Vasopressor_use, 0),
    RRT_use = replace_na(RRT_use, 0),
    Sedative_use = replace_na(Sedative_use, 0),
    baseline_Fluid = replace_na(baseline_Fluid, 0),

    fluid_input_24h = baseline_Fluid,
    WBC = wbc_max,
    Creatinine = creatinine_max,
    MAP = mbp_mean,
    Hemoglobin = hemoglobin_min,
    Platelet = platelets_min,
    PT = pt_max,
    APTT = ptt_max,
    BUN = bun_max,
    Glucose = glucose_max,
    Lactate = lactate_max,
    pH = ph_min,
    pCO2 = pco2_max,
    Heart_rate = heart_rate_mean,
    Temperature = temperature_mean,
    Respiratory_rate = resp_rate_mean,
    SpO2 = spo2_mean,
    SAPS_II = sapsii,
     APS_III = apsiii,
     SOFA_score = sofa_score,
     Diabetes = icd_diabetes,
     Renal = icd_renal,
     Liver = icd_liver,
    CAD = icd_cad,
    Stroke = icd_stroke,
    Malignancy = icd_malignancy,
    CHF = icd_chf
  ) %>%
  mutate(
    intime = as.POSIXct(intime),
    dod = as.POSIXct(dod),
    death_28_day = ifelse(!is.na(dod) & difftime(dod, intime, units="days") <= 28, 1, 0),
    time_28_day = case_when(
      death_28_day == 1 ~ as.numeric(difftime(dod, intime, units="days")),
      TRUE ~ 28
    ),
    time_28_day = ifelse(time_28_day <= 0, 0.1, time_28_day)
  ) %>%
  select(
    stay_id, subject_id, hadm_id, admission_age, gender, 
    los_icu, hospital_expire_flag, dod, death_28_day, time_28_day,
    baseline_AG, baseline_Troponin, fluid_input_24h,
    WBC, Creatinine, MAP, Hemoglobin, Platelet, PT, APTT, BUN, Glucose, Lactate, pH, pCO2,
    Heart_rate, Temperature, Respiratory_rate, SpO2,
    SOFA_score, SAPS_II, APS_III, charlson_comorbidity_index,
    Mechanical_ventilation_use, Vasopressor_use, RRT_use, Sedative_use,
    Diabetes, Renal, Liver, CAD, Stroke, Malignancy, CHF
  )

write_csv(final_df, "cohort_final_baseline.csv")

dbDisconnect(con)
