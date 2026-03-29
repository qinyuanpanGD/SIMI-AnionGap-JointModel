#   1. Sepsis 3.0 (Is_Sepsis == 1)
#   2. 年龄 >= 18
#   3. ICU 停留时间 >= 24h
#   4. 仅首次入 ICU
#   5. 排除既往心脏病史 (心梗, 心衰, 心肌病等) 和 COPD
#   6. SIMI 确诊: 入 ICU 前后 24 小时内 cTnI > 0.01
#   7. 具有明确的结局数据

rm(list = ls())
library(data.table)
library(dplyr)
library(stringr)

# 1. 设置路径
db_path <- "./"
out_path <- "./"

sepsis_df <- fread(paste0(db_path, "Derived_Sepsis.csv"))
baseline_df <- fread(paste0(db_path, "dtBaseline.csv"))
icd_df <- fread(paste0(db_path, "dtICD.csv"))
outcome_df <- fread(paste0(db_path, "dtOutCome.csv"))

system(paste0("head -n 1 ", db_path, "dtLab.csv > /tmp/zigong_troponin.csv"))
system(paste0("grep -i 'troponin' ", db_path, "dtLab.csv >> /tmp/zigong_troponin.csv"))
trop_df <- fread("/tmp/zigong_troponin.csv")

setnames(sepsis_df, "patient_id", "PATIENT_ID")
baseline_df[, ICU_discharge_time := NULL] # 移除重复列

# 2. 合并数据与初步 Sepsis 筛选
cohort <- merge(sepsis_df, baseline_df, by = c("PATIENT_ID", "INP_NO"), all.x = TRUE)

cohort <- cohort[Is_Sepsis == 1]
cohort <- cohort[Sepsis_Onset_Time >= (ICU_admission_time - 24) & Sepsis_Onset_Time <= (ICU_admission_time + 24)]

cohort <- cohort[Age >= 18]

cohort <- cohort[(ICU_discharge_time - ICU_admission_time) >= 24]

cohort <- cohort[order(PATIENT_ID, ICU_admission_time)]
cohort <- cohort[!duplicated(PATIENT_ID)]

# 3. 排除既往心脏病史与 COPD
exclusion_regex <- "^(I2[0-5]|I50|I4[0-3]|I3[3-9]|I46|J44|I51\\.4)"
desc_regex <- "(?i)(myocardial infarction|acute coronary|unstable angina|heart failure|cardiomyopathy|myocarditis|valvular|endocarditis|cardiac arrest|chronic obstructive pulmonary|COPD)"

excl_icd <- icd_df[str_detect(ICD_Code, exclusion_regex) | str_detect(ICD_desc, desc_regex)]
excl_inp_no <- unique(excl_icd$INP_NO)

cohort <- cohort[!INP_NO %in% excl_inp_no]

# 4. SIMI 确诊 (cTnI > 0.01 ng/mL, 入 ICU 前后 24 小时内)
trop_merged <- merge(trop_df, cohort[, .(INP_NO, ICU_admission_time)], by = "INP_NO", nomatch = NULL)
trop_window <- trop_merged[LabTime >= (ICU_admission_time - 24) & LabTime <= (ICU_admission_time + 24)]

trop_window[, num_value := as.numeric(gsub("[^0-9\\.]", "", LabValue))]
trop_max <- trop_window[, .(max_cTnI = max(num_value, na.rm = TRUE)), by = INP_NO]
simi_inp_no <- trop_max[max_cTnI > 0.01]$INP_NO

cohort_simi <- cohort[INP_NO %in% simi_inp_no]

# 5. 整合结局指标 (28天死亡率)
cohort_simi <- merge(cohort_simi, outcome_df[, .(PATIENT_ID, Death_Date, Follow_Vital, Follow_Date)], by = "PATIENT_ID", all.x = TRUE)

cohort_simi[, death_time_icu := Death_Date - ICU_admission_time]
cohort_simi[, follow_time_icu := Follow_Date - ICU_admission_time]

cohort_simi[, mortality_28d := ifelse(!is.na(death_time_icu) & death_time_icu <= 28 * 24, 1, 0)]
valid_outcome <- cohort_simi[!is.na(mortality_28d)]

    sprintf("(%.1f%%)\n", 100 * sum(valid_outcome$mortality_28d == 1) / nrow(valid_outcome)))

# 6. 保存最终队列数据
fwrite(valid_outcome, paste0(out_path, "zigong_simi_cohort.csv"))
