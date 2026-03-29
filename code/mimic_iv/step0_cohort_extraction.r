rm(list = ls())
library(tidyverse)
library(DBI)
library(odbc)
library(dbplyr)

# 1. 数据库连接 (复用配置)
db = "mimic4_v31"
con <- DBI::dbConnect(odbc::odbc(),
                      dsn = db,
                      database = db,
                      uid = "your_username",
                      pwd = "your_password",
                      server   = "localhost",
                      port     = 5432)

# 2. 定义排除标准的 ICD 代码 (来自 sql/population.sql)
# 1. 我们的 Troponin 筛选现在是 "ANY record" (n=77222)，而参考图是 "cTnT > 0.1"。

icd_acs <- c('41041','41011','41091','41092','41031','41042','41012','41051','I236','41081','I219','41082','I232','I235','41010','41032','41040','41090','I238','41052')
icd_ua <- c('I25110','I25700','I25790')
icd_cardiomyopathy <- c('4254','I255','I420','I422','I421','I426','4251','42518','4258','I427','4255','42511','4257','I425','O903','B3324','67454','67452','67451','4252')
icd_myocarditis <- c('I514','3643','I400','4290','D8685','42291','7423','42290','4220','I401','I408','I409','42293','42299','I090','J1082','42292','3980')
icd_valve <- c('I350','I340','Z952','V433','V422','Z953','I351','Z954','I361','T826XXA','99602','99671','T8203XA','I371','T8209XA','T82223A','T82228A','T8201XA','T8202XA','I088','I378','I089','I091','T826XXD','T82221A','T82222A','3971','I379','T826XXS','3979','I39','T8209XD','74609')
icd_endocarditis <- c('I330','4210','I38','42490','I339','4219','11281','B376','7422','3911','7420','42491','I011','M3211','I39')
icd_cardiac_arrest <- c('4275','I469','I462','I468','I97711','I97710','I97121','I97120')
icd_severe_arrhythmia <- c('I442','I490','I4901','I4902','I495','4260','42741','42742')
icd_copd <- c('J449','J441','J440','J440')
icd_chf <- c('4280','I5032','42832','I110','I130','I5033','I5022','42822','I5023','42823','42833','I132','42843','42842','I5043','I5042','40491','I5084','I1310','39891','I5082','40493','I0981','I1311','40490','I5040','40492','40401','40201','40411','40211','40410','40403','40400','40413')

all_codes <- unique(c(icd_acs, icd_ua, icd_cardiomyopathy, icd_myocarditis, icd_valve, icd_endocarditis, icd_cardiac_arrest, icd_severe_arrhythmia, icd_copd, icd_chf, '4294'))

# 3. 构建数据库表引用
sepsis_tbl <- tbl(con, I("mimiciv_derived.sepsis3"))
icu_detail_tbl <- tbl(con, I("mimiciv_derived.icustay_detail"))
diagnoses_tbl <- tbl(con, I("mimiciv_hosp.diagnoses_icd"))
lab_ag_tbl <- tbl(con, I("mimiciv_derived.first_day_lab"))
cardiac_marker_tbl <- tbl(con, I("mimiciv_derived.cardiac_marker"))

# 4. 获取需要排除的 HADM_ID
exclude_hadm <- diagnoses_tbl %>%
  filter(icd_code %in% all_codes | icd_code %like% 'I97%') %>%
  select(hadm_id) %>%
  distinct()

# 5. 构建核心 Sepsis 队列 (Sepsis3 + Demographics + Exclusions)
cohort_raw <- sepsis_tbl %>%
  select(subject_id, stay_id, sofa_score, suspected_infection_time) %>%
  inner_join(
    icu_detail_tbl %>% select(stay_id, hadm_id, admission_age, gender, dod, hospital_expire_flag, los_icu, icu_intime),
    by = "stay_id"
  ) %>%
  filter(admission_age >= 18) %>%
  collect()

cohort_raw$icu_intime <- as.POSIXct(cohort_raw$icu_intime)
cohort_raw$suspected_infection_time <- as.POSIXct(cohort_raw$suspected_infection_time)

cohort_filtered <- cohort_raw %>%
  filter(suspected_infection_time >= icu_intime - hours(24) & 
         suspected_infection_time <= icu_intime + hours(24))

cohort_first <- cohort_filtered %>%
  group_by(subject_id) %>%
  arrange(icu_intime) %>%
  mutate(seq = row_number()) %>%
  ungroup() %>%
  filter(seq == 1)

exclude_ids <- exclude_hadm %>% collect()
cohort_step1 <- cohort_first %>%
  anti_join(exclude_ids, by = "hadm_id")

# 6. 获取关键实验室指标 (从 Labevents 获取真实的 Troponin 数值)
labevents_tbl <- tbl(con, I("mimiciv_hosp.labevents"))

trop_lab_raw <- labevents_tbl %>%
  filter(itemid == 51003) %>%
  filter(!is.na(valuenum)) %>%
  select(subject_id, hadm_id, charttime, valuenum) %>%
  collect()

trop_matched <- trop_lab_raw %>%
  inner_join(cohort_step1 %>% select(subject_id, hadm_id, icu_intime), by = c("subject_id", "hadm_id")) %>%
  mutate(
    icu_intime = as.POSIXct(icu_intime),
    charttime = as.POSIXct(charttime)
  ) %>%
  filter(charttime >= icu_intime - hours(24) & charttime <= icu_intime + hours(24))

trop_patients <- trop_matched %>%
  group_by(subject_id, hadm_id) %>%
  summarise(max_trop = max(valuenum, na.rm = TRUE), .groups = "drop") %>%
  filter(max_trop > 0.01) %>%
  select(subject_id, hadm_id)

# 6.2 阴离子间隙 (AG) - 使用 First Day Lab 数据 (基线)
ag_data <- lab_ag_tbl %>%
  select(stay_id, aniongap_max) %>%
  collect()

# 7. 最终合并 (都在本地 R 环境中进行)
final_cohort_tbl <- cohort_step1 %>%
  inner_join(trop_patients, by = c("subject_id", "hadm_id")) %>%
  left_join(ag_data, by = "stay_id") %>%
  select(subject_id, stay_id, hadm_id, admission_age, gender, sofa_score, 
         aniongap_max, hospital_expire_flag, los_icu, dod)

# 8. 执行查询并下载数据
final_df <- final_cohort_tbl

# 9. 保存结果
write_csv(final_df, "cohort_test1.csv")

dbDisconnect(con)
