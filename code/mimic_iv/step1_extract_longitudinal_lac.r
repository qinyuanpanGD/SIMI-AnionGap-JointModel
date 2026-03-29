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

# 2. 读取基线队列 (确保读取的是最新的含28天结局的文件)
if(file.exists("../cohort_final_baseline.csv")) {
  cohort <- read_csv("cohort_final_baseline.csv")
} else if(file.exists("cohort_final_baseline.csv")) {
  cohort <- read_csv("cohort_final_baseline.csv")
} else {
  stop("找不到 cohort_final_baseline.csv，请先运行 step1_extract_baseline.r")
}

# 3. 准备时间映射表
icu_times <- tbl(con, I("mimiciv_icu.icustays")) %>%
  filter(stay_id %in% !!cohort$stay_id) %>%
  select(stay_id, intime) %>%
  collect() %>%
  mutate(
    intime = as.POSIXct(intime),
    d0_end = as.POSIXct(paste0(as.Date(intime) + 1, " 07:00:00"), tz = "UTC")
  )

# 4. 提取原始数据 (Raw Data with Timestamps)

# 4.1 Lactate (乳酸)

# 4.1.1 从 labevents 提取 (使用 subject_id 避免漏掉急诊/入ICU前没有 hadm_id 的记录)
item_lac_lab <- c(50813, 53154) # Blood Gas and Chemistry Lactate
lac_lab <- tbl(con, I("mimiciv_hosp.labevents")) %>%
  filter(subject_id %in% !!cohort$subject_id) %>%
  filter(itemid %in% item_lac_lab) %>%
  select(subject_id, charttime, valuenum) %>%
  collect() %>%
  rename(Lactate = valuenum) %>%
  inner_join(cohort %>% select(stay_id, subject_id), by = "subject_id") %>%
  select(stay_id, charttime, Lactate)

# 4.1.2 从 chartevents 提取 (补充 ICU 内的床旁血气/护理记录)
item_lac_chart <- c(225668) # Lactic Acid in chartevents
lac_chart <- tbl(con, I("mimiciv_icu.chartevents")) %>%
  filter(stay_id %in% !!cohort$stay_id) %>%
  filter(itemid %in% item_lac_chart) %>%
  select(stay_id, charttime, valuenum) %>%
  collect() %>%
  rename(Lactate = valuenum) %>%
  select(stay_id, charttime, Lactate)

# 4.1.3 合并两部分数据并去重质控
lac_long <- bind_rows(lac_lab, lac_chart) %>%
  filter(!is.na(Lactate)) %>%
  filter(Lactate >= 0 & Lactate <= 30) %>%
  group_by(stay_id, charttime) %>%
  summarise(Lactate = mean(Lactate, na.rm = TRUE), .groups = "drop")

# 4.2 Fluid Input
fluid_long <- tbl(con, I("mimiciv_icu.inputevents")) %>%
  filter(stay_id %in% !!cohort$stay_id) %>%
  filter(amountuom %in% c("ml", "mL", "L")) %>%
  select(stay_id, starttime, amount, amountuom) %>%
  collect() %>%
  mutate(amount_ml = ifelse(amountuom == "L", amount * 1000, amount)) %>%
  rename(charttime = starttime)

# 5. 映射到 Clinical Day (D0 - D28)
map_day <- function(df, time_col) {
  df %>%
    inner_join(icu_times, by = "stay_id") %>%
    mutate(
      time_val = !!sym(time_col),
      day_diff = as.numeric(difftime(time_val, d0_end, units = "days")),
      visit_day = case_when(
        time_val < intime ~ -1, 
        time_val >= intime & time_val < d0_end ~ 0,
        TRUE ~ floor(day_diff) + 1
      )
    ) %>%
    filter(visit_day >= 0 & visit_day <= 28)
}

lac_mapped <- map_day(lac_long, "charttime") %>%
  group_by(stay_id, visit_day) %>%
  summarise(lactate = max(Lactate, na.rm = TRUE), .groups = "drop") %>%
  mutate(lactate = ifelse(is.infinite(lactate), NA, lactate))

fluid_mapped <- map_day(fluid_long, "charttime") %>%
  group_by(stay_id, visit_day) %>%
  summarise(fluid_input_daily = sum(amount_ml, na.rm = TRUE), .groups = "drop")

# 6. 合并纵向数据 (扩展到 28 天)
los_info <- cohort %>% select(stay_id, los_icu)

full_grid <- expand_grid(
  stay_id = unique(cohort$stay_id),
  visit_day = 0:28
) %>%
  inner_join(los_info, by = "stay_id") %>%
  filter(visit_day <= ceiling(los_icu))

long_data <- full_grid %>%
  left_join(lac_mapped, by = c("stay_id", "visit_day")) %>%
  left_join(fluid_mapped, by = c("stay_id", "visit_day")) %>%
  mutate(
    fluid_input_daily = replace_na(fluid_input_daily, 0)
  )

# 7. 检查数据
stats <- long_data %>%
  group_by(stay_id) %>%
  summarise(
    n_lac = sum(!is.na(lactate)),
    n_fluid = sum(!is.na(fluid_input_daily))
  )

write_csv(long_data, "long_data_daily_lac.csv")

dbDisconnect(con)
