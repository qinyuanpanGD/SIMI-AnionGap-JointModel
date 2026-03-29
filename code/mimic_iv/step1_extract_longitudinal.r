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

# 2. 读取基线队列
if(file.exists("cohort_final_baseline.csv")) {
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
    intime_hour = hour(intime),
    d0_end = if_else(intime_hour < 7,
                     as.POSIXct(paste0(as.Date(intime), " 07:00:00"), tz = "UTC"),
                     as.POSIXct(paste0(as.Date(intime) + 1, " 07:00:00"), tz = "UTC"))
  )

# 4. 提取原始数据 (Raw Data with Timestamps)

# 4.1 AG Components
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

ag_long <- lab_raw %>%
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
  pivot_wider(names_from = label, values_from = valuenum, values_fn = mean) %>%

  # 1. 源头质控
  filter(
    (Na >= 110 & Na <= 180) &      # Na: 110-180
    (K >= 1 & K <= 10) &           # K: 1-10
    (Cl >= 80 & Cl <= 160) &       # Cl: 80-160
    (HCO3 >= 5 & HCO3 <= 50)       # HCO3: 5-50
  ) %>%

  filter(!is.na(Na) & !is.na(K) & !is.na(Cl) & !is.na(HCO3)) %>%

  # 2. 计算 AG (含钾公式)
  mutate(AG_strict = Na + K - Cl - HCO3) %>%

  # 3. 结果质控 (AG < 0 视为异常/缺失)
  filter(AG_strict >= 0 & AG_strict < 100) %>%

  select(hadm_id, charttime, AG_strict) %>%
  inner_join(cohort %>% select(stay_id, hadm_id), by = "hadm_id")

# 4.2 Troponin
item_trop <- c(51003)
trop_long <- tbl(con, I("mimiciv_hosp.labevents")) %>%
  filter(hadm_id %in% !!cohort$hadm_id) %>%
  filter(itemid %in% item_trop) %>%
  select(hadm_id, charttime, valuenum) %>%
  collect() %>%
  rename(Troponin = valuenum) %>%
  inner_join(cohort %>% select(stay_id, hadm_id), by = "hadm_id")

# 4.3 Fluid Input
fluid_long <- tbl(con, I("mimiciv_icu.inputevents")) %>%
  filter(stay_id %in% !!cohort$stay_id) %>%
  filter(amountuom %in% c("ml", "mL", "L")) %>%
  select(stay_id, starttime, amount, amountuom) %>%
  collect() %>%
  mutate(amount_ml = ifelse(amountuom == "L", amount * 1000, amount)) %>%
  rename(charttime = starttime) # Use starttime as reference

# 5. 映射到 Clinical Day (D0, D1, D2...)

map_day <- function(df, time_col) {
  df %>%
    inner_join(icu_times, by = "stay_id") %>%
    mutate(
      time_val = !!sym(time_col),
      day_diff = as.numeric(difftime(time_val, d0_end, units = "days")),
      visit_day = case_when(
        time_val < intime - hours(24) ~ -1, # Before admission (>24h)
        time_val >= intime - hours(24) & time_val < d0_end ~ 0, # Day 0 covers [intime-24h, d0_end]
        TRUE ~ floor(day_diff) + 1
      )
    ) %>%
    filter(visit_day >= 0 & visit_day <= 28) # Keep only 0-28 days
}

ag_mapped <- map_day(ag_long, "charttime") %>%
  group_by(stay_id, visit_day) %>%
  summarise(aniongap = max(AG_strict, na.rm = TRUE), .groups = "drop") %>%
  mutate(aniongap = ifelse(is.infinite(aniongap), NA, aniongap))

trop_mapped <- map_day(trop_long, "charttime") %>%
  group_by(stay_id, visit_day) %>%
  summarise(troponin_t = max(Troponin, na.rm = TRUE), .groups = "drop")

fluid_mapped <- map_day(fluid_long, "charttime") %>%
  group_by(stay_id, visit_day) %>%
  summarise(fluid_input_daily = sum(amount_ml, na.rm = TRUE), .groups = "drop")

# 6. 合并纵向数据

los_info <- cohort %>% select(stay_id, los_icu)

full_grid <- expand_grid(
  stay_id = unique(cohort$stay_id),
  visit_day = 0:28
) %>%
  inner_join(los_info, by = "stay_id") %>%
  filter(visit_day <= ceiling(los_icu)) # 限制在 ICU 住院期间 (向上取整)

long_data <- full_grid %>%
  left_join(ag_mapped, by = c("stay_id", "visit_day")) %>%
  left_join(trop_mapped, by = c("stay_id", "visit_day")) %>%
  left_join(fluid_mapped, by = c("stay_id", "visit_day")) %>%
  mutate(
    fluid_input_daily = replace_na(fluid_input_daily, 0)
  )

# 7. 检查数据丰度
stats <- long_data %>%
  group_by(stay_id) %>%
  summarise(
    n_ag = sum(!is.na(aniongap)),
    n_trop = sum(!is.na(troponin_t)),
    n_fluid = sum(!is.na(fluid_input_daily))
  )

write_csv(long_data, "long_data_daily.csv")

dbDisconnect(con)
