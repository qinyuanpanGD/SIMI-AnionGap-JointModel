library(tidyverse)

# 1. 读取数据
if(file.exists("cohort_final_baseline.csv")) {
  cohort <- read_csv("cohort_final_baseline.csv")
} else {
  stop("找不到 cohort_final_baseline.csv，请先运行 step1_extract_baseline.r")
}

# 2. 定义需要检查的变量
cont_vars <- colnames(cohort)

# 3. 计算缺失率
missing_summary <- cohort %>%
  summarise(across(everything(), ~ sum(is.na(.)))) %>%
  pivot_longer(everything(), names_to = "Variable", values_to = "Missing_Count") %>%
  mutate(
    Total_N = nrow(cohort),
    Missing_Rate = Missing_Count / Total_N,
    Missing_Percent = sprintf("%.2f%%", Missing_Rate * 100)
  ) %>%
  arrange(desc(Missing_Rate))

# 4. 打印结果

# 5. 保存结果
write_csv(missing_summary, "step1_missing.csv")

# 6. 简单的可视化 (可选)
high_missing <- missing_summary %>% filter(Missing_Rate > 0.2)
if(nrow(high_missing) > 0) {
} else {
}
