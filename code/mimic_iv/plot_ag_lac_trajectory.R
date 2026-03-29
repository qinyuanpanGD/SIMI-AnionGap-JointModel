# 1. 加载纵向AG和乳酸数据
# 2. 结合基线生存数据（28天死亡）
# 3. 绘制个体轨迹（细线）和拟合曲线（粗线）
# 4. 按生存状态（Survivors vs Non-survivors）分组展示

rm(list = ls()) 
gc() 

if (!require("pacman")) install.packages("pacman") 
pacman::p_load( 
  tidyverse,    # 数据处理和可视化 
  ggplot2,      # 绘图 
  dplyr,        # 数据操作 
  reshape2,     # 数据重塑 
  gridExtra,    # 图形排列 
  scales,       # 坐标轴格式化 
  here          # 路径管理 
) 

options(repr.plot.width = 16, repr.plot.height = 10) 

# 1. 数据加载 

base_path <- "./"

# 1. 加载基线数据 (包含生存状态)
cohort_file <- file.path(base_path, "cohort_imputed.csv")
cohort_data <- read.csv(cohort_file) %>%
  select(stay_id, death_28_day) %>%
  mutate(
    Group = factor(death_28_day, levels = c(0, 1), labels = c("Survivors", "Non-survivors"))
  )

# 2. 加载AG纵向数据
ag_file <- file.path(base_path, "long_data_daily.csv")
ag_data <- read.csv(ag_file) %>%
  select(stay_id, visit_day, aniongap) %>%
  rename(Value = aniongap) %>%
  mutate(Variable = "Anion Gap (mEq/L)")

# 3. 加载乳酸纵向数据
lac_file <- file.path(base_path, "script_lac/long_data_daily_lac.csv")
lac_data <- read.csv(lac_file) %>%
  select(stay_id, visit_day, lactate) %>%
  rename(Value = lactate) %>%
  mutate(Variable = "Lactate (mmol/L)")

# 2. 数据合并与预处理 

long_data <- bind_rows(ag_data, lac_data)

plot_data <- long_data %>%
  inner_join(cohort_data, by = "stay_id") %>%
  filter(!is.na(Value)) %>% # 移除缺失值
  filter(visit_day <= 28)   # 限制在28天内

set.seed(123)
sample_ids <- plot_data %>%
  select(stay_id, Group) %>%
  distinct() %>%
  group_by(Group) %>%
  slice_sample(n = 100) %>% # 每组抽100人，如果不足则全选
  pull(stay_id)

background_lines_data <- plot_data %>%
  filter(stay_id %in% sample_ids)

# 3. 绘图 

group_colors <- c("Survivors" = "#2E7D32", "Non-survivors" = "#D32F2F") # 绿 vs 红

p <- ggplot(data = plot_data, aes(x = visit_day, y = Value, color = Group, fill = Group)) +

  # 1. 绘制个体轨迹 (背景细线) - 只绘制抽样患者
  geom_line(data = background_lines_data, 
            aes(group = stay_id), 
            alpha = 0.15, linewidth = 0.3) +

  # 2. 绘制拟合平滑曲线 (LOESS 或 GAM) - 使用所有数据
  geom_smooth(method = "loess", se = TRUE, alpha = 0.2, linewidth = 1.2) +

  facet_wrap(~ Variable, scales = "free_y", ncol = 2) +

  scale_x_continuous(name = "Days after ICU Admission", breaks = seq(0, 28, 7)) +
  scale_y_continuous(name = "Concentration") +
  scale_color_manual(values = group_colors) +
  scale_fill_manual(values = group_colors) +

  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 12, color = "gray40"),
    strip.text = element_text(size = 14, face = "bold"),
    strip.background = element_rect(fill = "gray95", color = NA),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    panel.grid.minor = element_blank()
  ) +

  labs(
    title = "Longitudinal Trajectories of Anion Gap and Lactate",
    subtitle = "Comparison between 28-Day Survivors and Non-Survivors (Individual Lines + LOESS Fit)",
    caption = "Note: Thin lines represent a random sample of patients; Thick curves represent the population trend."
  )

# 4. 保存图片 

output_file <- file.path(base_path, "AG_Lactate_Trajectory_Plot.png")
ggsave(output_file, plot = p, width = 12, height = 8, dpi = 300, bg = "white")

