if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, forestploter, grid)

# 1. 读取数据
df <- read_csv("Table_S_Subgroup_Analysis.csv", show_col_types = FALSE)

# 2. 数据整理与提取
plot_df <- df %>%
  mutate(
    HR = as.numeric(str_extract(`HR (95% CI)`, "^[0-9.]+")),
    Lower = as.numeric(str_extract(`HR (95% CI)`, "(?<=\\()[0-9.]+")),
    Upper = as.numeric(str_extract(`HR (95% CI)`, "[0-9.]+(?=\\))")),
    ` ` = paste(rep(" ", 20), collapse = " ")
  )

structured_df <- data.frame(
  Subgroup = c(
    "Age", 
    "  <= 70", "  > 70",
    "Gender", 
    "  Male", "  Female",
    "SOFA score", 
    "  <= 3", "  > 3",
    "Vasopressor use", 
    "  Yes", "  No",
    "RRT use", 
    "  Yes", "  No",
    "Mechanical ventilation", 
    "  Yes", "  No",
    "SAPS II", 
    "  <= 44", "  > 44"
  ),
  stringsAsFactors = FALSE
)

structured_df$N <- NA
structured_df$Events <- NA
structured_df$`HR (95% CI)` <- NA
structured_df$`P Value` <- NA
structured_df$HR <- NA
structured_df$Lower <- NA
structured_df$Upper <- NA
structured_df$` ` <- paste(rep(" ", 20), collapse = " ")

map_data <- function(target_name, source_name) {
  idx_target <- which(structured_df$Subgroup == target_name)
  idx_source <- which(plot_df$Subgroup == source_name)
  if(length(idx_target) > 0 && length(idx_source) > 0) {
    structured_df$N[idx_target] <<- plot_df$N[idx_source]
    structured_df$Events[idx_target] <<- plot_df$Events[idx_source]
    structured_df$`HR (95% CI)`[idx_target] <<- plot_df$`HR (95% CI)`[idx_source]
    structured_df$`P Value`[idx_target] <<- plot_df$`P Value`[idx_source]
    structured_df$HR[idx_target] <<- plot_df$HR[idx_source]
    structured_df$Lower[idx_target] <<- plot_df$Lower[idx_source]
    structured_df$Upper[idx_target] <<- plot_df$Upper[idx_source]
  }
}

map_data("  <= 70", "Age <= 70")
map_data("  > 70", "Age > 70")
map_data("  Male", "Male")
map_data("  Female", "Female")
map_data("  <= 3", "SOFA <= 3")
map_data("  > 3", "SOFA > 3")
map_data("  Yes", "Vasopressor: Yes") # 注意Vasopressor下面有两个Yes,所以要按顺序或者修改映射逻辑。为了准确，我硬编码索引。

structured_df[2, 2:8] <- plot_df[1, 2:8] # Age <= 70
structured_df[3, 2:8] <- plot_df[2, 2:8] # Age > 70
structured_df[5, 2:8] <- plot_df[3, 2:8] # Male
structured_df[6, 2:8] <- plot_df[4, 2:8] # Female
structured_df[8, 2:8] <- plot_df[5, 2:8] # SOFA <= 3
structured_df[9, 2:8] <- plot_df[6, 2:8] # SOFA > 3
structured_df[11, 2:8] <- plot_df[7, 2:8] # Vasopressor: Yes
structured_df[12, 2:8] <- plot_df[8, 2:8] # Vasopressor: No
structured_df[14, 2:8] <- plot_df[9, 2:8] # RRT: Yes
structured_df[15, 2:8] <- plot_df[10, 2:8] # RRT: No
structured_df[17, 2:8] <- plot_df[11, 2:8] # MV: Yes
structured_df[18, 2:8] <- plot_df[12, 2:8] # MV: No
structured_df[20, 2:8] <- plot_df[13, 2:8] # SAPS II <= 44
structured_df[21, 2:8] <- plot_df[14, 2:8] # SAPS II > 44

structured_df[is.na(structured_df)] <- ""

# 3. 森林图主题设计 (参考用户提供的漂亮主题)
tm <- forest_theme( 
  base_size = 12, 
  ci_pch = 15, # 使用方形，更常见于学术论文
  ci_col = "#4575b4", 
  ci_fill = "#4575b4", 
  ci_lty = 1, 
  ci_lwd = 1.5, 
  refline_col = "red", 
  refline_lty = "dashed", 
  summary_fill = "#fc8d59", 
  xaxis_lwd = 1, 
  core = list( 
    bg_params = list(fill = c("white", "#f7f7f7")) 
  ) 
) 

data_to_plot <- structured_df %>%
  select(Subgroup, N, Events, ` `, `HR (95% CI)`, `P Value`)

structured_df$HR <- as.numeric(ifelse(structured_df$HR == "", NA, structured_df$HR))
structured_df$Lower <- as.numeric(ifelse(structured_df$Lower == "", NA, structured_df$Lower))
structured_df$Upper <- as.numeric(ifelse(structured_df$Upper == "", NA, structured_df$Upper))

# 4. 绘制森林图
p <- forest( 
  data = data_to_plot, 
  lower = structured_df$Lower, 
  upper = structured_df$Upper, 
  est = structured_df$HR, 
  ci_column = 4, # 绘图空列在第4列
  ref_line = 1, 
  theme = tm, 
  xlim = c(0.8, 1.25), # 根据实际数据范围调整 (最小0.919, 最大1.149)
  ticks_at = c(0.8, 1.0, 1.2) # 设置坐标轴刻度
) 

header_rows <- c(1, 4, 7, 10, 13, 16, 19)
for(i in header_rows) {
  p <- edit_plot(p, row = i, col = 1, gp = gpar(fontface = "bold"))
}

# 5. 保存图形
pdf("Figure_3_Subgroup_Forest_Plot.pdf", width = 10, height = 8)
dev.off()

