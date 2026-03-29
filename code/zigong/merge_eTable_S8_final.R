library(dplyr)
library(readr)
library(stringr)

# 1. Read the two CSV files
ag_data <- read_csv("Table_2_JMbayes_AG_Zigong.csv", show_col_types = FALSE)
lac_data <- read_csv("Table_2_JMbayes_Lac_Zigong.csv", show_col_types = FALSE)

# 2. Rename columns to indicate which model they belong to
colnames(ag_data) <- c("Variable", "AG_Model_HR (95% CI)", "AG_Model_P_Value")
colnames(lac_data) <- c("Variable", "Lac_Model_HR (95% CI)", "Lac_Model_P_Value")

# 3. Clean up variable names for better presentation
clean_var_names <- function(df) {
  df %>%
    mutate(Variable = case_when(
      Variable == "admission_age" ~ "Age (years)",
      Variable == "genderMale" ~ "Male gender",
      Variable == "SOFA_score" ~ "SOFA score",
      Variable == "charlson_comorbidity_index" ~ "Charlson comorbidity index",
      Variable == "fluid_input_24h" ~ "Initial fluid input (L/24h)",
      Variable == "WBC" ~ "WBC (x10^9/L)",
      Variable == "Vasopressor_use" ~ "Vasopressor use",
      Variable == "MAP" ~ "Mean arterial pressure (mmHg)",
      Variable == "Mechanical_ventilation_use" ~ "Mechanical ventilation",
      Variable == "RRT_use" ~ "Renal replacement therapy",
      Variable == "Creatinine" ~ "Creatinine (mg/dL)",
      Variable == "Sedative_use" ~ "Sedative use",
      Variable == "value(aniongap)" ~ "Longitudinal Anion Gap",
      Variable == "value(lactate)" ~ "Longitudinal Lactate",
      Variable == "value(fluid_input_daily)" ~ "Longitudinal Daily Fluid Input",
      TRUE ~ Variable
    ))
}

ag_data <- clean_var_names(ag_data)
lac_data <- clean_var_names(lac_data)

# 4. Merge the two datasets side by side
merged_data <- full_join(ag_data, lac_data, by = "Variable")

# 5. Reorder rows so that longitudinal variables are at the top/bottom logically
var_order <- c(
  "Longitudinal Anion Gap",
  "Longitudinal Lactate",
  "Longitudinal Daily Fluid Input",
  "Age (years)",
  "Male gender",
  "SOFA score",
  "Charlson comorbidity index",
  "Initial fluid input (L/24h)",
  "WBC (x10^9/L)",
  "Vasopressor use",
  "Mean arterial pressure (mmHg)",
  "Mechanical ventilation",
  "Renal replacement therapy",
  "Creatinine (mg/dL)",
  "Sedative use"
)

merged_data <- merged_data %>%
  arrange(match(Variable, var_order))

# 6. Replace NA with "-"
merged_data[is.na(merged_data)] <- "-"

# 7. Save the merged eTable S8
write_csv(merged_data, "/Users/qinyuanpan/Desktop/时间序列模拟2026-03-19算力平台备份最终版/study_ag_副本/script_zigong/eTable_S8_Joint_Model_Validation_Final.csv")

