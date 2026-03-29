rm(list = ls())
library(data.table)
library(dplyr)

out_path <- "./"


cohort <- fread(paste0(out_path, "zigong_simi_cohort.csv"))


baseline_all <- fread(paste0(out_path, "zigong_baseline_features.csv"))


baseline_628 <- baseline_all[INP_NO %in% cohort$INP_NO]


fwrite(baseline_628, paste0(out_path, "zigong_baseline_features_628.csv"))

fwrite(baseline_628, paste0(out_path, "zigong_baseline_features.csv"))
