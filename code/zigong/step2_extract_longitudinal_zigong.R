rm(list = ls())
library(data.table)
library(dplyr)

out_path <- "./"


cohort <- fread(paste0(out_path, "zigong_simi_cohort.csv"))


long_all <- fread(paste0(out_path, "zigong_longitudinal_features.csv"))


long_628 <- long_all[INP_NO %in% cohort$INP_NO]


fwrite(long_628, paste0(out_path, "zigong_longitudinal_features_628.csv"))

fwrite(long_628, paste0(out_path, "zigong_longitudinal_features.csv"))
