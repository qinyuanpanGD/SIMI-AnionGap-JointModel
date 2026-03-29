library(data.table)

mimic <- fread('/Users/qinyuanpan/Desktop/时间序列模拟2026-03-19算力平台备份最终版/study_ag_副本/script/cohort_imputed.csv')

if("Potassium" %in% colnames(mimic)) {
} else {
}
