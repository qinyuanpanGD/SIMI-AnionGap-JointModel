library(JMbayes2)

# 1. 检查AG模型
load("/Users/qinyuanpan/Desktop/时间序列模拟2026-03-19算力平台备份最终版/study_ag_副本/script_zigong/jmf_ag_zigong.RData")

pdf("/Users/qinyuanpan/Desktop/时间序列模拟2026-03-19算力平台备份最终版/study_ag_副本/script_zigong/Figure_S8_Diagnostics_AG_Zigong.pdf", width = 10, height = 8)
par(mfrow = c(2, 2))
plot(jmf_ag, which = "trace", parm = "alphas")
plot(jmf_ag, which = "density", parm = "alphas")
dev.off()

# 2. 检查Lac模型
load("/Users/qinyuanpan/Desktop/时间序列模拟2026-03-19算力平台备份最终版/study_ag_副本/script_zigong/jmf_lac_zigong.RData")

pdf("/Users/qinyuanpan/Desktop/时间序列模拟2026-03-19算力平台备份最终版/study_ag_副本/script_zigong/Figure_S8_Diagnostics_Lac_Zigong.pdf", width = 10, height = 8)
par(mfrow = c(2, 2))
plot(jmf_lac, which = "trace", parm = "alphas")
plot(jmf_lac, which = "density", parm = "alphas")
dev.off()

