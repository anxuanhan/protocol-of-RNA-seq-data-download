rm(list=ls())
load("GSE14520表达矩阵和分组信息.RData")

#确认分组的名字和矩阵的名字是一一对应的
identical(rownames(rt2),colnames(exp1))

group_list <- factor(rt2$group,levels = c("normal","tumor")) #表示tumor/normal

library(limma)
design <- model.matrix(~group_list)

#比较矩阵命名
design

#线性模型拟合
fit <- lmFit(exp1,design)

#贝叶斯检验
fit2 <- eBayes(fit)
deg <- topTable(fit2,coef = 2, number = Inf)

DEG <- na.omit(deg)#logFC负值表示在肿瘤中低表达，正值代表在肿瘤中高表达

#进行注释
logFC_cutoff <- 1
type1 = (DEG$P.Value<0.05)&(DEG$logFC < -logFC_cutoff)
type2 = (DEG$P.Value<0.05)&(DEG$logFC >  logFC_cutoff)

DEG$change <- ifelse(type1,"DOWN",ifelse(type2,"UP","NOT"))
write.csv(DEG,"GSE14520_DEGs.csv")

#火山图可视化
library(ggplot2)
library(cowplot)
library(extrafont)

ggplot(data = DEG,
       aes(x = logFC,
           y = -log10(P.Value))) +
  geom_point(alpha = 1, size = 1,
             aes(color = change)) +
  ylab("-log10(P.Value)") +
  scale_color_manual(values = c("#0072B5CC", "grey", "red")) +
  geom_vline(xintercept = c(-1, 1), lty = 5, col = "black", lwd = 0.8) +
  geom_hline(yintercept = -log10(0.05), lty = 4, col = "black", lwd = 0.8) +
  theme_half_open() +
  # 设置坐标轴粗细，宽度为0.8
  theme(axis.line = element_line(size = 0.8)) +
  # 设置所有文本字体为Arial，大小为12
  theme(text = element_text(family = "Arial", size = 12))

ggsave("TCGA_ssGSEA_DEG火山图.pdf", width = 7, height = 7, dpi = 300)


