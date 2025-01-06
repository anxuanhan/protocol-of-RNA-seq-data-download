#单因素cox回归

rm(list=ls())
library(survminer)
library(tidyverse)

load("TCGA-LIHC数据整理.RData")

genes <- read.csv("jVenn.csv")
colnames(genes) <- c("TCGA","GEO","comm")

index <- genes$comm
index <- index[index != ""]
index <- c("futime","fustat",index)

#bulk和sc-RNA的共有基因以及其os和os.time
dat1 <- surv[, index]

library(survival)
Coxoutput <- NULL # 初始化结果
for (i in 3:ncol(dat1)) {
  g <- colnames(dat1)[i]
  formula_str <- paste("Surv(futime, fustat) ~", colnames(dat1)[i])
  formula_obj <- as.formula(formula_str)
  cox <- coxph(formula_obj, data = dat1) # 单变量cox模型
  coxSummary <- summary(cox)
  
  # 取出模型的相应系数
  Coxoutput <- rbind.data.frame(Coxoutput,
                                data.frame(gene = g,
                                           HR = as.numeric(coxSummary$coefficients[,"exp(coef)"])[1],
                                           z = as.numeric(coxSummary$coefficients[,"z"])[1],
                                           pvalue = as.numeric(coxSummary$coefficients[,"Pr(>|z|)"])[1],
                                           lower = as.numeric(coxSummary$conf.int[,3][1]),
                                           upper = as.numeric(coxSummary$conf.int[,4][1]),
                                           stringsAsFactors = FALSE),
                                stringsAsFactors = FALSE)
}

Coxoutput1 <- Coxoutput[which(Coxoutput$pvalue < 0.05),]
gg <- Coxoutput1$gene

##单因素cox结果
write.csv(Coxoutput1,"单因素cox分析结果.csv")

