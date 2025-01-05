rm(list=ls())   
options(stringsAsFactors = F)

#########----------处理样本表达矩阵--------------###
#ESTIMATE包
library(utils)
library(data.table)
library(tidyverse)
##读取表达文件
log2FPKM_ens <- fread("./GDCdata/TCGA-LIHC.star_tpm.tsv",header = T)
library(clusterProfiler)
library(org.Hs.eg.db) #org.Mm.eg.db
## 有版本号，直接转不行的, 这句代码用于去除版本号
log2FPKM_ens$Ensembl_ID <- gsub("\\..*", "", log2FPKM_ens$Ensembl_ID)
ens.id <- log2FPKM_ens$Ensembl_ID
gene.symbol <- bitr(geneID = ens.id, 
                    fromType = "ENSEMBL",
                    toType = "SYMBOL",
                    OrgDb = org.Hs.eg.db)
log2FPKM <- log2FPKM_ens %>% 
  dplyr::rename(ENSEMBL=Ensembl_ID) %>% 
  inner_join(gene.symbol,by='ENSEMBL') %>% 
  drop_na() %>% 
  dplyr::select(SYMBOL,2:ncol(log2FPKM_ens)) %>% 
  distinct(SYMBOL,.keep_all = T) %>% 
  column_to_rownames(var = 'SYMBOL')
# #pheno <- read.table('TCGA-CESC.GDC_phenotype.tsv',header = T,row.names = 1,sep="\t")
# library(data.table)
# #data.table包提供了一个data.frame的高级版本，让你的程序做数据整型的运算速度大大的增加。
# pheno <- as.data.frame(fread("D:/ky/0.TCGA/1.CESC/TCGA-CESC.GDC_phenotype.tsv"))
# colnames(pheno)

#TCGA可以根据第14和第15位判断是癌组织还是癌旁组织。01表示癌症组织，11表示正常组织
#question 把它变成因子型，因子型是有顺序的。这里需要变换成因子型
sample <- colnames(log2FPKM)
group_list=factor(ifelse(as.numeric(substr(sample,14,15)) < 10,
                         'tumor','normal'))
#group_list记录了表达矩阵中患者ID的对应的组织类型，和表达矩阵的顺序一样。
names(group_list) <- sample
table(group_list)
group_list <- as.data.frame(group_list)
group_list$sample <- rownames(group_list)
tumor <- rownames(group_list[group_list$group_list=='tumor',])
normal <- rownames(group_list[group_list$group_list=='normal',])
tumor_logfpkm <- log2FPKM[,tumor]
normal_logfpkm <- log2FPKM[,normal]
save(log2FPKM,file='all_log2TPM.Rdata')
save(tumor_logfpkm,file='hcc_log2TPM.Rdata')
save(normal_logfpkm,file = 'normal_log2TPM.Rdata')




#########-------------------- estimate计算免疫得分---------
library(estimate)
filterCommonGenes(input.f = "CESC_log2FPKM.txt",   #输入文件名
                  output.f = "CESC_log2FPKM.gct",   #输出文件名
                  id = "GeneSymbol")   #行名为gene symbol
estimateScore("CESC_log2FPKM.gct",   
              "CESC_estimate_score.txt",  
              platform="affymetrix")   #默认平台

est <- read.table("CESC_estimate_score.txt",
                  sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
est <- est[,-1]   #移除第一列
colnames(est) <- est[1,]   #设置列名
est <- as.data.frame(t(est[-1,]))
#rownames(est) <- colnames(expr)
write.table(est, file = "CESC_est.txt",
            sep = "\t",row.names = T,col.names = NA,quote = F) 

###########---------------------Fig2A-C生存曲线-------------
rm(list=ls())   #清空环境变量
options(stringsAsFactors = F)
library(tidyverse)
##ESTIMATE评分
est <- read.table("CESC_est.txt",
                  sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
est$sample <- rownames(est)
#提取生存信息
clinical <- read_tsv("D:/过渡/zj/TCGA分析/原始数据/TCGA-KIRC.survival.tsv") 
clinical <- clinical %>% 
  dplyr::select(-"_PATIENT")

#--ImmuneScore生存分析
colnames(rt1) <- c("sample","score","group")
surv <- rt1 %>% inner_join(clinical,by = "sample") 
surv <- surv[!duplicated(surv$sample),]
#分组
library(pROC)
roc <- roc(surv$OS,surv$score)
plot(roc,
     legacy.axes = TRUE,
     print.auc=TRUE,
     main="ROC曲线最佳阈值点",
     thresholds="best", # 基于youden指数选择roc曲线最佳阈值点
     print.thres="best") # 在roc曲线上显示最佳阈值点
threshold <- coords(roc, "best")$threshold
threshold
surv$group <- ifelse(surv$score > threshold,"High","Low")
surv$group <- factor(surv$group, levels = c("Low","High")) 

table(surv$group)
## OS曲线
library(survival)
library(survminer)
library(ggplot2)
surfit <- survfit(Surv(OS.time,OS)~group,data=surv)
ggsurvplot(surfit,pval = TRUE,pval.method=TRUE,conf.int = TRUE,
           risk.table = TRUE,
           risk.table.col = 'strata',
           xlab='Time in months',
           linetype = 'strata',
           surv.median.line = 'hv',
           risk.table.y.text = FALSE, #在风险表图例中的文本注释中显示条形而不是名称
           palette="lancet",
           #ggtheme = theme_bw()
)
#--stromalScore生存分析
#分组
roc.stro <- roc(surv$OS,surv$StromalScore)
threshold_stro <- coords(roc.stro, "best")$threshold
threshold_stro
surv$stro.group <- ifelse(surv$StromalScore > threshold_stro,"High","Low")
surv$stro.group <- factor(surv$stro.group, levels = c("Low","High")) 
table(surv$stro.group)
## OS曲线
surfit_stro <- survfit(Surv(OS.time,OS)~stro.group,data=surv)
ggsurvplot(surfit_stro,pval = TRUE,pval.method=TRUE,conf.int = TRUE,
           risk.table = TRUE,
           risk.table.col = 'strata',
           xlab='Time in months',
           linetype = 'strata',
           surv.median.line = 'hv',
           risk.table.y.text = FALSE, #在风险表图例中的文本注释中显示条形而不是名称
           palette="lancet",
           #ggtheme = theme_bw()
)
#--estimateScore生存分析
#分组
roc.est <- roc(surv$OS,surv$ESTIMATEScore)
threshold_est <- coords(roc.est, "best")$threshold
threshold_est
surv$est.group <- ifelse(surv$ESTIMATEScore > threshold_est,"High","Low")
table(surv$est.group)
## OS曲线
surfit_est <- survfit(Surv(OS.time,OS)~est.group,data=surv)
ggsurvplot(surfit_est,pval = TRUE,pval.method=TRUE,conf.int = TRUE,
           risk.table = TRUE,
           risk.table.col = 'strata',
           xlab='Time in months',
           linetype = 'strata',
           surv.median.line = 'hv',
           risk.table.y.text = FALSE, #在风险表图例中的文本注释中显示条形而不是名称
           palette="lancet",
           #ggtheme = theme_bw()
)






# fitd <- survdiff(Surv(OS.time, OS) ~ immu.group,
#                  data      = surv,
#                  na.action = na.exclude)
# pValue <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
# pValue #查看P值
# fit <- survfit(Surv(OS.time, OS)~ immu.group, data = surv)
# summary(fit)
# p.lab <- paste0("P", ifelse(pValue < 0.001, " < 0.001", 
#                             paste0(" = ",round(pValue, 3))))


#
trans_3d_2d <- function(data, theta=135, phi=60) {
  pmat <- plot3D::perspbox(z=diag(2), plot=F, theta=theta, phi=phi)
  data <- as.data.frame(data)
  XY <- plot3D::trans3D(
    x = data$x,
    y = data$y,
    z = data$z,
    pmat = pmat) %>%
    data.frame()
  data$x <- XY$x
  data$y <- XY$y
  return(data[, c('x', 'y')])
}
require(onion)
data(bunny)
d <- trans_3d_2d(bunny)
require(ggplot2)
ggplot(d, aes(x, y, color=x)) + geom_point(alpha=.2)  + 
  coord_flip() + 
  theme(legend.position="none") + 
  theme_bw()+ #背景变为白色
  xlab('')+ylab('')+
  scale_y_reverse()
