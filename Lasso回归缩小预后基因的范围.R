##Lasso回归缩小预后基因的范围

rm(list = ls())

genes <- read.csv("单因素cox分析结果.csv")
index <- genes$gene
index <- index[index != ""]
index <- c("futime","fustate",index)

load(".Rdata")

dat1 <- surv2[, index]
colnames(dat1)[1] <- "futime"
colnames(dat1)[2] <- "fustat"
rt <- dat1

library(glmnet)
library(survival)

str(rt)
###3. 构建模型
set.seed(23)   #设定随机种子
x=as.matrix(rt[,c(3:ncol(rt))])
y=data.matrix(Surv(rt$futime,rt$fustat))
fit=glmnet(x, y, family = "cox", maxit = 1000)

plot(fit, xvar = "lambda", label = TRUE)
cvfit = cv.glmnet(x, y, family="cox", maxit = 1000)

plot(cvfit)
abline(v = log(c(cvfit$lambda.min,cvfit$lambda.1se)),lty="dashed")
coef=coef(fit, s = cvfit$lambda.min)

index=which(coef != 0)
actCoef=coef[index]
lassoGene=row.names(coef)[index]
geneCoef=cbind(Gene=lassoGene,Coef=actCoef)
geneCoef   #查看模型的相关系数

write.csv(geneCoef,"lasso回归的7个基因.csv")


#建模
FinalGeneExp = rt[,lassoGene]
myFun = function(x){crossprod(as.numeric(x),actCoef)}
riskScore = apply(FinalGeneExp,1,myFun)
outCol = c(lassoGene)
risk = as.vector(ifelse(riskScore > median(riskScore), "high", "low"))
trainRiskOut = cbind(rt[,c("futime","fustat",outCol)], riskScore=as.vector(riskScore), risk)


df <- trainRiskOut
df <- df[order(df$riskScore), ]







# 加载pROC包
library(pROC)      
library(ggplot2)  

gfit <- roc(trainRiskOut$fustat,trainRiskOut$riskScore)

plot(gfit,
     print.auc=TRUE,
     # 图像上输出AUC值,坐标为（x，y）
     auc.polygon=TRUE, 
     auc.polygon.col="white", # 设置ROC曲线下填充色
     smooth=F, # 绘制不平滑曲线
     col="#d86967",  # 曲线颜色
     legacy.axes=TRUE)   # 使横轴从0到1，表示为1-特异度

# 加载timeROC
library(timeROC)
trainRiskOut$futime <- trainRiskOut$futime/365

ROC3 <- timeROC(T=trainRiskOut$futime,   #结局时间
                delta=trainRiskOut$fustat,   #结局指标
                marker=trainRiskOut$riskScore,   #预测变量
                cause=1,   #阳性结局指标数值
                weighting="marginal",   #计算方法，默认为marginal
                times=c(1,3,5),   #时间点，选取3年，4年和5年的生存率
                iid=TRUE)
plot(ROC3,
     time=1,col="#ab4332",lwd=2)   #time是时间点，col是线条颜色
plot(ROC3,
     time=3,col="#426b97", add=TRUE,lwd=2)   #add指是否添加在上一张图中
plot(ROC3,
     time=5,col="#cc7f44", add=TRUE,lwd=2)
legend("bottomright",
       c("Year1", "Year3", "Year5"),
       col=c("#ab4332", "#426b97", "#cc7f44"),
       lwd=2)


#生存曲线
library(survival)
library(survminer)
diff=survdiff(Surv(futime, fustat) ~risk,data = trainRiskOut)
fit <- survfit(Surv(futime, fustat)~ risk, data =trainRiskOut )
ggsurvplot(fit,
           data = trainRiskOut,
           palette =c("red","blue"),
           pval=T)

ggsurvplot(fit,
           data = trainRiskOut,
           conf.int = TRUE, 
           palette =c("red","blue"),
           risk.table = T,
           pval=T)




