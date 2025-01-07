##SVM-REF
rm(list=ls())
genes <- read.csv("单因素cox分析结果.csv")
index <- genes$gene
index <- index[index != ""]
index <- c("futime","fustat",index)

load("TCGA-LIHC数据整理.Rdata")

dat1 <- surv[, index]

input <- dat1
input <- input[,-1]

write.csv(input,"SVM_input.csv")


library(tidyverse)
library(glmnet)
source('msvmRFE.R') 
BiocManager::install("sigFeature")
library(sigFeature)
library(e1071)
library(caret)
library(randomForest)

nfold=5#5倍交叉验证
nrows=nrow(input)
folds=rep(1:nfold,len=nrows)[sample(nrows)]
folds=lapply(1:nfold,function(x) which(folds == x))
results = lapply(folds, svmRFE.wrap,input,k=10,halve.above=100)#特征选择

top.features=WriteFeatures(results,input,save=FALSE)

#进行特征重要性扫描
featsweep=lapply(1:30,FeatSweep.wrap,results,input)

#计算最小信息率
no.info=min(prop.table(table(input[,1])))

#计算模型预测误差
errors=sapply(featsweep,function(x) ifelse(is.null(x),NA,x$error))

dev.new(width=4, height=4, bg='white')
pdf("B_svm-accuracy.pdf",width = 5,height = 5)
Plotaccuracy(1-errors,no.info=no.info) #查看准确率
dev.off()

PlotErrors(errors, no.info=no.info)














svmRFE(input, k = 5, halve.above = 100) #分割数据，分配随机数，

nfold = 5
nrows = nrow(input)
folds = rep(1:nfold, len=nrows)[sample(nrows)]
folds = lapply(1:nfold, function(x) which(folds == x))
results = lapply(folds, svmRFE.wrap, input, k=5, halve.above=100)


top.features = WriteFeatures(results, input, save=F) #查看主要变量
head(top.features)


#把SVM-REF找到的特征保存到文件
write.csv(top.features,"feature_svm.csv")

featsweep = lapply(1:25, FeatSweep.wrap, results, input)
save(featsweep,file = "featsweep.RData")
load("featsweep.RData")

no.info = min(prop.table(table(input[,1])))
errors = sapply(featsweep, function(x) ifelse(is.null(x), NA, x$error))

#dev.new(width=4, height=4, bg='white')
#pdf("B_svm-error.pdf",width = 5,height = 5)
PlotErrors(errors, no.info=no.info) #查看错误率
write.csv(errors,"SVMREFerrors.csv")

a <- data.frame(col1 = rep(NA, 12))
# 再给它添加num列并赋值
a$num <- c(1:12)
a$error <- errors[1:12]

#####画错误率和特征个数之间的折线图
p1 <- ggplot(data = a, aes(x = num, y = error)) +
  geom_line(show.legend = FALSE, col = "#5d7baf", lwd = 1.5) +
  geom_point(size = 2, color = ifelse(a$error == min(a$error), "#df584a", "black")) +  # 设置颜色条件
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank(), legend.position = c(0.2, 0.8),
        axis.text = element_text(color = "black")) +  # 将坐标轴上的文本颜色设置为黑色
  labs(x = "Number of Feature", y = "10× CV Error")

# 找到最低点的坐标
lowest_point <- a[which.min(a$error), ]

# 添加标签
p1 <- p1 + geom_text(aes(x = lowest_point$num, y = lowest_point$error, 
                         label = sprintf("%d %.3f", which.min(a$error), lowest_point$error)), 
                     vjust = -2, color = "#df584a")

p1

#####画准确率和特征个数之间的折线图
b <- a
b$accuracy <- 1-b$error

ggplot(data = b, aes(x = num, y = accuracy)) +
  geom_line(show.legend = FALSE, col = "#5d7baf", lwd = 1.5) +
  geom_point(size = 2, color = ifelse(b$accuracy == max(b$accuracy), "#df584a", "black")) +  # 设置颜色条件
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank(), legend.position = c(0.2, 0.8),
        axis.text = element_text(color = "black")) +  # 将坐标轴上的文本颜色设置为黑色
  labs(x = "Number of Feature", y = "5× CV Accuracy")+
  geom_text(aes(x = highest_point$num, y = highest_point$accuracy, 
                label = sprintf("%d, %.3f", which.max(b$accuracy), highest_point$accuracy)), 
            vjust = 1.5,hjust=0.2, color = "#df584a")+
  theme(
    panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
    text = element_text(family = "Arial", size = 12))


# 找到最高点的坐标
highest_point <- b[which.max(b$accuracy), ]

# 添加标签


p2



which.min(errors)
top<-top.features[1:which.min(errors), "FeatureName"]
write.csv(top,"SVMtop.csv")
