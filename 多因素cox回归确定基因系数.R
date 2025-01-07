gg <- read.csv("SVM_Lasso.csv")
gg <- gg$SVM.RFE.Lasso
gg <- gg[gg != ""]

##多因素cox回归
dat2 <- surv2 %>% dplyr::select(1,2,gg)
multiCox <- coxph(Surv(futime, fustat) ~ ., data = dat2)
multiCoxSum=summary(multiCox)

ggforest(multiCox,
         main = "Hazard ratio",
         cpositions = c(0.02,0.22, 0.4), 
         fontsize = 0.7, 
         refLabel = "reference", 
         noDigits = 2)

#输出模型相关信息
outMultiTab=data.frame()
outMultiTab=cbind(
  coef=round(multiCoxSum$coefficients[,"coef"],3),
  HR=round(multiCoxSum$conf.int[,"exp(coef)"],3),
  lower=round(multiCoxSum$conf.int[,"lower .95"],3),
  upper=round(multiCoxSum$conf.int[,"upper .95"],3),
  pvalue=round(multiCoxSum$coefficients[,"Pr(>|z|)"],4))
outMultiTab <- cbind(gene=rownames(outMultiTab), as.data.frame(outMultiTab))

Coxoutput2 <- outMultiTab[which(outMultiTab$pvalue < 0.05),] #p值改一下
write.csv(outMultiTab,"多因素cox分析结果.csv")

ltt <- cbind(c("Gene",Coxoutput2$gene),
             c("HR",Coxoutput2$HR),
             c("L CI",Coxoutput2$lower),
             c("U CI",Coxoutput2$upper),
             c("pvalue",Coxoutput2$pvalue))
##3.2 绘制森林图


library(forestplot)

forestplot(labeltext=ltt,
           mean=c(NA,as.numeric(Coxoutput2$HR)),
           lower=c(NA,as.numeric(Coxoutput2$lower)), 
           upper=c(NA,as.numeric(Coxoutput2$upper)),
           graph.pos=5,# 图在表中的列位置
           fn.ci_norm="fpDrawDiamondCI",# box类型选择钻石
           col=fpColors(box="red", lines="red"),# box颜色
           boxsize=0.4,# box大小固定
           zero=1,
           lwd.xaxis=1,
           xlab="Hazard ratios",
           txt_gp=fpTxtGp(label=gpar(cex=1),# 字体大小设置
                          ticks=gpar(cex=1),
                          xlab=gpar(cex=1),
                          title=gpar(cex=1)),
           hrzl_lines=list("1" = gpar(lwd=1.5, col="black"), 
                           "2" = gpar(lwd=1.5, col="black"), 
                           "5" = gpar(lwd=1.5, col="black")))



load("os and time 和基因表达矩阵.Rdata")
rt <- surv2
rt$sample <- rownames(rt)

gene <- read.csv("多因素cox分析结果.csv")
mulcoxgene <- gene$gene

actCoef <- gene$coef
#建模
FinalGeneExp = rt[,mulcoxgene]
myFun = function(x){crossprod(as.numeric(x),actCoef)}
riskScore = apply(FinalGeneExp,1,myFun)
outCol = c(mulcoxgene)
risk = as.vector(ifelse(riskScore > median(riskScore), "high", "low"))
trainRiskOut = cbind(rt[,c("sample","futime","fustat",outCol)], riskScore=as.vector(riskScore), risk)
rownames(trainRiskOut) <- trainRiskOut$sample
df <- trainRiskOut
df <- df[,-1]
df <- df[order(df$riskScore), ]




#绘制风险三联图
## ggrisk方法
library(ggrisk)
library(survival)
library(rms)
fit <- coxph(Surv(futime,fustat)~MPEG1+PLAUR+MIR155HG,df1)
ggforest(fit,
         main = "Hazard ratio",
         cpositions = c(0.02,0.22, 0.4), 
         fontsize = 0.7, 
         refLabel = "reference", 
         noDigits = 2)


fitSum=summary(fit)
fitTab=data.frame()
fitTab=cbind(
  coef=round(fitSum$coefficients[,"coef"],3),
  HR=round(fitSum$conf.int[,"exp(coef)"],3),
  lower=round(fitSum$conf.int[,"lower .95"],3),
  upper=round(fitSum$conf.int[,"upper .95"],3),
  pvalue=round(fitSum$coefficients[,"Pr(>|z|)"],4))
fitTab <- cbind(gene=rownames(fitTab), as.data.frame(fitTab))
write.csv(fitTab,"最终三个基因的多因素回归.csv")


actCoef <- fitTab$coef
FinalGeneExp = rt[,fitTab$gene]
myFun = function(x){crossprod(as.numeric(x),actCoef)}
riskScore = apply(FinalGeneExp,1,myFun)
outCol = c(fitTab$gene)
risk = as.vector(ifelse(riskScore > median(riskScore), "high", "low"))
trainRiskOut = cbind(rt[,c("sample","futime","fustat",outCol)], riskScore=as.vector(riskScore), risk)
rownames(trainRiskOut) <- trainRiskOut$sample
df <- trainRiskOut
df <- df[,-1]
df <- df[order(df$riskScore), ]
df3 <- df[,c(3,4,5,7)]

dd <- datadist(df3)
options(datadist="dd")

fit <- cph(Surv(futime,fustat)~MPEG1+PLAUR+MIR155HG,df1)
coef <- coef(fit)
riskScore_3 <- predict(fit,df1,type="risk")
df1$riskscore <- riskScore_3
df1$group <- ifelse(df1$riskscore > median(df1$riskscore), "high", "low")

dd <- datadist(df3)
options(datadist="dd")

ggrisk(fit,
       cutoff.value='median', #可选‘
       cutoff.x = 350,  #“cutoff”文本的水平位置
       cutoff.y = -1,  #“cutoff”文本的垂直位置
       code.highrisk = 'High PR Risk',#高风险标签，默认为 ’High’
       code.lowrisk = 'Low PR Risk', #低风险标签，默认为 ’Low’
       title.A.ylab='PRscore', #A图 y轴名称
       title.B.ylab='Survival Time(year)', #B图 y轴名称，注意区分year month day
       title.A.legend='Risk Group', #A图图例名称
       title.B.legend='Status',     #B图图例名称
       title.C.legend='Expression', #C图图例名称
       relative_heights=c(0.1,0.1,0.01,0.15), #A、B、热图注释和热图C的相对高度    
       color.A=c(low='#f1a540',high='#473f6f'),#A图中点的颜色
       color.B=c(code.0='#473f6f',code.1='#e53a32'), #B图中点的颜色
       color.C=c(low='#f1a540',median='white',high='#473f6f'), #C图中热图颜色
       vjust.A.ylab=1, #A图中y轴标签到y坐标轴的距离,默认是1
       vjust.B.ylab=2  #B图中y轴标签到y坐标轴的距离,默认是2
)

##手动方法
df2 <- df
df2$id <- 1:length(rownames(df))
df2$risk <- factor(df2$risk,levels = c('low','high')) #指定顺序
write.csv(df2,"三个基因加权求和计算pr riskscore.csv")

p1 <- ggplot(df2,aes(x = id,y = riskScore)) +
  geom_point(aes(col = risk)) +
  scale_colour_manual(values = c('#f1a540','#473f6f')) +
  geom_vline(xintercept = sum(df2$risk == "low"), colour="grey", linetype = "dashed", size = 0.8) +
  theme_bw()+
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        #panel.border = element_blank(),
        axis.title.x = element_blank())

p1

df2$status <- ifelse(df2$fustat==0,"Alive","Dead")
df2$status <- factor(df2$status)

p2 <- ggplot(df2,aes(x = id,y = futime)) +
  geom_point(aes(col = status)) +
  scale_colour_manual(values = c('#473f6f','#e53a32')) +
  theme_bw()+
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        #panel.border = element_blank(),
        axis.title.x = element_blank())
p2

plot_grid(p1,p2, nrow = 2,align = "v", axis = "tlbr")

library(ggplot2)
library(pheatmap)
library(ggplotify)
library(cowplot)

mycol <- colorRampPalette(c('#fbbf7a','#fffefe','#4f4775'))(100) # 自定义颜色

df3 <- df[, c(3, 4, 5)]
df3 <- t(df3)
annotation <- df2 %>% dplyr::select("risk")

identical(colnames(df3), rownames(annotation))

ann_colors <- list(risk = c('low' = '#f1a540',
                            'high' = '#473f6f'))


pheatmap(df3,
         col = mycol,
         cluster_rows = F,
         cluster_cols = F,
         show_colnames = F,
         annotation_col = annotation,
         annotation_colors = ann_colors,
         annotation_legend = F
)


df4 <- df[,c(3,4,5,7)]
write.csv(df4,"深度学习数据整理.csv")
