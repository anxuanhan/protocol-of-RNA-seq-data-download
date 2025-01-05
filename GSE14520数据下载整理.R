#GEO芯片数据下载处理通用代码
rm(list=ls())

#用GEOquery包下载芯片数据
library(GEOquery)
gset <- getGEO("GSE14520",destdir = ".", AnnotGPL = F, getGPL = F)
class(gset)
#获取gset中的基本信息
gset[[1]]

#### 1.提取gset中的临床信息，最重要的是分组信息

#提取gset中的临床信息
rt1 <- pData(gset[[1]])

#根据临床信息确定每一个样本属于normal还是tumor
library(tidyverse)
group_list <- ifelse(str_detect(rt1$characteristics_ch1,"Non"),"normal","tumor")

#单独创建个数据框rt2保存分组的信息
rt1$group <- group_list
rt2 <- rt1 %>% dplyr::select("group")

#### 2.提取gset中的表达矩阵，并进行标准化和归一化
exp <- exprs(gset[[1]])
boxplot(exp,outline=F,notch=T,las=2)
#数据的标准化
library(limma)
exp <- normalizeBetweenArrays(exp)
boxplot(exp,outline=F,notch=T,las=2)

range(exp)
#根据range(exp)选择是否进行归一化
#exp <- log2(exp+1)

identical(rownames(rt2),colnames(exp))

#### 3.id转换
gset[[1]]@annotation
BiocManager::install("hthgu133a.db")#注意根据GPL平台信息选择相对应的R包
library(hthgu133a.db)
ls("package:hthgu133a.db")
#ids是probe_id和gene_symbol相对应的一个表格
ids <- toTable(hthgu133aSYMBOL)
length(unique(ids$symbol))
table(ids$symbol)
table(sort(table(ids$symbol)))

library(tidyverse)
exp1 <- as.data.frame(exp)

#将exp1中的probe_id对应的symbol添加到exp1表格上
exp1 <- exp1 %>% 
  mutate(probe_id = rownames(exp1)) %>% 
  inner_join(ids,by = "probe_id") %>% 
  dplyr::select(probe_id,symbol,everything())

#对symbol进行去重，因为有多个probe_id对应同一个symbol的情况
exp1 <- exp1[!duplicated(exp1$symbol),]

#将列名转换成行名
rownames(exp1) <- exp1$symbol
exp1 <- exp1[,-(1:2)]

#保存表达矩阵和分组的信息
save(exp1,rt2,file = "GSE14520表达矩阵和分组信息.RData")
