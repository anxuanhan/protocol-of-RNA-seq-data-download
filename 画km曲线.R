##TCGA和GEO数据库根据焦亡基因评分高低画生存曲线

a <- read.csv("TCGA2个基因临床和pr riskscore.csv",row.names = 1)

## OS曲线
library(survival)
library(survminer)
library(ggplot2)
surfit <- survfit(Surv(futime,fustat)~risk,data=a)

ggsurvplot(surfit, pval = TRUE, pval.method = TRUE, conf.int = TRUE,
           risk.table = TRUE, risk.table.col = 'strata',
           xlab = 'Time in years', linetype = 'strata',
           surv.median.line = 'hv', risk.table.y.text = FALSE,
           palette = c('#473f6f','#f1a540'))




a <- read.csv("GSE296092个基因临床和pr riskscore.csv")


surfit <- survfit(Surv(futime,fustat)~risk,data=a)

ggsurvplot(surfit, pval = TRUE, pval.method = TRUE, conf.int = TRUE,
           risk.table = TRUE, risk.table.col = 'strata',
           xlab = 'Time in years', linetype = 'strata',
           surv.median.line = 'hv', risk.table.y.text = FALSE,
           palette = c('#473f6f','#f1a540'))

