library(here)
library(limma)
res <- readRDS(here("Results/TCGA paired data results May 2020", "results.TCGA.paired_brca.rds"))
names(res)
lfc <- readRDS(here("Results/TCGA paired data results May 2020", "lfc.TCGA.paired_brca.rds"))
names(lfc)
plot(lfc$lfc.mean.expHM, lfc$lfc.mean.lnHM, pch=20)
cor(lfc$lfc.mean.expHM, lfc$lfc.mean.lnHM) # 0.9977
plot(lfc$lfc.disp.expHM, lfc$lfc.disp.lnHM, pch=20)
cor(lfc$lfc.disp.expHM, lfc$lfc.disp.lnHM) # 0.9966
plot(res$p.mean.expHM.log, res$p.mean.lnHM.log, pch=20)
cor(res$p.mean.expHM.log, res$p.mean.lnHM.log) # 0.9812
cor(res$p.mean.expHM.log, res$p.mean.lnHM.log, method='spearman') # 0.9676
plot(res$p.disp.expHM.log, res$p.disp.lnHM.log, pch=20)
cor(res$p.disp.expHM.log, res$p.disp.lnHM.log) # 0.9807
cor(res$p.disp.expHM.log, res$p.disp.lnHM.log, method='spearman') # 0.9546


plot(log(rowMeans(res$counts)), lfc$lfc.mean.expHM, pch=20)
points(log(rowMeans(res$counts)), lfc$lfc.mean.lnHM, col='blue', pch=20)
