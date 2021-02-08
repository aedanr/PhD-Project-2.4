library(here)
library(dplyr)

## brca
voom.brca <- read.delim(
   here(
     "Results/GSEA results June 2020/jun23/20200623abslfcnlp.voom.brca.GseaPreranked.1592884198910", 
     "gsea_report_for_na_pos_1592884198910.xls"
   )
)
lnHMdisp.brca <- read.delim(
  here(
    "Results/GSEA results June 2020/jun23/20200623abslfcnlp.lnHMdisp.brca.GseaPreranked.1592883630665", 
    "gsea_report_for_na_pos_1592883630665.xls"
  )
)

voom.brca <- voom.brca[, c(1,4:8,10)]
lnHMdisp.brca <- lnHMdisp.brca[, c(1,4:8,10)]
dim(voom.brca) # 3286 7
dim(lnHMdisp.brca) # 1734 7

voom.brcafdr.05 <- voom.brca %>% 
   filter(FDR.q.val < 0.05)
lnHMdisp.brcafdr.05 <- lnHMdisp.brca %>% 
   filter(FDR.q.val < 0.05)
dim(voom.brcafdr.05) # 767 7
dim(lnHMdisp.brcafdr.05) # 42 7
sum(voom.brcafdr.05$NAME %in% lnHMdisp.brcafdr.05$NAME) # 4
# 763 GO terms identified by mean only, 38 by dispersion only with FDR < 0.05

voom.brcafdr.01 <- voom.brca %>% 
   filter(FDR.q.val < 0.01)
lnHMdisp.brcafdr.01 <- lnHMdisp.brca %>% 
   filter(FDR.q.val < 0.01)
dim(voom.brcafdr.01) # 410 7
dim(lnHMdisp.brcafdr.01) # 7 7
sum(voom.brcafdr.01$NAME %in% lnHMdisp.brcafdr.01$NAME) # 1
# 409 GO terms identified by mean only, 6 by dispersion only with FDR < 0.01


## kirc
voom.kirc <- read.delim(
   here(
      "Results/GSEA results June 2020/jun24/20200624abslfcnlp.voom.kirc.GseaPreranked.1592961675038", 
      "gsea_report_for_na_pos_1592961675038.xls"
   )
)
lnHMdisp.kirc <- read.delim(
   here(
      "Results/GSEA results June 2020/jun24/20200624abslfcnlp.lnHMdisp.kirc.GseaPreranked.1592961045354", 
      "gsea_report_for_na_pos_1592961045354.xls"
   )
)

voom.kirc <- voom.kirc[, c(1,4:8,10)]
lnHMdisp.kirc <- lnHMdisp.kirc[, c(1,4:8,10)]
dim(voom.kirc) # 3259 7
dim(lnHMdisp.kirc) # 1688 7

voom.kircfdr.05 <- voom.kirc %>% 
   filter(FDR.q.val < 0.05)
lnHMdisp.kircfdr.05 <- lnHMdisp.kirc %>% 
   filter(FDR.q.val < 0.05)
dim(voom.kircfdr.05) # 966 7
dim(lnHMdisp.kircfdr.05) # 127 7
sum(voom.kircfdr.05$NAME %in% lnHMdisp.kircfdr.05$NAME) # 8
# 958 GO terms identified by mean only, 119 by dispersion only with FDR < 0.05

voom.kircfdr.01 <- voom.kirc %>% 
   filter(FDR.q.val < 0.01)
lnHMdisp.kircfdr.01 <- lnHMdisp.kirc %>% 
   filter(FDR.q.val < 0.01)
dim(voom.kircfdr.01) # 586 7
dim(lnHMdisp.kircfdr.01) # 76 7
sum(voom.kircfdr.01$NAME %in% lnHMdisp.kircfdr.01$NAME) # 3
# 583 GO terms identified by mean only, 73 by dispersion only with FDR < 0.01


## lihc
voom.lihc <- read.delim(
   here(
      "Results/GSEA results June 2020/jun24/20200624abslfcnlp.voom.lihc.GseaPreranked.1592962868726", 
      "gsea_report_for_na_pos_1592962868726.xls"
   )
)
lnHMdisp.lihc <- read.delim(
   here(
      "Results/GSEA results June 2020/jun24/20200624abslfcnlp.lnHMdisp.lihc.GseaPreranked.1592962254943", 
      "gsea_report_for_na_pos_1592962254943.xls"
   )
)

voom.lihc <- voom.lihc[, c(1,4:8,10)]
lnHMdisp.lihc <- lnHMdisp.lihc[, c(1,4:8,10)]
dim(voom.lihc) # 3274 7
dim(lnHMdisp.lihc) # 2292 7

voom.lihcfdr.05 <- voom.lihc %>% 
   filter(FDR.q.val < 0.05)
lnHMdisp.lihcfdr.05 <- lnHMdisp.lihc %>% 
   filter(FDR.q.val < 0.05)
dim(voom.lihcfdr.05) # 642 7
dim(lnHMdisp.lihcfdr.05) # 281 7
sum(voom.lihcfdr.05$NAME %in% lnHMdisp.lihcfdr.05$NAME) # 153
# 489 GO terms identified by mean only, 128 by dispersion only with FDR < 0.05

voom.lihcfdr.01 <- voom.lihc %>% 
   filter(FDR.q.val < 0.01)
lnHMdisp.lihcfdr.01 <- lnHMdisp.lihc %>% 
   filter(FDR.q.val < 0.01)
dim(voom.lihcfdr.01) # 377 7
dim(lnHMdisp.lihcfdr.01) # 163 7
sum(voom.lihcfdr.01$NAME %in% lnHMdisp.lihcfdr.01$NAME) # 96
# 281 GO terms identified by mean only, 67 by dispersion only with FDR < 0.01

## thca
voom.thca <- read.delim(
   here(
      "Results/GSEA results June 2020/jun30/abslfcnlp.voom.thca.GseaPreranked.1593486735731", 
      "gsea_report_for_na_pos_1593486735731.xls"
   )
)
lnHMdisp.thca <- read.delim(
   here(
      "Results/GSEA results June 2020/jun30/abslfcnlp.lnHMdisp.thca.GseaPreranked.1593479849129", 
      "gsea_report_for_na_pos_1593479849129.xls"
   )
)

voom.thca <- voom.thca[, c(1,4:8,10)]
lnHMdisp.thca <- lnHMdisp.thca[, c(1,4:8,10)]
dim(voom.thca) # 3375 7
dim(lnHMdisp.thca) # 2836 7

voom.thcafdr.05 <- voom.thca %>% 
   filter(FDR.q.val < 0.05)
lnHMdisp.thcafdr.05 <- lnHMdisp.thca %>% 
   filter(FDR.q.val < 0.05)
dim(voom.thcafdr.05) # 796 7
dim(lnHMdisp.thcafdr.05) # 518 7
sum(voom.thcafdr.05$NAME %in% lnHMdisp.thcafdr.05$NAME) # 315
# 481 GO terms identified by mean only, 203 by dispersion only with FDR < 0.05

voom.thcafdr.01 <- voom.thca %>% 
   filter(FDR.q.val < 0.01)
lnHMdisp.thcafdr.01 <- lnHMdisp.thca %>% 
   filter(FDR.q.val < 0.01)
dim(voom.thcafdr.01) # 415 7
dim(lnHMdisp.thcafdr.01) # 247 7
sum(voom.thcafdr.01$NAME %in% lnHMdisp.thcafdr.01$NAME) # 147
# 268 GO terms identified by mean only, 100 by dispersion only with FDR < 0.01

## luad
voom.luad <- read.delim(
   here(
      "Results/GSEA results June 2020/jun30/abslfcnlp.voom.luad.GseaPreranked.1593482149447", 
      "gsea_report_for_na_pos_1593482149447.xls"
   )
)
lnHMdisp.luad <- read.delim(
   here(
      "Results/GSEA results June 2020/jun30/abslfcnlp.lnHMdisp.luad.GseaPreranked.1593476473505", 
      "gsea_report_for_na_pos_1593476473505.xls"
   )
)

voom.luad <- voom.luad[, c(1,4:8,10)]
lnHMdisp.luad <- lnHMdisp.luad[, c(1,4:8,10)]
dim(voom.luad) # 3373 7
dim(lnHMdisp.luad) # 2836 7

voom.luadfdr.05 <- voom.luad %>% 
   filter(FDR.q.val < 0.05)
lnHMdisp.luadfdr.05 <- lnHMdisp.luad %>% 
   filter(FDR.q.val < 0.05)
dim(voom.luadfdr.05) # 642 7
dim(lnHMdisp.luadfdr.05) # 153 7
sum(voom.luadfdr.05$NAME %in% lnHMdisp.luadfdr.05$NAME) # 22
# 620 GO terms identified by mean only, 131 by dispersion only with FDR < 0.05

voom.luadfdr.01 <- voom.luad %>% 
   filter(FDR.q.val < 0.01)
lnHMdisp.luadfdr.01 <- lnHMdisp.luad %>% 
   filter(FDR.q.val < 0.01)
dim(voom.luadfdr.01) # 269 7
dim(lnHMdisp.luadfdr.01) # 72 7
sum(voom.luadfdr.01$NAME %in% lnHMdisp.luadfdr.01$NAME) # 9
# 260 GO terms identified by mean only, 63 by dispersion only with FDR < 0.01


## lusc
voom.lusc <- read.delim(
   here(
      "Results/GSEA results June 2020/jun30/abslfcnlp.voom.lusc.GseaPreranked.1593483651876", 
      "gsea_report_for_na_pos_1593483651876.xls"
   )
)
lnHMdisp.lusc <- read.delim(
   here(
      "Results/GSEA results June 2020/jun30/abslfcnlp.lnHMdisp.lusc.GseaPreranked.1593477201503", 
      "gsea_report_for_na_pos_1593477201503.xls"
   )
)

voom.lusc <- voom.lusc[, c(1,4:8,10)]
lnHMdisp.lusc <- lnHMdisp.lusc[, c(1,4:8,10)]
dim(voom.lusc) # 3564 7
dim(lnHMdisp.lusc) # 3306 7

voom.luscfdr.05 <- voom.lusc %>% 
   filter(FDR.q.val < 0.05)
lnHMdisp.luscfdr.05 <- lnHMdisp.lusc %>% 
   filter(FDR.q.val < 0.05)
dim(voom.luscfdr.05) # 681 7
dim(lnHMdisp.luscfdr.05) # 256 7
sum(voom.luscfdr.05$NAME %in% lnHMdisp.luscfdr.05$NAME) # 61
# 620 GO terms identified by mean only, 195 by dispersion only with FDR < 0.05

voom.luscfdr.01 <- voom.lusc %>% 
   filter(FDR.q.val < 0.01)
lnHMdisp.luscfdr.01 <- lnHMdisp.lusc %>% 
   filter(FDR.q.val < 0.01)
dim(voom.luscfdr.01) # 314 7
dim(lnHMdisp.luscfdr.01) # 116 7
sum(voom.luscfdr.01$NAME %in% lnHMdisp.luscfdr.01$NAME) # 16
# 298 GO terms identified by mean only, 100 by dispersion only with FDR < 0.01


## prad
voom.prad <- read.delim(
   here(
      "Results/GSEA results June 2020/jun30/abslfcnlp.voom.prad.GseaPreranked.1593485427756", 
      "gsea_report_for_na_pos_1593485427756.xls"
   )
)
lnHMdisp.prad <- read.delim(
   here(
      "Results/GSEA results June 2020/jun30/abslfcnlp.lnHMdisp.prad.GseaPreranked.1593478835858", 
      "gsea_report_for_na_pos_1593478835858.xls"
   )
)

voom.prad <- voom.prad[, c(1,4:8,10)]
lnHMdisp.prad <- lnHMdisp.prad[, c(1,4:8,10)]
dim(voom.prad) # 3267 7
dim(lnHMdisp.prad) # 1755 7

voom.pradfdr.05 <- voom.prad %>% 
   filter(FDR.q.val < 0.05)
lnHMdisp.pradfdr.05 <- lnHMdisp.prad %>% 
   filter(FDR.q.val < 0.05)
dim(voom.pradfdr.05) # 775 7
dim(lnHMdisp.pradfdr.05) # 0 7
sum(voom.pradfdr.05$NAME %in% lnHMdisp.pradfdr.05$NAME) # 0
# 775 GO terms identified by mean only, 0 by dispersion only with FDR < 0.05

voom.pradfdr.01 <- voom.prad %>% 
   filter(FDR.q.val < 0.01)
lnHMdisp.pradfdr.01 <- lnHMdisp.prad %>% 
   filter(FDR.q.val < 0.01)
dim(voom.pradfdr.01) # 395 7
dim(lnHMdisp.pradfdr.01) # 0 7
sum(voom.pradfdr.01$NAME %in% lnHMdisp.pradfdr.01$NAME) # 0
# 395 GO terms identified by mean only, 0 by dispersion only with FDR < 0.01


## coad
voom.coad <- read.delim(
   here(
      "Results/GSEA results June 2020/jun30/abslfcnlp.voom.coad.GseaPreranked.1593480916111", 
      "gsea_report_for_na_pos_1593480916111.xls"
   )
)
lnHMdisp.coad <- read.delim(
   here(
      "Results/GSEA results June 2020/jun30/abslfcnlp.lnHMdisp.coad.GseaPreranked.1593475904625", 
      "gsea_report_for_na_pos_1593475904625.xls"
   )
)

voom.coad <- voom.coad[, c(1,4:8,10)]
lnHMdisp.coad <- lnHMdisp.coad[, c(1,4:8,10)]
dim(voom.coad) # 3041 7
dim(lnHMdisp.coad) # 1359 7

voom.coadfdr.05 <- voom.coad %>% 
   filter(FDR.q.val < 0.05)
lnHMdisp.coadfdr.05 <- lnHMdisp.coad %>% 
   filter(FDR.q.val < 0.05)
dim(voom.coadfdr.05) # 503 7
dim(lnHMdisp.coadfdr.05) # 144 7
sum(voom.coadfdr.05$NAME %in% lnHMdisp.coadfdr.05$NAME) # 3
# 500 GO terms identified by mean only, 141 by dispersion only with FDR < 0.05

voom.coadfdr.01 <- voom.coad %>% 
   filter(FDR.q.val < 0.01)
lnHMdisp.coadfdr.01 <- lnHMdisp.coad %>% 
   filter(FDR.q.val < 0.01)
dim(voom.coadfdr.01) # 259 7
dim(lnHMdisp.coadfdr.01) # 72 7
sum(voom.coadfdr.01$NAME %in% lnHMdisp.coadfdr.01$NAME) # 1
# 258 GO terms identified by mean only, 71 by dispersion only with FDR < 0.01


