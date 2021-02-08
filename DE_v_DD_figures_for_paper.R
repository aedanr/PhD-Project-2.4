library(here)
library(RColorBrewer)

qual_col_pals = brewer.pal.info[brewer.pal.info$category == "qual" & 
                                  brewer.pal.info$colorblind == T,]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_vector <- col_vector[c(1:3,9,11,26)]
rm(qual_col_pals)

## Import data ####
folder <- "Results/TCGA all results Sept 2020"
for (i in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {
  for (j in c("HM.HM", "voom.HM")) {
    assign(
      paste0(i, "_", j, "_mean_v_disp"), 
      read.csv(
        here(folder, paste0(i, "_mean_v_disp_", j, "_known.genes_lfc.tiebreaker.csv")), 
      )
    )
  }
  assign(
    paste0(i, "_MDSeq.MDSeq_mean_v_disp"), 
    read.csv(
      here(folder, paste0(i, "_mean_v_disp_MDSeq.MDSeq_known.genes.csv"))
    )
  )
}
rm(folder,i,j)


## Correlation and DE v DD plots for luad and prad, for p-value thresholds up to 0.1
par(mfcol=c(2,2), mar=c(2.5,3,4,0.5), mgp=c(3,1.5,0))
plot(luad_HM.HM_mean_v_disp$modp, luad_HM.HM_mean_v_disp$cor, 
     type="l", lwd=2, xlim=c(0,0.1),  cex.axis=2, xlab=NA, ylab=NA, yaxt='n', 
     main=paste0("Spearman correlation between differential expression\n", 
                 "and differential dispersion gene lists"), 
     cex.main=2)
axis(side=2, cex.axis=2, mgp=c(3,1,0))
mtext(side=3, line=-2, cex=1.5, text='Lung adenocarcinoma')
abline(h=0, lty=2)
plot(prad_HM.HM_mean_v_disp$modp, prad_HM.HM_mean_v_disp$cor, 
     type="l", lwd=2, xlim=c(0,0.1), cex.axis=2, ann=F, yaxt='n', 
     cex.main=2)
axis(side=2, cex.axis=2, mgp=c(3,1,0))
mtext(side=3, line=-2, cex=1.5, text='Prostate adenocarcinoma')
abline(h=0, lty=2)

plot(luad_voom.HM_mean_v_disp$pvals, luad_voom.HM_mean_v_disp$TP_DE, 
     type="l", lwd=2, col="red", xlim=c(0,0.1), ylim=c(0,280), cex.axis=2, xlab=NA, ylab=NA, yaxt='n', 
     main=paste0("Previously identified cancer genes identified by\n", 
                 "differences in mean and/or dispersion"), 
     cex.main=2)
lines(luad_voom.HM_mean_v_disp$pvals, luad_voom.HM_mean_v_disp$TP_DD, 
      col='blue', lwd=2)
lines(luad_voom.HM_mean_v_disp$pvals, luad_voom.HM_mean_v_disp$TP_both, 
      col='grey', lwd=2)
lines(luad_voom.HM_mean_v_disp$pvals, luad_voom.HM_mean_v_disp$TP_DEonly, 
      col='red', lty=2, lwd=2)
lines(luad_voom.HM_mean_v_disp$pvals, luad_voom.HM_mean_v_disp$TP_DDonly, 
      col='blue', lty=2, lwd=2)
axis(side=2, cex.axis=2, at=c(100, 200), mgp=c(3,1,0))
mtext(side=3, line=-2, cex=1.5, text='Lung adenocarcinoma')
legend(x=0.045, y=165, 
       col=c("red", "blue", "grey", "red", "blue"), lty=c(1,1,1,2,2), 
       lwd=2, bty='n', cex=1.5, ncol=2, 
       legend=c("Mean", 
                "Dispersion", 
                "Both", 
                "Mean only", 
                "Dispersion only"))
plot(prad_voom.HM_mean_v_disp$pvals, prad_voom.HM_mean_v_disp$TP_DE, 
     type="l", lwd=2, col="red", xlim=c(0,0.1), ylim=c(0,750), cex.axis=2, ann=F, yaxt='n')
lines(prad_voom.HM_mean_v_disp$pvals, prad_voom.HM_mean_v_disp$TP_DD, 
      col='blue', lwd=2)
lines(prad_voom.HM_mean_v_disp$pvals, prad_voom.HM_mean_v_disp$TP_both, 
      col='grey', lwd=2)
lines(prad_voom.HM_mean_v_disp$pvals, prad_voom.HM_mean_v_disp$TP_DEonly, 
      col='red', lty=2, lwd=2)
lines(prad_voom.HM_mean_v_disp$pvals, prad_voom.HM_mean_v_disp$TP_DDonly, 
      col='blue', lty=2, lwd=2)
axis(side=2, cex.axis=2, mgp=c(3,1,0))
mtext(side=3, line=-2, cex=1.5, text='Prostate adenocarcinoma')

## Alternative with sub-figures A and B
par(mfcol=c(2,2), mar=c(2.5,3.5,2.5,2), mgp=c(3,1.5,0))
plot(luad_HM.HM_mean_v_disp$modp, luad_HM.HM_mean_v_disp$cor, 
     type="l", lwd=2, xlim=c(0,0.1),  cex.axis=2, xlab=NA, ylab=NA, yaxt='n', 
     main="a", adj=0, cex.main=3)
axis(side=2, cex.axis=2, mgp=c(3,1,0))
mtext(side=3, line=-2, cex=1.5, text='Lung adenocarcinoma')
abline(h=0, lty=2)
plot(prad_HM.HM_mean_v_disp$modp, prad_HM.HM_mean_v_disp$cor, 
     type="l", lwd=2, xlim=c(0,0.1), cex.axis=2, ann=F, yaxt='n', 
     cex.main=2)
axis(side=2, cex.axis=2, mgp=c(3,1,0))
mtext(side=3, line=-2, cex=1.5, text='Prostate adenocarcinoma')
abline(h=0, lty=2)

plot(luad_voom.HM_mean_v_disp$pvals, luad_voom.HM_mean_v_disp$TP_DE, 
     type="l", lwd=2, col="red", xlim=c(0,0.1), ylim=c(0,280), cex.axis=2, xlab=NA, ylab=NA, yaxt='n', 
     main="b", adj=0, cex.main=3)
lines(luad_voom.HM_mean_v_disp$pvals, luad_voom.HM_mean_v_disp$TP_DD, 
      col='blue', lwd=2)
lines(luad_voom.HM_mean_v_disp$pvals, luad_voom.HM_mean_v_disp$TP_both, 
      col='grey', lwd=2)
lines(luad_voom.HM_mean_v_disp$pvals, luad_voom.HM_mean_v_disp$TP_DEonly, 
      col='red', lty=2, lwd=2)
lines(luad_voom.HM_mean_v_disp$pvals, luad_voom.HM_mean_v_disp$TP_DDonly, 
      col='blue', lty=2, lwd=2)
axis(side=2, cex.axis=2, at=c(100, 200), mgp=c(3,1,0))
mtext(side=3, line=-2, cex=1.5, text='Lung adenocarcinoma')
legend(x=0.045, y=165, 
       col=c("red", "blue", "grey", "red", "blue"), lty=c(1,1,1,2,2), 
       lwd=2, bty='n', cex=1.5, ncol=2, 
       legend=c("Mean", 
                "Dispersion", 
                "Both", 
                "Mean only", 
                "Dispersion only"))
plot(prad_voom.HM_mean_v_disp$pvals, prad_voom.HM_mean_v_disp$TP_DE, 
     type="l", lwd=2, col="red", xlim=c(0,0.1), ylim=c(0,750), cex.axis=2, ann=F, yaxt='n')
lines(prad_voom.HM_mean_v_disp$pvals, prad_voom.HM_mean_v_disp$TP_DD, 
      col='blue', lwd=2)
lines(prad_voom.HM_mean_v_disp$pvals, prad_voom.HM_mean_v_disp$TP_both, 
      col='grey', lwd=2)
lines(prad_voom.HM_mean_v_disp$pvals, prad_voom.HM_mean_v_disp$TP_DEonly, 
      col='red', lty=2, lwd=2)
lines(prad_voom.HM_mean_v_disp$pvals, prad_voom.HM_mean_v_disp$TP_DDonly, 
      col='blue', lty=2, lwd=2)
axis(side=2, cex.axis=2, mgp=c(3,1,0))
mtext(side=3, line=-2, cex=1.5, text='Prostate adenocarcinoma')


## All correlation plots for supplementary material ####
par(mfcol=c(3,8), mar=c(2.5,3.5,2.5,0.5), mgp=c(3,1.2,0))
for (i in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {
  dat <- get(paste0(i, "_HM.HM_mean_v_disp"))
  plot(dat$modp, 
       dat$cor, 
       type="l", lwd=2, xlab=NA, ylab=NA, 
       ylim=c(-1,1),
       cex.axis=2, 
       main=toupper(i), cex.main=3)
  abline(h=0, lty=2)
  plot(dat$modp, 
       dat$p.cor, 
       type="l", lwd=2, xlab=NA, ylab=NA, 
       # xlim=c(0,0.001),
       ylim=c(0,1),
       cex.axis=2)
  plot(dat$modp, 
       dat$cor, 
       type="l", lwd=2, xlab=NA, ylab=NA, 
       xlim=c(0,0.001),
       ylim=c(-1,1),
       cex.axis=2)
  abline(h=0, lty=2)
}


## All DE v DD TP plots for supplementary material ####
par(mfrow=c(2,4), mar=c(3.5,3.5,2.5,1), mgp=c(3,1.5,0))
for (i in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {
  dat <- get(paste0(i, "_voom.HM_mean_v_disp"))
plot(dat$pvals, dat$TP_DE, 
     type="l", lwd=2, col="red", xlim=c(0,1), cex.axis=2.5, xlab=NA, ylab=NA, 
     main=toupper(i), cex.main=3)
lines(dat$pvals, dat$TP_DD, 
      col='blue', lwd=2)
lines(dat$pvals, dat$TP_both, 
      col='grey', lwd=2)
lines(dat$pvals, dat$TP_DEonly, 
      col='red', lty=2, lwd=2)
lines(dat$pvals, dat$TP_DDonly, 
      col='blue', lty=2, lwd=2)
}


## Correlation plots and DE v DD together ####
par(mfcol=c(4,8), mar=c(2.5,3,2.5,0.5), mgp=c(3,1.5,0))
for (i in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {
  dat <- get(paste0(i, "_HM.HM_mean_v_disp"))
  plot(dat$modp, 
       dat$cor, 
       type="l", lwd=2, xlab=NA, ylab=NA, yaxt='n', xaxt='n', 
       ylim=c(-1,1),
       main=toupper(i), cex.main=3)
  axis(side=1, cex.axis=2.5, at=c(0,0.5,1), labels=c("0", "0.5", "1"))
  if(i == "brca") {
    axis(side=2, cex.axis=2.5, at=c(-0.8,0,0.8), mgp=c(3,1,0))
  } else {
    axis(side=2, cex.axis=2.5, at=c(-0.8,0,0.8), labels=F)
  }
  abline(h=0, lty=2)
  plot(dat$modp, 
       dat$p.cor, 
       type="l", lwd=2, ann=F, yaxt='n', xaxt='n', 
       ylim=c(0,1))
  axis(side=1, cex.axis=2.5, at=c(0,0.5,1), labels=c("0", "0.5", "1"))
  if(i == "brca") {
    axis(side=2, cex.axis=2.5, at=c(-0,1), mgp=c(3,1,0))
  } else {
    axis(side=2, cex.axis=2.5, at=c(0,1), labels=F)
  }
  plot(dat$modp, 
       dat$cor, 
       type="l", lwd=2, ann=F, yaxt='n', xaxt='n', 
       xlim=c(0,0.001),
       ylim=c(-1,1))
  axis(side=1, cex.axis=2.5, at=c(0,0.0005), labels=c("0", "0.0005"))
  if(i == "brca") {
    axis(side=2, cex.axis=2.5, at=c(-0.8,0,0.8), mgp=c(3,1,0))
  } else {
    axis(side=2, cex.axis=2.5, at=c(-0.8,0,0.8), labels=F)
  }
  abline(h=0, lty=2)
  dat <- get(paste0(i, "_voom.HM_mean_v_disp"))
  plot(dat$pvals, dat$TP_DE, 
       type="l", lwd=2, col="red", xlim=c(0,1), ann=F, yaxt='n', xaxt='n')
  lines(dat$pvals, dat$TP_DD, 
        col='blue', lwd=2)
  lines(dat$pvals, dat$TP_both, 
        col='grey', lwd=2)
  lines(dat$pvals, dat$TP_DEonly, 
        col='red', lty=2, lwd=2)
  lines(dat$pvals, dat$TP_DDonly, 
        col='blue', lty=2, lwd=2)
  axis(side=1, cex.axis=2.5, at=c(0,0.5,1), labels=c("0", "0.5", "1"))
  axis(side=2, cex.axis=2.5, mgp=c(3,1,0))
}
