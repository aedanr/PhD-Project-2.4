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


## Correlation plots ####
# Assess correlation plots for all cancers, by p-value and by number of discoveries (union)
par(mfrow=c(4,8), mar=c(2,2,2,0.5), mgp=c(3,0.7,0))
for (i in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {
  dat <- get(paste0(i, "_HM.HM_mean_v_disp"))
  plot(dat$modp, 
       dat$cor, 
       type="l", 
       # xlim=c(0,0.001),
       ylim=c(-1,1),
       col="red", main=i)
}
for (i in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {
  dat <- get(paste0(i, "_HM.HM_mean_v_disp"))
  plot(dat$modp, 
       dat$p.cor, 
       type="l", 
       # xlim=c(0,0.001),
       ylim=c(0,1),
       col="red", main=i)
}
for (i in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {
  dat <- get(paste0(i, "_HM.HM_mean_v_disp"))
  plot(dat$union, 
       dat$cor, 
       type="l", 
       # xlim=c(0,1000),
       ylim=c(-1,1),
       col="red", main=i)
}
for (i in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {
  dat <- get(paste0(i, "_HM.HM_mean_v_disp"))
  plot(dat$union, 
       dat$p.cor, 
       type="l", 
       # xlim=c(0,1000),
       ylim=c(0,1),
       col="red", main=i)
}

# brca only as shown in paper, but probably just show all together for thesis.
par(mfcol=c(3,1), mar=c(3,3,2.5,0.5), mgp=c(3,1,0))
plot(brca_HM.HM_mean_v_disp$modp, brca_HM.HM_mean_v_disp$cor, 
     type="l", lwd=2, xlim=c(0,1), ylim=c(-1,0.5), cex.axis=1.8, ylab=NA, xlab=NA, 
     main=paste0("Spearman correlation between differential expression and ", 
                                    "differential dispersion gene lists"), 
     cex.main=1.8)
abline(h=0, lty=2)
plot(brca_HM.HM_mean_v_disp$modp, brca_HM.HM_mean_v_disp$p.cor, 
     type="l", lwd=2, xlim=c(0,1), ylim=c(0,1), cex.axis=1.8, ylab=NA, xlab=NA,
     main="p-value for hypothesis test for negative correlation", 
     cex.main=1.8)
plot(brca_HM.HM_mean_v_disp$modp, brca_HM.HM_mean_v_disp$cor, 
     type="l", lwd=2, xlim=c(0,0.001), ylim=c(-1,0.5), cex.axis=1.8, ylab=NA, xlab=NA, 
     main=paste0("Spearman correlation for thresholds up to 0.001"), 
     cex.main=1.8)
abline(h=0, lty=2)


## Number of discoveries by DE and DD ####
par(mfcol=c(2,1), mar=c(2.5,3,2.5,0.5), mgp=c(3,1,0))
plot(luad_voom.HM_mean_v_disp$pvals, luad_voom.HM_mean_v_disp$TP_DE, 
     type="l", lwd=2, col="red", xlim=c(0,1), cex.axis=1.5, xlab=NA, ylab=NA, 
     main="Lung adenocarcinoma", cex.main=1.5)
lines(luad_voom.HM_mean_v_disp$pvals, luad_voom.HM_mean_v_disp$TP_DD, 
      col='blue', lwd=2)
lines(luad_voom.HM_mean_v_disp$pvals, luad_voom.HM_mean_v_disp$TP_both, 
      col='grey', lwd=2)
lines(luad_voom.HM_mean_v_disp$pvals, luad_voom.HM_mean_v_disp$TP_DEonly, 
      col='red', lty=2, lwd=2)
lines(luad_voom.HM_mean_v_disp$pvals, luad_voom.HM_mean_v_disp$TP_DDonly, 
      col='blue', lty=2, lwd=2)
legend("right", 
       col=c("red", "blue", "grey", "red", "blue"), lty=c(1,1,1,2,2), 
       lwd=2, bty='n', cex=1.2, ncol=1, 
       legend=c("Differential expression", 
                "Differential dispersion", 
                "Both", 
                "Differential expression only", 
                "Differential dispersion only"))
plot(prad_voom.HM_mean_v_disp$pvals, prad_voom.HM_mean_v_disp$TP_DE, 
     type="l", lwd=2, col="red", xlim=c(0,1), cex.axis=1.5, xlab=NA, ylab=NA, 
     main="Prostate adenocarcinoma", cex.main=1.5)
lines(prad_voom.HM_mean_v_disp$pvals, prad_voom.HM_mean_v_disp$TP_DD, 
      col='blue', lwd=2)
lines(prad_voom.HM_mean_v_disp$pvals, prad_voom.HM_mean_v_disp$TP_both, 
      col='grey', lwd=2)
lines(prad_voom.HM_mean_v_disp$pvals, prad_voom.HM_mean_v_disp$TP_DEonly, 
      col='red', lty=2, lwd=2)
lines(prad_voom.HM_mean_v_disp$pvals, prad_voom.HM_mean_v_disp$TP_DDonly, 
      col='blue', lty=2, lwd=2)

# Probably no need to include zoomed in plots
# plot(luad_voom.HM_mean_v_disp$pvals, luad_voom.HM_mean_v_disp$TP_DE, 
#      type="l", lwd=2, col="red", xlim=c(0,0.001), 
#      main="Lung adenocarcinoma, thresholds up to 0.001")
# lines(luad_voom.HM_mean_v_disp$pvals, luad_voom.HM_mean_v_disp$TP_DD, 
#       col='blue', lwd=2)
# lines(luad_voom.HM_mean_v_disp$pvals, luad_voom.HM_mean_v_disp$TP_both, 
#       col='grey', lwd=2)
# lines(luad_voom.HM_mean_v_disp$pvals, luad_voom.HM_mean_v_disp$TP_DEonly, 
#       col='red', lty=2, lwd=2)
# lines(luad_voom.HM_mean_v_disp$pvals, luad_voom.HM_mean_v_disp$TP_DDonly, 
#       col='blue', lty=2, lwd=2)
# plot(prad_voom.HM_mean_v_disp$pvals, prad_voom.HM_mean_v_disp$TP_DE, 
#      type="l", lwd=2, col="red", xlim=c(0,0.001), 
#      main="Prostate adenocarcinoma, thresholds up to 0.001")
# lines(prad_voom.HM_mean_v_disp$pvals, prad_voom.HM_mean_v_disp$TP_DD, 
#       col='blue', lwd=2)
# lines(prad_voom.HM_mean_v_disp$pvals, prad_voom.HM_mean_v_disp$TP_both, 
#       col='grey', lwd=2)
# lines(prad_voom.HM_mean_v_disp$pvals, prad_voom.HM_mean_v_disp$TP_DEonly, 
#       col='red', lty=2, lwd=2)
# lines(prad_voom.HM_mean_v_disp$pvals, prad_voom.HM_mean_v_disp$TP_DDonly, 
#       col='blue', lty=2, lwd=2)

# Probably no need to plot by discoveries
# par(mfrow=c(4,4), mar=c(0.2,1,0.2,0), mgp=c(2,0,0))
# # TP from DE (red), DD (blue) by modified p-value
# for (i in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {
#   dat <- get(paste0(i, "_voom.HM_mean_v_disp"))
#   plot(dat$pvals, dat$TP_DE, type="l", col="red", xaxt="n", ann=F, tcl=0)
#   lines(dat$pvals, dat$TP_DD, col='blue')
#   lines(dat$pvals, dat$TP_both, col='grey')
#   lines(dat$pvals, dat$TP_DEonly, col='red', lty=2)
#   lines(dat$pvals, dat$TP_DDonly, col='blue', lty=2)
#   legend('top', legend=i, bty="n")
# }
# # TP from DE (red), DD (blue) by number of discoveries
# for (i in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {
#   dat <- get(paste0(i, "_voom.HM_mean_v_disp"))
#   plot(dat$listDE, dat$TP_DE, type="l", col="red", xaxt="n", ann=F, tcl=0)
#   lines(dat$listDD, dat$TP_DD, col='blue')
#   legend('top', legend=i, bty="n")
# }
# rm(i,dat)


## All correlation plots ####
# Thought it would be better to also include p-values zoomed in to p_threshold < 0.001 
# as for correlations, but doesn't really add anything because p-values are all 
# indistinguishable from zero for thresholds between 10^-5 and 0.001, so stick with 
# zoomed in only for correlations. On the other hand, if I want to discuss the p-values 
# being small for the highest ranked genes, it's better to be able to show that, even 
# if it is with plots that don't look very nice.
par(mfcol=c(4,8), mar=c(2.5,3.5,2.5,0.5), mgp=c(3,1.2,0))
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
       ylim=c(0,1),
       cex.axis=2)
  plot(dat$modp, 
       dat$cor, 
       type="l", lwd=2, xlab=NA, ylab=NA, 
       xlim=c(0,0.001),
       ylim=c(-1,1),
       cex.axis=2)
  abline(h=0, lty=2)
  plot(dat$modp, 
       dat$p.cor, 
       type="l", lwd=2, xlab=NA, ylab=NA, 
       xlim=c(0,0.001),
       ylim=c(0,1),
       cex.axis=2)
}


## All DE v DD TP plots ####
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
