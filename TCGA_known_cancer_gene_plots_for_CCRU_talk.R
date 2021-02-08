library(here)
folder <- "Results/TCGA paired data results May 2020"

# Import overlap and correlation results
for (i in c("brca", "kirc", "thca", "luad", "lihc", 
            "lihc_long.chain", "lusc", "prad", "coad")) {
  for (j in c("HM", "HM_voom", "HMM_voom", "MDSeq")) {
    for (k in c("known.genes", "pathway.genes", "pathway.genes.minus.known.genes")) {
      assign(
        paste0(i, "_", j, "_mean_v_disp_", k), 
        read.csv(
          here(folder, paste0(i, "_", j, "_mean_v_disp_", k, ".csv")), 
        )
      )
    }
  }
}
rm(folder,i,j,k)

par(mfrow=c(1,1), mgp=c(3.5,1,0), mar=c(5,6,3,1))
plot(kirc_HM_voom_mean_v_disp_known.genes$pvals, kirc_HM_voom_mean_v_disp_known.genes$TPA, 
     type="l", col="red", lwd=3, main="Clear cell renal cell carcinoma", cex.main=3, 
     xlim=c(0,0.01), cex.axis=2, 
     xlab="p-value", ylab="Known cancer-related genes identified", cex.lab=2.5)
lines(kirc_HM_voom_mean_v_disp_known.genes$pvals, kirc_HM_voom_mean_v_disp_known.genes$TPB, col="blue", lwd=3)
lines(kirc_HM_voom_mean_v_disp_known.genes$pvals, kirc_HM_voom_mean_v_disp_known.genes$TPboth, col='grey', lwd=3)
lines(kirc_HM_voom_mean_v_disp_known.genes$pvals, kirc_HM_voom_mean_v_disp_known.genes$TPAonly, col='red', lty=2, lwd=3)
lines(kirc_HM_voom_mean_v_disp_known.genes$pvals, kirc_HM_voom_mean_v_disp_known.genes$TPBonly, col='blue', lty=2, lwd=3)
legend('topleft', bty='n', col=c("red", "blue", "grey", "red", "blue"), lwd=3, lty=c(1,1,1,2,2), cex=2, ncol=2, 
       legend=c(
         "Identified by differential expression", 
         "Identified by differential dispersion", 
         "Identified by both", 
         "Identified by differential expression only", 
         "Identified by differential dispersion only"
       ))


par(mfrow=c(9,12), mar=c(0.2,1,0.2,0), mgp=c(2,0,0))
for (i in c("brca", "kirc", "thca", "luad", "lihc", 
            "lihc_long.chain", "lusc", "prad", "coad")) {
  for (j in c("HM", "HM_voom", "HMM_voom", "MDSeq")) {
    for (k in c("known.genes", "pathway.genes", "pathway.genes.minus.known.genes")) {
      dat <- get(paste0(i, "_", j, "_mean_v_disp_", k))
      plot(dat$pvals, dat$TPA, type="l", col="red", xaxt="n", ann=F, tcl=0, xlim=c(0,0.01))
      lines(dat$pvals, dat$TPB, col='blue')
      lines(dat$pvals, dat$TPboth, col='grey')
      lines(dat$pvals, dat$TPAonly, col='red', lty=2)
      lines(dat$pvals, dat$TPBonly, col='blue', lty=2)
      legend('top', legend=paste0(i, " ", j), bty="n")
    }
  }
}
