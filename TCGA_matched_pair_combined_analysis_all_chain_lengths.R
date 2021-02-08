library(here)
library(RColorBrewer)

qual_col_pals = brewer.pal.info[brewer.pal.info$category == "qual" & 
                                  brewer.pal.info$colorblind == T,]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_vector <- col_vector[c(1:3,9,11,26)]
rm(qual_col_pals)

## Import data ####
folder <- "Results/TCGA all results Sept 2020"

# Fisher's test results for detection of cancer-related genes using threshold
for (i in c("fisher.tests_known.genes", 
            "fisher.tests_pathway.genes.minus.known.genes")) {
  assign(
    i, read.csv(
      here(folder, paste0(i, ".csv")), row.names=1
    )
  )
}

# Wilcoxon rank-sum test results for ranking of cancer-related genes
for (i in c("ranksum.tests_known.genes", 
            "ranksum.tests_pathway.genes.minus.known.genes")) {
  assign(
    i, read.csv(
      here(folder, paste0(i, ".csv")), row.names=1
    )
  )
}

# Overlap and correlation results
for (i in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {
  for (j in c("HM.HM", "MDSeq.MDSeq", "voom.HM")) {
    for (k in c("known.genes", "pathway.genes.minus.known.genes")) {
      assign(
        paste0(i, "_", j, "_mean_v_disp_", k), 
        read.csv(
          here(folder, paste0(i, "_mean_v_disp_", j, "_", k, ".csv")), 
        )
      )
    }
  }
}
rm(folder,i,j,k)


## Ranking test results ####
# Report lnHM log mean and disp, HMM, HM/HM, voom/HM, voom, maybe edgeR QL, MDSeq ZI disp
# rownames(ranksum.tests_known.genes)
ranksum.tests_known.genes[
  c("edgeR.ql", "voom", "mean.lnHM.log.20k", 
    "disp.MDSeq.zi", "disp.lnHM.log.20k", "lnHMM.20k", 
    "voom_lnHM.20k", "lnHM_lnHM.20k"), 
]
# Known cancer genes ranked significantly higher than other genes for all cancers for DE 
# (with edgeR and voom, but not for lusc with HM). For DD, cancer genes ranked significantly 
# higher than other genes for brca, thca, luad, lihc, lusc with MDSeq, and for thca, luad 
# and lusc with HM. For DEDD, cancer genes ranked significantly higher than other genes for 
# all cancers with voom/HM, for all but coad with HM/HM, and for luad, lihc and lusc with HMM.
ranksum.tests_pathway.genes.minus.known.genes[
  c("edgeR.ql", "voom", "mean.lnHM.log.20k", 
    "disp.MDSeq.zi", "disp.lnHM.log.20k", "lnHMM.20k", 
    "voom_lnHM.20k", "lnHM_lnHM.20k"), 
]
# Pathway genes with known genes excluded ranked significantly higher than other genes for 
# all but luad and prad with edgeR, all but prad with voom, and for kirc, thca, lihc and prad 
# with HM. For DD, pathway genes ranked significantly higher than other genes for thca, lihc 
# and lusc with MDSeq, and for thca, luad and lusc with HM. For DEDD, pathway genes ranked 
# significantly higher than other genes for all but prad and coad with voom/HM, for thca, luad, 
# lihc and lusc with HM/HM, and for thca, lihc and lusc with HMM.


## DE v DD correlation results ####
# Look at all possible plots first, with brca as an example
par(mfrow=c(2,2), mar=c(2,2,2,0.5), mgp=c(3,0.7,0))
plot(brca_HM.HM_mean_v_disp_known.genes$pvals, 
     brca_HM.HM_mean_v_disp_known.genes$cor, 
     type="l", 
     # xlim=c(0,0.01), 
     ylim=c(-1,1),
     col="red")
lines(brca_MDSeq.MDSeq_mean_v_disp_known.genes$pvals, 
      brca_MDSeq.MDSeq_mean_v_disp_known.genes$cor, 
     col="blue")
plot(brca_HM.HM_mean_v_disp_known.genes$pvals, 
     brca_HM.HM_mean_v_disp_known.genes$p.cor, 
     type="l", 
     # xlim=c(0,0.01), 
     ylim=c(0,1),
     col="red")
lines(brca_MDSeq.MDSeq_mean_v_disp_known.genes$pvals, 
      brca_MDSeq.MDSeq_mean_v_disp_known.genes$p.cor, 
      col="blue")
plot(brca_HM.HM_mean_v_disp_known.genes$pvals, 
     brca_HM.HM_mean_v_disp_known.genes$overlap / brca_HM.HM_mean_v_disp_known.genes$union, 
     type="l", 
     # xlim=c(0,0.01), 
     ylim=c(0,1),
     col="red")
lines(brca_MDSeq.MDSeq_mean_v_disp_known.genes$pvals, 
      brca_MDSeq.MDSeq_mean_v_disp_known.genes$overlap / brca_MDSeq.MDSeq_mean_v_disp_known.genes$union, 
      col="blue")
plot(brca_HM.HM_mean_v_disp_known.genes$pvals, 
     brca_HM.HM_mean_v_disp_known.genes$p.hyper, 
     type="l", 
     # xlim=c(0,0.01), 
     ylim=c(0,1),
     col="red")
lines(brca_MDSeq.MDSeq_mean_v_disp_known.genes$pvals, 
      brca_MDSeq.MDSeq_mean_v_disp_known.genes$p.hyper, 
      col="blue")
# Overlap or overlap/union is least convincing. Hypergeometric tests clearly show difference 
# between MDSeq and HM, but not really convincing in showing that DD identifies different 
# sets of genes from DE, because p-value is never very low, and is close to 1 for threshold 
# p-values up to nearly 0.1.

# See what each plot looks like for each cancer
par(mfrow=c(4,8), mar=c(2,2,2,0.5), mgp=c(3,0.7,0))
# par(mfrow=c(2,4), mar=c(2,2,2,0.5), mgp=c(3,0.7,0))
for (i in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {
  dat_HM <- get(paste0(i, "_HM.HM_mean_v_disp_known.genes"))
  dat_MD <- get(paste0(i, "_MDSeq.MDSeq_mean_v_disp_known.genes"))
  plot(dat_HM$pvals, 
       dat_HM$cor, 
       type="l", 
       xlim=c(0,0.0001),
       ylim=c(-1,1),
       col="red", main=i)
  lines(dat_MD$pvals, 
        dat_MD$cor, 
        col="blue")
}
# Correlation plots look pretty similar for all.
# Showing up to p = 1e-4, HM always starts at zero then quickly falls, while MDSeq looks to 
# do the same but on very small scale, so can't see without zooming in a lot more.
# par(mfrow=c(2,4), mar=c(2,2,2,0.5), mgp=c(3,0.7,0))
for (i in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {
  dat_HM <- get(paste0(i, "_HM.HM_mean_v_disp_known.genes"))
  dat_MD <- get(paste0(i, "_MDSeq.MDSeq_mean_v_disp_known.genes"))
  plot(dat_HM$pvals, 
       dat_HM$p.cor, 
       type="l", 
       xlim=c(0,0.0001),
       ylim=c(0,1),
       col="red", main=i)
  lines(dat_MD$pvals, 
        dat_MD$p.cor, 
        col="blue")
}
# Same pattern for correlation tests for all, but point at which HM p-val increases changes, 
# and for kirc and coad it never does.
# Showing up to p = 1e-4, for most cancers, MDSeq jumps up to 1 somewhere too close to zero 
# to be able to discern, but for some it increases more slowly (but still faster than HM), 
# and for prad and coad it stays at zero along with HM tp to p = 1e-4.
# par(mfrow=c(2,4), mar=c(2,2,2,0.5), mgp=c(3,0.7,0))
for (i in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {
  dat_HM <- get(paste0(i, "_HM.HM_mean_v_disp_known.genes"))
  dat_MD <- get(paste0(i, "_MDSeq.MDSeq_mean_v_disp_known.genes"))
  plot(dat_HM$pvals, 
       dat_HM$overlap / dat_HM$union, 
       type="l", 
       xlim=c(0,0.0001),
       ylim=c(0,1),
       col="red", main=i)
  lines(dat_MD$pvals, 
        dat_MD$overlap / dat_MD$union, 
        col="blue")
}
# Basically same pattern for overlap/union for all. Clearer difference for prad, and very 
# little difference for luad, lihc, lusc, coad.
# Showing up to p = 1e-4, HM always lower than MDSeq from where HM starts (which looks like 
# 5e-5), but MDSeq does start at zero, and increases too fast to be able to see where it 
# happens.
# par(mfrow=c(2,4), mar=c(2,2,2,0.5), mgp=c(3,0.7,0))
for (i in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {
  dat_HM <- get(paste0(i, "_HM.HM_mean_v_disp_known.genes"))
  dat_MD <- get(paste0(i, "_MDSeq.MDSeq_mean_v_disp_known.genes"))
  plot(dat_HM$pvals, 
       dat_HM$p.hyper, 
       type="l", 
       xlim=c(0,0.0001),
       ylim=c(0,1),
       col="red", main=i)
  lines(dat_MD$pvals, 
        dat_MD$p.hyper, 
        col="blue")
}
# Some big differences in hypergeometric tests between cancers, but p-values always lower for 
# HM than MD. Both stay near(ish) 1 the whole time for thca, luad, lusc and prad. HM very 
# close to 0 up to about 0.3 for kirc and coad, and MD gets close to 0 below about 0.05 for 
# coad.
# Showing up to p = 1e-4, HM and MDSeq both always start at 1, and stay there except for kirc 
# and coad where HM drops quickly to zero (at p = 5e-5).

# Difficult to make general inferences about differences between cancers in different 
# analyses. Probably best to stick to making a general case but using one cancer, and have 
# others in supplementary material. In that case, most general argument is using correlation 
# and correlation tests, and brca is as good as any to show plots for since it's chosen 
# without bias as the biggest sample size and the first that I looked at.
par(mfrow=c(2,1), mar=c(2,2,2,0.5), mgp=c(3,0.7,0))
plot(brca_HM.HM_mean_v_disp_known.genes$pvals, brca_HM.HM_mean_v_disp_known.genes$cor, 
     type="l", lwd=2, xlim=c(0.02,0.98), ylim=c(-1,0.8),
     col=col_vector[2], main=paste0("Spearman correlation between differential expression and ", 
                            "differential dispersion gene lists"))
lines(brca_MDSeq.MDSeq_mean_v_disp_known.genes$pvals, 
      brca_MDSeq.MDSeq_mean_v_disp_known.genes$cor, 
      lwd=2, col=col_vector[1])
legend("bottomright", fill=col_vector[1:2], bty='n', cex=1.5, ncol=1, 
       legend=c("MDSeq", "HM"))
plot(brca_HM.HM_mean_v_disp_known.genes$pvals, brca_HM.HM_mean_v_disp_known.genes$p.cor, 
     type="l", lwd=2, xlim=c(0.02,0.98), ylim=c(0,1),
     col=col_vector[2], main="p-value for hypothesis test for negative correlation")
lines(brca_MDSeq.MDSeq_mean_v_disp_known.genes$pvals, 
      brca_MDSeq.MDSeq_mean_v_disp_known.genes$p.cor, 
      lwd=2, col=col_vector[1])


## Cancer-related genes identified by DE, DD, DEDD ####
# Plot for all cancers, for known genes and pathway genes excluding known genes, for lnHM log DD 
# v voom, lnHM log DD v DE and MDSeq DD v DE
par(mfrow=c(8,6), mar=c(0.2,1,0.2,0), mgp=c(2,0,0))
for (i in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {
  for (j in c("HM.HM", "voom.HM", "MDSeq.MDSeq")) {
    for (k in c("known.genes", "pathway.genes.minus.known.genes")) {
      dat <- get(paste0(i, "_", j, "_mean_v_disp_", k))
      plot(dat$pvals, dat$TP_DE, type="l", col="red", xaxt="n", ann=F, tcl=0, xlim=c(0,0.01))
      lines(dat$pvals, dat$TP_DD, col='blue')
      lines(dat$pvals, dat$TP_both, col='grey')
      lines(dat$pvals, dat$TP_DEonly, col='red', lty=2)
      lines(dat$pvals, dat$TP_DDonly, col='blue', lty=2)
      legend('top', legend=paste0(i, " ", j), bty="n")
    }
  }
}
# Don't think there's any benefit in showing DE v DD for HM since it's clear that voom is better 
# at identifying DE. Also not really any need to show MDSeq, since results are similar and I've 
# already argued that MDSeq DD doesn't identify different sets of genes from DE as well as HM. 

# Include luad and prad in main paper, rest in supplementary material.
# Want to show that DD identifies known and potential cancer-related genes, but maybe only 
# include known genes in main paper since there's not much difference between known and pathway 
# genes for luad and prad, and only include pathway genes in supplementary material.

par(mfcol=c(2,2), mar=c(2,2,2,0.5), mgp=c(3,0.7,0))
plot(luad_voom.HM_mean_v_disp_known.genes$pvals, luad_voom.HM_mean_v_disp_known.genes$TP_DE, 
     type="l", lwd=2, col="red", xlim=c(0.03,0.97))
lines(luad_voom.HM_mean_v_disp_known.genes$pvals, luad_voom.HM_mean_v_disp_known.genes$TP_DD, 
      col='blue', lwd=2)
lines(luad_voom.HM_mean_v_disp_known.genes$pvals, luad_voom.HM_mean_v_disp_known.genes$TP_both, 
      col='grey', lwd=2)
lines(luad_voom.HM_mean_v_disp_known.genes$pvals, luad_voom.HM_mean_v_disp_known.genes$TP_DEonly, 
      col='red', lty=2, lwd=2)
lines(luad_voom.HM_mean_v_disp_known.genes$pvals, luad_voom.HM_mean_v_disp_known.genes$TP_DDonly, 
      col='blue', lty=2, lwd=2)
plot(prad_voom.HM_mean_v_disp_known.genes$pvals, prad_voom.HM_mean_v_disp_known.genes$TP_DE, 
     type="l", lwd=2, col="red", xlim=c(0.03,0.97))
lines(prad_voom.HM_mean_v_disp_known.genes$pvals, prad_voom.HM_mean_v_disp_known.genes$TP_DD, 
      col='blue', lwd=2)
lines(prad_voom.HM_mean_v_disp_known.genes$pvals, prad_voom.HM_mean_v_disp_known.genes$TP_both, 
      col='grey', lwd=2)
lines(prad_voom.HM_mean_v_disp_known.genes$pvals, prad_voom.HM_mean_v_disp_known.genes$TP_DEonly, 
      col='red', lty=2, lwd=2)
lines(prad_voom.HM_mean_v_disp_known.genes$pvals, prad_voom.HM_mean_v_disp_known.genes$TP_DDonly, 
      col='blue', lty=2, lwd=2)
plot(luad_voom.HM_mean_v_disp_known.genes$pvals, luad_voom.HM_mean_v_disp_known.genes$TP_DE, 
     type="l", lwd=2, col="red", xlim=c(0.00002,0.001), ylim=c(0,230))
lines(luad_voom.HM_mean_v_disp_known.genes$pvals, luad_voom.HM_mean_v_disp_known.genes$TP_DD, 
      col='blue', lwd=2)
lines(luad_voom.HM_mean_v_disp_known.genes$pvals, luad_voom.HM_mean_v_disp_known.genes$TP_both, 
      col='grey', lwd=2)
lines(luad_voom.HM_mean_v_disp_known.genes$pvals, luad_voom.HM_mean_v_disp_known.genes$TP_DEonly, 
      col='red', lty=2, lwd=2)
lines(luad_voom.HM_mean_v_disp_known.genes$pvals, luad_voom.HM_mean_v_disp_known.genes$TP_DDonly, 
      col='blue', lty=2, lwd=2)
plot(prad_voom.HM_mean_v_disp_known.genes$pvals, prad_voom.HM_mean_v_disp_known.genes$TP_DE, 
     type="l", lwd=2, col="red", xlim=c(0.00002,0.001), ylim=c(0,550))
lines(prad_voom.HM_mean_v_disp_known.genes$pvals, prad_voom.HM_mean_v_disp_known.genes$TP_DD, 
      col='blue', lwd=2)
lines(prad_voom.HM_mean_v_disp_known.genes$pvals, prad_voom.HM_mean_v_disp_known.genes$TP_both, 
      col='grey', lwd=2)
lines(prad_voom.HM_mean_v_disp_known.genes$pvals, prad_voom.HM_mean_v_disp_known.genes$TP_DEonly, 
      col='red', lty=2, lwd=2)
lines(prad_voom.HM_mean_v_disp_known.genes$pvals, prad_voom.HM_mean_v_disp_known.genes$TP_DDonly, 
      col='blue', lty=2, lwd=2)
# Not sure whether there's any benefit to zooming in on lower p-values. Trend is clearly 
# discernible from looking at whole range. Think I should still show p-values up to 1, rather 
# than only zoomed in, to give full view of number of cancer-related genes present. If zoom in 
# as far as p-values of 0 to 0.001, clearly shows that voom identifies more cancer-related genes 
# than HM DD at start, but that is just because of minimum p-value (and actually would not be 
# the case if I made a different arbitrary decision and called the minimal value p=0 instead of 
# 0.0001 as I have done), so probably not really a disadvantage to show that level of detail. 
# However, still don't really get any more information from zoomed-in version than showing full 
# range, so think I'll stick with just showing full range.

par(mfcol=c(2,1), mar=c(2,2,2,0.5), mgp=c(3,0.7,0))
plot(luad_voom.HM_mean_v_disp_known.genes$pvals, luad_voom.HM_mean_v_disp_known.genes$TP_DE, 
     type="l", lwd=2, col="red", xlim=c(0.03,0.97), 
     main="Lung adenocarcinoma")
lines(luad_voom.HM_mean_v_disp_known.genes$pvals, luad_voom.HM_mean_v_disp_known.genes$TP_DD, 
      col='blue', lwd=2)
lines(luad_voom.HM_mean_v_disp_known.genes$pvals, luad_voom.HM_mean_v_disp_known.genes$TP_both, 
      col='grey', lwd=2)
lines(luad_voom.HM_mean_v_disp_known.genes$pvals, luad_voom.HM_mean_v_disp_known.genes$TP_DEonly, 
      col='red', lty=2, lwd=2)
lines(luad_voom.HM_mean_v_disp_known.genes$pvals, luad_voom.HM_mean_v_disp_known.genes$TP_DDonly, 
      col='blue', lty=2, lwd=2)
legend("right", 
       col=c("red", "blue", "grey", "red", "blue"), lty=c(1,1,1,2,2), 
       lwd=2, bty='n', cex=1.2, ncol=1, 
       legend=c("Differential expression", 
                "Differential dispersion", 
                "Both", 
                "Differential expression only", 
                "Differential dispersion only"))
plot(prad_voom.HM_mean_v_disp_known.genes$pvals, prad_voom.HM_mean_v_disp_known.genes$TP_DE, 
     type="l", lwd=2, col="red", xlim=c(0.03,0.97), 
     main="Prostate adenocarcinoma")
lines(prad_voom.HM_mean_v_disp_known.genes$pvals, prad_voom.HM_mean_v_disp_known.genes$TP_DD, 
      col='blue', lwd=2)
lines(prad_voom.HM_mean_v_disp_known.genes$pvals, prad_voom.HM_mean_v_disp_known.genes$TP_both, 
      col='grey', lwd=2)
lines(prad_voom.HM_mean_v_disp_known.genes$pvals, prad_voom.HM_mean_v_disp_known.genes$TP_DEonly, 
      col='red', lty=2, lwd=2)
lines(prad_voom.HM_mean_v_disp_known.genes$pvals, prad_voom.HM_mean_v_disp_known.genes$TP_DDonly, 
      col='blue', lty=2, lwd=2)


# Would be best to also (or instead) plot by number of discoveries, and could then plot FDR 
# or specificity (or something like them since we don't know true positives, just some sort of 
# proxy) rather than just number of cancer-related genes identified.
# Even with 20k chains, for most cancers HM shows a big jump from zero to a large number of 
# genes identified, so these plots still aren't really useful. They could be if I can find a 
# way of ranking genes with the same p-value differently (like sorting by LFC).




## Overlap and correlation plots by number of genes instead of p-value ####
# Question is how to plot when there are different numbers of genes for DE and DD for given 
# cor, etc.May make sense to use union, i.e. number of genes identified by either method.
par(mfrow=c(4,8), mar=c(2,2,2,0.5), mgp=c(3,0.7,0))
# par(mfrow=c(2,4), mar=c(2,2,2,0.5), mgp=c(3,0.7,0))
for (i in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {
  dat_HM <- get(paste0(i, "_HM.HM_mean_v_disp_known.genes"))
  dat_MD <- get(paste0(i, "_MDSeq.MDSeq_mean_v_disp_known.genes"))
  plot(dat_HM$union, 
       dat_HM$cor, 
       type="l", 
       xlim=c(0,20000),
       ylim=c(-1,1),
       col="red", main=i)
  lines(dat_MD$union, 
        dat_MD$cor, 
        col="blue")
}
# Correlation always lower for HM than MDSeq where it exists, but looking at up to 1000, 
# correlation is well below zero for MDSeq for reasonable numbers of genes which HM can't 
# differentiate between.
# par(mfrow=c(2,4), mar=c(2,2,2,0.5), mgp=c(3,0.7,0))
for (i in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {
  dat_HM <- get(paste0(i, "_HM.HM_mean_v_disp_known.genes"))
  dat_MD <- get(paste0(i, "_MDSeq.MDSeq_mean_v_disp_known.genes"))
  plot(dat_HM$union, 
       dat_HM$p.cor, 
       type="l", 
       xlim=c(0,20000),
       ylim=c(0,1),
       col="red", main=i)
  lines(dat_MD$union, 
        dat_MD$p.cor, 
        col="blue")
}
# Correlation p-value is lower for HM than MDSeq wherever it exists and is zero or very 
# close to it for nearly the whole range, but looking at up to 1000 genes, it's also zero 
# or very close to it for MDSeq for reasonable numbers of genes.
# par(mfrow=c(2,4), mar=c(2,2,2,0.5), mgp=c(3,0.7,0))
for (i in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {
  dat_HM <- get(paste0(i, "_HM.HM_mean_v_disp_known.genes"))
  dat_MD <- get(paste0(i, "_MDSeq.MDSeq_mean_v_disp_known.genes"))
  plot(dat_HM$union, 
       dat_HM$overlap / dat_HM$union, 
       type="l", 
       xlim=c(0,20000),
       ylim=c(0,1),
       col="red", main=i)
  lines(dat_MD$union, 
        dat_MD$overlap / dat_MD$union, 
        col="blue")
}
# Overlap/union always lower for HM than MDSeq where it's defined, but very low for MDSeq 
# for up to at least 1000 genes, where it isn't defined for HM.
# par(mfrow=c(2,4), mar=c(2,2,2,0.5), mgp=c(3,0.7,0))
for (i in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {
  dat_HM <- get(paste0(i, "_HM.HM_mean_v_disp_known.genes"))
  dat_MD <- get(paste0(i, "_MDSeq.MDSeq_mean_v_disp_known.genes"))
  plot(dat_HM$union, 
       dat_HM$p.hyper, 
       type="l", 
       xlim=c(0,20000),
       ylim=c(0,1),
       col="red", main=i)
  lines(dat_MD$union, 
        dat_MD$p.hyper, 
        col="blue")
}
# Hypergeometric test p-values always lower for HM than MDSeq where they exist for HM, but 
# never below 1 for HM until at least around 17,000 genes, so nothing to differentiate 
# between the methods really.

# This shows that my conclusions from using p-values aren't valid for MDSeq. Need to find 
# a way to differentiate between high-ranked genes for HM to be able to make this comparison. 
# It's valid to use this data to show that DE and DD identify different sets of genes, but 
# not to say that HM is better than MDSeq really, because I can only say that for numbers of 
# genes that are never going to be used in practice.

