library(here)
library(RColorBrewer)

qual_col_pals = brewer.pal.info[brewer.pal.info$category == "qual" & 
                                  brewer.pal.info$colorblind == T,]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_vector <- col_vector[c(1:3,9,11,19)]
rm(qual_col_pals)

## Import data ####
folder <- "Results/TCGA paired data results May 2020"

# Fisher's test results for detection of cancer-related genes using threshold
for (i in c("fisher.tests.weak_known.genes", 
            "fisher.tests.strong_known.genes", 
            "fisher.tests.weak_pathway.genes", 
            "fisher.tests.strong_pathway.genes", 
            "fisher.tests_pathway.genes.minus.known.genes")) {
  assign(
    i, read.csv(
      here(folder, paste0(i, ".csv")), row.names=1
    )
  )
}

# Wilcoxon rank-sum test results for ranking of cancer-related genes
for (i in c("ranksum.tests.weak_known.genes", 
            "ranksum.tests.strong_known.genes", 
            "ranksum.tests.weak_pathway.genes", 
            "ranksum.tests.strong_pathway.genes", 
            "ranksum.tests_pathway.genes.minus.known.genes")) {
  assign(
    i, read.csv(
      here(folder, paste0(i, ".csv")), row.names=1
    )
  )
}

# Overlap and correlation results
for (i in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {
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


## Ranking test results ####
# Report lnHM log mean and disp, HMM, HM/HM, voom/HM, voom, maybe edgeR QL, MDSeq ZI disp
# rownames(ranksum.tests.weak_known.genes)
ranksum.tests.weak_known.genes[
  c("p.edgeR.ql", "p.voom", "p.mean.lnHM.log", 
    "p.disp.MDSeq.zi", "p.disp.lnHM.log", "prob.lnHMM", 
    "voom_lnHM.log", "lnHM.log_lnHM.log"), 
  c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")
]
# Known cancer genes ranked significantly higher than other genes for all cancers for DE 
# (with edgeR and voom, but not for lusc with HM). For DD, cancer genes ranked significantly 
# higher than other genes for brca, thca, luad, lihc, lusc with MDSeq, and for thca, luad 
# and lusc with HM. For DEDD, cancer genes ranked significantly higher than other genes for 
# all cancers with voom/HM, for brca, thca, luac, lihc, lusc and prad with HM/HM, and for 
# luad, lihc and lusc with HMM.
ranksum.tests_pathway.genes.minus.known.genes[
  c("p.edgeR.ql", "p.voom", "p.mean.lnHM.log", 
    "p.disp.MDSeq.zi", "p.disp.lnHM.log", "prob.lnHMM", 
    "voom_lnHM.log", "lnHM.log_lnHM.log"), 
  c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")
]
# Pathway genes with known genes excluded ranked significantly higher than other genes for 
# all but luad and prad with edgeR, all but prad with voom, and for kirc, thca and lihc with 
# HM. For DD, pathway genes ranked significantly higher than other genes for thca, lihc and 
# lusc with MDSeq, and for thca, luad and lusc with HM. For DEDD, pathway genes ranked 
# significantly higher than other genes for all but prad with voom/HM, for thca, luad, lihc 
# and lusc with HM/HM, and for thca and lusc with HMM.


## DE v DD correlation results ####
# Look at all possible plots first, with brca as an example
par(mfrow=c(2,2), mar=c(2,2,2,0.5), mgp=c(3,0.7,0))
plot(brca_HM_mean_v_disp_known.genes$pvals, 
     brca_HM_mean_v_disp_known.genes$cor, 
     type="l", 
     # xlim=c(0,0.01), 
     ylim=c(-1,1),
     col="red")
lines(brca_MDSeq_mean_v_disp_known.genes$pvals, 
      brca_MDSeq_mean_v_disp_known.genes$cor, 
     col="blue")
plot(brca_HM_mean_v_disp_known.genes$pvals, 
     brca_HM_mean_v_disp_known.genes$p.cor, 
     type="l", 
     # xlim=c(0,0.01), 
     ylim=c(0,1),
     col="red")
lines(brca_MDSeq_mean_v_disp_known.genes$pvals, 
      brca_MDSeq_mean_v_disp_known.genes$p.cor, 
      col="blue")
plot(brca_HM_mean_v_disp_known.genes$pvals, 
     brca_HM_mean_v_disp_known.genes$overlap / brca_HM_mean_v_disp_known.genes$union, 
     type="l", 
     # xlim=c(0,0.01), 
     ylim=c(0,1),
     col="red")
lines(brca_MDSeq_mean_v_disp_known.genes$pvals, 
      brca_MDSeq_mean_v_disp_known.genes$overlap / brca_MDSeq_mean_v_disp_known.genes$union, 
      col="blue")
plot(brca_HM_mean_v_disp_known.genes$pvals, 
     brca_HM_mean_v_disp_known.genes$p.hyper, 
     type="l", 
     # xlim=c(0,0.01), 
     ylim=c(0,1),
     col="red")
lines(brca_MDSeq_mean_v_disp_known.genes$pvals, 
      brca_MDSeq_mean_v_disp_known.genes$p.hyper, 
      col="blue")
# Overlap or overlap/union is least convincing. Hypergeometric tests clearly show difference 
# between MDSeq and HM, but not really convincing in showing that DD identifies different 
# sets of genes from DE, because p-value is never very low, and is close to 1 for threshold 
# p-values up to nearly 0.1.

# See what each plot looks like for each cancer
par(mfrow=c(4,8), mar=c(2,2,2,0.5), mgp=c(3,0.7,0))
# par(mfrow=c(2,4), mar=c(2,2,2,0.5), mgp=c(3,0.7,0))
for (i in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {
  dat_HM <- get(paste0(i, "_HM_mean_v_disp_known.genes"))
  dat_MD <- get(paste0(i, "_MDSeq_mean_v_disp_known.genes"))
  plot(dat_HM$pvals, 
       dat_HM$cor, 
       type="l", 
       # xlim=c(0,0.01), 
       ylim=c(-1,1),
       col="red", main=i)
  lines(dat_MD$pvals, 
        dat_MD$cor, 
        col="blue")
}
# Correlation plots look pretty similar for all.
# par(mfrow=c(2,4), mar=c(2,2,2,0.5), mgp=c(3,0.7,0))
for (i in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {
  dat_HM <- get(paste0(i, "_HM_mean_v_disp_known.genes"))
  dat_MD <- get(paste0(i, "_MDSeq_mean_v_disp_known.genes"))
  plot(dat_HM$pvals, 
       dat_HM$p.cor, 
       type="l", 
       # xlim=c(0,0.01), 
       ylim=c(0,1),
       col="red", main=i)
  lines(dat_MD$pvals, 
        dat_MD$p.cor, 
        col="blue")
}
# Same pattern for correlation tests for all, but point at which HM p-val increases changes, 
# and for kirc and coad it never does.
# par(mfrow=c(2,4), mar=c(2,2,2,0.5), mgp=c(3,0.7,0))
for (i in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {
  dat_HM <- get(paste0(i, "_HM_mean_v_disp_known.genes"))
  dat_MD <- get(paste0(i, "_MDSeq_mean_v_disp_known.genes"))
  plot(dat_HM$pvals, 
       dat_HM$overlap / dat_HM$union, 
       type="l", 
       # xlim=c(0,0.01), 
       ylim=c(0,1),
       col="red", main=i)
  lines(dat_MD$pvals, 
        dat_MD$overlap / dat_MD$union, 
        col="blue")
}
# Basically same pattern for overlap/union for all. Clearer difference for prad, and very 
# little difference for luad, lihc, lusc, coad.
# par(mfrow=c(2,4), mar=c(2,2,2,0.5), mgp=c(3,0.7,0))
for (i in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {
  dat_HM <- get(paste0(i, "_HM_mean_v_disp_known.genes"))
  dat_MD <- get(paste0(i, "_MDSeq_mean_v_disp_known.genes"))
  plot(dat_HM$pvals, 
       dat_HM$p.hyper, 
       type="l", 
       # xlim=c(0,0.01), 
       ylim=c(0,1),
       col="red", main=i)
  lines(dat_MD$pvals, 
        dat_MD$p.hyper, 
        col="blue")
}
# Some big differences in hypergeometric tests between cancers, but p-values always lower for 
# HM than MD. Both stay near(ish) 1 the whole time for thca, luad, lusc and prad. HM very 
# close to 0 up to about 0.3 for kirc and coad, and MD gets close to 0 below about 0.05 for 
# MD.

# Difficult to make general inferences about differences between cancers in different 
# analyses. Probably best to stick to making a general case but using one cancer, and have 
# others in supplementary material. In that case, most general argument is using correlation 
# and correlation tests, and brca is as good as any to show plots for since it's chosen 
# without bias as the biggest sample size and the first that I looked at.
par(mfrow=c(2,1), mar=c(2,2,2,0.5), mgp=c(3,0.7,0))
plot(brca_HM_mean_v_disp_known.genes$pvals, brca_HM_mean_v_disp_known.genes$cor, 
     type="l", lwd=2, xlim=c(0.02,0.98), ylim=c(-1,0.8),
     col=col_vector[2], main=paste0("Spearman correlation between differential expression and ", 
                            "differential dispersion gene lists"))
lines(brca_MDSeq_mean_v_disp_known.genes$pvals, 
      brca_MDSeq_mean_v_disp_known.genes$cor, 
      lwd=2, col=col_vector[1])
legend("bottomright", fill=col_vector[1:2], bty='n', cex=1.5, ncol=1, 
       legend=c("MDSeq", "HM"))
plot(brca_HM_mean_v_disp_known.genes$pvals, brca_HM_mean_v_disp_known.genes$p.cor, 
     type="l", lwd=2, xlim=c(0.02,0.98), ylim=c(0,1),
     col=col_vector[2], main="p-value for hypothesis test for negative correlation")
lines(brca_MDSeq_mean_v_disp_known.genes$pvals, 
      brca_MDSeq_mean_v_disp_known.genes$p.cor, 
      lwd=2, col=col_vector[1])


## Cancer-related genes identified by DE, DD, DEDD ####
# Plot for all cancers, for known genes and pathway genes excluding known genes, for lnHM log DD 
# v voom, lnHM log DD v DE, lnHMM v voom, and MDSeq DD v DE
par(mfrow=c(8,8), mar=c(0.2,1,0.2,0), mgp=c(2,0,0))
for (i in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {
  for (j in c("HM", "HM_voom", "HMM_voom", "MDSeq")) {
    for (k in c("known.genes", "pathway.genes.minus.known.genes")) {
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
# Don't think there's any benefit in showing DE v DD for HM since it's clear that voom is better 
# at identifying DE. Also not really any need to show MDSeq, since results are similar and I've 
# already argued that MDSeq DD doesn't identify different sets of genes from DE as well as HM. 
# Would like to include HMM to keep up the discussion about differential distribution - 
# hopefully to show that differential distribution identifies more cancer-related genes than 
# either DE or DD alone, except that that isn't the case generally - in most cases it identifies 
# nearly exactly the same number of genes as DE, but for thca and prad HMM identifies far fewer 
# than DE, and in all cases there are still genes that are identified by DE but not by HMM (as 
# well as genes that are identified by HMM but not DE). For now just include lnHM log DD v voom, 
# and think about adding HMM later.

# Include luad and prad in main paper, rest in supplementary material.
# Want to show that DD identifies known and potential cancer-related genes, but maybe only 
# include known genes in main paper since there's not much difference between known and pathway 
# genes for luad and prad, and only include pathway genes in supplementary material.

par(mfcol=c(2,2), mar=c(2,2,2,0.5), mgp=c(3,0.7,0))
plot(luad_HM_voom_mean_v_disp_known.genes$pvals, luad_HM_voom_mean_v_disp_known.genes$TPA, 
     type="l", lwd=2, col="red", xlim=c(0.03,0.97))
lines(luad_HM_voom_mean_v_disp_known.genes$pvals, luad_HM_voom_mean_v_disp_known.genes$TPB, 
      col='blue', lwd=2)
lines(luad_HM_voom_mean_v_disp_known.genes$pvals, luad_HM_voom_mean_v_disp_known.genes$TPboth, 
      col='grey', lwd=2)
lines(luad_HM_voom_mean_v_disp_known.genes$pvals, luad_HM_voom_mean_v_disp_known.genes$TPAonly, 
      col='red', lty=2, lwd=2)
lines(luad_HM_voom_mean_v_disp_known.genes$pvals, luad_HM_voom_mean_v_disp_known.genes$TPBonly, 
      col='blue', lty=2, lwd=2)
plot(prad_HM_voom_mean_v_disp_known.genes$pvals, prad_HM_voom_mean_v_disp_known.genes$TPA, 
     type="l", lwd=2, col="red", xlim=c(0.03,0.97))
lines(prad_HM_voom_mean_v_disp_known.genes$pvals, prad_HM_voom_mean_v_disp_known.genes$TPB, 
      col='blue', lwd=2)
lines(prad_HM_voom_mean_v_disp_known.genes$pvals, prad_HM_voom_mean_v_disp_known.genes$TPboth, 
      col='grey', lwd=2)
lines(prad_HM_voom_mean_v_disp_known.genes$pvals, prad_HM_voom_mean_v_disp_known.genes$TPAonly, 
      col='red', lty=2, lwd=2)
lines(prad_HM_voom_mean_v_disp_known.genes$pvals, prad_HM_voom_mean_v_disp_known.genes$TPBonly, 
      col='blue', lty=2, lwd=2)
plot(luad_HM_voom_mean_v_disp_known.genes$pvals, luad_HM_voom_mean_v_disp_known.genes$TPA, 
     type="l", lwd=2, col="red", xlim=c(0.0002,0.01), ylim=c(0,230))
lines(luad_HM_voom_mean_v_disp_known.genes$pvals, luad_HM_voom_mean_v_disp_known.genes$TPB, 
      col='blue', lwd=2)
lines(luad_HM_voom_mean_v_disp_known.genes$pvals, luad_HM_voom_mean_v_disp_known.genes$TPboth, 
      col='grey', lwd=2)
lines(luad_HM_voom_mean_v_disp_known.genes$pvals, luad_HM_voom_mean_v_disp_known.genes$TPAonly, 
      col='red', lty=2, lwd=2)
lines(luad_HM_voom_mean_v_disp_known.genes$pvals, luad_HM_voom_mean_v_disp_known.genes$TPBonly, 
      col='blue', lty=2, lwd=2)
plot(prad_HM_voom_mean_v_disp_known.genes$pvals, prad_HM_voom_mean_v_disp_known.genes$TPA, 
     type="l", lwd=2, col="red", xlim=c(0.0002,0.01), ylim=c(0,550))
lines(prad_HM_voom_mean_v_disp_known.genes$pvals, prad_HM_voom_mean_v_disp_known.genes$TPB, 
      col='blue', lwd=2)
lines(prad_HM_voom_mean_v_disp_known.genes$pvals, prad_HM_voom_mean_v_disp_known.genes$TPboth, 
      col='grey', lwd=2)
lines(prad_HM_voom_mean_v_disp_known.genes$pvals, prad_HM_voom_mean_v_disp_known.genes$TPAonly, 
      col='red', lty=2, lwd=2)
lines(prad_HM_voom_mean_v_disp_known.genes$pvals, prad_HM_voom_mean_v_disp_known.genes$TPBonly, 
      col='blue', lty=2, lwd=2)
# Not sure whether there's any benefit to zooming in on lower p-values. Trend is clearly 
# discernible from looking at whole range. Think I should still show p-values up to 1, rather 
# than only zoomed in, to give full view of number of cancer-related genes present. If zoom in 
# as far as p-values of 0 to 0.01, clearly shows that voom identifies more cancer-related genes 
# than HM DD at start, but that is just because of minimum p-value (and actually would not be 
# the case if I made a different arbitrary decision and called the minimal value p=0 instead of 
# 0.0001 as I have done), so probably not really a disadvantage to show that level of detail. 
# However, still don't really get any more information from zoomed-in version than showing full 
# range, so think I'll stick with just showing full range.

par(mfcol=c(2,1), mar=c(2,2,2,0.5), mgp=c(3,0.7,0))
plot(luad_HM_voom_mean_v_disp_known.genes$pvals, luad_HM_voom_mean_v_disp_known.genes$TPA, 
     type="l", lwd=2, col="red", xlim=c(0.03,0.97), 
     main="Lung adenocarcinoma")
lines(luad_HM_voom_mean_v_disp_known.genes$pvals, luad_HM_voom_mean_v_disp_known.genes$TPB, 
      col='blue', lwd=2)
lines(luad_HM_voom_mean_v_disp_known.genes$pvals, luad_HM_voom_mean_v_disp_known.genes$TPboth, 
      col='grey', lwd=2)
lines(luad_HM_voom_mean_v_disp_known.genes$pvals, luad_HM_voom_mean_v_disp_known.genes$TPAonly, 
      col='red', lty=2, lwd=2)
lines(luad_HM_voom_mean_v_disp_known.genes$pvals, luad_HM_voom_mean_v_disp_known.genes$TPBonly, 
      col='blue', lty=2, lwd=2)
legend("right", 
       col=c("red", "blue", "grey", "red", "blue"), lty=c(1,1,1,2,2), 
       lwd=2, bty='n', cex=1.2, ncol=1, 
       legend=c("Differential expression", 
                "Differential dispersion", 
                "Both", 
                "Differential expression only", 
                "Differential dispersion only"))
plot(prad_HM_voom_mean_v_disp_known.genes$pvals, prad_HM_voom_mean_v_disp_known.genes$TPA, 
     type="l", lwd=2, col="red", xlim=c(0.03,0.97), 
     main="Prostate adenocarcinoma")
lines(prad_HM_voom_mean_v_disp_known.genes$pvals, prad_HM_voom_mean_v_disp_known.genes$TPB, 
      col='blue', lwd=2)
lines(prad_HM_voom_mean_v_disp_known.genes$pvals, prad_HM_voom_mean_v_disp_known.genes$TPboth, 
      col='grey', lwd=2)
lines(prad_HM_voom_mean_v_disp_known.genes$pvals, prad_HM_voom_mean_v_disp_known.genes$TPAonly, 
      col='red', lty=2, lwd=2)
lines(prad_HM_voom_mean_v_disp_known.genes$pvals, prad_HM_voom_mean_v_disp_known.genes$TPBonly, 
      col='blue', lty=2, lwd=2)


# Would be best to also (or instead) plot by number of discoveries, and could then plot FDR 
# or specificity (or something like them since we don't know true positives, just some sort of 
# proxy) rather than just number of cancer-related genes identified.
par(mfrow=c(4,8), mar=c(0.2,1,0.2,0), mgp=c(2,0,0))
for (i in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {
  for (j in c("HM_voom", "HMM_voom")) {
    for (k in c("known.genes", "pathway.genes.minus.known.genes")) {
      if (k == "known.genes") {
        l <- "known"
      } else {
        l <- "pathway"
      }
      dat <- get(paste0(i, "_", j, "_mean_v_disp_", k))
      plot(dat$listA, dat$TPA, type="l", col="red", xaxt="n", ann=F, tcl=0)
      lines(dat$listB, dat$TPB, col='blue')
      # lines(dat$pvals, dat$TPboth, col='grey')
      # find a way of combining results from two methods for given number of disc
      # lines(dat$pvals, dat$TPAonly, col='red', lty=2)
      # lines(dat$pvals, dat$TPBonly, col='blue', lty=2)
      legend('top', legend=paste0(i, " ", j, " ", l), bty="n")
    }
  }
}
rm(i,j,k,l)
# This isn't really any use until/unless I can get HM(M) to return finer p-values (or use 
# smaller samples). But it does show that over the range where HMs go from 0 to some large  
# number of positives, voom generally isn't doing much better than picking randomly from among 
# the genes with minimal p-values from HM. This might end up being a good thing to include once 
# I have finer p-values from HM(M).





