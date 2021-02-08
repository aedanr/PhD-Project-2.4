# Experimenting with ways around limitation of having lots of equal minimum p-values for HM(M).
# Here, trying combining LFC with p-values by taking 10^-x where x is absolute value of LFC multiplied 
# by -log10 p-value, so measure used is p-value multiplied by 10^LFC.
# Another possibility is to use LFC purely as a tie-breaker: add a very small multiple of absolute LFC 
# to p-values - small enough that it won't change any ranking other than those with equal p-values.

library(here)

# Import p-value ranking and MDSeq results for comparison
folder <- "Results/TCGA all results Sept 2020"
for (i in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {
  for (j in c("HM.HM", "MDSeq.MDSeq", "voom.HM")) {
    assign(
      paste0(i, "_", j, "_mean_v_disp_known.genes"), 
      read.csv(
        here(folder, paste0(i, "_mean_v_disp_", j, "_known.genes.csv")), 
      )
    )
  }
}
rm(folder,i,j)

# Get p-values and calls
folder <- "Results/TCGA all results Sept 2020"
for (i in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {
  assign(paste0("calls.", i), 
         read.csv(here(folder, paste0("calls.", i, ".csv")), 
                  stringsAsFactors=F, header=T, row.names=1))
  assign(paste0("pvals.", i), 
         read.csv(here(folder, paste0("pvals.", i, ".csv")), 
                  stringsAsFactors=F, header=T, row.names=1))
}

# Get cancer-related genes
folder <- "Data sources/Cancer-related genes"
for (i in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {
  assign(paste0("genes.", i), 
         readRDS(here(folder, paste0(i, "_genes_info.rds"))))
  assign(paste0(i, "_related"), 
         rownames(get(paste0("pvals.", i))) %in% get(paste0("genes.", i))$ENSEMBL)
}
rm(i)

# Get LFCs for HM (lnHM long chain only; can do others if/when need to, but probably won't need to)
folder <- "Results/TCGA long chain results Aug 2020"
for (i in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {
  res <- readRDS(here(folder, paste0("results.TCGA_20k_iter_", i, "_lnHM.rds")))
  assign(paste0("lfc.lnHMmean_", i), res$lfc.mean)
  assign(paste0("lfc.lnHMdisp_", i), res$lfc.disp)
}
rm(i,res,folder)


## Rank genes by 10^-(|lfc| * (-log10(p))) = 10^-|lfc| * p for HM ####
# Keeps ranking variable in (0,1) but has potential to change rankings in a way that isn't easy to 
# predict, so almost certain to reverse some rankings compared to p-values alone.
# Only makes sense to do like-for-like comparisons, so if want to do this for voom-HM then need to 
# transform voom p-values in the same way. For now just do for HM because immediately interested in 
# comparing DE v DD for HM with DE v DD for MDSeq.
for (j in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {
  assign(paste0("p.explfc.", j, ".lnHMmean"), 
         get(paste0("pvals.", j))$mean.lnHM.log.20k * 10^-abs(get(paste0("lfc.lnHMmean_", j))))
  assign(paste0("p.explfc.", j, ".lnHMdisp"), 
         get(paste0("pvals.", j))$disp.lnHM.log.20k * 10^-abs(get(paste0("lfc.lnHMdisp_", j))))
  p.explfc <- sort(unique(c(get(paste0("p.explfc.", j, ".lnHMmean")), get(paste0("p.explfc.", j, ".lnHMdisp")))))
  res <- data.frame(
    p.explfc = p.explfc, 
    listDE = numeric(length(p.explfc)), 
    listDD = numeric(length(p.explfc)), 
    overlap = numeric(length(p.explfc)), 
    p.hyper = numeric(length(p.explfc)), 
    exp.overlap = numeric(length(p.explfc)), 
    union = numeric(length(p.explfc)), 
    p.cor = numeric(length(p.explfc)), 
    cor = numeric(length(p.explfc)), 
    TP_DE = numeric(length(p.explfc)), 
    TP_DD = numeric(length(p.explfc)), 
    TP_both = numeric(length(p.explfc)), 
    TP_DEonly = numeric(length(p.explfc)), 
    TP_DDonly = numeric(length(p.explfc))
  )
  for (i in seq_len(length(p.explfc))) {
    res$listDE[i] <- sum(get(paste0("p.explfc.", j, ".lnHMmean")) < res$p.explfc[i])
    res$listDD[i] <- sum(get(paste0("p.explfc.", j, ".lnHMdisp")) < res$p.explfc[i])
    res$overlap[i] <- sum(get(paste0("p.explfc.", j, ".lnHMmean")) < res$p.explfc[i] & 
                            get(paste0("p.explfc.", j, ".lnHMdisp")) < res$p.explfc[i])
    res$p.hyper[i] <- phyper(res$overlap[i], 
                             res$listDE[i], 
                             nrow(get(paste0("pvals.", j))) - res$listDE[i], 
                             res$listDD[i])
    res$exp.overlap[i] <- mean(get(paste0("p.explfc.", j, ".lnHMmean")) < res$p.explfc[i]) * 
      mean(get(paste0("p.explfc.", j, ".lnHMdisp")) < res$p.explfc[i]) * nrow(get(paste0("pvals.", j)))
    res$union[i] <- sum(get(paste0("p.explfc.", j, ".lnHMmean")) < res$p.explfc[i] | 
                          get(paste0("p.explfc.", j, ".lnHMdisp")) < res$p.explfc[i])
    if (sum(get(paste0("p.explfc.", j, ".lnHMmean")) < res$p.explfc[i] | 
            get(paste0("p.explfc.", j, ".lnHMdisp")) < res$p.explfc[i]) > 1) {
      cortest <- 
        cor.test(get(paste0("p.explfc.", j, ".lnHMmean"))[get(paste0("p.explfc.", j, ".lnHMmean")) < res$p.explfc[i] | 
                                                              get(paste0("p.explfc.", j, ".lnHMdisp")) < res$p.explfc[i]], 
                 get(paste0("p.explfc.", j, ".lnHMdisp"))[get(paste0("p.explfc.", j, ".lnHMmean")) < res$p.explfc[i] | 
                                                              get(paste0("p.explfc.", j, ".lnHMdisp")) < res$p.explfc[i]], 
                 method="spearman", alt="less")
      res$p.cor[i] <- cortest$p.val
      res$cor[i] <- cortest$estimate
      rm(cortest)
    }
    res$TP_DE[i] <- sum(get(paste0("p.explfc.", j, ".lnHMmean"))[get(paste0(j, "_related"))] < res$p.explfc[i])
    res$TP_DD[i] <- sum(get(paste0("p.explfc.", j, ".lnHMdisp"))[get(paste0(j, "_related"))] < res$p.explfc[i])
    res$TP_both[i] <- sum(get(paste0("p.explfc.", j, ".lnHMmean"))[get(paste0(j, "_related"))] < res$p.explfc[i] & 
                            get(paste0("p.explfc.", j, ".lnHMdisp"))[get(paste0(j, "_related"))] < res$p.explfc[i])
    res$TP_DEonly[i] <- sum(get(paste0("p.explfc.", j, ".lnHMmean"))[get(paste0(j, "_related"))] < res$p.explfc[i] & 
                              get(paste0("p.explfc.", j, ".lnHMdisp"))[get(paste0(j, "_related"))] >= res$p.explfc[i])
    res$TP_DDonly[i] <- sum(get(paste0("p.explfc.", j, ".lnHMmean"))[get(paste0(j, "_related"))] >= res$p.explfc[i] & 
                              get(paste0("p.explfc.", j, ".lnHMdisp"))[get(paste0(j, "_related"))] < res$p.explfc[i])
  }
  assign(paste0(j, ".HM_HM_p.explfc"), res)
  # write.csv(get(paste0(j, ".HM_HM_p.explfc")), 
  #           here("Results/TCGA all results Sept 2020", 
  #                paste0(j, "_mean_v_disp_HM.HM_known.genes_p.explfc.csv")), 
  #           row.names=F)
  rm(p.explfc, res)
}

# Compare HM using p-values (red) and p.explfc (blue)
par(mfrow=c(4,8), mar=c(2,2,2,0.5), mgp=c(3,0.7,0))
for (i in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {
  dat_pval <- get(paste0(i, "_HM.HM_mean_v_disp_known.genes"))
  dat_p.explfc <- get(paste0(i, ".HM_HM_p.explfc"))
  plot(dat_pval$pvals, 
       dat_pval$cor, 
       type="l", 
       xlim=c(0,1e-4),
       ylim=c(-1,1),
       col="red", main=i)
  lines(dat_p.explfc$p.explfc, 
        dat_p.explfc$cor, 
        col="blue")
}
# On full scale, correlation generally lower using p-values, especially near zero.
# Showing up to 1e-4, correlation using p.explfc is negative except extremely close to 
# zero where it doesn't exist using p-values.
for (i in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {
  dat_pval <- get(paste0(i, "_HM.HM_mean_v_disp_known.genes"))
  dat_p.explfc <- get(paste0(i, ".HM_HM_p.explfc"))
  plot(dat_pval$pvals, 
       dat_pval$p.cor, 
       type="l", 
       xlim=c(0,1e-4),
       ylim=c(0,1),
       col="red", main=i)
  lines(dat_p.explfc$p.explfc, 
        dat_p.explfc$p.cor, 
        col="blue")
}
# On full scale, correlation p-values close to zero for longer using p-values.
# Showing up to 1e-4, correlation p-values are at or very close to zero along the whole 
# scale using p.explfc, where they don't exist using p-values.
for (i in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {
  dat_pval <- get(paste0(i, "_HM.HM_mean_v_disp_known.genes"))
  dat_p.explfc <- get(paste0(i, ".HM_HM_p.explfc"))
  plot(dat_pval$pvals, 
       dat_pval$overlap / dat_pval$union, 
       type="l", 
       xlim=c(0,1e-4),
       ylim=c(0,1),
       col="red", main=i)
  lines(dat_p.explfc$p.explfc, 
        dat_p.explfc$overlap / dat_p.explfc$union, 
        col="blue")
}
# On full scale, less overlap using p-values.
# Showing up to 1e-4, overlap using p.explfc starts at zero but increases quite quickly, 
# and is always higher than using p-values from around 5e-5 where overlap can be defined 
# using p-values.
for (i in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {
  dat_pval <- get(paste0(i, "_HM.HM_mean_v_disp_known.genes"))
  dat_p.explfc <- get(paste0(i, ".HM_HM_p.explfc"))
  plot(dat_pval$pvals, 
       dat_pval$p.hyper, 
       type="l", 
       xlim=c(0,1e-4),
       ylim=c(0,1),
       col="red", main=i)
  lines(dat_p.explfc$p.explfc, 
        dat_p.explfc$p.hyper, 
        col="blue")
}
# On full scale, no consistent differences in hypergeometric p-values.
# Showing up to 1e-4, hypergeometric p-values are nearly always 1 across the range, but zero 
# for kirc and coad, for which they start at 1 and drop to zero at around 5e-5 using p-values.

# Compare MDSeq (red) with HM using p.explfc (blue)
par(mfrow=c(4,8), mar=c(2,2,2,0.5), mgp=c(3,0.7,0))
for (i in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {
  dat_pval <- get(paste0(i, "_MDSeq.MDSeq_mean_v_disp_known.genes"))
  dat_p.explfc <- get(paste0(i, ".HM_HM_p.explfc"))
  plot(dat_pval$pvals, 
       dat_pval$cor, 
       type="l", 
       xlim=c(0,1e-4),
       ylim=c(-1,1),
       col="red", main=i)
  lines(dat_p.explfc$p.explfc, 
        dat_p.explfc$cor, 
        col="blue")
}
# On full scale, correlation lower across whole range for HM except for thca, where they're 
# almost identical.
# Showing up to 1e-4, correlation still looks to be lower across the whole range for HM.
for (i in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {
  dat_pval <- get(paste0(i, "_MDSeq.MDSeq_mean_v_disp_known.genes"))
  dat_p.explfc <- get(paste0(i, ".HM_HM_p.explfc"))
  plot(dat_pval$pvals, 
       dat_pval$p.cor, 
       type="l", 
       xlim=c(0,1e-4),
       ylim=c(0,1),
       col="red", main=i)
  lines(dat_p.explfc$p.explfc, 
        dat_p.explfc$p.cor, 
        col="blue")
}
# On full scale, often no difference in correlation p-values, but lower for HM for small 
# p-values for brca and prad, and over the whole range for kirc.
# Showing up to 1e-4, correlation p-values lower across whole range for HM.
for (i in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {
  dat_pval <- get(paste0(i, "_MDSeq.MDSeq_mean_v_disp_known.genes"))
  dat_p.explfc <- get(paste0(i, ".HM_HM_p.explfc"))
  plot(dat_pval$pvals, 
       dat_pval$overlap / dat_pval$union, 
       type="l", 
       xlim=c(0,1e-4),
       ylim=c(0,1),
       col="red", main=i)
  lines(dat_p.explfc$p.explfc, 
        dat_p.explfc$overlap / dat_p.explfc$union, 
        col="blue")
}
# On full scale, overlap nearly always higher for HM.
# Showing up to 1e-4, not much difference in overlap, but where there's a difference it's 
# always lower for HM.
for (i in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {
  dat_pval <- get(paste0(i, "_MDSeq.MDSeq_mean_v_disp_known.genes"))
  dat_p.explfc <- get(paste0(i, ".HM_HM_p.explfc"))
  plot(dat_pval$pvals, 
       dat_pval$p.hyper, 
       type="l", 
       xlim=c(0,1e-4),
       ylim=c(0,1),
       col="red", main=i)
  lines(dat_p.explfc$p.explfc, 
        dat_p.explfc$p.hyper, 
        col="blue")
}
# On full scale, hypergeometric p-values look to always be lower for HM, which is very 
# strange given there's nearly always greater overlap between DE and DD lists for HM.
# Showing up to 1e-4, for most cancers there's no difference in hypergeometric p-values, 
# but zero across the whole range for HM for kirc and coad, where it's 1 for MDSeq, and 
# zero very close to zero for HM for brca where it's 1 for MDSeq.

# Compare HM using p-values (red) and p.explfc (blue) based on number of genes
par(mfrow=c(4,8), mar=c(2,2,2,0.5), mgp=c(3,0.7,0))
for (i in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {
  dat_pval <- get(paste0(i, "_HM.HM_mean_v_disp_known.genes"))
  dat_p.explfc <- get(paste0(i, ".HM_HM_p.explfc"))
  plot(dat_pval$union, 
       dat_pval$cor, 
       type="l", 
       xlim=c(0,1000),
       ylim=c(-1,1),
       col="red", main=i)
  lines(dat_p.explfc$union, 
        dat_p.explfc$cor, 
        col="blue")
}
# On full scale, correlation nearly always negative using p.explfc where it doesn't 
# exist using p-values, and crosses over to higher at a point where the ranking doesn't 
# really matter any more (at least 5000 genes).
# Showing up to 1000 genes, correlation is generally positive at start using p.explfc then 
# drops below zero between about 100 and 600 genes.
for (i in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {
  dat_pval <- get(paste0(i, "_HM.HM_mean_v_disp_known.genes"))
  dat_p.explfc <- get(paste0(i, ".HM_HM_p.explfc"))
  plot(dat_pval$union, 
       dat_pval$p.cor, 
       type="l", 
       xlim=c(0,1000),
       ylim=c(0,1),
       col="red", main=i)
  lines(dat_p.explfc$union, 
        dat_p.explfc$p.cor, 
        col="blue")
}
# On full scale, correlation p-values look to be high using p.explfc at the very low end 
# and for thca and lihc at the high end, and where there's a difference, they look to be 
# always higher using p.explfc than p-values, but that's only at the extreme high end.
# Showing up to 1000 genes, correlation p-values are always higher using p.explfc over 
# the entire low end except a very very small bit at the start, which could just be one or 
# two genes, where both are zero.
for (i in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {
  dat_pval <- get(paste0(i, "_HM.HM_mean_v_disp_known.genes"))
  dat_p.explfc <- get(paste0(i, ".HM_HM_p.explfc"))
  plot(dat_pval$union, 
       dat_pval$overlap / dat_pval$union, 
       type="l", 
       xlim=c(0,1000),
       ylim=c(0,1),
       col="red", main=i)
  lines(dat_p.explfc$union, 
        dat_p.explfc$overlap / dat_p.explfc$union, 
        col="blue")
}
# On full scale, overlap starts off low using p.explfc where it isn't defined using p-values, 
# and where it is defined using p-values they're very similar, but that's only for at least 
# about 10,000 genes.
# Showing up to 1000 genes, overlap is generally very low using p.explfc, and at or nearly at 
# zero at the start, except for thca. This is probably because there are no discoveries for 
# one test (probably DD since DE seems to go to much smaller values).
for (i in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {
  dat_pval <- get(paste0(i, "_HM.HM_mean_v_disp_known.genes"))
  dat_p.explfc <- get(paste0(i, ".HM_HM_p.explfc"))
  plot(dat_pval$union, 
       dat_pval$p.hyper, 
       type="l", 
       xlim=c(0,1000),
       ylim=c(0,1),
       col="red", main=i)
  lines(dat_p.explfc$union, 
        dat_p.explfc$p.hyper, 
        col="blue")
}
# On full scale, hypergeometric p-values look to start at 1 using both p.explfc and p-values, 
# and for most cancers stay around there, but for kirc and coad, both drop to zero but on a 
# finer scale using p.explfc, and for brca there is a point where the hypergeometric p-value 
# drops to nearly zero using p.explfc but stays at 1 using p-values, but there's no difference 
# for reasonable numbers of genes.
# Showing up to 1000 genes, hypergeometric p-value is 1 over the whole range except for tiny 
# bits for a couple of cancers. It's very strange that the p-value is so close to 1 when the 
# overlap/union is so small. It could just be the small numbers, but the numbers aren't really 
# small for up to 1000 genes, so need to check whether I'm doing the test properly. Looking at 
# all the numbers for some over the range where the overlap/union is low but p-value 1, it 
# looks like it's just because for these relatively small numbers of genes, the expected overlap 
# is extremely small, so looking at overlap/union alone doesn't really say anything about 
# whether there's less overlap than expected.

# Compare MDSeq (red) with HM using p.explfc (blue) based on number of genes
par(mfrow=c(4,8), mar=c(2,2,2,0.5), mgp=c(3,0.7,0))
for (i in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {
  dat_pval <- get(paste0(i, "_MDSeq.MDSeq_mean_v_disp_known.genes"))
  dat_p.explfc <- get(paste0(i, ".HM_HM_p.explfc"))
  plot(dat_pval$union, 
       dat_pval$cor, 
       type="l", 
       xlim=c(0,1000),
       ylim=c(-1,1),
       col="red", main=i)
  lines(dat_p.explfc$union, 
        dat_p.explfc$cor, 
        col="blue")
}
# On full scale, correlation lower across most of range for HM for most cancers, but for 
# most looks like it's lower for MDSeq at the very low end.
# Showing up to 1000 genes, correlation is lower at low end for all genes for MDSeq, and 
# generally similar or lower for HM after around 300-600 genes, except for luad and prad, for 
# which it's lower for MDSeq across the whole range. MDSeq seems to have a higher minimum 
# number of genes identified than HM when using p.explfc, so correlation jumps from zero to 
# something negative at around 100-400 genes for most cancers, whereas for HM it starts off 
# around zero or slightly positive, and more slowly drops off. Strangely though, looking at 
# kirc as an example, there is a correlation and a correlation p-value for MDSeq at a point 
# where there is no overlap between the DE and DD lists. This is because the correlation I've 
# used is between all genes in the union. I wonder if it makes more sense to use that or all 
# genes in the intersection (which would mean it's not defined for a longer period at the 
# beginning).
for (i in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {
  dat_pval <- get(paste0(i, "_MDSeq.MDSeq_mean_v_disp_known.genes"))
  dat_p.explfc <- get(paste0(i, ".HM_HM_p.explfc"))
  plot(dat_pval$union, 
       dat_pval$p.cor, 
       type="l", 
       xlim=c(0,1000),
       ylim=c(0,1),
       col="red", main=i)
  lines(dat_p.explfc$union, 
        dat_p.explfc$p.cor, 
        col="blue")
}
# On full scale, correlation p-values look to be always the same or lower for HM, but 
# possibly stays at zero at the very low end for MDSeq where it's high for HM.
# Showing up to 1000 genes, the correlation p-value is always lower for MDSeq at the low end, 
# up to from 100 to 600 genes.
for (i in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {
  dat_pval <- get(paste0(i, "_MDSeq.MDSeq_mean_v_disp_known.genes"))
  dat_p.explfc <- get(paste0(i, ".HM_HM_p.explfc"))
  plot(dat_pval$union, 
       dat_pval$overlap / dat_pval$union, 
       type="l", 
       xlim=c(0,1000),
       ylim=c(0,1),
       col="red", main=i)
  lines(dat_p.explfc$union, 
        dat_p.explfc$overlap / dat_p.explfc$union, 
        col="blue")
}
# On full scale, overlap generally lower for HM across most of range, but looks to 
# usually be lower near the low end for MDSeq, but low for both.
# Showing up to 1000 genes, overlap is generally similar for HM and MDSeq, but where it's 
# different it's usually lower for MDSeq.
for (i in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {
  dat_pval <- get(paste0(i, "_MDSeq.MDSeq_mean_v_disp_known.genes"))
  dat_p.explfc <- get(paste0(i, ".HM_HM_p.explfc"))
  plot(dat_pval$union, 
       dat_pval$p.hyper, 
       type="l", 
       xlim=c(0,1000),
       ylim=c(0,1),
       col="red", main=i)
  lines(dat_p.explfc$union, 
        dat_p.explfc$p.hyper, 
        col="blue")
}
# On full scale, hypergeometric p-values usually same, but more often lower for HM than 
# MDSeq, except at very low end seems to often be lower for MDSeq, but still close to 1.
# Showing up to 1000 genes, hypergeometric p-values are at or very close to 1 for both HM 
# and MDSeq across the whole range.


## Rank genes by p-value with lfc as tie-breaker for HM ####
# Easiest way to do it is to subtract a multiple of absolute lfc from p-values such that it's 
# guaranteed not to change the ranking of any pairs of genes with different p-values. That 
# means the factor subtracted needs to be guaranteed to be less than 1/21,000. Dividing 
# absolute lfc by 1e6 should do it but gives some negative numbers, so use 1e7.
for (j in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {
  assign(paste0("modp.", j, ".lnHMmean"), 
         get(paste0("pvals.", j))$mean.lnHM.log.20k - abs(get(paste0("lfc.lnHMmean_", j)))/1e7)
  assign(paste0("modp.", j, ".lnHMdisp"), 
         get(paste0("pvals.", j))$disp.lnHM.log.20k - abs(get(paste0("lfc.lnHMdisp_", j)))/1e7)
  modp <- sort(unique(c(get(paste0("modp.", j, ".lnHMmean")), get(paste0("modp.", j, ".lnHMdisp")))))
  res <- data.frame(
    modp = modp, 
    listDE = numeric(length(modp)), 
    listDD = numeric(length(modp)), 
    overlap = numeric(length(modp)), 
    p.hyper = numeric(length(modp)), 
    exp.overlap = numeric(length(modp)), 
    union = numeric(length(modp)), 
    p.cor = numeric(length(modp)), 
    cor = numeric(length(modp)), 
    TP_DE = numeric(length(modp)), 
    TP_DD = numeric(length(modp)), 
    TP_both = numeric(length(modp)), 
    TP_DEonly = numeric(length(modp)), 
    TP_DDonly = numeric(length(modp))
  )
  for (i in seq_len(length(modp))) {
    res$listDE[i] <- sum(get(paste0("modp.", j, ".lnHMmean")) < res$modp[i])
    res$listDD[i] <- sum(get(paste0("modp.", j, ".lnHMdisp")) < res$modp[i])
    res$overlap[i] <- sum(get(paste0("modp.", j, ".lnHMmean")) < res$modp[i] & 
                            get(paste0("modp.", j, ".lnHMdisp")) < res$modp[i])
    res$p.hyper[i] <- phyper(res$overlap[i], 
                             res$listDE[i], 
                             nrow(get(paste0("pvals.", j))) - res$listDE[i], 
                             res$listDD[i])
    res$exp.overlap[i] <- mean(get(paste0("modp.", j, ".lnHMmean")) < res$modp[i]) * 
      mean(get(paste0("modp.", j, ".lnHMdisp")) < res$modp[i]) * nrow(get(paste0("pvals.", j)))
    res$union[i] <- sum(get(paste0("modp.", j, ".lnHMmean")) < res$modp[i] | 
                          get(paste0("modp.", j, ".lnHMdisp")) < res$modp[i])
    if (sum(get(paste0("modp.", j, ".lnHMmean")) < res$modp[i] | 
            get(paste0("modp.", j, ".lnHMdisp")) < res$modp[i]) > 1) {
      cortest <- 
        cor.test(get(paste0("modp.", j, ".lnHMmean"))[get(paste0("modp.", j, ".lnHMmean")) < res$modp[i] | 
                                                            get(paste0("modp.", j, ".lnHMdisp")) < res$modp[i]], 
                 get(paste0("modp.", j, ".lnHMdisp"))[get(paste0("modp.", j, ".lnHMmean")) < res$modp[i] | 
                                                            get(paste0("modp.", j, ".lnHMdisp")) < res$modp[i]], 
                 method="spearman", alt="less")
      res$p.cor[i] <- cortest$p.val
      res$cor[i] <- cortest$estimate
      rm(cortest)
    }
    res$TP_DE[i] <- sum(get(paste0("modp.", j, ".lnHMmean"))[get(paste0(j, "_related"))] < res$modp[i])
    res$TP_DD[i] <- sum(get(paste0("modp.", j, ".lnHMdisp"))[get(paste0(j, "_related"))] < res$modp[i])
    res$TP_both[i] <- sum(get(paste0("modp.", j, ".lnHMmean"))[get(paste0(j, "_related"))] < res$modp[i] & 
                            get(paste0("modp.", j, ".lnHMdisp"))[get(paste0(j, "_related"))] < res$modp[i])
    res$TP_DEonly[i] <- sum(get(paste0("modp.", j, ".lnHMmean"))[get(paste0(j, "_related"))] < res$modp[i] & 
                              get(paste0("modp.", j, ".lnHMdisp"))[get(paste0(j, "_related"))] >= res$modp[i])
    res$TP_DDonly[i] <- sum(get(paste0("modp.", j, ".lnHMmean"))[get(paste0(j, "_related"))] >= res$modp[i] & 
                              get(paste0("modp.", j, ".lnHMdisp"))[get(paste0(j, "_related"))] < res$modp[i])
  }
  assign(paste0(j, ".HM_HM_modp"), res)
  # write.csv(get(paste0(j, ".HM_HM_modp")), 
  #           here("Results/TCGA all results Sept 2020", 
  #                paste0(j, "_mean_v_disp_HM.HM_known.genes_lfc.tiebreaker.csv")), 
  #           row.names=F)
  rm(modp, res)
}

# Compare HM using p-values (red) and modp (blue)
par(mfrow=c(4,8), mar=c(2,2,2,0.5), mgp=c(3,0.7,0))
for (i in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {
  dat_pval <- get(paste0(i, "_HM.HM_mean_v_disp_known.genes"))
  dat_modp <- get(paste0(i, ".HM_HM_modp"))
  plot(dat_pval$pvals, 
       dat_pval$cor, 
       type="l", 
       xlim=c(0,1),
       ylim=c(-1,1),
       col="red", main=i)
  lines(dat_modp$modp, 
        dat_modp$cor, 
        col="blue")
}
# On full scale, very similar using modp and p-values.
# Showing up to 1e-4, for most cancers, using modp correlation is low very close to the 
# start but higher than using p-values after about 5e-5.
for (i in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {
  dat_pval <- get(paste0(i, "_HM.HM_mean_v_disp_known.genes"))
  dat_modp <- get(paste0(i, ".HM_HM_modp"))
  plot(dat_pval$pvals, 
       dat_pval$p.cor, 
       type="l", 
       xlim=c(0,1),
       ylim=c(0,1),
       col="red", main=i)
  lines(dat_modp$modp, 
        dat_modp$p.cor, 
        col="blue")
}
# On full scale, correlation p-values generally stay close to zero for longer using p-values, 
# but at or very close to zero using both methods for low p-values.
# Showing up to 1e-4, correlation p-values using modp are zero nearly from the very start, 
# and can't see lines using p-values, so either same as using modp or don't exist.
for (i in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {
  dat_pval <- get(paste0(i, "_HM.HM_mean_v_disp_known.genes"))
  dat_modp <- get(paste0(i, ".HM_HM_modp"))
  plot(dat_pval$pvals, 
       dat_pval$overlap / dat_pval$union, 
       type="l", 
       xlim=c(0,1),
       ylim=c(0,1),
       col="red", main=i)
  lines(dat_modp$modp, 
        dat_modp$overlap / dat_modp$union, 
        col="blue")
}
# On full scale, overlap appears identical using modp and p-values.
# Showing up to 1e-4, overlap starts at zero using modp and quickly increases to same 
# point where it starts using p-values, and then almost exactly same.
for (i in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {
  dat_pval <- get(paste0(i, "_HM.HM_mean_v_disp_known.genes"))
  dat_modp <- get(paste0(i, ".HM_HM_modp"))
  plot(dat_pval$pvals, 
       dat_pval$p.hyper, 
       type="l", 
       xlim=c(0,1),
       ylim=c(0,1),
       col="red", main=i)
  lines(dat_modp$modp, 
        dat_modp$p.hyper, 
        col="blue")
}
# On full scale, hypergeometric p-values appear identical using modp and p-values.
# Showing up to 1e-4, hypergeometric p-values nearly always 1 using modp except for kirc 
# and coad, where they start at 1 and drop to zero immediately, whereas using p-values it 
# also drops from 1 to zero, but isn't defined between zero and about 5e-5.

# Compare MDSeq (red) with HM using modp (blue)
par(mfrow=c(4,8), mar=c(2,2,2,0.5), mgp=c(3,0.7,0))
for (i in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {
  dat_pval <- get(paste0(i, "_MDSeq.MDSeq_mean_v_disp_known.genes"))
  dat_modp <- get(paste0(i, ".HM_HM_modp"))
  plot(dat_pval$pvals, 
       dat_pval$cor, 
       type="l", 
       xlim=c(0,1),
       ylim=c(-1,1),
       col="red", main=i)
  lines(dat_modp$modp, 
        dat_modp$cor, 
        col="blue")
}
# On full scale, correlation always lower for HM.
# Showing up to 1e-4, where the correlation is defined, it's always lower for HM except at 
# the very beginning.
for (i in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {
  dat_pval <- get(paste0(i, "_MDSeq.MDSeq_mean_v_disp_known.genes"))
  dat_modp <- get(paste0(i, ".HM_HM_modp"))
  plot(dat_pval$pvals, 
       dat_pval$p.cor, 
       type="l", 
       xlim=c(0,1),
       ylim=c(0,1),
       col="red", main=i)
  lines(dat_modp$modp, 
        dat_modp$p.cor, 
        col="blue")
}
# On full scale, correlation p-values always same or lower for HM.
# Showing up to 1e-4, correlation p-values are nearly always lower for HM from where they 
# exist, starting at 1 and dropping immediately to zero.
for (i in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {
  dat_pval <- get(paste0(i, "_MDSeq.MDSeq_mean_v_disp_known.genes"))
  dat_modp <- get(paste0(i, ".HM_HM_modp"))
  plot(dat_pval$pvals, 
       dat_pval$overlap / dat_pval$union, 
       type="l", 
       xlim=c(0,1),
       ylim=c(0,1),
       col="red", main=i)
  lines(dat_modp$modp, 
        dat_modp$overlap / dat_modp$union, 
        col="blue")
}
# On full scale, overlap very similar but where it's different it's lower for HM
# Showing up to 1e-4, overlap is always lower for HM where it's defined, except at the very 
# start for thca.
for (i in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {
  dat_pval <- get(paste0(i, "_MDSeq.MDSeq_mean_v_disp_known.genes"))
  dat_modp <- get(paste0(i, ".HM_HM_modp"))
  plot(dat_pval$pvals, 
       dat_pval$p.hyper, 
       type="l", 
       xlim=c(0,1),
       ylim=c(0,1),
       col="red", main=i)
  lines(dat_modp$modp, 
        dat_modp$p.hyper, 
        col="blue")
}
# On full scale, hypergeometric p-values always lower for HM.
# Showing up to 1e-4, hypergeometric p-values are generally the same for HM and MDSeq, but 
# for kirc and coad, they're 1 across nearly the whole range for MDSeq and zero for HM.

# Compare HM using p-values (red) and modp (blue) based on number of genes
par(mfrow=c(4,8), mar=c(2,2,2,0.5), mgp=c(3,0.7,0))
for (i in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {
  dat_pval <- get(paste0(i, "_HM.HM_mean_v_disp_known.genes"))
  dat_modp <- get(paste0(i, ".HM_HM_modp"))
  plot(dat_pval$union, 
       dat_pval$cor, 
       type="l", 
       xlim=c(0,20000),
       ylim=c(-1,1),
       col="red", main=i)
  lines(dat_modp$union, 
        dat_modp$cor, 
        col="blue")
}
# On full scale, correlation using modp has positive peaks at the very start, but over most 
# of the range is lower than using p-values, and negative.
# Showing up to 1000 genes, correlation generally starts higher using modp, around zero or 
# positive, but drops to negative and below line using p-values before about 500 genes.
for (i in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {
  dat_pval <- get(paste0(i, "_HM.HM_mean_v_disp_known.genes"))
  dat_modp <- get(paste0(i, ".HM_HM_modp"))
  plot(dat_pval$union, 
       dat_pval$p.cor, 
       type="l", 
       xlim=c(0,20000),
       ylim=c(0,1),
       col="red", main=i)
  lines(dat_modp$union, 
        dat_modp$p.cor, 
        col="blue")
}
# On full scale, correlation p-values identical and at or very close to zero across most of 
# range using modp and p-values.
# Showing up to 1000 genes, correlation p-value starts of high using modp where it's on 
# a straight line from (0,0) using p-values. Drops to zero using modp between about 100 
# and 600 genes.
for (i in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {
  dat_pval <- get(paste0(i, "_HM.HM_mean_v_disp_known.genes"))
  dat_modp <- get(paste0(i, ".HM_HM_modp"))
  plot(dat_pval$union, 
       dat_pval$overlap / dat_pval$union, 
       type="l", 
       xlim=c(0,20000),
       ylim=c(0,1),
       col="red", main=i)
  lines(dat_modp$union, 
        dat_modp$overlap / dat_modp$union, 
        col="blue")
}
# On full scale, overlap identical using modp and p-values.
# Showing up to 1000 genes, overlap is low across the whole range using modp, where it's 
# not defined using p-values.
for (i in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {
  dat_pval <- get(paste0(i, "_HM.HM_mean_v_disp_known.genes"))
  dat_modp <- get(paste0(i, ".HM_HM_modp"))
  plot(dat_pval$union, 
       dat_pval$p.hyper, 
       type="l", 
       xlim=c(0,20000),
       ylim=c(0,1),
       col="red", main=i)
  lines(dat_modp$union, 
        dat_modp$p.hyper, 
        col="blue")
}
# On full scale, hypergeometric p-values usually the same using modp and p-values, but 
# sometimes lower using modp.
# Showing up to 1000 genes, hypergeometric p-value is 1 across nearly whole range using 
# modp, same as using p-values for most genes.

# Compare MDSeq (red) with HM using modp (blue) based on number of genes
par(mfrow=c(4,8), mar=c(2,2,2,0.5), mgp=c(3,0.7,0))
for (i in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {
  dat_pval <- get(paste0(i, "_MDSeq.MDSeq_mean_v_disp_known.genes"))
  dat_modp <- get(paste0(i, ".HM_HM_modp"))
  plot(dat_pval$union, 
       dat_pval$cor, 
       type="l", 
       xlim=c(0,20000),
       ylim=c(-1,1),
       col="red", main=i)
  lines(dat_modp$union, 
        dat_modp$cor, 
        col="blue")
}
# On full scale, correlations usually the same or lower for HM across most of the range, 
# but often lower for MDSeq near the start.
# Showing up to 1000 genes, correlation generally starts higher for HM than MDSeq but 
# for most genes lower for HM by 1000 genes, and generally similar. It's never lower 
# for HM than MDSeq at the start, but similar (and nearly always negative) for thca, 
# thca, lihc and coad.
for (i in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {
  dat_pval <- get(paste0(i, "_MDSeq.MDSeq_mean_v_disp_known.genes"))
  dat_modp <- get(paste0(i, ".HM_HM_modp"))
  plot(dat_pval$union, 
       dat_pval$p.cor, 
       type="l", 
       xlim=c(0,20000),
       ylim=c(0,1),
       col="red", main=i)
  lines(dat_modp$union, 
        dat_modp$p.cor, 
        col="blue")
}
# On full scale, correlation p-values nraly always lower for HM.
# Showing up to 1000 genes, correlation p-value is always lower for MDSeq than HM at 
# the start (up to 100-600 genes), then generally zero for both, but for brca and lusc 
# it goes down to zero for HM and increases for MDSeq.
for (i in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {
  dat_pval <- get(paste0(i, "_MDSeq.MDSeq_mean_v_disp_known.genes"))
  dat_modp <- get(paste0(i, ".HM_HM_modp"))
  plot(dat_pval$union, 
       dat_pval$overlap / dat_pval$union, 
       type="l", 
       xlim=c(0,20000),
       ylim=c(0,1),
       col="red", main=i)
  lines(dat_modp$union, 
        dat_modp$overlap / dat_modp$union, 
        col="blue")
}
# On full scale, overlap nearly always lower for HM, but for some cancers lower for 
# MDSeq near the start.
# Showing up to 1000 genes, overlap very similar for HM and MDSeq but where it's 
# different, generally lower for MDSeq.
for (i in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {
  dat_pval <- get(paste0(i, "_MDSeq.MDSeq_mean_v_disp_known.genes"))
  dat_modp <- get(paste0(i, ".HM_HM_modp"))
  plot(dat_pval$union, 
       dat_pval$p.hyper, 
       type="l", 
       xlim=c(0,20000),
       ylim=c(0,1),
       col="red", main=i)
  lines(dat_modp$union, 
        dat_modp$p.hyper, 
        col="blue")
}
# On full scale, hypergeometric p-values generally the same for MDSeq and HM, but where 
# they're different they're usually lower for HM.
# Showing up to 1000 genes, hypergeometric p-values very similar, and at or very close 
# to 1 across the whole range for both.


## HM DD v voom DE with lfc as tiebreaker for HM ####
for (j in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {
  assign(paste0("modp.", j, ".lnHMdisp"), 
         get(paste0("pvals.", j))$disp.lnHM.log.20k - abs(get(paste0("lfc.lnHMdisp_", j)))/1e7)
  pvals <- sort(unique(c(get(paste0("pvals.", j))$voom, get(paste0("modp.", j, ".lnHMdisp")))))
  res <- data.frame(
    pvals = pvals, 
    listDE = numeric(length(pvals)), 
    listDD = numeric(length(pvals)), 
    overlap = numeric(length(pvals)), 
    p.hyper = numeric(length(pvals)), 
    exp.overlap = numeric(length(pvals)), 
    union = numeric(length(pvals)), 
    p.cor = numeric(length(pvals)), 
    cor = numeric(length(pvals)), 
    TP_DE = numeric(length(pvals)), 
    TP_DD = numeric(length(pvals)), 
    TP_both = numeric(length(pvals)), 
    TP_DEonly = numeric(length(pvals)), 
    TP_DDonly = numeric(length(pvals))
  )
  for (i in seq_len(length(pvals))) {
    res$listDE[i] <- sum(get(paste0("pvals.", j))$voom < res$pvals[i])
    res$listDD[i] <- sum(get(paste0("modp.", j, ".lnHMdisp")) < res$pvals[i])
    res$overlap[i] <- sum(get(paste0("pvals.", j))$voom < res$pvals[i] & 
                            get(paste0("modp.", j, ".lnHMdisp")) < res$pvals[i])
    res$p.hyper[i] <- phyper(res$overlap[i], 
                             res$listDE[i], 
                             nrow(get(paste0("pvals.", j))) - res$listDE[i], 
                             res$listDD[i])
    res$exp.overlap[i] <- mean(get(paste0("pvals.", j))$voom < res$pvals[i]) * 
      mean(get(paste0("modp.", j, ".lnHMdisp")) < res$pvals[i]) * nrow(get(paste0("pvals.", j)))
    res$union[i] <- sum(get(paste0("pvals.", j))$voom < res$pvals[i] | 
                          get(paste0("modp.", j, ".lnHMdisp")) < res$pvals[i])
    if (sum(get(paste0("pvals.", j))$voom < res$pvals[i] | 
            get(paste0("modp.", j, ".lnHMdisp")) < res$pvals[i]) > 1) {
      cortest <- 
        cor.test(get(paste0("pvals.", j))$voom[get(paste0("pvals.", j))$voom < res$pvals[i] | 
                                                        get(paste0("modp.", j, ".lnHMdisp")) < res$pvals[i]], 
                 get(paste0("modp.", j, ".lnHMdisp"))[get(paste0("pvals.", j))$voom < res$pvals[i] | 
                                                        get(paste0("modp.", j, ".lnHMdisp")) < res$pvals[i]], 
                 method="spearman", alt="less")
      res$p.cor[i] <- cortest$p.val
      res$cor[i] <- cortest$estimate
      rm(cortest)
    }
    res$TP_DE[i] <- sum(get(paste0("pvals.", j))$voom[get(paste0(j, "_related"))] < res$pvals[i])
    res$TP_DD[i] <- sum(get(paste0("modp.", j, ".lnHMdisp"))[get(paste0(j, "_related"))] < res$pvals[i])
    res$TP_both[i] <- sum(get(paste0("pvals.", j))$voom[get(paste0(j, "_related"))] < res$pvals[i] & 
                            get(paste0("modp.", j, ".lnHMdisp"))[get(paste0(j, "_related"))] < res$pvals[i])
    res$TP_DEonly[i] <- sum(get(paste0("pvals.", j))$voom[get(paste0(j, "_related"))] < res$pvals[i] & 
                              get(paste0("modp.", j, ".lnHMdisp"))[get(paste0(j, "_related"))] >= res$pvals[i])
    res$TP_DDonly[i] <- sum(get(paste0("pvals.", j))$voom[get(paste0(j, "_related"))] >= res$pvals[i] & 
                              get(paste0("modp.", j, ".lnHMdisp"))[get(paste0(j, "_related"))] < res$pvals[i])
  }
  assign(paste0(j, ".voom_HM_lfc.tiebreaker"), res)
  # write.csv(get(paste0(j, ".voom_HM_lfc.tiebreaker")), 
  #           here("Results/TCGA all results Sept 2020", 
  #                paste0(j, "_mean_v_disp_voom.HM_known.genes_lfc.tiebreaker.csv")), 
  #           row.names=F)
  rm(pvals, res)
}



## Compare true positives between voom and lnHM disp using lfc as tiebreaker ####
# Not really any way of comparing unique TPs for each method on the same scale unless the 
# p-values from each method are on the same scale, since when they're very different, if 
# using the union or overlap as the common scale, the method with smaller p-values  idenitifies 
# hundreds genes before the method with larger p-values identifies any.
par(mfrow=c(4,4), mar=c(0.2,1,0.2,0), mgp=c(2,0,0))
# TP from DE (red), DD (blue) by modified p-value
for (i in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {
  dat <- get(paste0(i, ".voom_HM_lfc.tiebreaker"))
  plot(dat$pvals, dat$TP_DE, type="l", col="red", xaxt="n", ann=F, tcl=0)
  lines(dat$pvals, dat$TP_DD, col='blue')
  lines(dat$pvals, dat$TP_both, col='grey')
  lines(dat$pvals, dat$TP_DEonly, col='red', lty=2)
  lines(dat$pvals, dat$TP_DDonly, col='blue', lty=2)
  legend('top', legend=i, bty="n")
}
# TP from DE (red), DD (blue) by number of discoveries
for (i in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {
  dat <- get(paste0(i, ".voom_HM_lfc.tiebreaker"))
  plot(dat$listDE, dat$TP_DE, type="l", col="red", xaxt="n", ann=F, tcl=0)
  lines(dat$listDD, dat$TP_DD, col='blue')
  legend('top', legend=i, bty="n")
}
rm(i,dat)
# Look pretty similar to previous results on full scale, except low end by number of 
# discoveries, where there now isn't a big jump for HM. There are no cases where DD 
# identifies more TPs than DE consistently, but for thca it looks like DD identifies 
# more at the very low end and the two methods are very similar all the way along, and 
# for luad and lusc they're quite similar all the way along, and for lihc at the low end.



