library(here)

## Import original and long chain pvals and integrate ####
for (i in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {
  assign(paste0("calls.all_methods.", i), 
         read.csv(here("Results/TCGA paired data results May 2020", 
                       paste0("calls.", i, ".csv")), 
                  stringsAsFactors=F, 
                  header=T, 
                  row.names=1))
  assign(paste0("pvals.all_methods.", i), 
         read.csv(here("Results/TCGA paired data results May 2020", 
                       paste0("pvals.", i, ".csv")), 
                  stringsAsFactors=F, 
                  header=T, 
                  row.names=1))
  assign(paste0("calls.10k.", i), 
         read.csv(here("Results/TCGA long chain results Aug 2020", 
                       paste0("calls.", i, ".csv")), 
                  stringsAsFactors=F, 
                  header=T, 
                  row.names=1))
  assign(paste0("pvals.10k.", i), 
         read.csv(here("Results/TCGA long chain results Aug 2020", 
                       paste0("pvals.", i, ".csv")), 
                  stringsAsFactors=F, 
                  header=T, 
                  row.names=1))
  assign(paste0("calls.20k.", i), 
         read.csv(here("Results/TCGA long chain results Aug 2020", 
                       paste0("calls.", i, "_20k.csv")), 
                  stringsAsFactors=F, 
                  header=T, 
                  row.names=1))
  assign(paste0("pvals.20k.", i), 
         read.csv(here("Results/TCGA long chain results Aug 2020", 
                       paste0("pvals.", i, "_20k.csv")), 
                  stringsAsFactors=F, 
                  header=T, 
                  row.names=1))
}

for (i in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {
  temp <- get(paste0("pvals.all_methods.", i))
  names(temp)[13:25] <- paste0(names(pvals.all_methods.brca), "_short")[13:25]
  assign(paste0("pvals.all_methods.", i), temp)
  temp <- get(paste0("pvals.10k.", i))
  names(temp) <- paste0(names(pvals.10k.brca), "_10k")
  assign(paste0("pvals.10k.", i), temp)
  temp <- get(paste0("pvals.20k.", i))
  names(temp) <- paste0(names(pvals.10k.brca), "_20k")
  assign(paste0("pvals.20k.", i), temp)
}

for (i in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {
  assign(paste0("pvals.", i), 
         cbind(
    get(paste0("pvals.all_methods.", i)), 
    get(paste0("pvals.10k.", i)), 
    get(paste0("pvals.20k.", i))
  ))
}

for (i in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {
  temp <- get(paste0("pvals.", i))
  names(temp) <- c(
    "edgeR.ql", "edgeR.lr", "edgeR.et", "DESeq2.noif", "DESeq2.if", "voom", "DSS", "baySeq", 
    "mean.MDSeq.zi", "mean.MDSeq.nozi", "disp.MDSeq.zi", "disp.MDSeq.nozi", 
    "expHMM.short", "mean.expHM.untr.short", "mean.expHM.log.short", "disp.expHM.untr.short", "disp.expHM.log.short", 
    "lnHMM.short", "mean.lnHM.untr.short", "mean.lnHM.log.short", "disp.lnHM.untr.short", "disp.lnHM.log.short", 
    "voom_lnHM.short", "voom_MDSeq.zi", "lnHM_lnHM.short", 
    "expHMM.10k", "mean.expHM.untr.10k", "mean.expHM.log.10k", "disp.expHM.untr.10k", "disp.expHM.log.10k", 
    "lnHMM.10k", "mean.lnHM.untr.10k", "mean.lnHM.log.10k", "disp.lnHM.untr.10k", "disp.lnHM.log.10k", "lnHM_lnHM.10k", 
    "expHMM.20k", "mean.expHM.untr.20k", "mean.expHM.log.20k", "disp.expHM.untr.20k", "disp.expHM.log.20k", 
    "lnHMM.20k", "mean.lnHM.untr.20k", "mean.lnHM.log.20k", "disp.lnHM.untr.20k", "disp.lnHM.log.20k", "lnHM_lnHM.20k"
  )
  assign(paste0("pvals.", i), temp)
}


## Import original and long chain calls and integrate ####
for (i in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {
  assign(paste0("calls.all_methods.", i), 
         read.csv(here("Results/TCGA paired data results May 2020", 
                       paste0("calls.", i, ".csv")), 
                  stringsAsFactors=F, 
                  header=T, 
                  row.names=1))
  assign(paste0("calls.all_methods.", i), 
         read.csv(here("Results/TCGA paired data results May 2020", 
                       paste0("calls.", i, ".csv")), 
                  stringsAsFactors=F, 
                  header=T, 
                  row.names=1))
  assign(paste0("calls.10k.", i), 
         read.csv(here("Results/TCGA long chain results Aug 2020", 
                       paste0("calls.", i, ".csv")), 
                  stringsAsFactors=F, 
                  header=T, 
                  row.names=1))
  assign(paste0("calls.10k.", i), 
         read.csv(here("Results/TCGA long chain results Aug 2020", 
                       paste0("calls.", i, ".csv")), 
                  stringsAsFactors=F, 
                  header=T, 
                  row.names=1))
  assign(paste0("calls.20k.", i), 
         read.csv(here("Results/TCGA long chain results Aug 2020", 
                       paste0("calls.", i, "_20k.csv")), 
                  stringsAsFactors=F, 
                  header=T, 
                  row.names=1))
  assign(paste0("calls.20k.", i), 
         read.csv(here("Results/TCGA long chain results Aug 2020", 
                       paste0("calls.", i, "_20k.csv")), 
                  stringsAsFactors=F, 
                  header=T, 
                  row.names=1))
}

for (i in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {
  temp <- get(paste0("calls.all_methods.", i))
  names(temp)[12:28] <- paste0(names(calls.all_methods.brca), "_short")[12:28]
  assign(paste0("calls.all_methods.", i), temp)
  temp <- get(paste0("calls.10k.", i))
  names(temp) <- paste0(names(calls.10k.brca), "_10k")
  assign(paste0("calls.10k.", i), temp)
  temp <- get(paste0("calls.20k.", i))
  names(temp) <- paste0(names(calls.10k.brca), "_20k")
  assign(paste0("calls.20k.", i), temp)
}

for (i in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {
  assign(paste0("calls.", i), 
         cbind(
           get(paste0("calls.all_methods.", i)), 
           get(paste0("calls.10k.", i)), 
           get(paste0("calls.20k.", i))
         ))
}

for (i in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {
  temp <- get(paste0("calls.", i))
  names(temp) <- c(
    "edgeR.ql", "edgeR.lr", "edgeR.et", "DESeq2.noif", "DESeq2.if", "voom", "baySeq", 
    "mean.MDSeq.zi", "mean.MDSeq.nozi", "disp.MDSeq.zi", "disp.MDSeq.nozi", 
    "expHMM.short_bfdr", "mean.expHM.untr.short", "mean.expHM.log.short", "disp.expHM.untr.short", "disp.expHM.log.short", 
    "lnHMM.short_bfdr", "mean.lnHM.untr.short", "mean.lnHM.log.short", "disp.lnHM.untr.short", "disp.lnHM.log.short", 
    "expHMM.short_thr", "expHMM.short_point5", "lnHMM.short_thr", "lnHMM.short_point5", 
    "voom_lnHM.short", "voom_MDSeq.zi", "lnHM_lnHM.short", 
    "expHMM.10k_bfdr", "mean.expHM.untr.10k", "mean.expHM.log.10k", "disp.expHM.untr.10k", "disp.expHM.log.10k", 
    "lnHMM.10k_bfdr", "mean.lnHM.untr.10k", "mean.lnHM.log.10k", "disp.lnHM.untr.10k", "disp.lnHM.log.10k", 
    "expHMM.10k_thr", "expHMM.10k_point5", "lnHMM.10k_thr", "lnHMM.10k_point5", "lnHM_lnHM.10k", 
    "expHMM.20k_bfdr", "mean.expHM.untr.20k", "mean.expHM.log.20k", "disp.expHM.untr.20k", "disp.expHM.log.20k", 
    "lnHMM.20k_bfdr", "mean.lnHM.untr.20k", "mean.lnHM.log.20k", "disp.lnHM.untr.20k", "disp.lnHM.log.20k", 
    "expHMM.20k_thr", "expHMM.20k_point5", "lnHMM.20k_thr", "lnHMM.20k_point5", "lnHM_lnHM.20k"
  )
  assign(paste0("calls.", i), temp)
  rm(temp)
}


## Create missing hybrid results for long chains (i.e. combinations with other methods, i.e. with voom) ####
{
  pvals.brca$voom_lnHM.10k <- pmin(pvals.brca$voom, pvals.brca$disp.lnHM.log.10k)
  pvals.kirc$voom_lnHM.10k <- pmin(pvals.kirc$voom, pvals.kirc$disp.lnHM.log.10k)
  pvals.thca$voom_lnHM.10k <- pmin(pvals.thca$voom, pvals.thca$disp.lnHM.log.10k)
  pvals.luad$voom_lnHM.10k <- pmin(pvals.luad$voom, pvals.luad$disp.lnHM.log.10k)
  pvals.lihc$voom_lnHM.10k <- pmin(pvals.lihc$voom, pvals.lihc$disp.lnHM.log.10k)
  pvals.lusc$voom_lnHM.10k <- pmin(pvals.lusc$voom, pvals.lusc$disp.lnHM.log.10k)
  pvals.prad$voom_lnHM.10k <- pmin(pvals.prad$voom, pvals.prad$disp.lnHM.log.10k)
  pvals.coad$voom_lnHM.10k <- pmin(pvals.coad$voom, pvals.coad$disp.lnHM.log.10k)
  calls.brca$voom_lnHM.10k <- p.adjust(pvals.brca$voom_lnHM.10k, method="BH") < 0.05
  calls.kirc$voom_lnHM.10k <- p.adjust(pvals.kirc$voom_lnHM.10k, method="BH") < 0.05
  calls.thca$voom_lnHM.10k <- p.adjust(pvals.thca$voom_lnHM.10k, method="BH") < 0.05
  calls.luad$voom_lnHM.10k <- p.adjust(pvals.luad$voom_lnHM.10k, method="BH") < 0.05
  calls.lihc$voom_lnHM.10k <- p.adjust(pvals.lihc$voom_lnHM.10k, method="BH") < 0.05
  calls.lusc$voom_lnHM.10k <- p.adjust(pvals.lusc$voom_lnHM.10k, method="BH") < 0.05
  calls.prad$voom_lnHM.10k <- p.adjust(pvals.prad$voom_lnHM.10k, method="BH") < 0.05
  calls.coad$voom_lnHM.10k <- p.adjust(pvals.coad$voom_lnHM.10k, method="BH") < 0.05
  pvals.brca$voom_lnHM.20k <- pmin(pvals.brca$voom, pvals.brca$disp.lnHM.log.20k)
  pvals.kirc$voom_lnHM.20k <- pmin(pvals.kirc$voom, pvals.kirc$disp.lnHM.log.20k)
  pvals.thca$voom_lnHM.20k <- pmin(pvals.thca$voom, pvals.thca$disp.lnHM.log.20k)
  pvals.luad$voom_lnHM.20k <- pmin(pvals.luad$voom, pvals.luad$disp.lnHM.log.20k)
  pvals.lihc$voom_lnHM.20k <- pmin(pvals.lihc$voom, pvals.lihc$disp.lnHM.log.20k)
  pvals.lusc$voom_lnHM.20k <- pmin(pvals.lusc$voom, pvals.lusc$disp.lnHM.log.20k)
  pvals.prad$voom_lnHM.20k <- pmin(pvals.prad$voom, pvals.prad$disp.lnHM.log.20k)
  pvals.coad$voom_lnHM.20k <- pmin(pvals.coad$voom, pvals.coad$disp.lnHM.log.20k)
  calls.brca$voom_lnHM.20k <- p.adjust(pvals.brca$voom_lnHM.20k, method="BH") < 0.05
  calls.kirc$voom_lnHM.20k <- p.adjust(pvals.kirc$voom_lnHM.20k, method="BH") < 0.05
  calls.thca$voom_lnHM.20k <- p.adjust(pvals.thca$voom_lnHM.20k, method="BH") < 0.05
  calls.luad$voom_lnHM.20k <- p.adjust(pvals.luad$voom_lnHM.20k, method="BH") < 0.05
  calls.lihc$voom_lnHM.20k <- p.adjust(pvals.lihc$voom_lnHM.20k, method="BH") < 0.05
  calls.lusc$voom_lnHM.20k <- p.adjust(pvals.lusc$voom_lnHM.20k, method="BH") < 0.05
  calls.prad$voom_lnHM.20k <- p.adjust(pvals.prad$voom_lnHM.20k, method="BH") < 0.05
  calls.coad$voom_lnHM.20k <- p.adjust(pvals.coad$voom_lnHM.20k, method="BH") < 0.05
}


## Save full results - original, long chains, longer chains, other methods, hybrids ####
for (i in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {
  write.csv(get(paste0("calls.", i)), 
            here("Results/TCGA all results Sept 2020", 
                 paste0("calls.", i, ".csv")), 
            row.names=row.names(get(paste0("pvals.", i))))
  write.csv(get(paste0("pvals.", i)), 
            here("Results/TCGA all results Sept 2020", 
                 paste0("pvals.", i, ".csv")))
}


## Import cancer-related gene data ####
folder <- "Data sources/Cancer-related genes"
for (i in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {
  assign(paste0("genes.", i), 
         readRDS(here(folder, paste0(i, "_genes_info.rds"))))
  assign(paste0(i, "_related"), 
         rownames(get(paste0("pvals.", i))) %in% get(paste0("genes.", i))$ENSEMBL)
}
rm(i, folder)


## Test ability to identify cancer-related genes with threshold using one-sided Fisher's exact test ####
tests.fisher <- names(calls.brca)
for (j in tests.fisher) {
  for (i in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {
    assign(paste0("fisher_", i, "_", j), 
           fisher.test(get(paste0(i, "_related")), 
                       get(paste0("calls.", i))[[j]], 
                       alternative="greater")$p.val)
  }
}
rm(i,j)

fisher.tests <- data.frame(matrix(nrow=length(tests.fisher), ncol=8))
names(fisher.tests) <- c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")
rownames(fisher.tests) <- tests.fisher
for (j in seq_along(rownames(fisher.tests))) {
  for (i in names(fisher.tests)) {
    fisher.tests[[i]][j] <- get(paste0("fisher_", i, "_", tests.fisher[j]))
    rm(list=paste0("fisher_", i, "_", tests.fisher[j]))
  }
}
rm(i,j,tests.fisher)

write.csv(fisher.tests, 
          here("Results/TCGA all results Sept 2020", 
               "fisher.tests_known.genes.csv"), 
)


## Test ability to rank cancer-related genes above non-cancer-related genes using Wilcoxon rank-sum test ####
tests.ranksum <- names(pvals.brca)
for (j in tests.ranksum) {
  for (i in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {
    if(grepl("baySeq", j) | grepl("HMM", j)) {
      assign(paste0("ranksum_", i, "_", j), 
             wilcox.test(
               rank(1 - get(paste0("pvals.", i))[[j]])[get(paste0(i, "_related")) == T], 
               rank(1 - get(paste0("pvals.", i))[[j]])[get(paste0(i, "_related")) == F], 
               paired=F, 
               alternative="less")$p.val)
    } else {
      assign(paste0("ranksum_", i, "_", j), 
             wilcox.test(
               rank(get(paste0("pvals.", i))[[j]])[get(paste0(i, "_related")) == T], 
               rank(get(paste0("pvals.", i))[[j]])[get(paste0(i, "_related")) == F], 
               paired=F, 
               alternative="less")$p.val)
    }
  }
  rm(i)
}
rm(j)

ranksum.tests <- data.frame(matrix(nrow=length(tests.ranksum), ncol=8))
names(ranksum.tests) <- c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")
rownames(ranksum.tests) <- tests.ranksum
for (j in seq_along(rownames(ranksum.tests))) {
  for (i in names(ranksum.tests)) {
    ranksum.tests[[i]][j] <- get(paste0("ranksum_", i, "_", tests.ranksum[j]))
    rm(list=paste0("ranksum_", i, "_", tests.ranksum[j]))
  }
}
rm(i,j,tests.ranksum)

write.csv(ranksum.tests, 
          here("Results/TCGA all results Sept 2020", 
               "ranksum.tests_known.genes.csv"), 
)




## Assess whether DE, DD, DEDD identify different sets of cancer-related genes ####
# Informally assess whether DE and DD HMs identify different sets of cancer-related genes by plotting 
# hypergeometric and Spearman correlation test p-values with alternative hypothesis that overlap and 
# correlation, respectively, are lower than expected under null for sets of top-ranking genes from 1 
# to all genes (or possibly for p-value threshold from 0 to 1 for hypergeometric).
for (j in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {

  # DE from HM, DD from HM; longest chains
  pvals <- sort(unique(c(get(paste0("pvals.", j))$mean.lnHM.log.20k, get(paste0("pvals.", j))$disp.lnHM.log.20k)))
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
    res$listDE[i] <- sum(get(paste0("pvals.", j))$mean.lnHM.log.20k < res$pvals[i])
    res$listDD[i] <- sum(get(paste0("pvals.", j))$disp.lnHM.log.20k < res$pvals[i])
    res$overlap[i] <- sum(get(paste0("pvals.", j))$mean.lnHM.log.20k < res$pvals[i] & 
                            get(paste0("pvals.", j))$disp.lnHM.log.20k < res$pvals[i])
    res$p.hyper[i] <- phyper(res$overlap[i], 
                             res$listDE[i], 
                             nrow(get(paste0("pvals.", j))) - res$listDE[i], 
                             res$listDD[i])
    res$exp.overlap[i] <- mean(get(paste0("pvals.", j))$mean.lnHM.log.20k < res$pvals[i]) * 
      mean(get(paste0("pvals.", j))$disp.lnHM.log.20k < res$pvals[i]) * nrow(get(paste0("pvals.", j)))
    res$union[i] <- sum(get(paste0("pvals.", j))$mean.lnHM.log.20k < res$pvals[i] | 
                          get(paste0("pvals.", j))$disp.lnHM.log.20k < res$pvals[i])
    if (sum(get(paste0("pvals.", j))$mean.lnHM.log.20k < res$pvals[i] | 
            get(paste0("pvals.", j))$disp.lnHM.log.20k < res$pvals[i]) > 1) {
      cortest <- 
        cor.test(get(paste0("pvals.", j))$mean.lnHM.log.20k[get(paste0("pvals.", j))$mean.lnHM.log.20k < res$pvals[i] | 
                                                            get(paste0("pvals.", j))$disp.lnHM.log.20k < res$pvals[i]], 
                 get(paste0("pvals.", j))$disp.lnHM.log.20k[get(paste0("pvals.", j))$mean.lnHM.log.20k < res$pvals[i] | 
                                                            get(paste0("pvals.", j))$disp.lnHM.log.20k < res$pvals[i]], 
                 method="spearman", alt="less")
      res$p.cor[i] <- cortest$p.val
      res$cor[i] <- cortest$estimate
      rm(cortest)
    }
    res$TP_DE[i] <- sum(get(paste0("pvals.", j))$mean.lnHM.log.20k[get(paste0(j, "_related"))] < res$pvals[i])
    res$TP_DD[i] <- sum(get(paste0("pvals.", j))$disp.lnHM.log.20k[get(paste0(j, "_related"))] < res$pvals[i])
    res$TP_both[i] <- sum(get(paste0("pvals.", j))$mean.lnHM.log.20k[get(paste0(j, "_related"))] < res$pvals[i] & 
                           get(paste0("pvals.", j))$disp.lnHM.log.20k[get(paste0(j, "_related"))] < res$pvals[i])
    res$TP_DEonly[i] <- sum(get(paste0("pvals.", j))$mean.lnHM.log.20k[get(paste0(j, "_related"))] < res$pvals[i] & 
                            get(paste0("pvals.", j))$disp.lnHM.log.20k[get(paste0(j, "_related"))] >= res$pvals[i])
    res$TP_DDonly[i] <- sum(get(paste0("pvals.", j))$mean.lnHM.log.20k[get(paste0(j, "_related"))] >= res$pvals[i] & 
                            get(paste0("pvals.", j))$disp.lnHM.log.20k[get(paste0(j, "_related"))] < res$pvals[i])
  }
  assign(paste0(j, ".HM_HM"), res)
  write.csv(get(paste0(j, ".HM_HM")), 
            here("Results/TCGA paired data results May 2020", 
                 paste0(j, "_mean_v_disp_HM.HM_known.genes.csv")), 
            row.names=F)
  rm(pvals, res)
  
  # DE from MDSeq, DD from MDSeq
  pvals <- sort(unique(c(get(paste0("pvals.", j))$mean.MDSeq.zi, get(paste0("pvals.", j))$disp.MDSeq.zi)))
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
    res$listDE[i] <- sum(get(paste0("pvals.", j))$mean.MDSeq.zi < res$pvals[i], na.rm=T)
    res$listDD[i] <- sum(get(paste0("pvals.", j))$disp.MDSeq.zi < res$pvals[i], na.rm=T)
    res$overlap[i] <- sum(get(paste0("pvals.", j))$mean.MDSeq.zi < res$pvals[i] & 
                            get(paste0("pvals.", j))$disp.MDSeq.zi < res$pvals[i], na.rm=T)
    res$union[i] <- sum(get(paste0("pvals.", j))$mean.MDSeq.zi < res$pvals[i] | 
                          get(paste0("pvals.", j))$disp.MDSeq.zi < res$pvals[i], na.rm=T)
    res$exp.overlap[i] <- mean(get(paste0("pvals.", j))$mean.MDSeq.zi < res$pvals[i], na.rm=T) * 
      mean(get(paste0("pvals.", j))$disp.MDSeq.zi < res$pvals[i], na.rm=T) * nrow(get(paste0("pvals.", j)))
    res$p.hyper[i] <- phyper(res$overlap[i], 
                             res$listDE[i], 
                             nrow(get(paste0("pvals.", j))) - res$listDE[i], 
                             res$listDD[i])
    if (sum(get(paste0("pvals.", j))$mean.MDSeq.zi < res$pvals[i] | 
            get(paste0("pvals.", j))$disp.MDSeq.zi < res$pvals[i], na.rm=T) > 1) {
      cortest <- cor.test(get(paste0("pvals.", j))$mean.MDSeq.zi[get(paste0("pvals.", j))$mean.MDSeq.zi < res$pvals[i] | 
                                                                     get(paste0("pvals.", j))$disp.MDSeq.zi < res$pvals[i]], 
                          get(paste0("pvals.", j))$disp.MDSeq.zi[get(paste0("pvals.", j))$mean.MDSeq.zi < res$pvals[i] | 
                                                                     get(paste0("pvals.", j))$disp.MDSeq.zi < res$pvals[i]], 
                          method="spearman", alt="less")
      res$p.cor[i] <- cortest$p.val
      res$cor[i] <- cortest$estimate
      rm(cortest)
    }
    res$TP_DE[i] <- sum(get(paste0("pvals.", j))$mean.MDSeq.zi[get(paste0(j, "_related"))] < res$pvals[i], na.rm=T)
    res$TP_DD[i] <- sum(get(paste0("pvals.", j))$disp.MDSeq.zi[get(paste0(j, "_related"))] < res$pvals[i], na.rm=T)
    res$TP_both[i] <- sum(get(paste0("pvals.", j))$mean.MDSeq.zi[get(paste0(j, "_related"))] < res$pvals[i] & 
                           get(paste0("pvals.", j))$disp.MDSeq.zi[get(paste0(j, "_related"))] < res$pvals[i], na.rm=T)
    res$TP_DEonly[i] <- sum(get(paste0("pvals.", j))$mean.MDSeq.zi[get(paste0(j, "_related"))] < res$pvals[i] & 
                            get(paste0("pvals.", j))$disp.MDSeq.zi[get(paste0(j, "_related"))] >= res$pvals[i], na.rm=T)
    res$TP_DDonly[i] <- sum(get(paste0("pvals.", j))$mean.MDSeq.zi[get(paste0(j, "_related"))] >= res$pvals[i] & 
                            get(paste0("pvals.", j))$disp.MDSeq.zi[get(paste0(j, "_related"))] < res$pvals[i], na.rm=T)
  }
  assign(paste0(j, ".MD_MD"), res)
  write.csv(get(paste0(j, ".MD_MD")), 
            here("Results/TCGA paired data results May 2020", 
                 paste0(j, "_mean_v_disp_MDSeq.MDSeq_minus.known.genes.csv")), 
            row.names=F)
  rm(pvals, res)
  
  # DE from voom, DD from HM; longest chains
  pvals <- sort(unique(c(get(paste0("pvals.", j))$voom, get(paste0("pvals.", j))$disp.lnHM.log.20k)))
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
  for (i in seq_len(nrow(res))) {
    res$listDE[i] <- sum(get(paste0("pvals.", j))$voom < res$pvals[i])
    res$listDD[i] <- sum(get(paste0("pvals.", j))$disp.lnHM.log.20k < res$pvals[i])
    res$overlap[i] <- sum(get(paste0("pvals.", j))$voom < res$pvals[i] & 
                            get(paste0("pvals.", j))$disp.lnHM.log.20k < res$pvals[i])
    res$p.hyper[i] <- phyper(res$overlap[i], 
                             res$listDE[i], 
                             nrow(get(paste0("pvals.", j))) - res$listDE[i], 
                             res$listDD[i])
    res$exp.overlap[i] <- mean(get(paste0("pvals.", j))$voom < res$pvals[i]) * 
      mean(get(paste0("pvals.", j))$disp.lnHM.log.20k < res$pvals[i]) * nrow(get(paste0("pvals.", j)))
    res$union[i] <- sum(get(paste0("pvals.", j))$voom < res$pvals[i] | 
                          get(paste0("pvals.", j))$disp.lnHM.log.20k < res$pvals[i])
    if (sum(get(paste0("pvals.", j))$voom < res$pvals[i] | 
            get(paste0("pvals.", j))$disp.lnHM.log.20k < res$pvals[i]) > 1) {
      cortest <- cor.test(get(paste0("pvals.", j))$voom[get(paste0("pvals.", j))$voom < res$pvals[i] | 
                                                            get(paste0("pvals.", j))$disp.lnHM.log.20k < res$pvals[i]], 
                          get(paste0("pvals.", j))$disp.lnHM.log.20k[get(paste0("pvals.", j))$voom < res$pvals[i] | 
                                                                     get(paste0("pvals.", j))$disp.lnHM.log.20k < res$pvals[i]], 
                          method="spearman", alt="less")
      res$p.cor[i] <- cortest$p.val
      res$cor[i] <- cortest$estimate
      rm(cortest)
    }
    res$TP_DE[i] <- sum(get(paste0("pvals.", j))$voom[get(paste0(j, "_related"))] < res$pvals[i])
    res$TP_DD[i] <- sum(get(paste0("pvals.", j))$disp.lnHM.log.20k[get(paste0(j, "_related"))] < res$pvals[i])
    res$TP_both[i] <- sum(get(paste0("pvals.", j))$voom[get(paste0(j, "_related"))] < res$pvals[i] & 
                           get(paste0("pvals.", j))$disp.lnHM.log.20k[get(paste0(j, "_related"))] < res$pvals[i])
    res$TP_DEonly[i] <- sum(get(paste0("pvals.", j))$voom[get(paste0(j, "_related"))] < res$pvals[i] & 
                            get(paste0("pvals.", j))$disp.lnHM.log.20k[get(paste0(j, "_related"))] >= res$pvals[i])
    res$TP_DDonly[i] <- sum(get(paste0("pvals.", j))$voom[get(paste0(j, "_related"))] >= res$pvals[i] & 
                            get(paste0("pvals.", j))$disp.lnHM.log.20k[get(paste0(j, "_related"))] < res$pvals[i])
  }
  assign(paste0(j, ".voom_HM"), res)
  write.csv(get(paste0(j, ".voom_HM")), 
            here("Results/TCGA paired data results May 2020", 
                 paste0(j, "_mean_v_disp_voom_HM.known.genes.csv")), 
            row.names=F)
  rm(pvals, res)
}


