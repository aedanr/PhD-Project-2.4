library(here)

for (i in c("brca", "kirc", "thca", "luad", "lihc", 
            "lihc_long.chain", "lusc", "prad", "coad")) {
  assign(paste0("calls.", i), 
         read.csv(here("Results/TCGA paired data results May 2020", 
                       paste0("calls.", i, ".csv")), 
                  stringsAsFactors=F, 
                  header=T, 
                  row.names=1))
  assign(paste0("pvals.", i), 
         read.csv(here("Results/TCGA paired data results May 2020", 
                       paste0("pvals.", i, ".csv")), 
                  stringsAsFactors=F, 
                  header=T, 
                  row.names=1))
}

folder <- "Data sources/Cancer-related pathway genes"
for (i in c("brca", "kirc", "thca", "luad", "lihc", 
            "lusc", "prad", "coad")) {
  assign(paste0("genes.", i), 
         read.csv(here(folder, paste0(i, "_pathway_genes_recount.csv")), stringsAsFactors=F))
  assign(paste0("dup.genes.", i), 
         read.csv(here(folder, paste0(i, "_pathway_genes_recount_multiple_pathways.csv")), stringsAsFactors=F))
  assign(paste0(i, "_related_weak"), 
         rownames(get(paste0("pvals.", i))) %in% get(paste0("genes.", i))$ENSEMBL)
  assign(paste0(i, "_related_strong"), 
         rownames(get(paste0("pvals.", i))) %in% get(paste0("dup.genes.", i))$ENSEMBL)
}
rm(i, folder)
lihc_long.chain_related_weak <- lihc_related_weak
lihc_long.chain_related_strong <- lihc_related_strong


## Test ability to identify cancer-related genes with threshold using one-sided Fisher's exact test ####
tests.fisher <- names(calls.brca)
for (j in tests.fisher) {
  for (i in c("brca", "kirc", "thca", "luad", "lihc", 
              "lusc", "prad", "coad")) {
    assign(paste0("fisher.weak_", i, "_", j), 
           fisher.test(get(paste0(i, "_related_weak")), 
                       get(paste0("calls.", i))[[j]], 
                       alternative="greater")$p.val)
    assign(paste0("fisher.strong_", i, "_", j), 
           fisher.test(get(paste0(i, "_related_strong")), 
                       get(paste0("calls.", i))[[j]], 
                       alternative="greater")$p.val)
  }
  assign(paste0("fisher.weak_lihc_long.chain_", j), 
         fisher.test(lihc_related_weak, 
                     calls.lihc_long.chain[[j]], 
                     alternative="greater")$p.val)
  assign(paste0("fisher.strong_lihc_long.chain_", j), 
         fisher.test(lihc_related_strong, 
                     calls.lihc_long.chain[[j]], 
                     alternative="greater")$p.val)
}
rm(i,j)

fisher.tests.weak <- data.frame(matrix(nrow=length(tests.fisher), ncol=9))
names(fisher.tests.weak) <- c("brca", "kirc", "thca", "luad", "lihc", 
                              "lihc_long.chain",  "lusc", "prad", "coad")
rownames(fisher.tests.weak) <- tests.fisher
for (j in seq_along(rownames(fisher.tests.weak))) {
  for (i in names(fisher.tests.weak)) {
    fisher.tests.weak[[i]][j] <- get(paste0("fisher.weak_", i, "_", tests.fisher[j]))
    rm(list=paste0("fisher.weak_", i, "_", tests.fisher[j]))
  }
}

fisher.tests.strong <- data.frame(matrix(nrow=length(tests.fisher), ncol=9))
names(fisher.tests.strong) <- c("brca", "kirc", "thca", "luad", "lihc", 
                                "lihc_long.chain",  "lusc", "prad", "coad")
rownames(fisher.tests.strong) <- tests.fisher
for (j in seq_along(rownames(fisher.tests.strong))) {
  for (i in names(fisher.tests.strong)) {
    fisher.tests.strong[[i]][j] <- get(paste0("fisher.strong_", i, "_", tests.fisher[j]))
    rm(list=paste0("fisher.strong_", i, "_", tests.fisher[j]))
  }
}
rm(i,j,tests.fisher)

write.csv(fisher.tests.weak, 
          here("Results/TCGA paired data results May 2020", 
               "fisher.tests.weak_pathway.genes.csv"), 
)
write.csv(fisher.tests.strong, 
          here("Results/TCGA paired data results May 2020", 
               "fisher.tests.strong_pathway.genes.csv"), 
)


## Test ability to rank cancer-related genes above non-cancer-related genes using Wilcoxon rank-sum test ####
tests.ranksum <- names(pvals.brca)
for (j in tests.ranksum) {
  for (i in c("brca", "kirc", "thca", "luad", "lihc", 
              "lusc", "prad", "coad")) {
    if(grepl("prob", j)) {
      assign(paste0("ranksum.weak_", i, "_", j), 
             wilcox.test(
               rank(1 - get(paste0("pvals.", i))[[j]])[get(paste0(i, "_related_weak")) == T], 
               rank(1 - get(paste0("pvals.", i))[[j]])[get(paste0(i, "_related_weak")) == F], 
               paired=F, 
               alternative="less")$p.val)
      assign(paste0("ranksum.strong_", i, "_", j), 
             wilcox.test(
               rank(1 - get(paste0("pvals.", i))[[j]])[get(paste0(i, "_related_strong")) == T], 
               rank(1 - get(paste0("pvals.", i))[[j]])[get(paste0(i, "_related_strong")) == F], 
               paired=F, 
               alternative="less")$p.val)
    } else {
      assign(paste0("ranksum.weak_", i, "_", j), 
             wilcox.test(
               rank(get(paste0("pvals.", i))[[j]])[get(paste0(i, "_related_weak")) == T], 
               rank(get(paste0("pvals.", i))[[j]])[get(paste0(i, "_related_weak")) == F], 
               paired=F, 
               alternative="less")$p.val)
      assign(paste0("ranksum.strong_", i, "_", j), 
             wilcox.test(
               rank(get(paste0("pvals.", i))[[j]])[get(paste0(i, "_related_strong")) == T], 
               rank(get(paste0("pvals.", i))[[j]])[get(paste0(i, "_related_strong")) == F], 
               paired=F, 
               alternative="less")$p.val)
    }
  }
  rm(i)
  if (grepl("prob", j)) {
    assign(paste0("ranksum.weak_lihc_long.chain_", j), 
           wilcox.test(
             rank(1 - pvals.lihc_long.chain[[j]])[lihc_related_weak == T], 
             rank(1 - pvals.lihc_long.chain[[j]])[lihc_related_weak == F], 
             paired=F, 
             alternative="less")$p.val)
    assign(paste0("ranksum.strong_lihc_long.chain_", j), 
           wilcox.test(
             rank(1 - pvals.lihc_long.chain[[j]])[lihc_related_strong == T], 
             rank(1 - pvals.lihc_long.chain[[j]])[lihc_related_strong == F], 
             paired=F, 
             alternative="less")$p.val)
  } else {
    assign(paste0("ranksum.weak_lihc_long.chain_", j), 
           wilcox.test(
             rank(pvals.lihc_long.chain[[j]])[lihc_related_weak == T], 
             rank(pvals.lihc_long.chain[[j]])[lihc_related_weak == F], 
             paired=F, 
             alternative="less")$p.val)
    assign(paste0("ranksum.strong_lihc_long.chain_", j), 
           wilcox.test(
             rank(pvals.lihc_long.chain[[j]])[lihc_related_strong == T], 
             rank(pvals.lihc_long.chain[[j]])[lihc_related_strong == F], 
             paired=F, 
             alternative="less")$p.val)
  }
}
rm(j)

ranksum.tests.weak <- data.frame(matrix(nrow=length(tests.ranksum), ncol=9))
names(ranksum.tests.weak) <- c("brca", "kirc", "thca", "luad", "lihc", 
                               "lihc_long.chain",  "lusc", "prad", "coad")
rownames(ranksum.tests.weak) <- tests.ranksum
for (j in seq_along(rownames(ranksum.tests.weak))) {
  for (i in names(ranksum.tests.weak)) {
    ranksum.tests.weak[[i]][j] <- get(paste0("ranksum.weak_", i, "_", tests.ranksum[j]))
    rm(list=paste0("ranksum.weak_", i, "_", tests.ranksum[j]))
  }
}

ranksum.tests.strong <- data.frame(matrix(nrow=length(tests.ranksum), ncol=9))
names(ranksum.tests.strong) <- c("brca", "kirc", "thca", "luad", "lihc", 
                                 "lihc_long.chain",  "lusc", "prad", "coad")
rownames(ranksum.tests.strong) <- tests.ranksum
for (j in seq_along(rownames(ranksum.tests.strong))) {
  for (i in names(ranksum.tests.strong)) {
    ranksum.tests.strong[[i]][j] <- get(paste0("ranksum.strong_", i, "_", tests.ranksum[j]))
    rm(list=paste0("ranksum.strong_", i, "_", tests.ranksum[j]))
  }
}
rm(i,j,tests.ranksum)

write.csv(ranksum.tests.weak, 
          here("Results/TCGA paired data results May 2020", 
               "ranksum.tests.weak_pathway.genes.csv"), 
)
write.csv(ranksum.tests.strong, 
          here("Results/TCGA paired data results May 2020", 
               "ranksum.tests.strong_pathway.genes.csv"), 
)


## Compare methods for ability to identify cancer-related genes with threshold using McNemar's test ####
c(mcnemar.test(calls.brca$q.voom == brca_related_weak, 
               calls.brca$q.mean.lnHM.log == brca_related_weak)$p.val, 
  mean(calls.brca$q.voom == brca_related_weak), 
  mean(calls.brca$q.mean.lnHM.log == brca_related_weak)) # voom better
c(mcnemar.test(calls.kirc$q.voom == kirc_related_weak, 
               calls.kirc$q.mean.lnHM.log == kirc_related_weak)$p.val, 
  mean(calls.kirc$q.voom == kirc_related_weak), 
  mean(calls.kirc$q.mean.lnHM.log == kirc_related_weak)) # voom better
c(mcnemar.test(calls.thca$q.voom == thca_related_weak, 
               calls.thca$q.mean.lnHM.log == thca_related_weak)$p.val, 
  mean(calls.thca$q.voom == thca_related_weak), 
  mean(calls.thca$q.mean.lnHM.log == thca_related_weak)) # HM better
c(mcnemar.test(calls.luad$q.voom == luad_related_weak, 
               calls.luad$q.mean.lnHM.log == luad_related_weak)$p.val, 
  mean(calls.luad$q.voom == luad_related_weak), 
  mean(calls.luad$q.mean.lnHM.log == luad_related_weak)) # voom better
c(mcnemar.test(calls.lihc$q.voom == lihc_related_weak, 
               calls.lihc$q.mean.lnHM.log == lihc_related_weak)$p.val, 
  mean(calls.lihc$q.voom == lihc_related_weak), 
  mean(calls.lihc$q.mean.lnHM.log == lihc_related_weak)) # voom better
c(mcnemar.test(calls.lihc_long.chain$q.voom == lihc_related_weak, 
               calls.lihc_long.chain$q.mean.lnHM.log == lihc_related_weak)$p.val, 
  mean(calls.lihc_long.chain$q.voom == lihc_related_weak), 
  mean(calls.lihc_long.chain$q.mean.lnHM.log == lihc_related_weak)) # voom better
c(mcnemar.test(calls.lusc$q.voom == lusc_related_weak, 
               calls.lusc$q.mean.lnHM.log == lusc_related_weak)$p.val, 
  mean(calls.lusc$q.voom == lusc_related_weak), 
  mean(calls.lusc$q.mean.lnHM.log == lusc_related_weak)) # voom better
c(mcnemar.test(calls.prad$q.voom == prad_related_weak, 
               calls.prad$q.mean.lnHM.log == prad_related_weak)$p.val, 
  mean(calls.prad$q.voom == prad_related_weak), 
  mean(calls.prad$q.mean.lnHM.log == prad_related_weak)) # voom better
c(mcnemar.test(calls.coad$q.voom == coad_related_weak, 
               calls.coad$q.mean.lnHM.log == coad_related_weak)$p.val, 
  mean(calls.coad$q.voom == coad_related_weak), 
  mean(calls.coad$q.mean.lnHM.log == coad_related_weak)) # voom better

c(mcnemar.test(calls.brca$q.voom == brca_related_weak, 
               calls.brca$q.disp.lnHM.log == brca_related_weak)$p.val, 
  mean(calls.brca$q.voom == brca_related_weak), 
  mean(calls.brca$q.disp.lnHM.log == brca_related_weak)) # voom better
c(mcnemar.test(calls.kirc$q.voom == kirc_related_weak, 
               calls.kirc$q.disp.lnHM.log == kirc_related_weak)$p.val, 
  mean(calls.kirc$q.voom == kirc_related_weak), 
  mean(calls.kirc$q.disp.lnHM.log == kirc_related_weak)) # HM better
c(mcnemar.test(calls.thca$q.voom == thca_related_weak, 
               calls.thca$q.disp.lnHM.log == thca_related_weak)$p.val, 
  mean(calls.thca$q.voom == thca_related_weak), 
  mean(calls.thca$q.disp.lnHM.log == thca_related_weak)) # HM better
c(mcnemar.test(calls.luad$q.voom == luad_related_weak, 
               calls.luad$q.disp.lnHM.log == luad_related_weak)$p.val, 
  mean(calls.luad$q.voom == luad_related_weak), 
  mean(calls.luad$q.disp.lnHM.log == luad_related_weak)) # voom better
c(mcnemar.test(calls.lihc$q.voom == lihc_related_weak, 
               calls.lihc$q.disp.lnHM.log == lihc_related_weak)$p.val, 
  mean(calls.lihc$q.voom == lihc_related_weak), 
  mean(calls.lihc$q.disp.lnHM.log == lihc_related_weak)) # voom better
c(mcnemar.test(calls.lihc_long.chain$q.voom == lihc_related_weak, 
               calls.lihc_long.chain$q.disp.lnHM.log == lihc_related_weak)$p.val, 
  mean(calls.lihc_long.chain$q.voom == lihc_related_weak), 
  mean(calls.lihc_long.chain$q.disp.lnHM.log == lihc_related_weak)) # voom better
c(mcnemar.test(calls.lusc$q.voom == lusc_related_weak, 
               calls.lusc$q.disp.lnHM.log == lusc_related_weak)$p.val, 
  mean(calls.lusc$q.voom == lusc_related_weak), 
  mean(calls.lusc$q.disp.lnHM.log == lusc_related_weak)) # voom better
c(mcnemar.test(calls.prad$q.voom == prad_related_weak, 
               calls.prad$q.disp.lnHM.log == prad_related_weak)$p.val, 
  mean(calls.prad$q.voom == prad_related_weak), 
  mean(calls.prad$q.disp.lnHM.log == prad_related_weak)) # HM better
c(mcnemar.test(calls.coad$q.voom == coad_related_weak, 
               calls.coad$q.disp.lnHM.log == coad_related_weak)$p.val, 
  mean(calls.coad$q.voom == coad_related_weak), 
  mean(calls.coad$q.disp.lnHM.log == coad_related_weak)) # voom better

c(mcnemar.test(calls.brca$q.disp.MDSeq.zi == brca_related_weak, 
               calls.brca$q.disp.lnHM.log == brca_related_weak)$p.val, 
  mean(calls.brca$q.disp.MDSeq.zi == brca_related_weak, na.rm=T), 
  mean(calls.brca$q.disp.lnHM.log == brca_related_weak)) # MDSeq better
c(mcnemar.test(calls.kirc$q.disp.MDSeq.zi == kirc_related_weak, 
               calls.kirc$q.disp.lnHM.log == kirc_related_weak)$p.val, 
  mean(calls.kirc$q.disp.MDSeq.zi == kirc_related_weak, na.rm=T), 
  mean(calls.kirc$q.disp.lnHM.log == kirc_related_weak)) # MDSeq better
c(mcnemar.test(calls.thca$q.disp.MDSeq.zi == thca_related_weak, 
               calls.thca$q.disp.lnHM.log == thca_related_weak)$p.val, 
  mean(calls.thca$q.disp.MDSeq.zi == thca_related_weak, na.rm=T), 
  mean(calls.thca$q.disp.lnHM.log == thca_related_weak)) # MDSeq better
c(mcnemar.test(calls.luad$q.disp.MDSeq.zi == luad_related_weak, 
               calls.luad$q.disp.lnHM.log == luad_related_weak)$p.val, 
  mean(calls.luad$q.disp.MDSeq.zi == luad_related_weak, na.rm=T), 
  mean(calls.luad$q.disp.lnHM.log == luad_related_weak)) # MDSeq better
c(mcnemar.test(calls.lihc$q.disp.MDSeq.zi == lihc_related_weak, 
               calls.lihc$q.disp.lnHM.log == lihc_related_weak)$p.val, 
  mean(calls.lihc$q.disp.MDSeq.zi == lihc_related_weak, na.rm=T), 
  mean(calls.lihc$q.disp.lnHM.log == lihc_related_weak)) # MDSeq better
c(mcnemar.test(calls.lihc_long.chain$q.disp.MDSeq.zi == lihc_related_weak, 
               calls.lihc_long.chain$q.disp.lnHM.log == lihc_related_weak)$p.val, 
  mean(calls.lihc_long.chain$q.disp.MDSeq.zi == lihc_related_weak, na.rm=T), 
  mean(calls.lihc_long.chain$q.disp.lnHM.log == lihc_related_weak)) # MDSeq better
c(mcnemar.test(calls.lusc$q.disp.MDSeq.zi == lusc_related_weak, 
               calls.lusc$q.disp.lnHM.log == lusc_related_weak)$p.val, 
  mean(calls.lusc$q.disp.MDSeq.zi == lusc_related_weak, na.rm=T), 
  mean(calls.lusc$q.disp.lnHM.log == lusc_related_weak)) # MDSeq better
c(mcnemar.test(calls.prad$q.disp.MDSeq.zi == prad_related_weak, 
               calls.prad$q.disp.lnHM.log == prad_related_weak)$p.val, 
  mean(calls.prad$q.disp.MDSeq.zi == prad_related_weak, na.rm=T), 
  mean(calls.prad$q.disp.lnHM.log == prad_related_weak)) # No difference
c(mcnemar.test(calls.coad$q.disp.MDSeq.zi == coad_related_weak, 
               calls.coad$q.disp.lnHM.log == coad_related_weak)$p.val, 
  mean(calls.coad$q.disp.MDSeq.zi == coad_related_weak, na.rm=T), 
  mean(calls.coad$q.disp.lnHM.log == coad_related_weak)) # MDSeq better


## Compare methods on ranking of cancer-related genes using Wilcoxon signed-rank test ####
c(wilcox.test(rank(pvals.brca$p.voom)[brca_related_weak == T], 
              rank(pvals.brca$p.mean.lnHM.log)[brca_related_weak == T], paired=T)$p.val, 
  median((rank(pvals.brca$p.voom) - rank(pvals.brca$p.mean.lnHM.log))[brca_related_weak == T])) # voom better
c(wilcox.test(rank(pvals.kirc$p.voom)[kirc_related_weak == T], 
              rank(pvals.kirc$p.mean.lnHM.log)[kirc_related_weak == T], paired=T)$p.val, 
  median((rank(pvals.kirc$p.voom) - rank(pvals.kirc$p.mean.lnHM.log))[kirc_related_weak == T])) # voom better
c(wilcox.test(rank(pvals.thca$p.voom)[thca_related_weak == T], 
              rank(pvals.thca$p.mean.lnHM.log)[thca_related_weak == T], paired=T)$p.val, 
  median((rank(pvals.thca$p.voom) - rank(pvals.thca$p.mean.lnHM.log))[thca_related_weak == T])) # voom better
c(wilcox.test(rank(pvals.luad$p.voom)[luad_related_weak == T], 
              rank(pvals.luad$p.mean.lnHM.log)[luad_related_weak == T], paired=T)$p.val, 
  median((rank(pvals.luad$p.voom) - rank(pvals.luad$p.mean.lnHM.log))[luad_related_weak == T])) # voom better
c(wilcox.test(rank(pvals.lihc$p.voom)[lihc_related_weak == T], 
              rank(pvals.lihc$p.mean.lnHM.log)[lihc_related_weak == T], paired=T)$p.val, 
  median((rank(pvals.lihc$p.voom) - rank(pvals.lihc$p.mean.lnHM.log))[lihc_related_weak == T])) # voom better
c(wilcox.test(rank(pvals.lihc_long.chain$p.voom)[lihc_related_weak == T], 
              rank(pvals.lihc_long.chain$p.mean.lnHM.log)[lihc_related_weak == T], paired=T)$p.val, 
  median((rank(pvals.lihc_long.chain$p.voom) - 
            rank(pvals.lihc_long.chain$p.mean.lnHM.log))[lihc_related_weak == T])) # voom better
c(wilcox.test(rank(pvals.lusc$p.voom)[lusc_related_weak == T], 
              rank(pvals.lusc$p.mean.lnHM.log)[lusc_related_weak == T], paired=T)$p.val, 
  median((rank(pvals.lusc$p.voom) - rank(pvals.lusc$p.mean.lnHM.log))[lusc_related_weak == T])) # voom better
c(wilcox.test(rank(pvals.prad$p.voom)[prad_related_weak == T], 
              rank(pvals.prad$p.mean.lnHM.log)[prad_related_weak == T], paired=T)$p.val, 
  median((rank(pvals.prad$p.voom) - rank(pvals.prad$p.mean.lnHM.log))[prad_related_weak == T])) # voom better
c(wilcox.test(rank(pvals.coad$p.voom)[coad_related_weak == T], 
              rank(pvals.coad$p.mean.lnHM.log)[coad_related_weak == T], paired=T)$p.val, 
  median((rank(pvals.coad$p.voom) - rank(pvals.coad$p.mean.lnHM.log))[coad_related_weak == T])) # voom better

c(wilcox.test(rank(pvals.brca$p.voom)[brca_related_weak == T], 
              rank(pvals.brca$p.disp.lnHM.log)[brca_related_weak == T], paired=T)$p.val, 
  median((rank(pvals.brca$p.voom) - rank(pvals.brca$p.disp.lnHM.log))[brca_related_weak == T])) # voom better
c(wilcox.test(rank(pvals.kirc$p.voom)[kirc_related_weak == T], 
              rank(pvals.kirc$p.disp.lnHM.log)[kirc_related_weak == T], paired=T)$p.val, 
  median((rank(pvals.kirc$p.voom) - rank(pvals.kirc$p.disp.lnHM.log))[kirc_related_weak == T])) # voom better
c(wilcox.test(rank(pvals.thca$p.voom)[thca_related_weak == T], 
              rank(pvals.thca$p.disp.lnHM.log)[thca_related_weak == T], paired=T)$p.val, 
  median((rank(pvals.thca$p.voom) - rank(pvals.thca$p.disp.lnHM.log))[thca_related_weak == T])) # No difference
c(wilcox.test(rank(pvals.luad$p.voom)[luad_related_weak == T], 
              rank(pvals.luad$p.disp.lnHM.log)[luad_related_weak == T], paired=T)$p.val, 
  median((rank(pvals.luad$p.voom) - rank(pvals.luad$p.disp.lnHM.log))[luad_related_weak == T])) # No difference
c(wilcox.test(rank(pvals.lihc$p.voom)[lihc_related_weak == T], 
              rank(pvals.lihc$p.disp.lnHM.log)[lihc_related_weak == T], paired=T)$p.val, 
  median((rank(pvals.lihc$p.voom) - rank(pvals.lihc$p.disp.lnHM.log))[lihc_related_weak == T])) # voom better
c(wilcox.test(rank(pvals.lihc_long.chain$p.voom)[lihc_related_weak == T], 
              rank(pvals.lihc_long.chain$p.disp.lnHM.log)[lihc_related_weak == T], paired=T)$p.val, 
  median((rank(pvals.lihc_long.chain$p.voom) - 
            rank(pvals.lihc_long.chain$p.disp.lnHM.log))[lihc_related_weak == T])) # voom better
c(wilcox.test(rank(pvals.lusc$p.voom)[lusc_related_weak == T], 
              rank(pvals.lusc$p.disp.lnHM.log)[lusc_related_weak == T], paired=T)$p.val, 
  median((rank(pvals.lusc$p.voom) - rank(pvals.lusc$p.disp.lnHM.log))[lusc_related_weak == T])) # No difference
c(wilcox.test(rank(pvals.prad$p.voom)[prad_related_weak == T], 
              rank(pvals.prad$p.disp.lnHM.log)[prad_related_weak == T], paired=T)$p.val, 
  median((rank(pvals.prad$p.voom) - rank(pvals.prad$p.disp.lnHM.log))[prad_related_weak == T])) # voom better
c(wilcox.test(rank(pvals.coad$p.voom)[coad_related_weak == T], 
              rank(pvals.coad$p.disp.lnHM.log)[coad_related_weak == T], paired=T)$p.val, 
  median((rank(pvals.coad$p.voom) - rank(pvals.coad$p.disp.lnHM.log))[coad_related_weak == T])) # voom better

c(wilcox.test(rank(pvals.brca$p.disp.MDSeq.zi)[brca_related_weak == T], 
              rank(pvals.brca$p.disp.lnHM.log)[brca_related_weak == T], paired=T)$p.val, 
  median((rank(pvals.brca$p.disp.MDSeq.zi) - rank(pvals.brca$p.disp.lnHM.log))[brca_related_weak == T])) # MD better
c(wilcox.test(rank(pvals.kirc$p.disp.MDSeq.zi)[kirc_related_weak == T], 
              rank(pvals.kirc$p.disp.lnHM.log)[kirc_related_weak == T], paired=T)$p.val, 
  median((rank(pvals.kirc$p.disp.MDSeq.zi) - rank(pvals.kirc$p.disp.lnHM.log))[kirc_related_weak == T])) # No difference
c(wilcox.test(rank(pvals.thca$p.disp.MDSeq.zi)[thca_related_weak == T], 
              rank(pvals.thca$p.disp.lnHM.log)[thca_related_weak == T], paired=T)$p.val, 
  median((rank(pvals.thca$p.disp.MDSeq.zi) - rank(pvals.thca$p.disp.lnHM.log))[thca_related_weak == T])) # No difference
c(wilcox.test(rank(pvals.luad$p.disp.MDSeq.zi)[luad_related_weak == T], 
              rank(pvals.luad$p.disp.lnHM.log)[luad_related_weak == T], paired=T)$p.val, 
  median((rank(pvals.luad$p.disp.MDSeq.zi) - rank(pvals.luad$p.disp.lnHM.log))[luad_related_weak == T])) # No difference
c(wilcox.test(rank(pvals.lihc$p.disp.MDSeq.zi)[lihc_related_weak == T], 
              rank(pvals.lihc$p.disp.lnHM.log)[lihc_related_weak == T], paired=T)$p.val, 
  median((rank(pvals.lihc$p.disp.MDSeq.zi) - rank(pvals.lihc$p.disp.lnHM.log))[lihc_related_weak == T])) # MD better
c(wilcox.test(rank(pvals.lihc_long.chain$p.disp.MDSeq.zi)[lihc_related_weak == T], 
              rank(pvals.lihc_long.chain$p.disp.lnHM.log)[lihc_related_weak == T], paired=T)$p.val, 
  median((rank(pvals.lihc_long.chain$p.disp.MDSeq.zi) - 
            rank(pvals.lihc_long.chain$p.disp.lnHM.log))[lihc_related_weak == T])) # MD better
c(wilcox.test(rank(pvals.lusc$p.disp.MDSeq.zi)[lusc_related_weak == T], 
              rank(pvals.lusc$p.disp.lnHM.log)[lusc_related_weak == T], paired=T)$p.val, 
  median((rank(pvals.lusc$p.disp.MDSeq.zi) - rank(pvals.lusc$p.disp.lnHM.log))[lusc_related_weak == T])) # No difference
c(wilcox.test(rank(pvals.prad$p.disp.MDSeq.zi)[prad_related_weak == T], 
              rank(pvals.prad$p.disp.lnHM.log)[prad_related_weak == T], paired=T)$p.val, 
  median((rank(pvals.prad$p.disp.MDSeq.zi) - rank(pvals.prad$p.disp.lnHM.log))[prad_related_weak == T])) # No difference
c(wilcox.test(rank(pvals.coad$p.disp.MDSeq.zi)[coad_related_weak == T], 
              rank(pvals.coad$p.disp.lnHM.log)[coad_related_weak == T], paired=T)$p.val, 
  median((rank(pvals.coad$p.disp.MDSeq.zi) - rank(pvals.coad$p.disp.lnHM.log))[coad_related_weak == T])) # MD better

c(wilcox.test(rank(pvals.brca$p.voom)[brca_related_weak == T], 
              rank(pvals.brca$lnHM.log_lnHM.log)[brca_related_weak == T], paired=T)$p.val, 
  median((rank(pvals.brca$p.voom) - rank(pvals.brca$lnHM.log_lnHM.log))[brca_related_weak == T])) # voom better
c(wilcox.test(rank(pvals.kirc$p.voom)[kirc_related_weak == T], 
              rank(pvals.kirc$lnHM.log_lnHM.log)[kirc_related_weak == T], paired=T)$p.val, 
  median((rank(pvals.kirc$p.voom) - rank(pvals.kirc$lnHM.log_lnHM.log))[kirc_related_weak == T])) # voom better
c(wilcox.test(rank(pvals.thca$p.voom)[thca_related_weak == T], 
              rank(pvals.thca$lnHM.log_lnHM.log)[thca_related_weak == T], paired=T)$p.val, 
  median((rank(pvals.thca$p.voom) - rank(pvals.thca$lnHM.log_lnHM.log))[thca_related_weak == T])) # voom better
c(wilcox.test(rank(pvals.luad$p.voom)[luad_related_weak == T], 
              rank(pvals.luad$lnHM.log_lnHM.log)[luad_related_weak == T], paired=T)$p.val, 
  median((rank(pvals.luad$p.voom) - rank(pvals.luad$lnHM.log_lnHM.log))[luad_related_weak == T])) # voom better
c(wilcox.test(rank(pvals.lihc$p.voom)[lihc_related_weak == T], 
              rank(pvals.lihc$lnHM.log_lnHM.log)[lihc_related_weak == T], paired=T)$p.val, 
  median((rank(pvals.lihc$p.voom) - rank(pvals.lihc$lnHM.log_lnHM.log))[lihc_related_weak == T])) # voom better
c(wilcox.test(rank(pvals.lihc_long.chain$p.voom)[lihc_related_weak == T], 
              rank(pvals.lihc_long.chain$lnHM.log_lnHM.log)[lihc_related_weak == T], paired=T)$p.val, 
  median((rank(pvals.lihc_long.chain$p.voom) - 
            rank(pvals.lihc_long.chain$lnHM.log_lnHM.log))[lihc_related_weak == T])) # voom better
c(wilcox.test(rank(pvals.lusc$p.voom)[lusc_related_weak == T], 
              rank(pvals.lusc$lnHM.log_lnHM.log)[lusc_related_weak == T], paired=T)$p.val, 
  median((rank(pvals.lusc$p.voom) - rank(pvals.lusc$lnHM.log_lnHM.log))[lusc_related_weak == T])) # voom better
c(wilcox.test(rank(pvals.prad$p.voom)[prad_related_weak == T], 
              rank(pvals.prad$lnHM.log_lnHM.log)[prad_related_weak == T], paired=T)$p.val, 
  median((rank(pvals.prad$p.voom) - rank(pvals.prad$lnHM.log_lnHM.log))[prad_related_weak == T])) # voom better
c(wilcox.test(rank(pvals.coad$p.voom)[coad_related_weak == T], 
              rank(pvals.coad$lnHM.log_lnHM.log)[coad_related_weak == T], paired=T)$p.val, 
  median((rank(pvals.coad$p.voom) - rank(pvals.coad$lnHM.log_lnHM.log))[coad_related_weak == T])) # voom better

c(wilcox.test(rank(pvals.brca$p.voom)[brca_related_weak == T], 
              rank(1 - pvals.brca$prob.lnHM)[brca_related_weak == T], paired=T)$p.val, 
  median((rank(pvals.brca$p.voom) - rank(1 - pvals.brca$prob.lnHM))[brca_related_weak == T])) # voom better
c(wilcox.test(rank(pvals.kirc$p.voom)[kirc_related_weak == T], 
              rank(1 - pvals.kirc$prob.lnHM)[kirc_related_weak == T], paired=T)$p.val, 
  median((rank(pvals.kirc$p.voom) - rank(1 - pvals.kirc$prob.lnHM))[kirc_related_weak == T])) # voom better
c(wilcox.test(rank(pvals.thca$p.voom)[thca_related_weak == T], 
              rank(1 - pvals.thca$prob.lnHM)[thca_related_weak == T], paired=T)$p.val, 
  median((rank(pvals.thca$p.voom) - rank(1 - pvals.thca$prob.lnHM))[thca_related_weak == T])) # voom better
c(wilcox.test(rank(pvals.luad$p.voom)[luad_related_weak == T], 
              rank(1 - pvals.luad$prob.lnHM)[luad_related_weak == T], paired=T)$p.val, 
  median((rank(pvals.luad$p.voom) - rank(1 - pvals.luad$prob.lnHM))[luad_related_weak == T])) # voom better
c(wilcox.test(rank(pvals.lihc$p.voom)[lihc_related_weak == T], 
              rank(1 - pvals.lihc$prob.lnHM)[lihc_related_weak == T], paired=T)$p.val, 
  median((rank(pvals.lihc$p.voom) - rank(1 - pvals.lihc$prob.lnHM))[lihc_related_weak == T])) # voom better
c(wilcox.test(rank(pvals.lihc_long.chain$p.voom)[lihc_related_weak == T], 
              rank(1 - pvals.lihc_long.chain$prob.lnHM)[lihc_related_weak == T], paired=T)$p.val, 
  median((rank(pvals.lihc_long.chain$p.voom) - 
            rank(1 - pvals.lihc_long.chain$prob.lnHM))[lihc_related_weak == T])) # voom better
c(wilcox.test(rank(pvals.lusc$p.voom)[lusc_related_weak == T], 
              rank(1 - pvals.lusc$prob.lnHM)[lusc_related_weak == T], paired=T)$p.val, 
  median((rank(pvals.lusc$p.voom) - rank(1 - pvals.lusc$prob.lnHM))[lusc_related_weak == T])) # voom better
c(wilcox.test(rank(pvals.prad$p.voom)[prad_related_weak == T], 
              rank(1 - pvals.prad$prob.lnHM)[prad_related_weak == T], paired=T)$p.val, 
  median((rank(pvals.prad$p.voom) - rank(1 - pvals.prad$prob.lnHM))[prad_related_weak == T])) # voom better
c(wilcox.test(rank(pvals.coad$p.voom)[coad_related_weak == T], 
              rank(1 - pvals.coad$prob.lnHM)[coad_related_weak == T], paired=T)$p.val, 
  median((rank(pvals.coad$p.voom) - rank(1 - pvals.coad$prob.lnHM))[coad_related_weak == T])) # voom better


## Assess whether DE, DD, DEDD identify different sets of cancer-related genes ####
## Assess whether DE, DD, DEDD identify different sets of cancer-related genes ####
# Informally assess whether DE and DD HMs identify different sets of cancer-related genes by plotting 
# hypergeometric and Spearman correlation test p-values with alternative hypothesis that overlap and 
# correlation, respectively, are lower than expected under null for sets of top-ranking genes from 1 
# to all genes (or possibly for p-value threshold from 0 to 1 for hypergeometric).
for (j in c("brca", "kirc", "thca", "luad", "lihc", 
            "lihc_long.chain",  "lusc", "prad", "coad")) {
  pvals <- sort(unique(c(get(paste0("pvals.", j))$p.mean.lnHM.log, get(paste0("pvals.", j))$p.disp.lnHM.log)))
  res <- data.frame(
    pvals = pvals, 
    listA = numeric(length(pvals)), 
    listB = numeric(length(pvals)), 
    overlap = numeric(length(pvals)), 
    p.hyper = numeric(length(pvals)), 
    exp.overlap = numeric(length(pvals)), 
    union. = numeric(length(pvals)), 
    p.cor = numeric(length(pvals)), 
    cor = numeric(length(pvals)), 
    TPA = numeric(length(pvals)), 
    TPB = numeric(length(pvals)), 
    TPboth = numeric(length(pvals)), 
    TPAonly = numeric(length(pvals)), 
    TPBonly = numeric(length(pvals))
  )
  for (i in seq_len(length(pvals))) {
    res$listA[i] <- sum(get(paste0("pvals.", j))$p.mean.lnHM.log < res$pvals[i])
    res$listB[i] <- sum(get(paste0("pvals.", j))$p.disp.lnHM.log < res$pvals[i])
    res$overlap[i] <- sum(get(paste0("pvals.", j))$p.mean.lnHM.log < res$pvals[i] & 
                            get(paste0("pvals.", j))$p.disp.lnHM.log < res$pvals[i])
    res$p.hyper[i] <- phyper(res$overlap[i], 
                             res$listA[i], 
                             nrow(get(paste0("pvals.", j))) - res$listA[i], 
                             res$listB[i])
    res$exp.overlap[i] <- mean(get(paste0("pvals.", j))$p.mean.lnHM.log < res$pvals[i]) * 
      mean(get(paste0("pvals.", j))$p.disp.lnHM.log < res$pvals[i]) * nrow(get(paste0("pvals.", j)))
    res$union[i] <- sum(get(paste0("pvals.", j))$p.mean.lnHM.log < res$pvals[i] | 
                          get(paste0("pvals.", j))$p.disp.lnHM.log < res$pvals[i])
    if (sum(get(paste0("pvals.", j))$p.mean.lnHM.log < res$pvals[i] | 
            get(paste0("pvals.", j))$p.disp.lnHM.log < res$pvals[i]) > 1) {
      cortest <- 
        cor.test(get(paste0("pvals.", j))$p.mean.lnHM.log[get(paste0("pvals.", j))$p.mean.lnHM.log < res$pvals[i] | 
                                                            get(paste0("pvals.", j))$p.disp.lnHM.log < res$pvals[i]], 
                 get(paste0("pvals.", j))$p.disp.lnHM.log[get(paste0("pvals.", j))$p.mean.lnHM.log < res$pvals[i] | 
                                                            get(paste0("pvals.", j))$p.disp.lnHM.log < res$pvals[i]], 
                 method="spearman", alt="less")
      res$p.cor[i] <- cortest$p.val
      res$cor[i] <- cortest$estimate
      rm(cortest)
    }
    res$TPA[i] <- sum(get(paste0("pvals.", j))$p.mean.lnHM.log[get(paste0(j, "_related_weak"))] < res$pvals[i])
    res$TPB[i] <- sum(get(paste0("pvals.", j))$p.disp.lnHM.log[get(paste0(j, "_related_weak"))] < res$pvals[i])
    res$TPboth[i] <- sum(get(paste0("pvals.", j))$p.mean.lnHM.log[get(paste0(j, "_related_weak"))] < res$pvals[i] & 
                           get(paste0("pvals.", j))$p.disp.lnHM.log[get(paste0(j, "_related_weak"))] < res$pvals[i])
    res$TPAonly[i] <- sum(get(paste0("pvals.", j))$p.mean.lnHM.log[get(paste0(j, "_related_weak"))] < res$pvals[i] & 
                            get(paste0("pvals.", j))$p.disp.lnHM.log[get(paste0(j, "_related_weak"))] >= res$pvals[i])
    res$TPBonly[i] <- sum(get(paste0("pvals.", j))$p.mean.lnHM.log[get(paste0(j, "_related_weak"))] >= res$pvals[i] & 
                            get(paste0("pvals.", j))$p.disp.lnHM.log[get(paste0(j, "_related_weak"))] < res$pvals[i])
  }
  assign(paste0(j, ".HM"), res)
  write.csv(get(paste0(j, ".HM")), 
            here("Results/TCGA paired data results May 2020", 
                 paste0(j, "_HM_mean_v_disp_pathway.genes.csv")), 
            row.names=F)
  rm(pvals, res)
  
  pvals <- sort(unique(c(get(paste0("pvals.", j))$p.mean.MDSeq.zi, get(paste0("pvals.", j))$p.disp.MDSeq.zi)))
  res <- data.frame(
    pvals = pvals, 
    listA = numeric(length(pvals)), 
    listB = numeric(length(pvals)), 
    overlap = numeric(length(pvals)), 
    p.hyper = numeric(length(pvals)), 
    exp.overlap = numeric(length(pvals)), 
    union. = numeric(length(pvals)), 
    p.cor = numeric(length(pvals)), 
    cor = numeric(length(pvals)), 
    TPA = numeric(length(pvals)), 
    TPB = numeric(length(pvals)), 
    TPboth = numeric(length(pvals)), 
    TPAonly = numeric(length(pvals)), 
    TPBonly = numeric(length(pvals))
  )
  for (i in seq_len(length(pvals))) {
    res$listA[i] <- sum(get(paste0("pvals.", j))$p.mean.MDSeq.zi < res$pvals[i], na.rm=T)
    res$listB[i] <- sum(get(paste0("pvals.", j))$p.disp.MDSeq.zi < res$pvals[i], na.rm=T)
    res$overlap[i] <- sum(get(paste0("pvals.", j))$p.mean.MDSeq.zi < res$pvals[i] & 
                            get(paste0("pvals.", j))$p.disp.MDSeq.zi < res$pvals[i], na.rm=T)
    res$union[i] <- sum(get(paste0("pvals.", j))$p.mean.MDSeq.zi < res$pvals[i] | 
                          get(paste0("pvals.", j))$p.disp.MDSeq.zi < res$pvals[i], na.rm=T)
    res$exp.overlap[i] <- mean(get(paste0("pvals.", j))$p.mean.MDSeq.zi < res$pvals[i], na.rm=T) * 
      mean(get(paste0("pvals.", j))$p.disp.MDSeq.zi < res$pvals[i], na.rm=T) * nrow(get(paste0("pvals.", j)))
    res$p.hyper[i] <- phyper(res$overlap[i], 
                             res$listA[i], 
                             nrow(get(paste0("pvals.", j))) - res$listA[i], 
                             res$listB[i])
    if (sum(get(paste0("pvals.", j))$p.mean.MDSeq.zi < res$pvals[i] | 
            get(paste0("pvals.", j))$p.disp.MDSeq.zi < res$pvals[i], na.rm=T) > 1) {
      cortest <- cor.test(get(paste0("pvals.", j))$p.mean.MDSeq.zi[get(paste0("pvals.", j))$p.mean.MDSeq.zi < res$pvals[i] | 
                                                                     get(paste0("pvals.", j))$p.disp.MDSeq.zi < res$pvals[i]], 
                          get(paste0("pvals.", j))$p.disp.MDSeq.zi[get(paste0("pvals.", j))$p.mean.MDSeq.zi < res$pvals[i] | 
                                                                     get(paste0("pvals.", j))$p.disp.MDSeq.zi < res$pvals[i]], 
                          method="spearman", alt="less")
      res$p.cor[i] <- cortest$p.val
      res$cor[i] <- cortest$estimate
      rm(cortest)
    }
    res$TPA[i] <- sum(get(paste0("pvals.", j))$p.mean.MDSeq.zi[get(paste0(j, "_related_weak"))] < res$pvals[i], na.rm=T)
    res$TPB[i] <- sum(get(paste0("pvals.", j))$p.disp.MDSeq.zi[get(paste0(j, "_related_weak"))] < res$pvals[i], na.rm=T)
    res$TPboth[i] <- sum(get(paste0("pvals.", j))$p.mean.MDSeq.zi[get(paste0(j, "_related_weak"))] < res$pvals[i] & 
                           get(paste0("pvals.", j))$p.disp.MDSeq.zi[get(paste0(j, "_related_weak"))] < res$pvals[i], na.rm=T)
    res$TPAonly[i] <- sum(get(paste0("pvals.", j))$p.mean.MDSeq.zi[get(paste0(j, "_related_weak"))] < res$pvals[i] & 
                            get(paste0("pvals.", j))$p.disp.MDSeq.zi[get(paste0(j, "_related_weak"))] >= res$pvals[i], na.rm=T)
    res$TPBonly[i] <- sum(get(paste0("pvals.", j))$p.mean.MDSeq.zi[get(paste0(j, "_related_weak"))] >= res$pvals[i] & 
                            get(paste0("pvals.", j))$p.disp.MDSeq.zi[get(paste0(j, "_related_weak"))] < res$pvals[i], na.rm=T)
  }
  assign(paste0(j, ".MD"), res)
  write.csv(get(paste0(j, ".MD")), 
            here("Results/TCGA paired data results May 2020", 
                 paste0(j, "_MDSeq_mean_v_disp_pathway.genes.csv")), 
            row.names=F)
  rm(pvals, res)
  
  pvals <- sort(unique(c(get(paste0("pvals.", j))$p.voom, get(paste0("pvals.", j))$p.disp.lnHM.log)))
  res <- data.frame(
    pvals = pvals, 
    listA = numeric(length(pvals)), 
    listB = numeric(length(pvals)), 
    overlap = numeric(length(pvals)), 
    p.hyper = numeric(length(pvals)), 
    exp.overlap = numeric(length(pvals)), 
    union. = numeric(length(pvals)), 
    p.cor = numeric(length(pvals)), 
    cor = numeric(length(pvals)), 
    TPA = numeric(length(pvals)), 
    TPB = numeric(length(pvals)), 
    TPboth = numeric(length(pvals)), 
    TPAonly = numeric(length(pvals)), 
    TPBonly = numeric(length(pvals))
  )
  for (i in seq_len(nrow(res))) {
    res$listA[i] <- sum(get(paste0("pvals.", j))$p.voom < res$pvals[i])
    res$listB[i] <- sum(get(paste0("pvals.", j))$p.disp.lnHM.log < res$pvals[i])
    res$overlap[i] <- sum(get(paste0("pvals.", j))$p.voom < res$pvals[i] & 
                            get(paste0("pvals.", j))$p.disp.lnHM.log < res$pvals[i])
    res$p.hyper[i] <- phyper(res$overlap[i], 
                             res$listA[i], 
                             nrow(get(paste0("pvals.", j))) - res$listA[i], 
                             res$listB[i])
    res$exp.overlap[i] <- mean(get(paste0("pvals.", j))$p.voom < res$pvals[i]) * 
      mean(get(paste0("pvals.", j))$p.disp.lnHM.log < res$pvals[i]) * nrow(get(paste0("pvals.", j)))
    res$union[i] <- sum(get(paste0("pvals.", j))$p.voom < res$pvals[i] | 
                          get(paste0("pvals.", j))$p.disp.lnHM.log < res$pvals[i])
    if (sum(get(paste0("pvals.", j))$p.voom < res$pvals[i] | 
            get(paste0("pvals.", j))$p.disp.lnHM.log < res$pvals[i]) > 1) {
      cortest <- cor.test(get(paste0("pvals.", j))$p.voom[get(paste0("pvals.", j))$p.voom < res$pvals[i] | 
                                                            get(paste0("pvals.", j))$p.disp.lnHM.log < res$pvals[i]], 
                          get(paste0("pvals.", j))$p.disp.lnHM.log[get(paste0("pvals.", j))$p.voom < res$pvals[i] | 
                                                                     get(paste0("pvals.", j))$p.disp.lnHM.log < res$pvals[i]], 
                          method="spearman", alt="less")
      res$p.cor[i] <- cortest$p.val
      res$cor[i] <- cortest$estimate
      rm(cortest)
    }
    res$TPA[i] <- sum(get(paste0("pvals.", j))$p.voom[get(paste0(j, "_related_weak"))] < res$pvals[i])
    res$TPB[i] <- sum(get(paste0("pvals.", j))$p.disp.lnHM.log[get(paste0(j, "_related_weak"))] < res$pvals[i])
    res$TPboth[i] <- sum(get(paste0("pvals.", j))$p.voom[get(paste0(j, "_related_weak"))] < res$pvals[i] & 
                           get(paste0("pvals.", j))$p.disp.lnHM.log[get(paste0(j, "_related_weak"))] < res$pvals[i])
    res$TPAonly[i] <- sum(get(paste0("pvals.", j))$p.voom[get(paste0(j, "_related_weak"))] < res$pvals[i] & 
                            get(paste0("pvals.", j))$p.disp.lnHM.log[get(paste0(j, "_related_weak"))] >= res$pvals[i])
    res$TPBonly[i] <- sum(get(paste0("pvals.", j))$p.voom[get(paste0(j, "_related_weak"))] >= res$pvals[i] & 
                            get(paste0("pvals.", j))$p.disp.lnHM.log[get(paste0(j, "_related_weak"))] < res$pvals[i])
  }
  assign(paste0(j, ".HMvoom"), res)
  write.csv(get(paste0(j, ".HMvoom")), 
            here("Results/TCGA paired data results May 2020", 
                 paste0(j, "_HM_voom_mean_v_disp_pathway.genes.csv")), 
            row.names=F)
  rm(pvals, res)
  
  pvals <- sort(unique(c(get(paste0("pvals.", j))$p.voom, 1 - get(paste0("pvals.", j))$prob.lnHMM)))
  res <- data.frame(
    pvals = pvals, 
    listA = numeric(length(pvals)), 
    listB = numeric(length(pvals)), 
    overlap = numeric(length(pvals)), 
    p.hyper = numeric(length(pvals)), 
    exp.overlap = numeric(length(pvals)), 
    union. = numeric(length(pvals)), 
    p.cor = numeric(length(pvals)), 
    cor = numeric(length(pvals)), 
    TPA = numeric(length(pvals)), 
    TPB = numeric(length(pvals)), 
    TPboth = numeric(length(pvals)), 
    TPAonly = numeric(length(pvals)), 
    TPBonly = numeric(length(pvals))
  )
  for (i in seq_len(nrow(res))) {
    res$listA[i] <- sum(get(paste0("pvals.", j))$p.voom < res$pvals[i])
    res$listB[i] <- sum(1 - get(paste0("pvals.", j))$prob.lnHMM < res$pvals[i])
    res$overlap[i] <- sum(get(paste0("pvals.", j))$p.voom < res$pvals[i] & 
                            1 - get(paste0("pvals.", j))$prob.lnHMM < res$pvals[i])
    res$p.hyper[i] <- phyper(res$overlap[i], 
                             res$listA[i], 
                             nrow(get(paste0("pvals.", j))) - res$listA[i], 
                             res$listB[i])
    res$exp.overlap[i] <- mean(get(paste0("pvals.", j))$p.voom < res$pvals[i]) * 
      mean(1 - get(paste0("pvals.", j))$prob.lnHMM < res$pvals[i]) * nrow(get(paste0("pvals.", j)))
    res$union[i] <- sum(get(paste0("pvals.", j))$p.voom < res$pvals[i] | 
                          1 - get(paste0("pvals.", j))$prob.lnHMM < res$pvals[i])
    if (sum(get(paste0("pvals.", j))$p.voom < res$pvals[i] | 
            1 - get(paste0("pvals.", j))$prob.lnHMM < res$pvals[i]) > 1) {
      cortest <- cor.test(get(paste0("pvals.", j))$p.voom[get(paste0("pvals.", j))$p.voom < res$pvals[i] | 
                                                            1 - get(paste0("pvals.", j))$prob.lnHMM < res$pvals[i]], 
                          1 - get(paste0("pvals.", j))$prob.lnHMM[get(paste0("pvals.", j))$p.voom < res$pvals[i] | 
                                                                    1 - get(paste0("pvals.", j))$prob.lnHMM < res$pvals[i]], 
                          method="spearman", alt="less")
      res$p.cor[i] <- cortest$p.val
      res$cor[i] <- cortest$estimate
      rm(cortest)
    }
    res$TPA[i] <- sum(get(paste0("pvals.", j))$p.voom[get(paste0(j, "_related_weak"))] < res$pvals[i])
    res$TPB[i] <- sum(1 - get(paste0("pvals.", j))$prob.lnHMM[get(paste0(j, "_related_weak"))] < res$pvals[i])
    res$TPboth[i] <- sum(get(paste0("pvals.", j))$p.voom[get(paste0(j, "_related_weak"))] < res$pvals[i] & 
                           1 - get(paste0("pvals.", j))$prob.lnHMM[get(paste0(j, "_related_weak"))] < res$pvals[i])
    res$TPAonly[i] <- sum(get(paste0("pvals.", j))$p.voom[get(paste0(j, "_related_weak"))] < res$pvals[i] & 
                            1 - get(paste0("pvals.", j))$prob.lnHMM[get(paste0(j, "_related_weak"))] >= res$pvals[i])
    res$TPBonly[i] <- sum(get(paste0("pvals.", j))$p.voom[get(paste0(j, "_related_weak"))] >= res$pvals[i] & 
                            1 - get(paste0("pvals.", j))$prob.lnHMM[get(paste0(j, "_related_weak"))] < res$pvals[i])
  }
  assign(paste0(j, ".HMMvoom"), res)
  write.csv(get(paste0(j, ".HMMvoom")), 
            here("Results/TCGA paired data results May 2020", 
                 paste0(j, "_HMM_voom_mean_v_disp_pathway.genes.csv")), 
            row.names=F)
  rm(pvals, res)
}
rm(i,j)


par(mfrow=c(2,2), mar=c(4,2,2,1), mgp=c(1.5,0.5,0))
plot(luad.HM$pvals, luad.HM$cor, type='l', ylim=c(-1,1), col="red", 
     main="Correlation between ranks of genes identified by differential expression and dispersion", 
     xlab="p-value threshold", ylab="")
lines(luad.MD$pvals, luad.MD$cor, col="blue")
plot(luad.HM$pvals, luad.HM$p.cor, type='l', ylim=c(0,1), col="red", 
     main="p-value of correlation test", 
     xlab="p-value threshold for differential expression or dispersion", ylab="")
lines(luad.MD$pvals, luad.MD$p.cor, col="blue")
plot(luad.HM$pvals, luad.HM$overlap / luad.HM$union, type='l', ylim=c(0,1), col="red", 
     main="Union / Intersection", 
     xlab="p-value threshold", ylab="")
lines(luad.MD$pvals, luad.MD$overlap / luad.MD$union, col='blue', 
      xlab="p-value threshold", ylab="")
plot(luad.HM$pvals, luad.HM$p.hyper, type='l', ylim=c(0,1), col="red", 
     main="p-value of test of overlap", 
     xlab="p-value threshold", ylab="")
lines(luad.MD$pvals, luad.MD$p.hyper, col='blue')


# genes identified by diff disp/dist missed by DE
par(mfrow=c(2,2), mar=c(4,2,2,1), mgp=c(1.5,0.5,0))
plot(luad.HM$pvals, luad.HM$TPA, type='l', col='red')
lines(luad.HM$pvals, luad.HM$TPB, col='blue')
lines(luad.HM$pvals, luad.HM$TPboth, col='grey')
lines(luad.HM$pvals, luad.HM$TPAonly, col='red', lty=2)
lines(luad.HM$pvals, luad.HM$TPBonly, col='blue', lty=2)
plot(luad.MD$pvals, luad.MD$TPA, type='l', col='red')
lines(luad.MD$pvals, luad.MD$TPB, col='blue')
lines(luad.MD$pvals, luad.MD$TPboth, col='grey')
lines(luad.MD$pvals, luad.MD$TPAonly, col='red', lty=2)
lines(luad.MD$pvals, luad.MD$TPBonly, col='blue', lty=2)
plot(luad.HMvoom$pvals, luad.HMvoom$TPA, type='l', col='red', 
     main="Genes identified by differential expression (red) and differential dispersion (blue)", 
     ylab="", xlab="p-value threshold")
lines(luad.HMvoom$pvals, luad.HMvoom$TPB, col='blue')
lines(luad.HMvoom$pvals, luad.HMvoom$TPboth, col='grey')
lines(luad.HMvoom$pvals, luad.HMvoom$TPAonly, col='red', lty=2)
lines(luad.HMvoom$pvals, luad.HMvoom$TPBonly, col='blue', lty=2)
plot(luad.HMMvoom$pvals, luad.HMMvoom$TPA, type='l', col='red', 
     main="Genes identified by differential expression (red) and differential distribution (blue)", 
     ylab="", xlab="p-value threshold")
lines(luad.HMMvoom$pvals, luad.HMMvoom$TPB, col='blue')
lines(luad.HMMvoom$pvals, luad.HMMvoom$TPboth, col='grey')
lines(luad.HMMvoom$pvals, luad.HMMvoom$TPAonly, col='red', lty=2)
lines(luad.HMMvoom$pvals, luad.HMMvoom$TPBonly, col='blue', lty=2)


par(mfrow=c(2,1), mar=c(4,2,2,1), mgp=c(1.5,0.5,0))
# plot(pvals.HM, TP.mean.HM, type='l', col='red')
# lines(pvals.HM, TP.disp.HM, col='blue')
# lines(pvals.HM, TP.both.HM, col='grey')
# lines(pvals.HM, TP.mean.only.HM, col='red', lty=2)
# lines(pvals.HM, TP.disp.only.HM, col='blue', lty=2)
# plot(pvals.MD, TP.mean.MD, type='l', col='red')
# lines(pvals.MD, TP.disp.MD, col='blue')
# lines(pvals.MD, TP.both.MD, col='grey')
# lines(pvals.MD, TP.mean.only.MD, col='red', lty=2)
# lines(pvals.MD, TP.disp.only.MD, col='blue', lty=2)
plot(pvals.HMvoom, TP.mean.HMvoom, type='l', col='red', 
     main="Genes identified by differential expression (red) and differential dispersion (blue)", 
     ylab="", xlab="p-value threshold")
lines(pvals.HMvoom, TP.disp.HMvoom, col='blue')
lines(pvals.HMvoom, TP.both.HMvoom, col='grey')
lines(pvals.HMvoom, TP.mean.only.HMvoom, col='red', lty=2)
lines(pvals.HMvoom, TP.disp.only.HMvoom, col='blue', lty=2)
plot(pvals.HMMvoom, TP.mean.HMMvoom, type='l', col='red', 
     main="Genes identified by differential expression (red) and differential distribution (blue)", 
     ylab="", xlab="p-value threshold")
lines(pvals.HMMvoom, TP.disp.HMMvoom, col='blue')
lines(pvals.HMMvoom, TP.both.HMMvoom, col='grey')
lines(pvals.HMMvoom, TP.mean.only.HMMvoom, col='red', lty=2)
lines(pvals.HMMvoom, TP.disp.only.HMMvoom, col='blue', lty=2)


plot(P.mean.HM, TP.mean.HM, type='l', col='red')
lines(P.disp.HM, TP.disp.HM, col='blue')
# lines(pvals.HM, TP.both.HM, col='grey')
lines(P.mean.HM, TP.mean.only.HM, col='red', lty=2)
lines(P.disp.HM, TP.disp.only.HM, col='blue', lty=2)
# plotting mean and disp discoveries against number of genes found by each method rather 
# than by p-value threshold, but what to do for overlapping genes?
# straight line from 0 to about 13k because no numbers of genes discovered in between.

plot(P.mean.HM, TP.mean.HM / P.mean.HM, type='l', col='red', ylim=c(0,0.3))
lines(P.disp.HM, TP.disp.HM / P.disp.HM, col='blue')
# lines(pvals.HM, TP.both.HM, col='grey')
lines(P.mean.HM, TP.mean.only.HM / P.mean.HM, col='red', lty=2)
lines(P.disp.HM, TP.disp.only.HM / P.disp.HM, col='blue', lty=2)
# Same but using proportion of known cancer genes discovered to head off criticism that 
# number of related genes alone doesn't mean much if there are a lot of false positives.


