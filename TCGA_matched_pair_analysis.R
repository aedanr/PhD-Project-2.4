library(here)

folder <- paste0("Results/TCGA paired data results May 2020")
for (i in c("brca", "kirc", "thca", "luad", "lihc", 
            "lihc_long.chain", "lusc", "prad", "coad")) {
  assign(paste0("res.", i), 
         readRDS(here(folder, paste0("results.TCGA.paired_", i, ".rds"))))
  assign(paste0("pvals.", i), 
         as.data.frame(
           get(paste0("res.", i))[c(2,3,4,8,9,12,14,17,19,20,23,24,27,31,32,35,36,39,43,44,47,48)]
         ))
  assign(paste0("qvals.", i), 
         as.data.frame(
           get(paste0("res.", i))[c(5,6,7,10,11,13,18,21,22,25,26,30,33,34,37,38,42,45,46,49,50)]
         ))
}

{
  rownames(pvals.brca) <- gsub('\\..*', '', rownames(res.brca$counts))
  rownames(pvals.kirc) <- gsub('\\..*', '', rownames(res.kirc$counts))
  rownames(pvals.thca) <- gsub('\\..*', '', rownames(res.thca$counts))
  rownames(pvals.luad) <- gsub('\\..*', '', rownames(res.luad$counts))
  rownames(pvals.lihc) <- gsub('\\..*', '', rownames(res.lihc$counts))
  rownames(pvals.lihc_long.chain) <- gsub('\\..*', '', rownames(res.lihc_long.chain$counts))
  rownames(pvals.lusc) <- gsub('\\..*', '', rownames(res.lusc$counts))
  rownames(pvals.prad) <- gsub('\\..*', '', rownames(res.prad$counts))
  rownames(pvals.coad) <- gsub('\\..*', '', rownames(res.coad$counts))
  rownames(qvals.brca) <- rownames(pvals.brca$counts)
  rownames(qvals.kirc) <- rownames(pvals.kirc$counts)
  rownames(qvals.thca) <- rownames(pvals.thca$counts)
  rownames(qvals.luad) <- rownames(pvals.luad$counts)
  rownames(qvals.lihc) <- rownames(pvals.lihc$counts)
  rownames(qvals.lihc_long.chain) <- rownames(pvals.lihc_long.chain$counts)
  rownames(qvals.lusc) <- rownames(pvals.lusc$counts)
  rownames(qvals.prad) <- rownames(pvals.prad$counts)
  rownames(qvals.coad) <- rownames(pvals.coad$counts)
}

for (i in c("brca", "kirc", "thca", "luad", "lihc", 
            "lihc_long.chain", "lusc", "prad", "coad")) {
  assign(paste0("calls.", i), 
         data.frame(
           get(paste0("qvals.", i)) < 0.05
         ))
  rm(list=paste0("qvals.", i))
}

{
  pvals.brca$voom_lnHM.log <- pmin(pvals.brca$p.voom, pvals.brca$p.disp.lnHM.log)
  pvals.brca$voom_MDSeq.zi <- pmin(pvals.brca$p.voom, pvals.brca$p.disp.MDSeq.zi)
  pvals.brca$lnHM.log_lnHM.log <- pmin(pvals.brca$p.mean.lnHM.log, pvals.brca$p.disp.lnHM.log)
  pvals.kirc$voom_lnHM.log <- pmin(pvals.kirc$p.voom, pvals.kirc$p.disp.lnHM.log)
  pvals.kirc$voom_MDSeq.zi <- pmin(pvals.kirc$p.voom, pvals.kirc$p.disp.MDSeq.zi)
  pvals.kirc$lnHM.log_lnHM.log <- pmin(pvals.kirc$p.mean.lnHM.log, pvals.kirc$p.disp.lnHM.log)
  pvals.thca$voom_lnHM.log <- pmin(pvals.thca$p.voom, pvals.thca$p.disp.lnHM.log)
  pvals.thca$voom_MDSeq.zi <- pmin(pvals.thca$p.voom, pvals.thca$p.disp.MDSeq.zi)
  pvals.thca$lnHM.log_lnHM.log <- pmin(pvals.thca$p.mean.lnHM.log, pvals.thca$p.disp.lnHM.log)
  pvals.luad$voom_lnHM.log <- pmin(pvals.luad$p.voom, pvals.luad$p.disp.lnHM.log)
  pvals.luad$voom_MDSeq.zi <- pmin(pvals.luad$p.voom, pvals.luad$p.disp.MDSeq.zi)
  pvals.luad$lnHM.log_lnHM.log <- pmin(pvals.luad$p.mean.lnHM.log, pvals.luad$p.disp.lnHM.log)
  pvals.lihc$voom_lnHM.log <- pmin(pvals.lihc$p.voom, pvals.lihc$p.disp.lnHM.log)
  pvals.lihc$voom_MDSeq.zi <- pmin(pvals.lihc$p.voom, pvals.lihc$p.disp.MDSeq.zi)
  pvals.lihc$lnHM.log_lnHM.log <- pmin(pvals.lihc$p.mean.lnHM.log, pvals.lihc$p.disp.lnHM.log)
  pvals.lihc_long.chain$voom_lnHM.log <- pmin(pvals.lihc_long.chain$p.voom, pvals.lihc_long.chain$p.disp.lnHM.log)
  pvals.lihc_long.chain$voom_MDSeq.zi <- pmin(pvals.lihc_long.chain$p.voom, pvals.lihc_long.chain$p.disp.MDSeq.zi)
  pvals.lihc_long.chain$lnHM.log_lnHM.log <- pmin(pvals.lihc_long.chain$p.mean.lnHM.log, pvals.lihc_long.chain$p.disp.lnHM.log)
  pvals.lusc$voom_lnHM.log <- pmin(pvals.lusc$p.voom, pvals.lusc$p.disp.lnHM.log)
  pvals.lusc$voom_MDSeq.zi <- pmin(pvals.lusc$p.voom, pvals.lusc$p.disp.MDSeq.zi)
  pvals.lusc$lnHM.log_lnHM.log <- pmin(pvals.lusc$p.mean.lnHM.log, pvals.lusc$p.disp.lnHM.log)
  pvals.prad$voom_lnHM.log <- pmin(pvals.prad$p.voom, pvals.prad$p.disp.lnHM.log)
  pvals.prad$voom_MDSeq.zi <- pmin(pvals.prad$p.voom, pvals.prad$p.disp.MDSeq.zi)
  pvals.prad$lnHM.log_lnHM.log <- pmin(pvals.prad$p.mean.lnHM.log, pvals.prad$p.disp.lnHM.log)
  pvals.coad$voom_lnHM.log <- pmin(pvals.coad$p.voom, pvals.coad$p.disp.lnHM.log)
  pvals.coad$voom_MDSeq.zi <- pmin(pvals.coad$p.voom, pvals.coad$p.disp.MDSeq.zi)
  pvals.coad$lnHM.log_lnHM.log <- pmin(pvals.coad$p.mean.lnHM.log, pvals.coad$p.disp.lnHM.log)
  calls.brca$thr.expHMM <- pvals.brca$prob.expHMM > res.brca$thr.expHMM
  calls.brca$point5.expHMM <- pvals.brca$prob.expHMM > 0.5
  calls.brca$thr.lnHMM <- pvals.brca$prob.lnHMM > res.brca$thr.lnHMM
  calls.brca$point5.lnHMM <- pvals.brca$prob.lnHMM > 0.5
  calls.kirc$thr.expHMM <- pvals.kirc$prob.expHMM > res.kirc$thr.expHMM
  calls.kirc$point5.expHMM <- pvals.kirc$prob.expHMM > 0.5
  calls.kirc$thr.lnHMM <- pvals.kirc$prob.lnHMM > res.kirc$thr.lnHMM
  calls.kirc$point5.lnHMM <- pvals.kirc$prob.lnHMM > 0.5
  calls.thca$thr.expHMM <- pvals.thca$prob.expHMM > res.thca$thr.expHMM
  calls.thca$point5.expHMM <- pvals.thca$prob.expHMM > 0.5
  calls.thca$thr.lnHMM <- pvals.thca$prob.lnHMM > res.thca$thr.lnHMM
  calls.thca$point5.lnHMM <- pvals.thca$prob.lnHMM > 0.5
  calls.luad$thr.expHMM <- pvals.luad$prob.expHMM > res.luad$thr.expHMM
  calls.luad$point5.expHMM <- pvals.luad$prob.expHMM > 0.5
  calls.luad$thr.lnHMM <- pvals.luad$prob.lnHMM > res.luad$thr.lnHMM
  calls.luad$point5.lnHMM <- pvals.luad$prob.lnHMM > 0.5
  calls.lihc$thr.expHMM <- pvals.lihc$prob.expHMM > res.lihc$thr.expHMM
  calls.lihc$point5.expHMM <- pvals.lihc$prob.expHMM > 0.5
  calls.lihc$thr.lnHMM <- pvals.lihc$prob.lnHMM > res.lihc$thr.lnHMM
  calls.lihc$point5.lnHMM <- pvals.lihc$prob.lnHMM > 0.5
  calls.lihc_long.chain$thr.expHMM <- pvals.lihc_long.chain$prob.expHMM > res.lihc_long.chain$thr.expHMM
  calls.lihc_long.chain$point5.expHMM <- pvals.lihc_long.chain$prob.expHMM > 0.5
  calls.lihc_long.chain$thr.lnHMM <- pvals.lihc_long.chain$prob.lnHMM > res.lihc_long.chain$thr.lnHMM
  calls.lihc_long.chain$point5.lnHMM <- pvals.lihc_long.chain$prob.lnHMM > 0.5
  calls.lusc$thr.expHMM <- pvals.lusc$prob.expHMM > res.lusc$thr.expHMM
  calls.lusc$point5.expHMM <- pvals.lusc$prob.expHMM > 0.5
  calls.lusc$thr.lnHMM <- pvals.lusc$prob.lnHMM > res.lusc$thr.lnHMM
  calls.lusc$point5.lnHMM <- pvals.lusc$prob.lnHMM > 0.5
  calls.prad$thr.expHMM <- pvals.prad$prob.expHMM > res.prad$thr.expHMM
  calls.prad$point5.expHMM <- pvals.prad$prob.expHMM > 0.5
  calls.prad$thr.lnHMM <- pvals.prad$prob.lnHMM > res.prad$thr.lnHMM
  calls.prad$point5.lnHMM <- pvals.prad$prob.lnHMM > 0.5
  calls.coad$thr.expHMM <- pvals.coad$prob.expHMM > res.coad$thr.expHMM
  calls.coad$point5.expHMM <- pvals.coad$prob.expHMM > 0.5
  calls.coad$thr.lnHMM <- pvals.coad$prob.lnHMM > res.coad$thr.lnHMM
  calls.coad$point5.lnHMM <- pvals.coad$prob.lnHMM > 0.5
  calls.brca$voom_lnHM.log <- p.adjust(pvals.brca$voom_lnHM.log, method="BH") < 0.05
  calls.brca$voom_MDSeq.zi <- p.adjust(pvals.brca$voom_MDSeq.zi, method="BH") < 0.05
  calls.brca$lnHM.log_lnHM.log <- p.adjust(pvals.brca$lnHM.log_lnHM.log, method="BH") < 0.05
  calls.kirc$voom_lnHM.log <- p.adjust(pvals.kirc$voom_lnHM.log, method="BH") < 0.05
  calls.kirc$voom_MDSeq.zi <- p.adjust(pvals.kirc$voom_MDSeq.zi, method="BH") < 0.05
  calls.kirc$lnHM.log_lnHM.log <- p.adjust(pvals.kirc$lnHM.log_lnHM.log, method="BH") < 0.05
  calls.thca$voom_lnHM.log <- p.adjust(pvals.thca$voom_lnHM.log, method="BH") < 0.05
  calls.thca$voom_MDSeq.zi <- p.adjust(pvals.thca$voom_MDSeq.zi, method="BH") < 0.05
  calls.thca$lnHM.log_lnHM.log <- p.adjust(pvals.thca$lnHM.log_lnHM.log, method="BH") < 0.05
  calls.luad$voom_lnHM.log <- p.adjust(pvals.luad$voom_lnHM.log, method="BH") < 0.05
  calls.luad$voom_MDSeq.zi <- p.adjust(pvals.luad$voom_MDSeq.zi, method="BH") < 0.05
  calls.luad$lnHM.log_lnHM.log <- p.adjust(pvals.luad$lnHM.log_lnHM.log, method="BH") < 0.05
  calls.lihc$voom_lnHM.log <- p.adjust(pvals.lihc$voom_lnHM.log, method="BH") < 0.05
  calls.lihc$voom_MDSeq.zi <- p.adjust(pvals.lihc$voom_MDSeq.zi, method="BH") < 0.05
  calls.lihc$lnHM.log_lnHM.log <- p.adjust(pvals.lihc$lnHM.log_lnHM.log, method="BH") < 0.05
  calls.lihc_long.chain$voom_lnHM.log <- p.adjust(pvals.lihc_long.chain$voom_lnHM.log, method="BH") < 0.05
  calls.lihc_long.chain$voom_MDSeq.zi <- p.adjust(pvals.lihc_long.chain$voom_MDSeq.zi, method="BH") < 0.05
  calls.lihc_long.chain$lnHM.log_lnHM.log <- p.adjust(pvals.lihc_long.chain$lnHM.log_lnHM.log, method="BH") < 0.05
  calls.lusc$voom_lnHM.log <- p.adjust(pvals.lusc$voom_lnHM.log, method="BH") < 0.05
  calls.lusc$voom_MDSeq.zi <- p.adjust(pvals.lusc$voom_MDSeq.zi, method="BH") < 0.05
  calls.lusc$lnHM.log_lnHM.log <- p.adjust(pvals.lusc$lnHM.log_lnHM.log, method="BH") < 0.05
  calls.prad$voom_lnHM.log <- p.adjust(pvals.prad$voom_lnHM.log, method="BH") < 0.05
  calls.prad$voom_MDSeq.zi <- p.adjust(pvals.prad$voom_MDSeq.zi, method="BH") < 0.05
  calls.prad$lnHM.log_lnHM.log <- p.adjust(pvals.prad$lnHM.log_lnHM.log, method="BH") < 0.05
  calls.coad$voom_lnHM.log <- p.adjust(pvals.coad$voom_lnHM.log, method="BH") < 0.05
  calls.coad$voom_MDSeq.zi <- p.adjust(pvals.coad$voom_MDSeq.zi, method="BH") < 0.05
  calls.coad$lnHM.log_lnHM.log <- p.adjust(pvals.coad$lnHM.log_lnHM.log, method="BH") < 0.05
}
for (i in c("brca", "kirc", "thca", "luad", "lihc", 
            "lihc_long.chain", "lusc", "prad", "coad")) {
  rm(list=c(paste0("res.", i)))
  write.csv(get(paste0("calls.", i)), 
            here("Results/TCGA paired data results May 2020", 
                 paste0("calls.", i, ".csv")), 
            row.names=row.names(get(paste0("pvals.", i))))
  write.csv(get(paste0("pvals.", i)), 
            here("Results/TCGA paired data results May 2020", 
                 paste0("pvals.", i, ".csv")))
}

folder <- "Data sources/Cancer-related genes"
for (i in c("brca", "kirc", "thca", "luad", "lihc", 
            "lusc", "prad", "coad")) {
  assign(paste0("genes.", i), 
         readRDS(here(folder, paste0(i, "_genes_info.rds"))))
  assign(paste0("dup.genes.", i), 
         readRDS(here(folder, paste0(i, "_genes_at_least_two_dbs_info.rds"))))
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
               "fisher.tests.weak_known.genes.csv"), 
)
write.csv(fisher.tests.strong, 
          here("Results/TCGA paired data results May 2020", 
               "fisher.tests.strong_known.genes.csv"), 
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
               "ranksum.tests.weak_known.genes.csv"), 
)
write.csv(ranksum.tests.strong, 
          here("Results/TCGA paired data results May 2020", 
               "ranksum.tests.strong_known.genes.csv"), 
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
  median((rank(pvals.kirc$p.disp.MDSeq.zi) - rank(pvals.kirc$p.disp.lnHM.log))[kirc_related_weak == T])) # MD better
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
  median((rank(pvals.luad$p.voom) - rank(pvals.luad$lnHM.log_lnHM.log))[luad_related_weak == T])) # No difference
c(wilcox.test(rank(pvals.lihc$p.voom)[lihc_related_weak == T], 
              rank(pvals.lihc$lnHM.log_lnHM.log)[lihc_related_weak == T], paired=T)$p.val, 
  median((rank(pvals.lihc$p.voom) - rank(pvals.lihc$lnHM.log_lnHM.log))[lihc_related_weak == T])) # voom better
c(wilcox.test(rank(pvals.lihc_long.chain$p.voom)[lihc_related_weak == T], 
              rank(pvals.lihc_long.chain$lnHM.log_lnHM.log)[lihc_related_weak == T], paired=T)$p.val, 
  median((rank(pvals.lihc_long.chain$p.voom) - 
            rank(pvals.lihc_long.chain$lnHM.log_lnHM.log))[lihc_related_weak == T])) # voom better
c(wilcox.test(rank(pvals.lusc$p.voom)[lusc_related_weak == T], 
              rank(pvals.lusc$lnHM.log_lnHM.log)[lusc_related_weak == T], paired=T)$p.val, 
  median((rank(pvals.lusc$p.voom) - rank(pvals.lusc$lnHM.log_lnHM.log))[lusc_related_weak == T])) # No difference
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
  median((rank(pvals.luad$p.voom) - rank(1 - pvals.luad$prob.lnHM))[luad_related_weak == T])) # No difference
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
                 paste0(j, "_HM_mean_v_disp_known.genes.csv")), 
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
                 paste0(j, "_MDSeq_mean_v_disp_known.genes.csv")), 
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
                 paste0(j, "_HM_voom_mean_v_disp_known.genes.csv")), 
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
                 paste0(j, "_HMM_voom_mean_v_disp_known.genes.csv")), 
            row.names=F)
  rm(pvals, res)
}


par(mfrow=c(2,2), mar=c(4,2,2,1), mgp=c(1.5,0.5,0))
plot(lihc.HM$pvals, lihc.HM$cor, type='l', ylim=c(-1,1), col="red", 
     main="Correlation between ranks of genes identified by differential expression and dispersion", 
     xlab="p-value threshold", ylab="")
lines(lihc.MD$pvals, lihc.MD$cor, col="blue")
plot(lihc.HM$pvals, lihc.HM$p.cor, type='l', ylim=c(0,1), col="red", 
     main="p-value of correlation test", 
     xlab="p-value threshold for differential expression or dispersion", ylab="")
lines(lihc.MD$pvals, lihc.MD$p.cor, col="blue")
plot(lihc.HM$pvals, lihc.HM$overlap / lihc.HM$union, type='l', ylim=c(0,1), col="red", 
     main="Union / Intersection", 
     xlab="p-value threshold", ylab="")
lines(lihc.MD$pvals, lihc.MD$overlap / lihc.MD$union, col='blue', 
      xlab="p-value threshold", ylab="")
plot(lihc.HM$pvals, lihc.HM$p.hyper, type='l', ylim=c(0,1), col="red", 
     main="p-value of test of overlap", 
     xlab="p-value threshold", ylab="")
lines(lihc.MD$pvals, lihc.MD$p.hyper, col='blue')


# genes identified by diff disp/dist missed by DE
par(mfrow=c(2,2), mar=c(4,2,2,1), mgp=c(1.5,0.5,0))
plot(lihc.HM$pvals, lihc.HM$TPA, type='l', col='red')
lines(lihc.HM$pvals, lihc.HM$TPB, col='blue')
lines(lihc.HM$pvals, lihc.HM$TPboth, col='grey')
lines(lihc.HM$pvals, lihc.HM$TPAonly, col='red', lty=2)
lines(lihc.HM$pvals, lihc.HM$TPBonly, col='blue', lty=2)
plot(lihc.MD$pvals, lihc.MD$TPA, type='l', col='red')
lines(lihc.MD$pvals, lihc.MD$TPB, col='blue')
lines(lihc.MD$pvals, lihc.MD$TPboth, col='grey')
lines(lihc.MD$pvals, lihc.MD$TPAonly, col='red', lty=2)
lines(lihc.MD$pvals, lihc.MD$TPBonly, col='blue', lty=2)
plot(lihc.HMvoom$pvals, lihc.HMvoom$TPA, type='l', col='red', 
     main="Genes identified by differential expression (red) and differential dispersion (blue)", 
     ylab="", xlab="p-value threshold")
lines(lihc.HMvoom$pvals, lihc.HMvoom$TPB, col='blue')
lines(lihc.HMvoom$pvals, lihc.HMvoom$TPboth, col='grey')
lines(lihc.HMvoom$pvals, lihc.HMvoom$TPAonly, col='red', lty=2)
lines(lihc.HMvoom$pvals, lihc.HMvoom$TPBonly, col='blue', lty=2)
plot(lihc.HMMvoom$pvals, lihc.HMMvoom$TPA, type='l', col='red', 
     main="Genes identified by differential expression (red) and differential distribution (blue)", 
     ylab="", xlab="p-value threshold")
lines(lihc.HMMvoom$pvals, lihc.HMMvoom$TPB, col='blue')
lines(lihc.HMMvoom$pvals, lihc.HMMvoom$TPboth, col='grey')
lines(lihc.HMMvoom$pvals, lihc.HMMvoom$TPAonly, col='red', lty=2)
lines(lihc.HMMvoom$pvals, lihc.HMMvoom$TPBonly, col='blue', lty=2)





