library(here)

folder <- "Results/TCGA long chain results Aug 2020"
for (i in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {
  for (j in c("exp", "ln")) {
    assign(paste0("res.", i, "_", j), 
           readRDS(here(folder, paste0("results.TCGA_20k_iter_", i, "_", j, "HM.rds"))))
    assign(paste0("pvals.", i, "_", j), 
           as.data.frame(
             get(paste0("res.", i, "_", j))[c(1,5,6,10,11)]
           ))
    assign(paste0("qvals.", i, "_", j), 
           as.data.frame(
             get(paste0("res.", i, "_", j))[c(4,7,8,12,13)]
           ))
  }
}
rm(j)

names(pvals.brca_exp) <- paste0(names(pvals.brca_exp), ".expHM")
names(pvals.kirc_exp) <- paste0(names(pvals.kirc_exp), ".expHM")
names(pvals.thca_exp) <- paste0(names(pvals.thca_exp), ".expHM")
names(pvals.luad_exp) <- paste0(names(pvals.luad_exp), ".expHM")
names(pvals.lihc_exp) <- paste0(names(pvals.lihc_exp), ".expHM")
names(pvals.lusc_exp) <- paste0(names(pvals.lusc_exp), ".expHM")
names(pvals.prad_exp) <- paste0(names(pvals.prad_exp), ".expHM")
names(pvals.coad_exp) <- paste0(names(pvals.coad_exp), ".expHM")
names(pvals.brca_ln) <- paste0(names(pvals.brca_ln), ".lnHM")
names(pvals.kirc_ln) <- paste0(names(pvals.kirc_ln), ".lnHM")
names(pvals.thca_ln) <- paste0(names(pvals.thca_ln), ".lnHM")
names(pvals.luad_ln) <- paste0(names(pvals.luad_ln), ".lnHM")
names(pvals.lihc_ln) <- paste0(names(pvals.lihc_ln), ".lnHM")
names(pvals.lusc_ln) <- paste0(names(pvals.lusc_ln), ".lnHM")
names(pvals.prad_ln) <- paste0(names(pvals.prad_ln), ".lnHM")
names(pvals.coad_ln) <- paste0(names(pvals.coad_ln), ".lnHM")
pvals.brca <- cbind(pvals.brca_exp, pvals.brca_ln)
pvals.kirc <- cbind(pvals.kirc_exp, pvals.kirc_ln)
pvals.thca <- cbind(pvals.thca_exp, pvals.thca_ln)
pvals.luad <- cbind(pvals.luad_exp, pvals.luad_ln)
pvals.lihc <- cbind(pvals.lihc_exp, pvals.lihc_ln)
pvals.lusc <- cbind(pvals.lusc_exp, pvals.lusc_ln)
pvals.prad <- cbind(pvals.prad_exp, pvals.prad_ln)
pvals.coad <- cbind(pvals.coad_exp, pvals.coad_ln)
rm(pvals.brca_exp, pvals.kirc_exp, pvals.thca_exp, pvals.luad_exp, 
   pvals.lihc_exp, pvals.lusc_exp, pvals.prad_exp, pvals.coad_exp, 
   pvals.brca_ln, pvals.kirc_ln, pvals.thca_ln, pvals.luad_ln, 
   pvals.lihc_ln, pvals.lusc_ln, pvals.prad_ln, pvals.coad_ln)

names(qvals.brca_exp) <- paste0(names(qvals.brca_exp), ".expHM")
names(qvals.kirc_exp) <- paste0(names(qvals.kirc_exp), ".expHM")
names(qvals.thca_exp) <- paste0(names(qvals.thca_exp), ".expHM")
names(qvals.luad_exp) <- paste0(names(qvals.luad_exp), ".expHM")
names(qvals.lihc_exp) <- paste0(names(qvals.lihc_exp), ".expHM")
names(qvals.lusc_exp) <- paste0(names(qvals.lusc_exp), ".expHM")
names(qvals.prad_exp) <- paste0(names(qvals.prad_exp), ".expHM")
names(qvals.coad_exp) <- paste0(names(qvals.coad_exp), ".expHM")
names(qvals.brca_ln) <- paste0(names(qvals.brca_ln), ".lnHM")
names(qvals.kirc_ln) <- paste0(names(qvals.kirc_ln), ".lnHM")
names(qvals.thca_ln) <- paste0(names(qvals.thca_ln), ".lnHM")
names(qvals.luad_ln) <- paste0(names(qvals.luad_ln), ".lnHM")
names(qvals.lihc_ln) <- paste0(names(qvals.lihc_ln), ".lnHM")
names(qvals.lusc_ln) <- paste0(names(qvals.lusc_ln), ".lnHM")
names(qvals.prad_ln) <- paste0(names(qvals.prad_ln), ".lnHM")
names(qvals.coad_ln) <- paste0(names(qvals.coad_ln), ".lnHM")
qvals.brca <- cbind(qvals.brca_exp, qvals.brca_ln)
qvals.kirc <- cbind(qvals.kirc_exp, qvals.kirc_ln)
qvals.thca <- cbind(qvals.thca_exp, qvals.thca_ln)
qvals.luad <- cbind(qvals.luad_exp, qvals.luad_ln)
qvals.lihc <- cbind(qvals.lihc_exp, qvals.lihc_ln)
qvals.lusc <- cbind(qvals.lusc_exp, qvals.lusc_ln)
qvals.prad <- cbind(qvals.prad_exp, qvals.prad_ln)
qvals.coad <- cbind(qvals.coad_exp, qvals.coad_ln)
rm(qvals.brca_exp, qvals.kirc_exp, qvals.thca_exp, qvals.luad_exp, 
   qvals.lihc_exp, qvals.lusc_exp, qvals.prad_exp, qvals.coad_exp, 
   qvals.brca_ln, qvals.kirc_ln, qvals.thca_ln, qvals.luad_ln, 
   qvals.lihc_ln, qvals.lusc_ln, qvals.prad_ln, qvals.coad_ln)


for (i in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {
  assign(paste0("genes.", i), 
         rownames(
           readRDS(here("Results/TCGA paired data results May 2020", 
                        paste0("results.TCGA.paired_", i, ".rds")))$counts)
  )
}
rm(i)
rownames(pvals.brca) <- gsub('\\..*', '', genes.brca)
rownames(pvals.kirc) <- gsub('\\..*', '', genes.kirc)
rownames(pvals.thca) <- gsub('\\..*', '', genes.thca)
rownames(pvals.luad) <- gsub('\\..*', '', genes.luad)
rownames(pvals.lihc) <- gsub('\\..*', '', genes.lihc)
rownames(pvals.lusc) <- gsub('\\..*', '', genes.lusc)
rownames(pvals.prad) <- gsub('\\..*', '', genes.prad)
rownames(pvals.coad) <- gsub('\\..*', '', genes.coad)
rownames(qvals.brca) <- rownames(pvals.brca)
rownames(qvals.kirc) <- rownames(pvals.kirc)
rownames(qvals.thca) <- rownames(pvals.thca)
rownames(qvals.luad) <- rownames(pvals.luad)
rownames(qvals.lihc) <- rownames(pvals.lihc)
rownames(qvals.lusc) <- rownames(pvals.lusc)
rownames(qvals.prad) <- rownames(pvals.prad)
rownames(qvals.coad) <- rownames(pvals.coad)
rm(genes.brca, genes.kirc, genes.thca, genes.luad, genes.lihc, genes.lusc, genes.prad, genes.coad)

for (i in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {
  assign(paste0("calls.", i), 
         data.frame(
           get(paste0("qvals.", i)) < 0.05
         ))
  rm(list=paste0("qvals.", i))
}
rm(i)

# Will take a bit of messing around to import voom results to combine voom DE with HM DD, so for now just 
# combine HM with HM, and add in voom-HM later if I decide I need it.
{
  pvals.brca$log.lnHM_log.lnHM <- pmin(pvals.brca$p.mean.log.lnHM, pvals.brca$p.disp.log.lnHM)
  pvals.kirc$log.lnHM_log.lnHM <- pmin(pvals.kirc$p.mean.log.lnHM, pvals.kirc$p.disp.log.lnHM)
  pvals.thca$log.lnHM_log.lnHM <- pmin(pvals.thca$p.mean.log.lnHM, pvals.thca$p.disp.log.lnHM)
  pvals.luad$log.lnHM_log.lnHM <- pmin(pvals.luad$p.mean.log.lnHM, pvals.luad$p.disp.log.lnHM)
  pvals.lihc$log.lnHM_log.lnHM <- pmin(pvals.lihc$p.mean.log.lnHM, pvals.lihc$p.disp.log.lnHM)
  pvals.lusc$log.lnHM_log.lnHM <- pmin(pvals.lusc$p.mean.log.lnHM, pvals.lusc$p.disp.log.lnHM)
  pvals.prad$log.lnHM_log.lnHM <- pmin(pvals.prad$p.mean.log.lnHM, pvals.prad$p.disp.log.lnHM)
  pvals.coad$log.lnHM_log.lnHM <- pmin(pvals.coad$p.mean.log.lnHM, pvals.coad$p.disp.log.lnHM)
  calls.brca$thr.expHM <- pvals.brca$prob.expHM > res.brca_exp$thr
  calls.brca$point5.expHM <- pvals.brca$prob.expHM > 0.5
  calls.brca$thr.lnHM <- pvals.brca$prob.lnHM > res.brca_ln$thr
  calls.brca$point5.lnHM <- pvals.brca$prob.lnHM > 0.5
  calls.kirc$thr.expHM <- pvals.kirc$prob.expHM > res.kirc_exp$thr
  calls.kirc$point5.expHM <- pvals.kirc$prob.expHM > 0.5
  calls.kirc$thr.lnHM <- pvals.kirc$prob.lnHM > res.kirc_ln$thr
  calls.kirc$point5.lnHM <- pvals.kirc$prob.lnHM > 0.5
  calls.thca$thr.expHM <- pvals.thca$prob.expHM > res.thca_exp$thr
  calls.thca$point5.expHM <- pvals.thca$prob.expHM > 0.5
  calls.thca$thr.lnHM <- pvals.thca$prob.lnHM > res.thca_ln$thr
  calls.thca$point5.lnHM <- pvals.thca$prob.lnHM > 0.5
  calls.luad$thr.expHM <- pvals.luad$prob.expHM > res.luad_exp$thr
  calls.luad$point5.expHM <- pvals.luad$prob.expHM > 0.5
  calls.luad$thr.lnHM <- pvals.luad$prob.lnHM > res.luad_ln$thr
  calls.luad$point5.lnHM <- pvals.luad$prob.lnHM > 0.5
  calls.lihc$thr.expHM <- pvals.lihc$prob.expHM > res.lihc_exp$thr
  calls.lihc$point5.expHM <- pvals.lihc$prob.expHM > 0.5
  calls.lihc$thr.lnHM <- pvals.lihc$prob.lnHM > res.lihc_ln$thr
  calls.lihc$point5.lnHM <- pvals.lihc$prob.lnHM > 0.5
  calls.lusc$thr.expHM <- pvals.lusc$prob.expHM > res.lusc_exp$thr
  calls.lusc$point5.expHM <- pvals.lusc$prob.expHM > 0.5
  calls.lusc$thr.lnHM <- pvals.lusc$prob.lnHM > res.lusc_ln$thr
  calls.lusc$point5.lnHM <- pvals.lusc$prob.lnHM > 0.5
  calls.prad$thr.expHM <- pvals.prad$prob.expHM > res.prad_exp$thr
  calls.prad$point5.expHM <- pvals.prad$prob.expHM > 0.5
  calls.prad$thr.lnHM <- pvals.prad$prob.lnHM > res.prad_ln$thr
  calls.prad$point5.lnHM <- pvals.prad$prob.lnHM > 0.5
  calls.coad$thr.expHM <- pvals.coad$prob.expHM > res.coad_exp$thr
  calls.coad$point5.expHM <- pvals.coad$prob.expHM > 0.5
  calls.coad$thr.lnHM <- pvals.coad$prob.lnHM > res.coad_ln$thr
  calls.coad$point5.lnHM <- pvals.coad$prob.lnHM > 0.5
  calls.brca$log.lnHM_log.lnHM <- p.adjust(pvals.brca$log.lnHM_log.lnHM, method="BH") < 0.05
  calls.kirc$log.lnHM_log.lnHM <- p.adjust(pvals.kirc$log.lnHM_log.lnHM, method="BH") < 0.05
  calls.thca$log.lnHM_log.lnHM <- p.adjust(pvals.thca$log.lnHM_log.lnHM, method="BH") < 0.05
  calls.luad$log.lnHM_log.lnHM <- p.adjust(pvals.luad$log.lnHM_log.lnHM, method="BH") < 0.05
  calls.lihc$log.lnHM_log.lnHM <- p.adjust(pvals.lihc$log.lnHM_log.lnHM, method="BH") < 0.05
  calls.lusc$log.lnHM_log.lnHM <- p.adjust(pvals.lusc$log.lnHM_log.lnHM, method="BH") < 0.05
  calls.prad$log.lnHM_log.lnHM <- p.adjust(pvals.prad$log.lnHM_log.lnHM, method="BH") < 0.05
  calls.coad$log.lnHM_log.lnHM <- p.adjust(pvals.coad$log.lnHM_log.lnHM, method="BH") < 0.05
}
for (i in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {
  rm(list=c(paste0("res.", i, "_exp")))
  rm(list=c(paste0("res.", i, "_ln")))
  write.csv(get(paste0("calls.", i)), 
            here(folder, paste0("calls.", i, "_20k.csv")), 
            row.names=row.names(get(paste0("pvals.", i))))
  write.csv(get(paste0("pvals.", i)), 
            here(folder, paste0("pvals.", i, "_20k.csv")))
}
rm(i)

folder <- "Data sources/Cancer-related genes"
for (i in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {
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


## Test ability to identify cancer-related genes with threshold using one-sided Fisher's exact test ####
tests.fisher <- names(calls.brca)
for (j in tests.fisher) {
  for (i in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {
    assign(paste0("fisher.weak_", i, "_", j), 
           fisher.test(get(paste0(i, "_related_weak")), 
                       get(paste0("calls.", i))[[j]], 
                       alternative="greater")$p.val)
    assign(paste0("fisher.strong_", i, "_", j), 
           fisher.test(get(paste0(i, "_related_strong")), 
                       get(paste0("calls.", i))[[j]], 
                       alternative="greater")$p.val)
  }
}
rm(i,j)

fisher.tests.weak <- data.frame(matrix(nrow=length(tests.fisher), ncol=8))
names(fisher.tests.weak) <- c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")
rownames(fisher.tests.weak) <- tests.fisher
for (j in seq_along(rownames(fisher.tests.weak))) {
  for (i in names(fisher.tests.weak)) {
    fisher.tests.weak[[i]][j] <- get(paste0("fisher.weak_", i, "_", tests.fisher[j]))
    rm(list=paste0("fisher.weak_", i, "_", tests.fisher[j]))
  }
}

fisher.tests.strong <- data.frame(matrix(nrow=length(tests.fisher), ncol=8))
names(fisher.tests.strong) <- c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")
rownames(fisher.tests.strong) <- tests.fisher
for (j in seq_along(rownames(fisher.tests.strong))) {
  for (i in names(fisher.tests.strong)) {
    fisher.tests.strong[[i]][j] <- get(paste0("fisher.strong_", i, "_", tests.fisher[j]))
    rm(list=paste0("fisher.strong_", i, "_", tests.fisher[j]))
  }
}
rm(i,j,tests.fisher)

folder <- "Results/TCGA long chain results Aug 2020"
write.csv(fisher.tests.weak, here(folder, "fisher.tests.weak_known.genes_20k.csv"))
write.csv(fisher.tests.strong, here(folder, "fisher.tests.strong_known.genes_20k.csv"))


## Test ability to rank cancer-related genes above non-cancer-related genes using Wilcoxon rank-sum test ####
tests.ranksum <- names(pvals.brca)
for (j in tests.ranksum) {
  for (i in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {
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
}
rm(j)

ranksum.tests.weak <- data.frame(matrix(nrow=length(tests.ranksum), ncol=8))
names(ranksum.tests.weak) <- c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")
rownames(ranksum.tests.weak) <- tests.ranksum
for (j in seq_along(rownames(ranksum.tests.weak))) {
  for (i in names(ranksum.tests.weak)) {
    ranksum.tests.weak[[i]][j] <- get(paste0("ranksum.weak_", i, "_", tests.ranksum[j]))
    rm(list=paste0("ranksum.weak_", i, "_", tests.ranksum[j]))
  }
}

ranksum.tests.strong <- data.frame(matrix(nrow=length(tests.ranksum), ncol=8))
names(ranksum.tests.strong) <- c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")
rownames(ranksum.tests.strong) <- tests.ranksum
for (j in seq_along(rownames(ranksum.tests.strong))) {
  for (i in names(ranksum.tests.strong)) {
    ranksum.tests.strong[[i]][j] <- get(paste0("ranksum.strong_", i, "_", tests.ranksum[j]))
    rm(list=paste0("ranksum.strong_", i, "_", tests.ranksum[j]))
  }
}
rm(i,j,tests.ranksum)

write.csv(ranksum.tests.weak, here(folder, "ranksum.tests.weak_known.genes_20k.csv"))
write.csv(ranksum.tests.strong, here(folder, "ranksum.tests.strong_known.genes_20k.csv"))


