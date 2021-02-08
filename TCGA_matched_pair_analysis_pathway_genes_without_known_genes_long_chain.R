library(here)

#### Original long chain (3500 iterations) ####
for (i in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {
  assign(paste0("calls.", i), 
         read.csv(here("Results/TCGA long chain results Aug 2020", 
                       paste0("calls.", i, ".csv")), 
                  stringsAsFactors=F, 
                  header=T, 
                  row.names=1))
  assign(paste0("pvals.", i), 
         read.csv(here("Results/TCGA long chain results Aug 2020", 
                       paste0("pvals.", i, ".csv")), 
                  stringsAsFactors=F, 
                  header=T, 
                  row.names=1))
}

folder <- "Data sources/Cancer-related pathway genes"
for (i in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {
  assign(paste0("genes.", i), 
         read.csv(here(folder, paste0(i, "_pathway_genes_recount_without_known_genes.csv")), stringsAsFactors=F))
  assign(paste0(i, "_related"), 
         rownames(get(paste0("pvals.", i))) %in% get(paste0("genes.", i))$ENSEMBL)
}
rm(i)

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
          here("Results/TCGA long chain results Aug 2020", 
               "fisher.tests_pathway.genes.minus.known.genes.csv"), 
)

## Test ability to rank cancer-related genes above non-cancer-related genes using Wilcoxon rank-sum test ####
tests.ranksum <- names(pvals.brca)
for (j in tests.ranksum) {
  for (i in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {
    if(grepl("prob", j)) {
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
          here("Results/TCGA long chain results Aug 2020", 
               "ranksum.tests_pathway.genes.minus.known.genes.csv"), 
)


#### Longer chain (two times 3500 iterations) ####
for (i in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {
  assign(paste0("calls.", i), 
         read.csv(here("Results/TCGA long chain results Aug 2020", 
                       paste0("calls.", i, "_20k.csv")), 
                  stringsAsFactors=F, 
                  header=T, 
                  row.names=1))
  assign(paste0("pvals.", i), 
         read.csv(here("Results/TCGA long chain results Aug 2020", 
                       paste0("pvals.", i, "_20k.csv")), 
                  stringsAsFactors=F, 
                  header=T, 
                  row.names=1))
}

folder <- "Data sources/Cancer-related pathway genes"
for (i in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {
  assign(paste0("genes.", i), 
         read.csv(here(folder, paste0(i, "_pathway_genes_recount_without_known_genes.csv")), stringsAsFactors=F))
  assign(paste0(i, "_related"), 
         rownames(get(paste0("pvals.", i))) %in% get(paste0("genes.", i))$ENSEMBL)
}
rm(i)

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
          here("Results/TCGA long chain results Aug 2020", 
               "fisher.tests_pathway.genes.minus.known.genes_20k.csv"), 
)

## Test ability to rank cancer-related genes above non-cancer-related genes using Wilcoxon rank-sum test ####
tests.ranksum <- names(pvals.brca)
for (j in tests.ranksum) {
  for (i in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {
    if(grepl("prob", j)) {
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
          here("Results/TCGA long chain results Aug 2020", 
               "ranksum.tests_pathway.genes.minus.known.genes_20k.csv"), 
)
