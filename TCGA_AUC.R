library(here)
library(ROCR)

for (i in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {
  assign(paste0("pvals.", i), 
         read.csv(here("Results/TCGA paired data results May 2020", 
                       paste0("pvals.", i, ".csv")), row.names=1))
}

for (i in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {
  assign(paste0("genes.", i), 
         readRDS(here("Data sources/Cancer-related genes", paste0(i, "_genes_info.rds"))))
  assign(paste0(i, "_related"), 
         rownames(get(paste0("pvals.", i))) %in% get(paste0("genes.", i))$ENSEMBL)
  rm(list=paste0("genes.", i))
}

aucs <- as.data.frame(matrix(nrow=25, ncol=8), row.names=names(pvals.brca))
names(aucs) <- c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")

for (i in names(aucs)) {
  for (j in rownames(aucs)) {
    NAs <- which(is.na(get(paste0("pvals.", i))[[j]]))
    if (length(NAs) > 0) {
      pred <- prediction(1 - get(paste0("pvals.", i))[[j]][-NAs], get(paste0(i, "_related"))[-NAs])
    }
    else {
      pred <- prediction(1 - get(paste0("pvals.", i))[[j]], get(paste0(i, "_related")))
    }
    if (grepl("prob", j)) {
      aucs[i][j,] <- 1 - performance(pred, "auc")@y.values[[1]]
    }
    else {
      aucs[i][j,] <- performance(pred, "auc")@y.values[[1]]
    }
  }
}
boxplot(t(aucs))
rowMeans(aucs)
write.csv(aucs, file=here("Results/TCGA paired data results May 2020", "AUCs.weak_known_genes.csv"))


# Long chains (original ones, 3500 runs, August 2020)
for (i in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {
  assign(paste0("pvals.", i), 
         read.csv(here("Results/TCGA long chain results Aug 2020", 
                       paste0("pvals.", i, ".csv")), row.names=1))
}

for (i in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {
  assign(paste0("genes.", i), 
         readRDS(here("Data sources/Cancer-related genes", paste0(i, "_genes_info.rds"))))
  assign(paste0(i, "_related"), 
         rownames(get(paste0("pvals.", i))) %in% get(paste0("genes.", i))$ENSEMBL)
  rm(list=paste0("genes.", i))
}

aucs <- as.data.frame(matrix(nrow=11, ncol=8), row.names=names(pvals.brca))
names(aucs) <- c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")

for (i in names(aucs)) {
  for (j in rownames(aucs)) {
    NAs <- which(is.na(get(paste0("pvals.", i))[[j]]))
    if (length(NAs) > 0) {
      pred <- prediction(1 - get(paste0("pvals.", i))[[j]][-NAs], get(paste0(i, "_related"))[-NAs])
    }
    else {
      pred <- prediction(1 - get(paste0("pvals.", i))[[j]], get(paste0(i, "_related")))
    }
    if (grepl("prob", j)) {
      aucs[i][j,] <- 1 - performance(pred, "auc")@y.values[[1]]
    }
    else {
      aucs[i][j,] <- performance(pred, "auc")@y.values[[1]]
    }
  }
}
boxplot(t(aucs))
rowMeans(aucs)
write.csv(aucs, file=here("Results/TCGA long chain results Aug 2020", "AUCs.weak_known_genes.csv"))


# Longer chains (two times 3500 runs, September 2020)
for (i in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {
  assign(paste0("pvals.", i), 
         read.csv(here("Results/TCGA long chain results Aug 2020", 
                       paste0("pvals.", i, "_20k.csv")), row.names=1))
}

for (i in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {
  assign(paste0("genes.", i), 
         readRDS(here("Data sources/Cancer-related genes", paste0(i, "_genes_info.rds"))))
  assign(paste0(i, "_related"), 
         rownames(get(paste0("pvals.", i))) %in% get(paste0("genes.", i))$ENSEMBL)
  rm(list=paste0("genes.", i))
}

aucs <- as.data.frame(matrix(nrow=11, ncol=8), row.names=names(pvals.brca))
names(aucs) <- c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")

for (i in names(aucs)) {
  for (j in rownames(aucs)) {
    NAs <- which(is.na(get(paste0("pvals.", i))[[j]]))
    if (length(NAs) > 0) {
      pred <- prediction(1 - get(paste0("pvals.", i))[[j]][-NAs], get(paste0(i, "_related"))[-NAs])
    }
    else {
      pred <- prediction(1 - get(paste0("pvals.", i))[[j]], get(paste0(i, "_related")))
    }
    if (grepl("prob", j)) {
      aucs[i][j,] <- 1 - performance(pred, "auc")@y.values[[1]]
    }
    else {
      aucs[i][j,] <- performance(pred, "auc")@y.values[[1]]
    }
  }
}
boxplot(t(aucs))
rowMeans(aucs)
write.csv(aucs, file=here("Results/TCGA long chain results Aug 2020", "AUCs.weak_known_genes_20k.csv"))

