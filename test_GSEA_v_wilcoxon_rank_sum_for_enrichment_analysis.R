# Initially used WRS to test for each method's ability to rank ca-related genes above non-related genes, 
# using p-values. Then used GSEA to test for enrichment of gene sets (defined by GO terms), using 
# absolute value of log fold change multiplied by -log10 p-values. But these two tests are essentially 
# looking to do the same thing, the only real difference being that ca-related genes are testing a 
# single, specifically created gene set. So, should I use the same method for both? And if so, which 
# one?
# 
# WRS can handle ties (when using normal approximation for null, which is default unless less than 50 
# samples or exact=T), but GSEA just takes ranked list in order, so to deal with ties would need to 
# permute order of tied genes and run many times, so using p-values isn't really an option for GSEA for 
# HMs.


## Related but partly separate - look at correlation among cancer-related genes
# Mean pairwise correlation for each set of cancer-related genes, across whole datasets and by group.
library(here)
library(recount)
library(dplyr)

# First get counts and groups
brca <- readRDS(here("recount data/TCGA", "brca_matched_pairs.rds"))
brca <- brca[, brca$cgc_case_histological_diagnosis == "Infiltrating Ductal Carcinoma"]
brca <- brca[which(rowMeans(apply(assay(brca), 2, function(x) x*1e6/sum(x))) > 0.5), ]
kirc <- readRDS(here("recount data/TCGA", "kirc_matched_pairs.rds"))
kirc <- kirc[which(rowMeans(apply(assay(kirc), 2, function(x) x*1e6/sum(x))) > 0.5), ]
thca <- readRDS(here("recount data/TCGA", "thca_matched_pairs.rds"))
thca <- thca[, thca$cgc_case_histological_diagnosis == "Thyroid Papillary Carcinoma - Classical/usual"]
thca <- thca[which(rowMeans(apply(assay(thca), 2, function(x) x*1e6/sum(x))) > 0.5), ]
luad <- readRDS(here("recount data/TCGA", "luad_matched_pairs.rds"))
luad <- luad[, luad$cgc_case_histological_diagnosis == "Lung Adenocarcinoma- Not Otherwise Specified (NOS)"]
luad <- luad[which(rowMeans(apply(assay(luad), 2, function(x) x*1e6/sum(x))) > 0.5), ]
lihc <- readRDS(here("recount data/TCGA", "lihc_matched_pairs.rds"))
lihc <- lihc[, lihc$cgc_case_histological_diagnosis == "Hepatocellular Carcinoma"]
lihc <- lihc[which(rowMeans(apply(assay(lihc), 2, function(x) x*1e6/sum(x))) > 0.5), ]
lusc <- readRDS(here("recount data/TCGA", "lusc_matched_pairs.rds"))
lusc <- lusc[, lusc$cgc_case_histological_diagnosis == "Lung Squamous Cell Carcinoma- Not Otherwise Specified (NOS)" |
             is.na(lusc$cgc_case_histological_diagnosis)]
lusc <- lusc[which(rowMeans(apply(assay(lusc), 2, function(x) x*1e6/sum(x))) > 0.5), ]
prad <- readRDS(here("recount data/TCGA", "prad_matched_pairs.rds"))
prad <- prad[which(rowMeans(apply(assay(prad), 2, function(x) x*1e6/sum(x))) > 0.5), ]
coad <- readRDS(here("recount data/TCGA", "coad_matched_pairs.rds"))
coad <- coad[, coad$cgc_case_histological_diagnosis == "Colon Adenocarcinoma"] # coad
coad <- coad[which(rowMeans(apply(assay(coad), 2, function(x) x*1e6/sum(x))) > 0.5), ]
for (i in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {
  group <- get(i)$gdc_cases.samples.sample_type
  group[group == "Solid Tissue Normal"] <- 1
  group[group == "Primary Tumor"] <- 2
  group <- factor(group)
  assign(paste0("group_", i), group)
  assign(i, assay(get(i)))
}

# Next get gene lists
for (i in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {
  folder <- "Data sources/Cancer-related genes"
  assign(paste0("genes_", i), 
         readRDS(here(folder, paste0(i, "_genes_info.rds"))))
  assign(paste0(i, "_related"), 
         gsub('\\..*', '', rownames(get(i))) %in% get(paste0("genes_", i))$ENSEMBL)
  
  folder <- "Data sources/Cancer-related pathway genes"
  assign(paste0("pathway_genes_", i), 
         read.csv(here(folder, paste0(i, "_pathway_genes_recount_without_known_genes.csv")), stringsAsFactors=F))
  assign(paste0(i, "_pathway"), 
         gsub('\\..*', '', rownames(get(i))) %in% get(paste0("pathway_genes_", i))$ENSEMBL)
}

# Finally get average correlations for related and pathway genes, overall and per group
correlations <- matrix(rep(0, 6*8), nrow=6, ncol=8)
rownames(correlations) <- c("overall_related", "normal_related", "tumour_related", 
                            "overall_pathway", "normal_pathway", "tumour_pathway")
colnames(correlations) <- c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")
for (i in 1:8) {
  ca <- colnames(correlations)[i]
  overall_related <- get(ca)[get(paste0(ca, "_related")), ]
  normal_related <- overall_related[, get(paste0("group_", ca)) == 1]
  tumour_related <- overall_related[, get(paste0("group_", ca)) == 2]
  overall_pathway <- get(ca)[get(paste0(ca, "_pathway")), ]
  normal_pathway <- overall_pathway[, get(paste0("group_", ca)) == 1]
  tumour_pathway <- overall_pathway[, get(paste0("group_", ca)) == 2]
  temp_overall_related <- 0
  temp_normal_related <- 0
  temp_tumour_related <- 0
  for (j in 1:nrow(overall_related)) {
    for (k in 1:nrow(overall_related)) {
      if (j != k) {
        temp_overall_related <- temp_overall_related + cor(overall_related[i, ], overall_related[j, ]) / 
          (nrow(overall_related) * (nrow(overall_related) - 1))
        temp_normal_related <- temp_normal_related + cor(normal_related[i, ], normal_related[j, ]) / 
          (nrow(overall_related) * (nrow(overall_related) - 1))
        temp_tumour_related <- temp_tumour_related + cor(tumour_related[i, ], tumour_related[j, ]) / 
          (nrow(overall_related) * (nrow(overall_related) - 1))
      }
    }
  }
  temp_overall_pathway <- 0
  temp_normal_pathway <- 0
  temp_tumour_pathway <- 0
  for (j in 1:nrow(overall_pathway)) {
    for (k in 1:nrow(overall_pathway)) {
      if (j != k) {
        temp_overall_pathway <- temp_overall_pathway + cor(overall_pathway[i, ], overall_pathway[j, ]) / 
          (nrow(overall_pathway) * (nrow(overall_pathway) - 1))
        temp_normal_pathway <- temp_normal_pathway + cor(normal_pathway[i, ], normal_pathway[j, ]) / 
          (nrow(overall_pathway) * (nrow(overall_pathway) - 1))
        temp_tumour_pathway <- temp_tumour_pathway + cor(tumour_pathway[i, ], tumour_pathway[j, ]) / 
          (nrow(overall_pathway) * (nrow(overall_pathway) - 1))
      }
    }
  }
  correlations[, i] <- c(
    temp_overall_related, 
    temp_normal_related, 
    temp_tumour_related, 
    temp_overall_pathway, 
    temp_normal_pathway, 
    temp_tumour_pathway
  )
}
#                       brca       kirc       thca        luad      lihc        lusc       prad       coad
# overall_related 0.06952284 0.13103635 0.10651580 -0.05753162 0.2899658  0.10820625 0.25203261 0.12442384
# normal_related  0.15138350 0.26916280 0.16166441 -0.07700190 0.3011495  0.11147362 0.27860088 0.15235122
# tumour_related  0.06904193 0.06707197 0.07077850  0.04183525 0.2154612  0.09469193 0.27861215 0.17775270
# overall_pathway 0.04499471 0.15358180 0.06327812  0.09032268 0.3054245  0.03754718 0.11098752 0.13492643
# normal_pathway  0.14898788 0.23865890 0.09345190 -0.15938012 0.2772315 -0.02435031 0.13829530 0.06106758
# tumour_pathway  0.03551066 0.14032549 0.06592996  0.19531951 0.2289660  0.07478431 0.07394432 0.22505980
# Generally not much correlation, but quite high for lihc, and for prad and coad for related genes, and 
# kirc for pathway genes.


## Repeat WRS tests with genes ranked by absolute value of log fold change multiplied by -log10 p-value 
# For now, only for voom and lnHMdisp since that's all I have abslfcnlp for so far.
for (i in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {
  assign(paste0("abslfcnlp.voom_", i), 
         read.table(here("GSEA data", paste0("abslfcnlp.voom.", i, ".rnk")), header=F, row.names=1))
  assign(paste0("abslfcnlp.lnHMdisp_", i), 
         read.table(here("GSEA data", paste0("abslfcnlp.lnHMdisp.", i, ".rnk")), header=F, row.names=1))
}

for (i in c("brca", "kirc", "thca", "luad", "lihc", 
            "lusc", "prad", "coad")) {
  assign(paste0("ranksum.voom_", i, "_related"), 
         wilcox.test(
           rank(get(paste0("abslfcnlp.voom_", i)))[get(paste0(i, "_related")) == T], 
           rank(get(paste0("abslfcnlp.voom_", i)))[get(paste0(i, "_related")) == F], 
           paired=F, 
           alternative="greater")$p.val)
  assign(paste0("ranksum.lnHMdisp_", i, "_related"), 
         wilcox.test(
           rank(get(paste0("abslfcnlp.lnHMdisp_", i)))[get(paste0(i, "_related")) == T], 
           rank(get(paste0("abslfcnlp.lnHMdisp_", i)))[get(paste0(i, "_related")) == F], 
           paired=F, 
           alternative="greater")$p.val)
  assign(paste0("ranksum.voom_", i, "_pathway"), 
         wilcox.test(
           rank(get(paste0("abslfcnlp.voom_", i)))[get(paste0(i, "_pathway")) == T], 
           rank(get(paste0("abslfcnlp.voom_", i)))[get(paste0(i, "_pathway")) == F], 
           paired=F, 
           alternative="greater")$p.val)
  assign(paste0("ranksum.lnHMdisp_", i, "_pathway"), 
         wilcox.test(
           rank(get(paste0("abslfcnlp.lnHMdisp_", i)))[get(paste0(i, "_pathway")) == T], 
           rank(get(paste0("abslfcnlp.lnHMdisp_", i)))[get(paste0(i, "_pathway")) == F], 
           paired=F, 
           alternative="greater")$p.val)
}
ranksum.tests <- data.frame(matrix(nrow=4, ncol=8))
names(ranksum.tests) <- c("brca", "kirc", "thca", "luad", "lihc",  "lusc", "prad", "coad")
rownames(ranksum.tests) <- c("voom_related", "lnHMdisp_related", "voom_pathway", "lnHMdisp_pathway")
for (i in names(ranksum.tests)) {
  ranksum.tests[[i]][1] <- get(paste0("ranksum.voom_", i, "_related"))
  ranksum.tests[[i]][2] <- get(paste0("ranksum.lnHMdisp_", i, "_related"))
  ranksum.tests[[i]][3] <- get(paste0("ranksum.voom_", i, "_pathway"))
  ranksum.tests[[i]][4] <- get(paste0("ranksum.lnHMdisp_", i, "_pathway"))
  # rm(list=paste0("ranksum_", i, "_", tests.ranksum[j]))
}

for (i in c("ranksum.tests.weak_known.genes", 
            "ranksum.tests_pathway.genes.minus.known.genes")) {
  assign(
    i, read.csv(
      here("Results/TCGA paired data results May 2020", paste0(i, ".csv")), row.names=1
    )
  )
}
ranksum.tests_related <- ranksum.tests.weak_known.genes[c("p.voom", "p.disp.lnHM.log"), ] %>% 
  select(-lihc_long.chain)
ranksum.tests_pathway <- ranksum.tests_pathway.genes.minus.known.genes[c("p.voom", "p.disp.lnHM.log"), ] %>% 
  select(-lihc_long.chain)
ranksum.tests_pval <- rbind(ranksum.tests_related, ranksum.tests_pathway)
rownames(ranksum.tests_pval) <- rownames(ranksum.tests)

ranksum.tests > ranksum.tests_pval
ranksum.tests < 0.05
ranksum.tests_pval < 0.05
# In most cases, WRS test gives smaller p-values when genes are ranked by p-value than |LFC log p| for 
# voom, and the opposite for lnHMdisp - i.e. if add extra information and ensure that there are no ties, 
# voom more strongly differentiates cancer-related genes, and lnHMdisp less strongly differentiates 
# cancer-related genes. In most cases the decision at 0.05 level doesn't change, but there is a very clear 
# pattern. The question is whether this really says anything or not. Do the ties in p-values with HM 
# disguise a lack of discriminative ability? Is that possible?



