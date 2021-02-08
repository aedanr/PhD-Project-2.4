library(here)

for (i in c("brca", "kirc", "thca", "luad", "lihc", 
            "lihc_long.chain", "lusc", "prad", "coad")) {
  assign(paste0("pvals.", i), 
         read.csv(here("Results/TCGA paired data results May 2020", 
                       paste0("pvals.", i, ".csv")), 
                  stringsAsFactors=F, 
                  header=T, 
                  row.names=1))
}
rank.brca.voom <- rank(pvals.brca$p.voom)
names(rank.brca.voom) <- rownames(pvals.brca)
pvals.brca.voom <- pvals.brca$p.voom
names(pvals.brca.voom) <- rownames(pvals.brca)
nlogp.brca.voom <- -log10(pvals.brca$p.voom)
names(nlogp.brca.voom) <- rownames(pvals.brca)
write.table(rank.brca.voom, here("rank.brca.voom.rnk"), quote=F, col.names=F, sep="\t")
write.table(pvals.brca.voom, here("pvals.brca.voom.rnk"), quote=F, col.names=F, sep="\t")
write.table(nlogp.brca.voom, here("nlogp.brca.voom.rnk"), quote=F, col.names=F, sep="\t")
rank.brca.mean.lnHM.log <- rank(pvals.brca$p.mean.lnHM.log)
names(rank.brca.mean.lnHM.log) <- rownames(pvals.brca)
pvals.brca.mean.lnHM.log <- pvals.brca$p.mean.lnHM.log
names(pvals.brca.mean.lnHM.log) <- rownames(pvals.brca)
nlogp.brca.mean.lnHM.log <- -log10(pvals.brca$p.mean.lnHM.log)
names(nlogp.brca.mean.lnHM.log) <- rownames(pvals.brca)
write.table(rank.brca.mean.lnHM.log, here("rank.brca.mean.lnHM.log.rnk"), quote=F, col.names=F, sep="\t")
write.table(pvals.brca.mean.lnHM.log, here("pvals.brca.mean.lnHM.log.rnk"), quote=F, col.names=F, sep="\t")
write.table(nlogp.brca.mean.lnHM.log, here("nlogp.brca.mean.lnHM.log.rnk"), quote=F, col.names=F, sep="\t")


# GSEA recommends to avoid using GSEA to collapse ranked list to gene symbols for GSEAPreranked
# Need to avoid duplicate symbols and Ensembl IDs:
# If one Ensembl ID maps to more than one symbol, need to choose which symbol to use. If using 
# one of GSEA's mapping (.chip) files, this hopefully won't be an issue.
# If multiple Ensembl IDs map to the same symbol, not sure what to do. Again hopefully won't be 
# an issue with GSEA's mapping, but if it is, probably need to arbitrarily choose one to keep.

# Using brca as an example to work out how to convert Ensembl IDs to gene symbols
# Have seen people average whatever metric is being used across multiple IDs with same symbol.
# Probably really doesn't matter much since it's going to be a very small number of genes that are 
# affected by duplicates in either direction.
dim(pvals.brca)
# brca list has 23211 genes
length(unique(rownames(pvals.brca)))
# All 23211 ensembl IDs are unique

library(recount)
library(org.Hs.eg.db)
gencode <- gsub('\\..*', '', names(recount_genes))
gene_info <- unique(select(org.Hs.eg.db, gencode, c('ENTREZID', 'SYMBOL','ENSEMBL'), 'ENSEMBL'))
dim(gene_info)
dim(unique(gene_info))
names(gene_info)
length(unique(gene_info$ENTREZID))
length(unique(gene_info$SYMBOL))
dim(unique(gene_info[, 2:3]))
length(unique(gene_info$ENSEMBL))
# 58191 unique combinations of ensembl ID, Entrez ID and gene symbol in recount, including 25665 
# unique Entrez IDs and gene symbols, which have a one-to-one correspondence, and 57992 unique 
# ensembl IDs, meaning at most only 199 ensembl IDs mapping to more than one symbol, but a lot of 
# symbols mapping to more than one ensembl ID.
sum(rownames(pvals.brca) %in% gene_info$ENSEMBL)
# All 23211 ensembl IDs in brca list are in recount (which must be true since it's recount data).
gene_info_brca <- gene_info[gene_info$ENSEMBL %in% rownames(pvals.brca), ]
dim(gene_info_brca)
# 23343 ensembl IDs in recount in brca list, so up to 132 ensembl IDs mapping to more than one symbol.
length(unique(gene_info_brca$SYMBOL))
# 17344 symbols in recount have an ensembl ID in brca list, so 5867 ensembl IDs in brca list either 
# don't map to a symbol or map to a symbol that also maps to other ensembl IDs.
sum(is.na(gene_info_brca$SYMBOL))
sum(is.na(gene_info_brca$ENTREZID))
# 5942 ensembl IDs in recount that are in brca list don't have a symbol.
gene_info_brca <- gene_info_brca[-which(is.na(gene_info_brca$SYMBOL)), ]
dim(gene_info_brca)
length(unique(gene_info_brca$ENSEMBL))
length(unique(gene_info_brca$SYMBOL))
# 17401 ensembl IDs in recount with a symbol in brca list, of which 17269 are unique, so up to 132 
# ensembl IDs mapping to more than one symbol, and 17343 unique symbols, so up to 58 symbols mapping 
# to more than one ensembl ID.
table(table(gene_info_brca$ENSEMBL))
# 125 ensembl IDs map to two symbols, 2 to three symbols, 1 to four symbols. For these, need to 
# either choose one symbol to use or exclude the genes.
table(table(gene_info_brca$SYMBOL))
# 58 symbols map to two ensembl IDs. For these, need to either somehow combine the ranking metrics 
# of each pair of ensemble IDs or exclude these genes.
dup_ensembl <- gene_info_brca$ENSEMBL[duplicated(gene_info_brca$ENSEMBL)]
dup_symbol <- gene_info_brca$SYMBOL[duplicated(gene_info_brca$SYMBOL)]
gene_info_brca[which(gene_info_brca$ENSEMBL %in% dup_ensembl), ]
gene_info_brca[which(gene_info_brca$SYMBOL %in% dup_symbol), ]
# Easiest and probably safest is to exclude all duplicates. Duplicate symbols is probably likely to 
# be splice variants, but in the first example I looked at (ENSG00000099974 and ENSG00000099977, 
# DDTL), one of the symbols looks to be wrong - ENSG00000099977 corresponds to DDT on ensembl.org, 
# although the locations overlap (but strands are opposite).
dim(gene_info_brca[-which(
  gene_info_brca$ENSEMBL %in% dup_ensembl | gene_info_brca$SYMBOL %in% dup_symbol
), ])
# Excluding all duplicates leaves 17058 genes.

Human_ENSEMBL_Gene_MSigDB.v7.1.chip <- read.delim("Human_ENSEMBL_Gene_MSigDB.v7.1.chip", stringsAsFactors=F)
dim(Human_ENSEMBL_Gene_MSigDB.v7.1.chip)
dim(unique(Human_ENSEMBL_Gene_MSigDB.v7.1.chip))
names(Human_ENSEMBL_Gene_MSigDB.v7.1.chip)
length(unique(Human_ENSEMBL_Gene_MSigDB.v7.1.chip$Gene.Symbol))
length(unique(Human_ENSEMBL_Gene_MSigDB.v7.1.chip$Gene.Title))
dim(unique(Human_ENSEMBL_Gene_MSigDB.v7.1.chip[, 2:3]))
length(unique(Human_ENSEMBL_Gene_MSigDB.v7.1.chip$Probe.Set.ID))
# 42959 unique combinations of ensembl ID, gene symbol and gene name in MSigDB 7.1 ensembl list, 
# including 38404 unique symbols and 38275 unique names, with all names mapping to one symbol, and 
# 42959 unique ensembl IDs, so no ensembl IDs map to more than one symbol, but up to 4555 symbols 
# map to more than one ensembl ID.
sum(rownames(pvals.brca) %in% Human_ENSEMBL_Gene_MSigDB.v7.1.chip$Probe.Set.ID)
# 19379 ensembl IDs in brca list are in MSigDB 7.1 ensembl list, which presumably means that even 
# if I had symbols for the other 3832 genes, they wouldn't be in any gene sets.
msig_brca <- Human_ENSEMBL_Gene_MSigDB.v7.1.chip[
  Human_ENSEMBL_Gene_MSigDB.v7.1.chip$Probe.Set.ID %in% rownames(pvals.brca), 
  ]
dim(msig_brca)
# 19379 ensembl IDs in MSigDB in brca list, same as number in brca list in MSigDB, so no ensembl 
# IDs mapping to more than one symbol, as already knew.
length(unique(msig_brca$Gene.Symbol))
sum(is.na(msig_brca$Gene.Symbol))
# 19367 symbols in MSigDB ensembl list have an ensembl ID in brca list, and no NA, so 12 ensembl 
# IDs in brca list map to a symbol that also maps to other ensembl IDs.
table(table(msig_brca$Probe.Set.ID))
# 19379 Ensembl IDs map to one symbol, as already knew
table(table(msig_brca$Gene.Symbol))
# 12 symbols map to two ensembl IDs.
dup_symbol <- msig_brca$Gene.Symbol[duplicated(msig_brca$Gene.Symbol)]
msig_brca[which(msig_brca$Gene.Symbol %in% dup_symbol), ]
# check out duplicates when ensembl.org working again but looks like this will definitely be a 
# better way to do it than using recount list anyway.
pvals.brca[msig_brca[which(msig_brca$Gene.Symbol %in% dup_symbol), ]$Probe.Set.ID, ]
# p-values pretty similar for most pairs, so probably ok to average metrics.


# GSEA says to make sure its weighted scoring scheme - incrementing a running sum statistic by the 
# absolute value of the ranking metric when a gene belongs to a set - applies to the ranking 
# statistic used - i.e. that the magnitude of the ranking statistic is biologically meaningful. 
# This clearly means that using ranks is no good. I think -log10 p-value should be fine (or actually 
# just log p-value since it uses magnitude). The magnitude of the statistic should increase with 
# (evidence for) biological relevance as for t-statistic. Could possibly use something like 1-p. 
# Not sure if there's any good reason for one or the other. If I think the difference between say 
# 0.1 and 0.01 is the same as the difference between 0.01 and 0.001, then should use log. If I think 
# the difference between 0.1 and 0.07 is the same as the difference between 0.07 and 0.04 then I 
# should use 1-p. GSEA recommend to use Enrichment statistic = 'classic' "when in doubt", i.e. if 
# not sure whether magnitude of ranking statistic is biologically meaningful. This sets the weight 
# to 0 instead of the default, 1. Some people use sign of FC multiplied by 1/p. This, like t-stat, 
# puts relevant genes at top and bottom of list. But since GSEA talks about the magnitude of the 
# ranking statistic being biologically meaningful, I'm not sure if this matters. But maybe it does, 
# in which case I can just do the same (with whatever statistic I use) - multiply by the sign of 
# the change in mean or dispersion. (I wouldn't be able to do this with diff dist though.) Could 
# also just use log2FC (Gordon Smyth recommended shrunken log2FC values in a message board post, 
# not sure if the shrunken part means that evidence for change is accounted for where it isn't if 
# just use normal log2FC though). Another issue is that people talk about gene sets being 
# upregulated or downregulated, which I think means there is an overall increase or decrease in 
# expression for genes in the set, but I don't necessarily expect a consistent change in mean or 
# dispersion across genes in a set (although this also really applies to a traditional DE 
# analysis as far as I can see - there's no reason to think that all genes in a set should change 
# expression in the same direction). From a blog post, a positive enrichment score means that 
# genes in a set are over-represented at the top of the list, and a negative enrichment score 
# means that genes in a set are over-represented at the bottom of a list. Since I don't think 
# there's any reason to expect dispersion (or mean) to be consistently changed in one direction 
# for genes in a given set, I could use a ranking statistic that doesn't distinguish by direction 
# and only consider gene sets with a positive enrichment score. A negative enrichment score would 
# I suppose in a way mean that a set is really not enriched. The same blog post suggests 
# signed FC * -log10p as ranking statistic (but I think this actually means just sign of FC). Maybe 
# I could use magnitude of FC * -10logp, which would also get around the issue of identical ranks. 
# Using only the magnitude of change without direction was suggested in a comment on the blog post 
# that talked about enrichment scores for the same reason I'm thinking of, and there was a reply that 
# this would be valid. Another possibility (recommended by Kevin Blighe in a Biostars thread) is to set 
# an adjusted p-value threshold and then rank based on absolute log2FC. But this was explicitly in 
# reference to gene set enrichment analysis as a general concept rather than the GSEA software itself; 
# when asked about that, he said he would rank by fold change, but didn't elaborate. Other people say 
# to include all genes for GSEA. Seem to be a lot of different ranking options for normal GSEA (i.e. 
# not pre-ranked), so really probably is pretty arbitrary. Even though I don't think I've seen 
# anyone else use it, strongly consider multiplying log2FC by -log10p. Surely this combines ranking 
# by evidence for difference and strength of difference and so shouldn't be invalid at all since 
# either method alone is frequently used and ranking using product of these should be in between the 
# two extremes that either alone would give.


# Check concordance of different metrics using voom results ####
lfc.brca <- readRDS(file=here("Results/TCGA paired data results May 2020", paste0("lfc.TCGA.paired_brca.rds")))
lfc.mean.expHM <- lfc.brca$lfc.mean.expHM
lfc.disp.expHM <- lfc.brca$lfc.disp.expHM
lfc.mean.lnHM <- lfc.brca$lfc.mean.lnHM
lfc.disp.lnHM <- lfc.brca$lfc.disp.lnHM
rm(lfc.brca)
dat <- readRDS(here("recount data/TCGA", "brca_matched_pairs.rds"))
dat <- dat[, dat$cgc_case_histological_diagnosis == "Infiltrating Ductal Carcinoma"]
dat <- dat[which(rowMeans(apply(assay(dat), 2, function(x) x*1e6/sum(x))) > 0.5), ]
group <- dat$gdc_cases.samples.sample_type
group[group == "Solid Tissue Normal"] <- 1
group[group == "Primary Tumor"] <- 2
group <- factor(group)
design <- model.matrix(~group)
counts <- assay(dat)
rm(dat)
library(limma)
library(edgeR)
nf <- calcNormFactors(counts, method="TMM")
dat.edgeR <- estimateDisp(DGEList(counts=counts, norm.factors=nf, group=group), design)
lfc.voom <- topTable(eBayes(lmFit(voom(dat.edgeR), design)), number=Inf, sort.by='none')$logFC
rm(group, design, counts, nf, dat.edgeR)

nlogp.voom <- -log10(pvals.brca$p.voom)
names(nlogp.voom) <- rownames(pvals.brca)
names(lfc.voom) <- rownames(pvals.brca)
abslfc.voom <- abs(lfc.voom)
lfcnlogp.voom <- lfc.voom * nlogp.voom
abslfcnlogp.voom <- abslfc.voom * nlogp.voom
signlfcnlogp.voom <- sign(lfc.voom) * nlogp.voom
write.table(nlogp.voom, here("nlogp.voom.rnk"), quote=F, col.names=F, sep="\t")
write.table(lfc.voom, here("lfc.voom.rnk"), quote=F, col.names=F, sep="\t")
write.table(abslfc.voom, here("abslfc.voom.rnk"), quote=F, col.names=F, sep="\t")
write.table(lfcnlogp.voom, here("lfcnlogp.voom.rnk"), quote=F, col.names=F, sep="\t")
write.table(abslfcnlogp.voom, here("abslfcnlogp.voom.rnk"), quote=F, col.names=F, sep="\t")
write.table(signlfcnlogp.voom, here("signlfcnlogp.voom.rnk"), quote=F, col.names=F, sep="\t")
# When collapse using MSigDB Ensembl 7.1 in GSEA, end up with 19367 genes, same as through my manual 
# sorting, so that's reassuring and means I can just do it in GSEA.


## Compare ranking metrics on brca voom data with hallmark gene sets ####
folder <- "Results/GSEA results June 2020/Testing ranking metrics"
abslfc <- read.csv(here(folder, "brca_voom_hallmarks_abslfc.csv"), stringsAsFactors=F)
abslfcnlogp <- read.csv(here(folder, "brca_voom_hallmarks_abslfcnlogp.csv"), stringsAsFactors=F)
lfc <- read.csv(here(folder, "brca_voom_hallmarks_lfc.csv"), stringsAsFactors=F)
lfcnlogp <- read.csv(here(folder, "brca_voom_hallmarks_lfcnlogp.csv"), stringsAsFactors=F)
nlogp <- read.csv(here(folder, "brca_voom_hallmarks_nlogp.csv"), stringsAsFactors=F)
signlfcnlogp <- read.csv(here(folder, "brca_voom_hallmarks_signlfcnlogp.csv"), stringsAsFactors=F)

nlogp <- nlogp[order(nlogp$GS), ]
lfc <- lfc[order(lfc$GS), ]
abslfc <- abslfc[order(abslfc$GS), ]
lfcnlogp <- lfcnlogp[order(lfcnlogp$GS), ]
abslfcnlogp <- abslfcnlogp[order(abslfcnlogp$GS), ]
signlfcnlogp <- signlfcnlogp[order(signlfcnlogp$GS), ]
NES <- cbind(nlogp$NES, lfc$NES, signlfcnlogp$NES, abslfc$NES, lfcnlogp$NES, abslfcnlogp$NES)
colnames(NES) <- c("nlogp", "lfc", "signlfcnlogp", "abslfc", "lfcnlogp", "abslfcnlogp")
pairs(NES, upper.panel=NULL, pch=20)
# nlogp, lfc and signlfcnlogp are established metrics.
cor(nlogp$NES, lfc$NES) # 0.2645732
cor(nlogp$NES, signlfcnlogp$NES) # 0.3745789
cor(lfc$NES, signlfcnlogp$NES) # 0.9745734
# Not much correlation between nlogp and lfc or nlogp and signlfcnlogp. Not surprising since nlogp is 
# unidirectional and the others are bidirectional. High correlation between lfc and signlfcnlogp, which 
# implies that the differences between nlogp and lfc are mainly because of directionality, and the 
# magnitudes of nlogp and lfc give comparable metrics.
cor(nlogp$NES, lfcnlogp$NES) # 0.3413462
cor(lfc$NES, lfcnlogp$NES) # 0.9891913
cor(signlfcnlogp$NES, lfcnlogp$NES) # 0.9938079
# Correlation between nlogp and lfcnlogp is similar to correlation between nlogp and signlfcnlogp and 
# between nlogp and lfc. The scatterplots of NES look very similar in all three cases, so in each case 
# can safely attribute most of the differences in results to directionality.
# Very high correlation between lfc and lfcnlogp, and signlfcnlogp and lfcnlogp. This suggests that 
# combining the magnitude of lfc with nlogp doesn't affect the ranking much when compared to lfc alone 
# or directional nlogp, and that multiplying lfc by nlogp is valid at least for bidirectional analysis. 
cor(lfc$NES, abslfc$NES) # -0.2332679
cor(lfcnlogp$NES, abslfcnlogp$NES) # 0.05618524
# Very low or negative correlation between unidirection and bidirectional versions of lfc and lfcnlogp, 
# but no difference in validity - difference just reflects difference in directionality, although it's 
# interesting that the scatterplots look quite different between these pairs than the do between the 
# unidirectional and birectional versions of nlogp.
cor(nlogp$NES, abslfc$NES) # 0.3563862
cor(nlogp$NES, abslfcnlogp$NES) # 0.8182463
cor(abslfc$NES, abslfcnlogp$NES) # 0.7367668
# Not much correlation between the two basic unidirectional metrics, nlogp and abslfc, but high 
# correlation between each one and the combination. For all three pairs, scatterplots show higher 
# correlation for higher NES values, which is further reassuring since it's the genesets with the 
# highest enrichment scores that we're interested in - negative NES for unidirection measures doesn't 
# really mean anything. Since nlogp and abslfc are either established or easily justifiable metrics for 
# unidirectional analysis, this means that abslfcnlogp can be easily justified as well.

pairs(cbind(nlogp.voom, lfc.voom, abslfc.voom, lfcnlogp.voom, abslfcnlogp.voom, signlfcnlogp.voom), 
      upper.panel=NULL, pch=20)
# Looking at correlation of the metrics themselves, abslfcnlogp is more strongly correlated with abslfc and 
# nlogp than abslfc and nlogp are with each other. lfcnlogp is about as strongly correlated with lfc and 
# signlfcnlogp as lfc and signlfcnlogp are with each other. Since lfc and signlfcnlogp are both established 
# metrics, and using absolute value to test for ranking without requiring consistent direction of change, this 
# means that abslfcnlogp is easily justifiable as a ranking metric.
cor(abslfcnlogp.voom, abslfc.voom) # 0.8721041
cor(abslfcnlogp.voom, nlogp.voom) # 0.8302976
cor(abslfc.voom, nlogp.voom) # 0.759716
cor(lfcnlogp.voom, lfc.voom) # 0.8639055
cor(lfcnlogp.voom, signlfcnlogp.voom) # 0.838368
cor(lfc.voom, signlfcnlogp.voom) # 0.8652331


## Test run GSEA ####
# For inital run, compare voom and lnHM disp for brca, kirc and lihc
library(here)
library(recount)
library(edgeR)
library(limma)
for (i in c("brca", "kirc", "lihc")) {
  assign(paste0("pvals.", i), 
         read.csv(here("Results/TCGA paired data results May 2020", 
                       paste0("pvals.", i, ".csv")), 
                  stringsAsFactors=F, 
                  header=T, 
                  row.names=1))
}
for (i in c("brca", "kirc", "lihc")) {
  assign(paste0("lfc.", i), 
         readRDS(here("Results/TCGA paired data results May 2020", 
                       paste0("lfc.TCGA.paired_", i, ".rds"))))
}

dat <- readRDS(here("recount data/TCGA", "brca_matched_pairs.rds"))
dat <- dat[, dat$cgc_case_histological_diagnosis == "Infiltrating Ductal Carcinoma"]
dat <- dat[which(rowMeans(apply(assay(dat), 2, function(x) x*1e6/sum(x))) > 0.5), ]
group <- dat$gdc_cases.samples.sample_type
group[group == "Solid Tissue Normal"] <- 1
group[group == "Primary Tumor"] <- 2
group <- factor(group)
design <- model.matrix(~group)
counts <- assay(dat)
rm(dat)
nf <- calcNormFactors(counts, method="TMM")
dat.edgeR <- estimateDisp(DGEList(counts=counts, norm.factors=nf, group=group), design)
lfc.brca.voom <- topTable(eBayes(lmFit(voom(dat.edgeR), design)), number=Inf, sort.by='none')$logFC
rm(group, design, counts, nf, dat.edgeR)

dat <- readRDS(here("recount data/TCGA", "kirc_matched_pairs.rds"))
dat <- dat[which(rowMeans(apply(assay(dat), 2, function(x) x*1e6/sum(x))) > 0.5), ]
group <- dat$gdc_cases.samples.sample_type
group[group == "Solid Tissue Normal"] <- 1
group[group == "Primary Tumor"] <- 2
group <- factor(group)
design <- model.matrix(~group)
counts <- assay(dat)
rm(dat)
nf <- calcNormFactors(counts, method="TMM")
dat.edgeR <- estimateDisp(DGEList(counts=counts, norm.factors=nf, group=group), design)
lfc.kirc.voom <- topTable(eBayes(lmFit(voom(dat.edgeR), design)), number=Inf, sort.by='none')$logFC
rm(group, design, counts, nf, dat.edgeR)

dat <- readRDS(here("recount data/TCGA", "lihc_matched_pairs.rds"))
dat <- dat[, dat$cgc_case_histological_diagnosis == "Hepatocellular Carcinoma"]
dat <- dat[which(rowMeans(apply(assay(dat), 2, function(x) x*1e6/sum(x))) > 0.5), ]
group <- dat$gdc_cases.samples.sample_type
group[group == "Solid Tissue Normal"] <- 1
group[group == "Primary Tumor"] <- 2
group <- factor(group)
design <- model.matrix(~group)
counts <- assay(dat)
rm(dat)
nf <- calcNormFactors(counts, method="TMM")
dat.edgeR <- estimateDisp(DGEList(counts=counts, norm.factors=nf, group=group), design)
lfc.lihc.voom <- topTable(eBayes(lmFit(voom(dat.edgeR), design)), number=Inf, sort.by='none')$logFC
rm(group, design, counts, nf, dat.edgeR)

for (i in c("brca", "kirc", "lihc")) {
  assign(paste0("abslfcnlp.voom.", i), 
         abs(get(paste0("lfc.", i, ".voom"))) * -log10(get(paste0("pvals.", i))$p.voom))
  assign(paste0("abslfcnlp.lnHMdisp.", i), 
         abs(get(paste0("lfc.", i))$lfc.disp.lnHM) * -log10(get(paste0("pvals.", i))$p.disp.lnHM.log))
}
rm(i)
names(abslfcnlp.voom.brca) <- rownames(pvals.brca)
names(abslfcnlp.lnHMdisp.brca) <- rownames(pvals.brca)
names(abslfcnlp.voom.kirc) <- rownames(pvals.kirc)
names(abslfcnlp.lnHMdisp.kirc) <- rownames(pvals.kirc)
names(abslfcnlp.voom.lihc) <- rownames(pvals.lihc)
names(abslfcnlp.lnHMdisp.lihc) <- rownames(pvals.lihc)
write.table(abslfcnlp.voom.brca, here("GSEA data", "abslfcnlp.voom.brca.rnk"), 
            quote=F, col.names=F, sep="\t")
write.table(abslfcnlp.lnHMdisp.brca, here("GSEA data", "abslfcnlp.lnHMdisp.brca.rnk"), 
            quote=F, col.names=F, sep="\t")
write.table(abslfcnlp.voom.kirc, here("GSEA data", "abslfcnlp.voom.kirc.rnk"), 
            quote=F, col.names=F, sep="\t")
write.table(abslfcnlp.lnHMdisp.kirc, here("GSEA data", "abslfcnlp.lnHMdisp.kirc.rnk"), 
            quote=F, col.names=F, sep="\t")
write.table(abslfcnlp.voom.lihc, here("GSEA data", "abslfcnlp.voom.lihc.rnk"), 
            quote=F, col.names=F, sep="\t")
write.table(abslfcnlp.lnHMdisp.lihc, here("GSEA data", "abslfcnlp.lnHMdisp.lihc.rnk"), 
            quote=F, col.names=F, sep="\t")

## Repeat for other cancers
library(here)
library(recount)
library(edgeR)
library(limma)
for (i in c("thca", "luad", "lusc", "prad", "coad")) {
  assign(paste0("pvals.", i), 
         read.csv(here("Results/TCGA paired data results May 2020", 
                       paste0("pvals.", i, ".csv")), 
                  stringsAsFactors=F, 
                  header=T, 
                  row.names=1))
}
for (i in c("thca", "luad", "lusc", "prad", "coad")) {
  assign(paste0("lfc.", i), 
         readRDS(here("Results/TCGA paired data results May 2020", 
                      paste0("lfc.TCGA.paired_", i, ".rds"))))
}

dat <- readRDS(here("recount data/TCGA", "thca_matched_pairs.rds"))
dat <- dat[, dat$cgc_case_histological_diagnosis == "Thyroid Papillary Carcinoma - Classical/usual"] # thca
dat <- dat[which(rowMeans(apply(assay(dat), 2, function(x) x*1e6/sum(x))) > 0.5), ]
group <- dat$gdc_cases.samples.sample_type
group[group == "Solid Tissue Normal"] <- 1
group[group == "Primary Tumor"] <- 2
group <- factor(group)
design <- model.matrix(~group)
counts <- assay(dat)
rm(dat)
nf <- calcNormFactors(counts, method="TMM")
dat.edgeR <- estimateDisp(DGEList(counts=counts, norm.factors=nf, group=group), design)
lfc.thca.voom <- topTable(eBayes(lmFit(voom(dat.edgeR), design)), number=Inf, sort.by='none')$logFC
rm(group, design, counts, nf, dat.edgeR)

dat <- readRDS(here("recount data/TCGA", "luad_matched_pairs.rds"))
dat <- dat[, dat$cgc_case_histological_diagnosis == "Lung Adenocarcinoma- Not Otherwise Specified (NOS)"] # luad
dat <- dat[which(rowMeans(apply(assay(dat), 2, function(x) x*1e6/sum(x))) > 0.5), ]
group <- dat$gdc_cases.samples.sample_type
group[group == "Solid Tissue Normal"] <- 1
group[group == "Primary Tumor"] <- 2
group <- factor(group)
design <- model.matrix(~group)
counts <- assay(dat)
rm(dat)
nf <- calcNormFactors(counts, method="TMM")
dat.edgeR <- estimateDisp(DGEList(counts=counts, norm.factors=nf, group=group), design)
lfc.luad.voom <- topTable(eBayes(lmFit(voom(dat.edgeR), design)), number=Inf, sort.by='none')$logFC
rm(group, design, counts, nf, dat.edgeR)

dat <- readRDS(here("recount data/TCGA", "lusc_matched_pairs.rds"))
dat <- dat[, dat$cgc_case_histological_diagnosis == "Lung Squamous Cell Carcinoma- Not Otherwise Specified (NOS)" |
             is.na(dat$cgc_case_histological_diagnosis)] # lusc
dat <- dat[which(rowMeans(apply(assay(dat), 2, function(x) x*1e6/sum(x))) > 0.5), ]
group <- dat$gdc_cases.samples.sample_type
group[group == "Solid Tissue Normal"] <- 1
group[group == "Primary Tumor"] <- 2
group <- factor(group)
design <- model.matrix(~group)
counts <- assay(dat)
rm(dat)
nf <- calcNormFactors(counts, method="TMM")
dat.edgeR <- estimateDisp(DGEList(counts=counts, norm.factors=nf, group=group), design)
lfc.lusc.voom <- topTable(eBayes(lmFit(voom(dat.edgeR), design)), number=Inf, sort.by='none')$logFC
rm(group, design, counts, nf, dat.edgeR)

dat <- readRDS(here("recount data/TCGA", "prad_matched_pairs.rds"))
dat <- dat[which(rowMeans(apply(assay(dat), 2, function(x) x*1e6/sum(x))) > 0.5), ]
group <- dat$gdc_cases.samples.sample_type
group[group == "Solid Tissue Normal"] <- 1
group[group == "Primary Tumor"] <- 2
group <- factor(group)
design <- model.matrix(~group)
counts <- assay(dat)
rm(dat)
nf <- calcNormFactors(counts, method="TMM")
dat.edgeR <- estimateDisp(DGEList(counts=counts, norm.factors=nf, group=group), design)
lfc.prad.voom <- topTable(eBayes(lmFit(voom(dat.edgeR), design)), number=Inf, sort.by='none')$logFC
rm(group, design, counts, nf, dat.edgeR)

dat <- readRDS(here("recount data/TCGA", "coad_matched_pairs.rds"))
dat <- dat[, dat$cgc_case_histological_diagnosis == "Colon Adenocarcinoma"] # coad
dat <- dat[which(rowMeans(apply(assay(dat), 2, function(x) x*1e6/sum(x))) > 0.5), ]
group <- dat$gdc_cases.samples.sample_type
group[group == "Solid Tissue Normal"] <- 1
group[group == "Primary Tumor"] <- 2
group <- factor(group)
design <- model.matrix(~group)
counts <- assay(dat)
rm(dat)
nf <- calcNormFactors(counts, method="TMM")
dat.edgeR <- estimateDisp(DGEList(counts=counts, norm.factors=nf, group=group), design)
lfc.coad.voom <- topTable(eBayes(lmFit(voom(dat.edgeR), design)), number=Inf, sort.by='none')$logFC
rm(group, design, counts, nf, dat.edgeR)

for (i in c("thca", "luad", "lusc", "prad", "coad")) {
  assign(paste0("abslfcnlp.voom.", i), 
         abs(get(paste0("lfc.", i, ".voom"))) * -log10(get(paste0("pvals.", i))$p.voom))
  assign(paste0("abslfcnlp.lnHMdisp.", i), 
         abs(get(paste0("lfc.", i))$lfc.disp.lnHM) * -log10(get(paste0("pvals.", i))$p.disp.lnHM.log))
}
rm(i)
names(abslfcnlp.voom.thca) <- rownames(pvals.thca)
names(abslfcnlp.lnHMdisp.thca) <- rownames(pvals.thca)
names(abslfcnlp.voom.luad) <- rownames(pvals.luad)
names(abslfcnlp.lnHMdisp.luad) <- rownames(pvals.luad)
names(abslfcnlp.voom.lusc) <- rownames(pvals.lusc)
names(abslfcnlp.lnHMdisp.lusc) <- rownames(pvals.lusc)
names(abslfcnlp.voom.prad) <- rownames(pvals.prad)
names(abslfcnlp.lnHMdisp.prad) <- rownames(pvals.prad)
names(abslfcnlp.voom.coad) <- rownames(pvals.coad)
names(abslfcnlp.lnHMdisp.coad) <- rownames(pvals.coad)
write.table(abslfcnlp.voom.thca, here("GSEA data", "abslfcnlp.voom.thca.rnk"), 
            quote=F, col.names=F, sep="\t")
write.table(abslfcnlp.lnHMdisp.thca, here("GSEA data", "abslfcnlp.lnHMdisp.thca.rnk"), 
            quote=F, col.names=F, sep="\t")
write.table(abslfcnlp.voom.luad, here("GSEA data", "abslfcnlp.voom.luad.rnk"), 
            quote=F, col.names=F, sep="\t")
write.table(abslfcnlp.lnHMdisp.luad, here("GSEA data", "abslfcnlp.lnHMdisp.luad.rnk"), 
            quote=F, col.names=F, sep="\t")
write.table(abslfcnlp.voom.lusc, here("GSEA data", "abslfcnlp.voom.lusc.rnk"), 
            quote=F, col.names=F, sep="\t")
write.table(abslfcnlp.lnHMdisp.lusc, here("GSEA data", "abslfcnlp.lnHMdisp.lusc.rnk"), 
            quote=F, col.names=F, sep="\t")
write.table(abslfcnlp.voom.prad, here("GSEA data", "abslfcnlp.voom.prad.rnk"), 
            quote=F, col.names=F, sep="\t")
write.table(abslfcnlp.lnHMdisp.prad, here("GSEA data", "abslfcnlp.lnHMdisp.prad.rnk"), 
            quote=F, col.names=F, sep="\t")
write.table(abslfcnlp.voom.coad, here("GSEA data", "abslfcnlp.voom.coad.rnk"), 
            quote=F, col.names=F, sep="\t")
write.table(abslfcnlp.lnHMdisp.coad, here("GSEA data", "abslfcnlp.lnHMdisp.coad.rnk"), 
            quote=F, col.names=F, sep="\t")

par(mfcol=c(2,5))
for (i in c("thca", "luad", "lusc", "prad", "coad")) {
  for (j in c("voom", "lnHMdisp")) {
    plot(get(paste0("abslfcnlp.", j, ".", i)), pch=20, main=paste0(i, " ", j))
  }
}

for (i in c("thca", "luad", "lusc", "prad", "coad")) {
    plot(density(
      get(paste0("lfc.", i))$lfc.disp.lnHM
    ), main=paste0("abs lfc ", i))
  plot(density(
    -log10(get(paste0("pvals.", i))$p.disp.lnHM.log)
  ), main=paste0("-log10 p ", i))
}

mean(lfc.thca$lfc.disp.lnHM < 0)
mean(lfc.luad$lfc.disp.lnHM < 0)
mean(lfc.lusc$lfc.disp.lnHM < 0)
mean(lfc.prad$lfc.disp.lnHM < 0)
mean(lfc.coad$lfc.disp.lnHM < 0)

for (i in c("thca", "luad", "lusc", "prad", "coad")) {
  plot(get(paste0("lfc.", i))$lfc.disp.lnHM, 
    -log10(get(paste0("pvals.", i))$p.disp.lnHM.log), pch=20, 
    main=paste0("-log10 p vs lfc, ", i))
}





