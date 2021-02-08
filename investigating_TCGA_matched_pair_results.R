# Check that results are all there and look like they should.

library(here)
folder <- paste0("Results/TCGA paired data results May 2020")
res.brca <- readRDS(here(folder, "results.TCGA.paired_brca.rds"))
res.kirc <- readRDS(here(folder, "results.TCGA.paired_kirc.rds"))
res.thca <- readRDS(here(folder, "results.TCGA.paired_thca.rds"))
res.lihc <- readRDS(here(folder, "results.TCGA.paired_lihc.rds"))
names(res.brca)
identical(names(res.brca), names(res.kirc)) # TRUE
identical(names(res.brca), names(res.thca)) # TRUE
identical(names(res.brca), names(res.lihc)) # TRUE
si.brca <- readRDS(here(folder, "sessionInfo.TCGA.paired_brca.rds"))
si.kirc <- readRDS(here(folder, "sessionInfo.TCGA.paired_kirc.rds"))
si.thca <- readRDS(here(folder, "sessionInfo.TCGA.paired_thca.rds"))
si.lihc <- readRDS(here(folder, "sessionInfo.TCGA.paired_lihc.rds"))
si.brca
si.kirc
identical(si.kirc, si.thca) # TRUE
identical(si.kirc, si.lihc) # TRUE
# kirc, thca, lihc run on iHPC, brca run on my laptop
length(names(res.brca)) # 50
names(res.brca) # [[1]] is counts
par(mfrow=c(5,10), mar=c(2,2,1,1), mgp=c(2,0.5,0))
for (i in 2:50) {
  plot(res.brca[[i]], pch=20, cex=0.5)
}
# q.DSS and lfdr.DSS are all 1
# bfdr.expHM increases from left to right at top of range
par(mfrow=c(5,10), mar=c(2,2,1,1), mgp=c(2,0.5,0))
for (i in 2:50) {
  plot(res.kirc[[i]], pch=20, cex=0.5)
}
# Vertical bands in all DESeq2 plots
# q.DSS and lfdr.DSS are all 1
# bfdr.expHM increases from left to right at top of range
par(mfrow=c(5,10), mar=c(2,2,1,1), mgp=c(2,0.5,0))
for (i in 2:50) {
  plot(res.thca[[i]], pch=20, cex=0.5)
}
# Vertical bands in all DESeq2 plots
# q.DSS and lfdr.DSS are all 1
# bfdr.expHM and lnHM increase from left to right at top of range
par(mfrow=c(5,10), mar=c(2,2,1,1), mgp=c(2,0.5,0))
for (i in 2:50) {
  plot(res.lihc[[i]], pch=20, cex=0.5)
}
# Vertical bands in all DESeq2 plots
# Discrete horizontal bands for q.DSS
# bfdr.expHM increases from left to right at top of range


# Test DSS with other datasets ####
dat <- readRDS("C:/Users/aedan/OneDrive - UTS/PhD/Project 2/03 Comparison on simulated and recount data/recount data/GTEx/muscle_50_samples_per_group/muscle_50_set1_DE.rds")
counts <- dat$counts
group <- factor(c(rep(1, 50), rep(2, 50)))
libsizes <- colSums(counts)
nf <- calcNormFactors(counts, method="TMM")
els <- nf * libsizes
sf <- els / exp(mean(log(libsizes)))
dat.DSS <- newSeqCountSet(counts=matrix(counts, ncol=length(group)),
                          designs=as.numeric(group), normalizationFactor=sf)
dat.DSS <- estDispersion(dat.DSS)
res.DSS <- waldTest(dat.DSS, 1, 2)[order(waldTest(dat.DSS, 1, 2)$geneIndex),]
p.DSS <- res.DSS$pval
q.DSS <- res.DSS$fdr
lfdr.DSS <- res.DSS$local.fdr
plot(p.DSS)
plot(q.DSS)
plot(lfdr.DSS)
plot(p.adjust(p.DSS, method="BH"))
# No issues with fdr or lfdr here, so problem seems to be specifically with the TCGA datasets.

i <- "brca"
folder <- paste0("recount data/TCGA")
dat <- readRDS(here(folder, paste0(i, "_matched_pairs.rds")))
dat <- dat[, dat$cgc_case_histological_diagnosis == "Infiltrating Ductal Carcinoma"]
dat <- dat[which(rowMeans(apply(assay(dat), 2, function(x) x*1e6/sum(x))) > 0.5), ]
counts <- assay(dat)
group <- dat$gdc_cases.samples.sample_type
# group <- dat$gdc_cases.samples.sample_type[sample(ncol(dat))]
group[group == "Solid Tissue Normal"] <- 1
group[group == "Primary Tumor"] <- 2
group <- factor(group)
design <- model.matrix(~group)
libsizes <- colSums(counts)
nf <- calcNormFactors(counts, method="TMM")
els <- nf * libsizes
sf <- els / exp(mean(log(libsizes)))
dat.DSS <- newSeqCountSet(counts=matrix(counts, ncol=length(group)),
                          designs=as.numeric(group), normalizationFactor=sf)
dat.DSS <- estDispersion(dat.DSS)
res.DSS <- waldTest(dat.DSS, 1, 2)[order(waldTest(dat.DSS, 1, 2)$geneIndex),]
p.DSS <- res.DSS$pval
q.DSS <- res.DSS$fdr
lfdr.DSS <- res.DSS$local.fdr
plot(p.DSS)
plot(q.DSS)
plot(lfdr.DSS)
plot(p.adjust(p.DSS, method="BH"))
mean(p.DSS < 0.05)
mean(q.DSS < 0.05)
mean(lfdr.DSS < 0.05)
mean(p.adjust(p.DSS, method="BH") < 0.05)
# No problems when group labels randomly permuted.
# Looks suspiciously like there are systematic differences between tumour and normal samples 
# unrelated to any true differences. This would explain the very high levels of differential 
# expression found by all methods (around 75% genes differentially expressed for all datasets).

plot(density(rowMeans(log(counts))))
lines(density(rowMeans(log(counts[, group == 1]))), col='red')
lines(density(rowMeans(log(counts[, group == 2]))), col='blue')
# Doesn't look to be any real overall difference in distributions of means between groups.

plot(rowMeans(log(counts[, group == 1])), rowMeans(log(counts[, group == 2])), pch=20, cex=0.2)
lines(c(1,13), c(1,13))
# No obvious systematic differences in means between groups.

dat.edgeR <- DGEList(counts=counts, norm.factors=nf, group=group)
dat.edgeR <- estimateDisp(dat.edgeR, design)
fit.edgeR.lr <- glmFit(dat.edgeR, design)
test.edgeR.lr <- glmLRT(fit.edgeR.lr)
plotMD(test.edgeR.lr)
# Large LFCs (range -10 to 10; example in edgeR user guide has range -8 to 4).
# Significant differences for very small LFCs.
summary(decideTests(test.edgeR.lr))

dat_small <- dat[, 1:6]
group_small <- group[1:6]
counts_small <- counts[, 1:6]
design_small <- model.matrix(~group_small)
libsizes <- colSums(counts)
nf_small <- calcNormFactors(counts_small, method="TMM")
dat.edgeR_small <- DGEList(counts=counts_small, norm.factors=nf_small, group=group_small)
dat.edgeR_small <- estimateDisp(dat.edgeR_small, design_small)
colnames(dat.edgeR_small) <- paste0(substr(dat_small$gdc_cases.case_id, 1, 4), "-", 
                                    as.character(group_small))
fit.edgeR.lr_small <- glmFit(dat.edgeR_small, design_small)
test.edgeR.lr_small <- glmLRT(fit.edgeR.lr_small)
plotMD(test.edgeR.lr_small)
summary(decideTests(test.edgeR.lr_small))
# Still large range of LFCs, but significant differences only for LFC > 1.5
plotMDS(dat.edgeR_small)
plotBCV(dat.edgeR_small)

colnames(dat.edgeR) <- paste0(substr(dat$gdc_cases.case_id, 1, 4), "-", as.character(group))
plotMDS(dat.edgeR)
plotBCV(dat.edgeR)
plotMD(test.edgeR.lr)

# Nothing seems really off about the data. Seems to just be the very large sample size that is 
# causing a huge number of genes to be called differentially expressed. This isn't a problem in 
# itself for analysis, but the huge number of genes with the same minimum p-value for the HMs 
# will probably be a problem for them (around 13k genes for brca data). See how results look, but 
# might need to do comparisons on smaller sample sizes to see benefit from HMs.


