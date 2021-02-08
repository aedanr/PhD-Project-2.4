library(here)
library(coda)
library(HDInterval)
library(limma)
library(edgeR)
library(DESeq2)
library(DSS)
library(baySeq)
library(MDSeq)
library(recount)

cores <- detectCores() - 1
if(require("parallel")) cl <- makeCluster(cores) else cl <- NULL

cluster <- F

if (cluster) {
  source(here('Data/scripts', '2019-04-03_exponential_hmm_adaptive_proposals_three_chains_function.R'))
  source(here('Data/scripts', '2019-03-27_lognormal_hmm_adaptive_proposals_three_chains_function.R'))
  source(here('Data/scripts', '2020-05-15_hpd_tail_prob_function.R'))
  source(here('Data/scripts', '2019-05-03_bfdr_function.R'))
  folder <- "Data/recount data/TCGA"
} else {
  source(here('scripts', '2019-04-03_exponential_hmm_adaptive_proposals_three_chains_function.R'))
  source(here('scripts', '2019-03-27_lognormal_hmm_adaptive_proposals_three_chains_function.R'))
  source(here('scripts', '2020-05-15_hpd_tail_prob_function.R'))
  source(here('scripts', '2019-05-03_bfdr_function.R'))
  folder <- "recount data/TCGA"
}

# Load and filter data
i <- "lihc"
dat <- readRDS(here(folder, paste0(i, "_matched_pairs.rds")))
# dat <- dat[, dat$cgc_case_histological_diagnosis == "Infiltrating Ductal Carcinoma"] # brca
# dat <- dat[, dat$cgc_case_histological_diagnosis == "Hepatocellular Carcinoma"] # lihc
dat <- dat[which(rowMeans(apply(assay(dat), 2, function(x) x*1e6/sum(x))) > 0.5), ]

# Data
group <- dat$gdc_cases.samples.sample_type
group[group == "Solid Tissue Normal"] <- 1
group[group == "Primary Tumor"] <- 2
group <- factor(group)
design <- model.matrix(~group)
contrasts <- get.model.matrix(group)
counts <- assay(dat)
rm(dat)

# Normalise
libsizes <- colSums(counts)
nf <- calcNormFactors(counts, method="TMM")
els <- nf * libsizes
sf <- els / exp(mean(log(libsizes)))
norm.counts <- t(t(counts) / sf)

# Create data objects
dat.edgeR <- DGEList(counts=counts, norm.factors=nf, group=group)
dat.edgeR <- estimateDisp(dat.edgeR, design)
dat.DESeq2 <- DESeqDataSetFromMatrix(countData=counts, colData=data.frame(group), design=~group)
sizeFactors(dat.DESeq2) <- sf
dat.DSS <- newSeqCountSet(counts=matrix(counts, ncol=length(group)),
                          designs=as.numeric(group), normalizationFactor=sf)
dat.DSS <- estDispersion(dat.DSS)
dat.baySeq <- new('countData', data=counts, replicates=group,
                  groups=list(NDE=rep(1,length(group)), DE=as.numeric(group)))
dat.baySeq@annotation <- data.frame(name = 1:nrow(dat.baySeq@data))
libsizes(dat.baySeq) <- els

# edgeR
fit.edgeR.ql <- glmQLFit(dat.edgeR, design)
test.edgeR.ql <- glmQLFTest(fit.edgeR.ql)
fit.edgeR.lr <- glmFit(dat.edgeR, design)
test.edgeR.lr <- glmLRT(fit.edgeR.lr)
edgeR.et <- exactTest(dat.edgeR)
p.edgeR.ql <- test.edgeR.ql$table$PValue
p.edgeR.lr <- test.edgeR.lr$table$PValue
p.edgeR.et <- edgeR.et$table$PValue
q.edgeR.ql <- topTags(test.edgeR.ql, n=nrow(dat.edgeR$counts), sort='none')$table$FDR
q.edgeR.lr <- topTags(test.edgeR.lr, n=nrow(dat.edgeR$counts), sort='none')$table$FDR
q.edgeR.et <- topTags(edgeR.et, n=nrow(dat.edgeR$counts), sort='none')$table$FDR
rm(fit.edgeR.ql, test.edgeR.ql, fit.edgeR.lr, test.edgeR.lr, edgeR.et)

# DESeq2
fit.DESeq2 <- DESeq(dat.DESeq2, minReplicatesForReplace=Inf)
res.DESeq2.noif <- results(fit.DESeq2, independentFiltering=F, cooksCutoff=F)
res.DESeq2.if <- results(fit.DESeq2, cooksCutoff=F, alpha=0.05)
p.DESeq2.noif <- res.DESeq2.noif$pvalue
p.DESeq2.if <- res.DESeq2.if$pvalue
q.DESeq2.noif <- res.DESeq2.noif$padj
q.DESeq2.if <- res.DESeq2.if$padj
rm(dat.DESeq2, fit.DESeq2, res.DESeq2.noif, res.DESeq2.if)

# limma-voom
dat.voom <- voom(dat.edgeR)
fit.voom <- lmFit(dat.voom, design)
res.voom <- eBayes(fit.voom)
p.voom <- topTable(res.voom, number=Inf, sort.by='none')$P.Value
q.voom <- topTable(res.voom, number=Inf, sort.by='none')$adj.P.Val
rm(dat.edgeR, dat.voom, fit.voom, res.voom)

# DSS
res.DSS <- waldTest(dat.DSS, 1, 2)[order(waldTest(dat.DSS, 1, 2)$geneIndex),]
p.DSS <- res.DSS$pval
q.DSS <- res.DSS$fdr
lfdr.DSS <- res.DSS$local.fdr
rm(dat.DSS, res.DSS)

# baySeq
dat.baySeq <- getPriors.NB(dat.baySeq, samplesize=10000, cl=cl)
dat.baySeq <- getLikelihoods(dat.baySeq, cl=cl) # error on cluster using cl here
# dat.baySeq <- getLikelihoods(dat.baySeq)
res.baySeq <- topCounts(dat.baySeq, group='DE', number=Inf)[
  order(topCounts(dat.baySeq, group='DE', number=Inf)$name), ]
prob.baySeq <- res.baySeq$likes
q.baySeq <- res.baySeq$FDR.DE
rm(dat.baySeq, res.baySeq)

# MDSeq
fit.MDSeq.zi <- MDSeq(counts, offsets=sf, contrast=contrasts, mc.cores=cores) # hangs on Phoenix
# fit.MDSeq.zi <- MDSeq(counts, offsets=sf, contrast=contrasts)
res.MDSeq.zi <- extract.ZIMD(fit.MDSeq.zi, compare=list(A="1",B="2"))
fit.MDSeq.nozi <- MDSeq(counts, offsets=sf, contrast=contrasts, test.ZI=F, mc.cores=cores) # hangs on Phoenix
# fit.MDSeq.nozi <- MDSeq(counts, offsets=sf, contrast=contrasts, test.ZI=F)
res.MDSeq.nozi <- extract.ZIMD(fit.MDSeq.nozi, compare=list(A="1",B="2"))
p.mean.MDSeq.zi <- res.MDSeq.zi$Pvalue.mean
q.mean.MDSeq.zi <- res.MDSeq.zi$FDR.mean
p.mean.MDSeq.nozi <- res.MDSeq.nozi$Pvalue.mean
q.mean.MDSeq.nozi <- res.MDSeq.nozi$FDR.mean
p.disp.MDSeq.zi <- res.MDSeq.zi$Pvalue.disp
q.disp.MDSeq.zi <- res.MDSeq.zi$FDR.disp
p.disp.MDSeq.nozi <- res.MDSeq.nozi$Pvalue.disp
q.disp.MDSeq.nozi <- res.MDSeq.nozi$FDR.disp
rm(fit.MDSeq.zi, res.MDSeq.zi, fit.MDSeq.nozi, res.MDSeq.nozi)

# expHM
expHM <- exp_hmm_adapt_3_chains(counts=t(norm.counts), groups=group)
# expHM <- exp_hmm_adapt_3_chains(counts=t(norm.counts), groups=group, chain.length=5000)
prob.expHMM <- unname(colMeans(as.matrix(expHM$indicators)))
post.prop.expHMM <- mean(as.matrix(expHM$proportion))
thr.expHMM <- sort(prob.expHMM, decreasing=T)[round(nrow(counts) * post.prop.expHMM)]
bfdr.expHMM <- bfdr(prob.expHMM)
mean.diff.expHM.untr <- unname(as.matrix(expHM$means1) - as.matrix(expHM$means2))
p.mean.expHM.untr <- apply(mean.diff.expHM.untr,2,hpd.pval)
rm(mean.diff.expHM.untr); gc()
mean.diff.expHM.log <- unname(log(as.matrix(expHM$means1)) - log(as.matrix(expHM$means2)))
p.mean.expHM.log <- apply(mean.diff.expHM.log,2,hpd.pval)
rm(mean.diff.expHM.log); gc()
q.mean.expHM.untr <- p.adjust(p.mean.expHM.untr, method='BH')
q.mean.expHM.log <- p.adjust(p.mean.expHM.log, method='BH')
disp.diff.expHM.untr <- unname(as.matrix(expHM$disps1) - as.matrix(expHM$disps2))
p.disp.expHM.untr <- apply(disp.diff.expHM.untr,2,hpd.pval)
rm(disp.diff.expHM.untr); gc()
disp.diff.expHM.log <- unname(log(as.matrix(expHM$disps1)) - log(as.matrix(expHM$disps2)))
rm(expHM); gc()
p.disp.expHM.log <- apply(disp.diff.expHM.log,2,hpd.pval)
rm(disp.diff.expHM.log); gc()
q.disp.expHM.untr <- p.adjust(p.disp.expHM.untr, method='BH')
q.disp.expHM.log <- p.adjust(p.disp.expHM.log, method='BH')

# lnHM
lnHM <- ln_hmm_adapt_3_chains(counts=t(norm.counts), groups=group)
# lnHM <- ln_hmm_adapt_3_chains(counts=t(norm.counts), groups=group, chain.length=5000)
prob.lnHMM <- unname(colMeans(as.matrix(lnHM$indicators)))
post.prop.lnHMM <- mean(as.matrix(lnHM$proportion))
thr.lnHMM <- sort(prob.lnHMM, decreasing=T)[round(nrow(counts) * post.prop.lnHMM)]
bfdr.lnHMM <- bfdr(prob.lnHMM)
mean.diff.lnHM.untr <- unname(as.matrix(lnHM$means1) - as.matrix(lnHM$means2))
p.mean.lnHM.untr <- apply(mean.diff.lnHM.untr,2,hpd.pval)
rm(mean.diff.lnHM.untr); gc()
mean.diff.lnHM.log <- unname(log(as.matrix(lnHM$means1)) - log(as.matrix(lnHM$means2)))
p.mean.lnHM.log <- apply(mean.diff.lnHM.log,2,hpd.pval)
rm(mean.diff.lnHM.log); gc()
q.mean.lnHM.untr <- p.adjust(p.mean.lnHM.untr, method='BH')
q.mean.lnHM.log <- p.adjust(p.mean.lnHM.log, method='BH')
disp.diff.lnHM.untr <- unname(as.matrix(lnHM$disps1) - as.matrix(lnHM$disps2))
p.disp.lnHM.untr <- apply(disp.diff.lnHM.untr,2,hpd.pval)
rm(disp.diff.lnHM.untr); gc()
disp.diff.lnHM.log <- unname(log(as.matrix(lnHM$disps1)) - log(as.matrix(lnHM$disps2)))
rm(lnHM); gc()
p.disp.lnHM.log <- apply(disp.diff.lnHM.log,2,hpd.pval)
rm(disp.diff.lnHM.log); gc()
q.disp.lnHM.untr <- p.adjust(p.disp.lnHM.untr, method='BH')
q.disp.lnHM.log <- p.adjust(p.disp.lnHM.log, method='BH')

results <- list(counts = counts, 
                p.edgeR.ql = p.edgeR.ql, 
                p.edgeR.lr = p.edgeR.lr, 
                p.edgeR.et = p.edgeR.et, 
                q.edgeR.ql = q.edgeR.ql, 
                q.edgeR.lr = q.edgeR.lr, 
                q.edgeR.et = q.edgeR.et, 
                p.DESeq2.noif = p.DESeq2.noif, 
                p.DESeq2.if = p.DESeq2.if, 
                q.DESeq2.noif = q.DESeq2.noif, 
                q.DESeq2.if = q.DESeq2.if, 
                p.voom = p.voom,
                q.voom = q.voom, 
                p.DSS = p.DSS,
                q.DSS = q.DSS,
                lfdr.DSS = lfdr.DSS,
                prob.baySeq = prob.baySeq,
                q.baySeq = q.baySeq,
                p.mean.MDSeq.zi = p.mean.MDSeq.zi,
                p.mean.MDSeq.nozi = p.mean.MDSeq.nozi,
                q.mean.MDSeq.zi = q.mean.MDSeq.zi,
                q.mean.MDSeq.nozi = q.mean.MDSeq.nozi,
                p.disp.MDSeq.zi = p.disp.MDSeq.zi,
                p.disp.MDSeq.nozi = p.disp.MDSeq.nozi,
                q.disp.MDSeq.zi = q.disp.MDSeq.zi,
                q.disp.MDSeq.nozi = q.disp.MDSeq.nozi,
                prob.expHMM = prob.expHMM,
                prop.expHMM = post.prop.expHMM,
                thr.expHMM = thr.expHMM,
                bfdr.expHMM = bfdr.expHMM,
                p.mean.expHM.untr = p.mean.expHM.untr,
                p.mean.expHM.log = p.mean.expHM.log,
                q.mean.expHM.untr = q.mean.expHM.untr,
                q.mean.expHM.log = q.mean.expHM.log,
                p.disp.expHM.untr = p.disp.expHM.untr,
                p.disp.expHM.log = p.disp.expHM.log,
                q.disp.expHM.untr = q.disp.expHM.untr,
                q.disp.expHM.log = q.disp.expHM.log,
                prob.lnHMM = prob.lnHMM,
                prop.lnHMM = post.prop.lnHMM,
                thr.lnHMM = thr.lnHMM,
                bfdr.lnHMM = bfdr.lnHMM,
                p.mean.lnHM.untr = p.mean.lnHM.untr,
                p.mean.lnHM.log = p.mean.lnHM.log,
                q.mean.lnHM.untr = q.mean.lnHM.untr,
                q.mean.lnHM.log = q.mean.lnHM.log,
                p.disp.lnHM.untr = p.disp.lnHM.untr,
                p.disp.lnHM.log = p.disp.lnHM.log,
                q.disp.lnHM.untr = q.disp.lnHM.untr,
                q.disp.lnHM.log = q.disp.lnHM.log)

if (cluster) {
  folder <- paste0("Data/TCGA paired data results May 2020")
} else {
  folder <- paste0("Results/TCGA paired data results May 2020")
}
saveRDS(results, file=here(folder, paste0("results.TCGA.paired_", i, ".rds")))
saveRDS(sessionInfo(), file=here(folder, paste0("sessionInfo.TCGA.paired_", i, ".rds")))

# saveRDS(results, file=here(folder, paste0("results.TCGA.paired_", i, "_long.chain.rds")))
# saveRDS(sessionInfo(), file=here(folder, paste0("sessionInfo.TCGA.paired_", i, "_long.chain.rds")))


