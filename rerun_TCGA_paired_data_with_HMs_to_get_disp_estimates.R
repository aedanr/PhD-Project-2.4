library(here)
library(coda)
library(HDInterval)
library(edgeR)
library(recount)

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
i <- "coad"
dat <- readRDS(here(folder, paste0(i, "_matched_pairs.rds")))
# dat <- dat[, dat$cgc_case_histological_diagnosis == "Infiltrating Ductal Carcinoma"] # brca
# dat <- dat[, dat$cgc_case_histological_diagnosis == "Thyroid Papillary Carcinoma - Classical/usual"] # thca
# dat <- dat[, dat$cgc_case_histological_diagnosis == "Lung Adenocarcinoma- Not Otherwise Specified (NOS)"] # luad
# dat <- dat[, dat$cgc_case_histological_diagnosis == "Hepatocellular Carcinoma"] # lihc
# dat <- dat[, dat$cgc_case_histological_diagnosis == "Lung Squamous Cell Carcinoma- Not Otherwise Specified (NOS)" |
#              is.na(dat$cgc_case_histological_diagnosis)] # lusc
# dat <- dat[, dat$cgc_case_histological_diagnosis == "Colon Adenocarcinoma"] # coad
dat <- dat[which(rowMeans(apply(assay(dat), 2, function(x) x*1e6/sum(x))) > 0.5), ]

# Data
group <- dat$gdc_cases.samples.sample_type
group[group == "Solid Tissue Normal"] <- 1
group[group == "Primary Tumor"] <- 2
group <- factor(group)
counts <- assay(dat)
rm(dat)

# Normalise
libsizes <- colSums(counts)
nf <- calcNormFactors(counts, method="TMM")
els <- nf * libsizes
sf <- els / exp(mean(log(libsizes)))
norm.counts <- t(t(counts) / sf)

# expHM
expHM <- exp_hmm_adapt_3_chains(counts=t(norm.counts), groups=group)
lfc.mean.expHM <- unname(log2(colMeans(as.matrix(expHM$means2))) - log2(colMeans(as.matrix(expHM$means1))))
lfc.disp.expHM <- unname(log2(colMeans(as.matrix(expHM$disps2))) - log2(colMeans(as.matrix(expHM$disps1))))
rm(expHM)

# lnHM
lnHM <- ln_hmm_adapt_3_chains(counts=t(norm.counts), groups=group)
lfc.mean.lnHM <- unname(log2(colMeans(as.matrix(lnHM$means2))) - log2(colMeans(as.matrix(lnHM$means1))))
lfc.disp.lnHM <- unname(log2(colMeans(as.matrix(lnHM$disps2))) - log2(colMeans(as.matrix(lnHM$disps1))))
rm(lnHM)

results <- list(counts = counts, 
                lfc.mean.expHM = lfc.mean.expHM, 
                lfc.disp.expHM = lfc.disp.expHM, 
                lfc.mean.lnHM = lfc.mean.lnHM, 
                lfc.disp.lnHM = lfc.disp.lnHM)

if (cluster) {
  folder <- paste0("Data/TCGA paired data results May 2020")
} else {
  folder <- paste0("Results/TCGA paired data results May 2020")
}
saveRDS(results, file=here(folder, paste0("lfc.TCGA.paired_", i, ".rds")))
saveRDS(sessionInfo(), file=here(folder, paste0("sessionInfo.lfc.TCGA.paired_", i, ".rds")))

