# Testing different versions of HMM function to try to maximise size of posterior sample.
# Looks like best is, regardless of whether use mcmc.list(), to save raw results and import one by one to 
# process.

library(coda)
library(here)
library(recount)
library(edgeR)

## Run analyses and save raw results ####
# Run two (maybe 3 later) runs, save separately and then combine. May not have memory to combine 3 without summarising first.
dat <- readRDS(here("recount data/TCGA", "coad_matched_pairs.rds"))

# dat <- dat[, dat$cgc_case_histological_diagnosis == "Infiltrating Ductal Carcinoma"] # brca
# dat <- dat[, dat$cgc_case_histological_diagnosis == "Thyroid Papillary Carcinoma - Classical/usual"] # thca
# dat <- dat[, dat$cgc_case_histological_diagnosis == "Lung Adenocarcinoma- Not Otherwise Specified (NOS)"] # luad
# dat <- dat[, dat$cgc_case_histological_diagnosis == "Hepatocellular Carcinoma"] # lihc
# dat <- dat[, dat$cgc_case_histological_diagnosis == "Lung Squamous Cell Carcinoma- Not Otherwise Specified (NOS)" |
#              is.na(dat$cgc_case_histological_diagnosis)] # lusc
# dat <- dat[, dat$cgc_case_histological_diagnosis == "Colon Adenocarcinoma"] # coad

dat <- dat[which(rowMeans(apply(assay(dat), 2, function(x) x*1e6/sum(x))) > 0.5), ]
group <- dat$gdc_cases.samples.sample_type
group[group == "Solid Tissue Normal"] <- 1
group[group == "Primary Tumor"] <- 2
group <- factor(group)
counts <- assay(dat)
rm(dat)
libsizes <- colSums(counts)
nf <- calcNormFactors(counts, method="TMM")
els <- nf * libsizes
sf <- els / exp(mean(log(libsizes)))
norm.counts <- t(t(counts) / sf)
rm(counts, libsizes, nf, els, sf)

source(here("scripts", "2020-07-17_conditional_posterior_functions_exponential_hmm.R"))
source(here("scripts", "2020-07-23_exponential_hmm_one_chain_function.R"))
source(here("scripts", "2020-07-23_exponential_hmm_three_chains_function.R"))
source(here("scripts", "2020-07-25_exponential_hmm_adaptive_proposals_three_chains_function.R"))
source(here("scripts", "2020-09-02_run_expHMM.R"))
exp.res <- expHMM(t(norm.counts), groups=group, chain.length=3500, return.raw.results=T)
saveRDS(exp.res, here("Results/TCGA long chain results Aug 2020", "coad_raw_expHM_2.rds"))
rm(exp.res)
gc()
source(here("scripts", "2020-07-17_conditional_posterior_functions_lognormal_hmm.R"))
source(here("scripts", "2020-07-23_lognormal_hmm_one_chain_function.R"))
source(here("scripts", "2020-07-23_lognormal_hmm_three_chains_function.R"))
source(here("scripts", "2020-07-25_lognormal_hmm_adaptive_proposals_three_chains_function.R"))
source(here("scripts", "2020-09-02_run_lnHMM.R"))
ln.res <- lnHMM(t(norm.counts), groups=group, chain.length=3500, return.raw.results=T)
saveRDS(ln.res, here("Results/TCGA long chain results Aug 2020", "coad_raw_lnHM_2.rds"))
rm(ln.res)
gc()


## Import results and combine runs to make overall summaries ####
library(coda)
library(here)
folder <- "Results/TCGA long chain results Aug 2020"

i <- "brca"
for (j in c("expHM", "lnHM")) {
  means1_run1 <- readRDS(here(folder, paste0(i, "_raw_", j, ".rds")))$means1
  means1_run2 <- readRDS(here(folder, paste0(i, "_raw_", j, "_2.rds")))$means1
  means1 <- rbind(as.matrix(means1_run1), as.matrix(means1_run2))
  rm(means1_run1, means1_run2)
  saveRDS(means1, here(
    paste0(folder, "/Combined chain intermediate results"), 
    paste0(i, "_means1_", j, ".rds")
  ))
  rm(means1)
  means2_run1 <- readRDS(here(folder, paste0(i, "_raw_", j, ".rds")))$means2
  means2_run2 <- readRDS(here(folder, paste0(i, "_raw_", j, "_2.rds")))$means2
  means2 <- rbind(as.matrix(means2_run1), as.matrix(means2_run2))
  rm(means2_run1, means2_run2)
  saveRDS(means2, here(
    paste0(folder, "/Combined chain intermediate results"), 
    paste0(i, "_means2_", j, ".rds")
  ))
  rm(means2)
  disps1_run1 <- readRDS(here(folder, paste0(i, "_raw_", j, ".rds")))$disps1
  disps1_run2 <- readRDS(here(folder, paste0(i, "_raw_", j, "_2.rds")))$disps1
  disps1 <- rbind(as.matrix(disps1_run1), as.matrix(disps1_run2))
  rm(disps1_run1, disps1_run2)
  saveRDS(disps1, here(
    paste0(folder, "/Combined chain intermediate results"), 
    paste0(i, "_disps1_", j, ".rds")
  ))
  rm(disps1)
  disps2_run1 <- readRDS(here(folder, paste0(i, "_raw_", j, ".rds")))$disps2
  disps2_run2 <- readRDS(here(folder, paste0(i, "_raw_", j, "_2.rds")))$disps2
  disps2 <- rbind(as.matrix(disps2_run1), as.matrix(disps2_run2))
  rm(disps2_run1, disps2_run2)
  saveRDS(disps2, here(
    paste0(folder, "/Combined chain intermediate results"), 
    paste0(i, "_disps2_", j, ".rds")
  ))
  rm(disps2)
  raw_run1 <- readRDS(here(folder, paste0(i, "_raw_", j, ".rds")))
  prob_run1 <- as.matrix(raw_run1$indicators)
  post.prop_run1 <- as.matrix(raw_run1$proportion)
  rm(raw_run1)
  raw_run2 <- readRDS(here(folder, paste0(i, "_raw_", j, "_2.rds")))
  prob_run2 <- as.matrix(raw_run2$indicators)
  post.prop_run2 <- as.matrix(raw_run2$proportion)
  rm(raw_run2)
  prob <- unname(colMeans(rbind(prob_run1, prob_run2)))
  rm(prob_run1, prob_run2)
  saveRDS(prob, here(
    paste0(folder, "/Combined chain intermediate results"), 
    paste0(i, "_prob_", j, ".rds")
  ))
  rm(prob)
  post.prop <- mean(c(as.numeric(post.prop_run1), as.numeric(post.prop_run2)))
  rm(post.prop_run1, post.prop_run2)
  saveRDS(post.prop, here(
    paste0(folder, "/Combined chain intermediate results"), 
    paste0(i, "_post.prop_", j, ".rds")
  ))
  rm(post.prop)
}

## Import and process intermediate results and save separately ####
library(here)
folder <- "Results/TCGA long chain results Aug 2020"

i <- "lusc"
for (j in c("expHM", "lnHM")) {
  means1 <- readRDS(
    here(
      paste0(folder, "/Combined chain intermediate results"), 
      paste0(i, "_means1_", j, ".rds")
    )
  )
  means2 <- readRDS(
    here(
      paste0(folder, "/Combined chain intermediate results"), 
      paste0(i, "_means2_", j, ".rds")
    )
  )
  lfc.mean <- unname(log2(colMeans(means2)) - log2(colMeans(means1)))
  saveRDS(lfc.mean, here(
    paste0(folder, "/Combined chain intermediate results"), 
    paste0(i, "_lfc.mean_", j, ".rds")
  ))
  rm(lfc.mean)
  mean.diff.untr <- unname(means1 - means2)
  saveRDS(mean.diff.untr, here(
    paste0(folder, "/Combined chain intermediate results"), 
    paste0(i, "_mean.diff.untr_", j, ".rds")
  ))
  rm(mean.diff.untr)
  mean.diff.log <- unname(log(means1) - log(means2))
  saveRDS(mean.diff.log, here(
    paste0(folder, "/Combined chain intermediate results"), 
    paste0(i, "_mean.diff.log_", j, ".rds")
  ))
  rm(mean.diff.log)
  rm(means1, means2)
  
  disps1 <- readRDS(
    here(
      paste0(folder, "/Combined chain intermediate results"), 
      paste0(i, "_disps1_", j, ".rds")
    )
  )
  disps2 <- readRDS(
    here(
      paste0(folder, "/Combined chain intermediate results"), 
      paste0(i, "_disps2_", j, ".rds")
    )
  )
  lfc.disp <- unname(log2(colMeans(disps2)) - log2(colMeans(disps1)))
  saveRDS(lfc.disp, here(
    paste0(folder, "/Combined chain intermediate results"), 
    paste0(i, "_lfc.disp_", j, ".rds")
  ))
  rm(lfc.disp)
  disp.diff.untr <- unname(disps1 - disps2)
  saveRDS(disp.diff.untr, here(
    paste0(folder, "/Combined chain intermediate results"), 
    paste0(i, "_disp.diff.untr_", j, ".rds")
  ))
  rm(disp.diff.untr)
  disp.diff.log <- unname(log(disps1) - log(disps2))
  saveRDS(disp.diff.log, here(
    paste0(folder, "/Combined chain intermediate results"), 
    paste0(i, "_disp.diff.log_", j, ".rds")
  ))
  rm(disp.diff.log)
  rm(disps1, disps2)
}


## Import separate results and save in usual format for analysis ####
library(here)
source(here("scripts", "2020-07-23_hpd_tail_prob_function.R"))
source(here("scripts", "2019-05-03_bfdr_function.R"))
folder <- "Results/TCGA long chain results Aug 2020"

i <- "lusc"
for (j in c("expHM", "lnHM")) {
  mean.diff.untr <- readRDS(
    here(
      paste0(folder, "/Combined chain intermediate results"), 
      paste0(i, "_mean.diff.untr_", j, ".rds")
    )
  )
  p.mean.untr <- apply(mean.diff.untr,2,hpd.pval)
  rm(mean.diff.untr)
  q.mean.untr <- p.adjust(p.mean.untr, method='BH')
  
  mean.diff.log <- readRDS(
    here(
      paste0(folder, "/Combined chain intermediate results"), 
      paste0(i, "_mean.diff.log_", j, ".rds")
    )
  )
  p.mean.log <- apply(mean.diff.log,2,hpd.pval)
  rm(mean.diff.log)
  q.mean.log <- p.adjust(p.mean.log, method='BH')
  
  disp.diff.untr <- readRDS(
    here(
      paste0(folder, "/Combined chain intermediate results"), 
      paste0(i, "_disp.diff.untr_", j, ".rds")
    )
  )
  p.disp.untr <- apply(disp.diff.untr,2,hpd.pval)
  rm(disp.diff.untr)
  q.disp.untr <- p.adjust(p.disp.untr, method='BH')
  
  disp.diff.log <- readRDS(
    here(
      paste0(folder, "/Combined chain intermediate results"), 
      paste0(i, "_disp.diff.log_", j, ".rds")
    )
  )
  p.disp.log <- apply(disp.diff.log,2,hpd.pval)
  rm(disp.diff.log)
  q.disp.log <- p.adjust(p.disp.log, method='BH')
  
  lfc.mean <- readRDS(
    here(
      paste0(folder, "/Combined chain intermediate results"), 
      paste0(i, "_lfc.mean_", j, ".rds")
    )
  )
  
  lfc.disp <- readRDS(
    here(
      paste0(folder, "/Combined chain intermediate results"), 
      paste0(i, "_lfc.disp_", j, ".rds")
    )
  )
  
  prob <- readRDS(
    here(
      paste0(folder, "/Combined chain intermediate results"), 
      paste0(i, "_prob_", j, ".rds")
    )
  )
  
  post.prop <- readRDS(
    here(
      paste0(folder, "/Combined chain intermediate results"), 
      paste0(i, "_post.prop_", j, ".rds")
    )
  )
  
  thr <- sort(prob, decreasing=T)[round(length(prob) * post.prop)]
  BFDR <- bfdr(prob)
  
  res <- list(prob = prob, 
              prop = post.prop, 
              thr = thr, 
              bfdr = BFDR, 
              p.mean.untr = p.mean.untr, 
              p.mean.log = p.mean.log, 
              q.mean.untr = q.mean.untr, 
              q.mean.log = q.mean.log, 
              lfc.mean = lfc.mean, 
              p.disp.untr = p.disp.untr, 
              p.disp.log = p.disp.log, 
              q.disp.untr = q.disp.untr, 
              q.disp.log = q.disp.log, 
              lfc.disp = lfc.disp)
  
  saveRDS(res, file=here(folder, paste0("results.TCGA_20k_iter_", i, "_", j, ".rds")))
  
  rm(lfc.mean, p.mean.untr, p.mean.log, q.mean.untr, q.mean.log, 
     lfc.disp, p.disp.untr, p.disp.log, q.disp.untr, q.disp.log, 
     prob, post.prop, thr, BFDR, res)
}




