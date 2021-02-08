# check dimensions of TCGA data used to make sure processing consistently
library(here)
library(recount)

for (i in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "coad")) {
  assign(paste0("res.", i), 
         readRDS(here("Results/TCGA paired data results May 2020", 
                       paste0("results.TCGA.paired_", i, ".rds"))))
}
dim(res.brca$counts) # 23211 182
dim(res.kirc$counts) # 22507 144
dim(res.thca$counts) # 22406 98
dim(res.luad$counts) # 23416 86
dim(res.lihc$counts) # 19987 98
dim(res.lusc$counts) # 23283 92
dim(res.prad$counts) # 22717 92
dim(res.coad$counts) # 22559 78


dat <- readRDS(here("recount data/TCGA", "brca_matched_pairs.rds"))
dat <- dat[, dat$cgc_case_histological_diagnosis == "Infiltrating Ductal Carcinoma"]
dat <- dat[which(rowMeans(apply(assay(dat), 2, function(x) x*1e6/sum(x))) > 0.5), ]
dim(dat) # 23211 182
dat <- readRDS(here("recount data/TCGA", "kirc_matched_pairs.rds"))
dat <- dat[which(rowMeans(apply(assay(dat), 2, function(x) x*1e6/sum(x))) > 0.5), ]
dim(dat) # 22507 144
dat <- readRDS(here("recount data/TCGA", "thca_matched_pairs.rds"))
dat <- dat[, dat$cgc_case_histological_diagnosis == "Thyroid Papillary Carcinoma - Classical/usual"]
dat <- dat[which(rowMeans(apply(assay(dat), 2, function(x) x*1e6/sum(x))) > 0.5), ]
dim(dat) # 22406 98
dat <- readRDS(here("recount data/TCGA", "luad_matched_pairs.rds"))
dat <- dat[, dat$cgc_case_histological_diagnosis == "Lung Adenocarcinoma- Not Otherwise Specified (NOS)"]
dat <- dat[which(rowMeans(apply(assay(dat), 2, function(x) x*1e6/sum(x))) > 0.5), ]
dim(dat) # 23416 86
dat <- readRDS(here("recount data/TCGA", "lihc_matched_pairs.rds"))
dat <- dat[, dat$cgc_case_histological_diagnosis == "Hepatocellular Carcinoma"]
dat <- dat[which(rowMeans(apply(assay(dat), 2, function(x) x*1e6/sum(x))) > 0.5), ]
dim(dat) # 19987 98
dat <- readRDS(here("recount data/TCGA", "lusc_matched_pairs.rds"))
dat <- dat[, dat$cgc_case_histological_diagnosis == "Lung Squamous Cell Carcinoma- Not Otherwise Specified (NOS)" |
             is.na(dat$cgc_case_histological_diagnosis)]
dat <- dat[which(rowMeans(apply(assay(dat), 2, function(x) x*1e6/sum(x))) > 0.5), ]
dim(dat) # 23283 92
dat <- readRDS(here("recount data/TCGA", "prad_matched_pairs.rds"))
dat <- dat[which(rowMeans(apply(assay(dat), 2, function(x) x*1e6/sum(x))) > 0.5), ]
dim(dat) # 22717 92
dat <- readRDS(here("recount data/TCGA", "coad_matched_pairs.rds"))
dat <- dat[, dat$cgc_case_histological_diagnosis == "Colon Adenocarcinoma"]
dat <- dat[which(rowMeans(apply(assay(dat), 2, function(x) x*1e6/sum(x))) > 0.5), ]
dim(dat) # 22559 78

