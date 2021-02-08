library(here)
library(dplyr)
library(recount)
library(org.Hs.eg.db)

folder <- "Data sources/Cancer-related pathway genes"

# Recount genes following 
# http://bioconductor.org/packages/release/bioc/vignettes/recount/inst/doc/recount-quickstart.html
gencode <- gsub('\\..*', '', names(recount_genes))
gene_info <- select(org.Hs.eg.db, 
                    gencode, 
                    c('ENSEMBL', 'ENTREZID', 'SYMBOL', 'GENENAME', 'ALIAS'),
                    'ENSEMBL')
# dim(gene_info) # 115974 5


# Breast
pathway_genes_breast <- read.csv(here("Data sources/Cancer-related pathway genes", 
                                      "KEGG_pathway_genes_Breast_cancer.csv"), 
                                 stringsAsFactors=F)
dim(pathway_genes_breast) # 1236 3
all.genes.brca <- unique(pathway_genes_breast)
dim(all.genes.brca) # 974 3
# Identify and remove genes in brca list missing from recount list
missing <- all.genes.brca[-which(all.genes.brca$EntrezID %in% gene_info$ENTREZID), ]
dim(missing) # 3 3
all.genes.brca <- all.genes.brca[which(all.genes.brca$EntrezID %in% gene_info$ENTREZID), ]
dim(all.genes.brca) # 971 3
rm(missing)
gene_info.brca <- unique(gene_info[which(gene_info$ENTREZID %in% all.genes.brca$EntrezID), 1:4])
dim(gene_info.brca) # 974 4
# 3 duplicates
gene_info.brca[which(gene_info.brca$ENTREZID %in% gene_info.brca$ENTREZID[duplicated(gene_info.brca$ENTREZID)]), ]
# 3 genes with two Ensemble IDs. Other fields all match so assume no error and keep all.
dup.genes.brca <- unique(pathway_genes_breast$EntrezID[duplicated(pathway_genes_breast$EntrezID)])[
  which(unique(pathway_genes_breast$EntrezID[duplicated(pathway_genes_breast$EntrezID)]) %in% gene_info.brca$ENTREZID)
  ]
length(dup.genes.brca) # 211
gene_info.brca.dup <- gene_info.brca[which(gene_info.brca$ENTREZID %in% dup.genes.brca), ]
dim(gene_info.brca.dup) # 211 4
write.csv(gene_info.brca, here(folder, "brca_pathway_genes_recount.csv"), row.names=F)
write.csv(gene_info.brca.dup, here(folder, "brca_pathway_genes_recount_multiple_pathways.csv"), row.names=F)
rm(gene_info.brca, gene_info.brca.dup, dup.genes.brca, all.genes.brca, pathway_genes_breast)

# Renal
pathway_genes_renal <- read.csv(here("Data sources/Cancer-related pathway genes", 
                                      "KEGG_pathway_genes_Renal_cell_carcinoma.csv"), 
                                 stringsAsFactors=F)
dim(pathway_genes_renal) # 726 3
all.genes.kirc <- unique(pathway_genes_renal)
dim(all.genes.kirc) # 622 3
# Identify and remove genes in kirc list missing from recount list
missing <- all.genes.kirc[-which(all.genes.kirc$EntrezID %in% gene_info$ENTREZID), ]
dim(missing) # 0 3
all.genes.kirc <- all.genes.kirc[which(all.genes.kirc$EntrezID %in% gene_info$ENTREZID), ]
dim(all.genes.kirc) # 622 3
rm(missing)
gene_info.kirc <- unique(gene_info[which(gene_info$ENTREZID %in% all.genes.kirc$EntrezID), 1:4])
dim(gene_info.kirc) # 624 4
# 2 duplicates
gene_info.kirc[which(gene_info.kirc$ENTREZID %in% gene_info.kirc$ENTREZID[duplicated(gene_info.kirc$ENTREZID)]), ]
# 2 genes with two Ensemble IDs. Other fields all match so assume no error and keep all.
dup.genes.kirc <- unique(pathway_genes_renal$EntrezID[duplicated(pathway_genes_renal$EntrezID)])[
  which(unique(pathway_genes_renal$EntrezID[duplicated(pathway_genes_renal$EntrezID)]) %in% gene_info.kirc$ENTREZID)
  ]
length(dup.genes.kirc) # 90
gene_info.kirc.dup <- gene_info.kirc[which(gene_info.kirc$ENTREZID %in% dup.genes.kirc), ]
dim(gene_info.kirc.dup) # 90 4
write.csv(gene_info.kirc, here(folder, "kirc_pathway_genes_recount.csv"), row.names=F)
write.csv(gene_info.kirc.dup, here(folder, "kirc_pathway_genes_recount_multiple_pathways.csv"), row.names=F)
rm(gene_info.kirc, gene_info.kirc.dup, dup.genes.kirc, all.genes.kirc, pathway_genes_renal)

# Lung adenocarcinoma
pathway_genes_lungad <- read.csv(here("Data sources/Cancer-related pathway genes", 
                                      "KEGG_pathway_genes_Non_small_cell_lung_cancer.csv"), 
                                 stringsAsFactors=F)
dim(pathway_genes_lungad) # 1354 3
all.genes.luad <- unique(pathway_genes_lungad)
dim(all.genes.luad) # 922 3
# Identify and remove genes in luad list missing from recount list
missing <- all.genes.luad[-which(all.genes.luad$EntrezID %in% gene_info$ENTREZID), ]
dim(missing) # 1 3
all.genes.luad <- all.genes.luad[which(all.genes.luad$EntrezID %in% gene_info$ENTREZID), ]
dim(all.genes.luad) # 921 3
rm(missing)
gene_info.luad <- unique(gene_info[which(gene_info$ENTREZID %in% all.genes.luad$EntrezID), 1:4])
dim(gene_info.luad) # 924 4
# 3 duplicates
gene_info.luad[which(gene_info.luad$ENTREZID %in% gene_info.luad$ENTREZID[duplicated(gene_info.luad$ENTREZID)]), ]
# 3 genes with two Ensemble IDs. Other fields all match so assume no error and keep all.
dup.genes.luad <- unique(pathway_genes_lungad$EntrezID[duplicated(pathway_genes_lungad$EntrezID)])[
  which(unique(pathway_genes_lungad$EntrezID[duplicated(pathway_genes_lungad$EntrezID)]) %in% gene_info.luad$ENTREZID)
  ]
length(dup.genes.luad) # 264
gene_info.luad.dup <- gene_info.luad[which(gene_info.luad$ENTREZID %in% dup.genes.luad), ]
dim(gene_info.luad.dup) # 265 4
gene_info.luad.dup[which(gene_info.luad.dup$ENTREZID %in% gene_info.luad.dup$ENTREZID[duplicated(gene_info.luad.dup$ENTREZID)]), ]
# 1 gene with two Ensemble IDs, already considered above.
write.csv(gene_info.luad, here(folder, "luad_pathway_genes_recount.csv"), row.names=F)
write.csv(gene_info.luad.dup, here(folder, "luad_pathway_genes_recount_multiple_pathways.csv"), row.names=F)
rm(gene_info.luad, gene_info.luad.dup, dup.genes.luad, all.genes.luad, pathway_genes_lungad)

# Hepatocellular
pathway_genes_hepatocellular <- read.csv(here("Data sources/Cancer-related pathway genes", 
                                      "KEGG_pathway_genes_Hepatocellular_carcinoma.csv"), 
                                 stringsAsFactors=F)
dim(pathway_genes_hepatocellular) # 1757 3
all.genes.lihc <- unique(pathway_genes_hepatocellular)
dim(all.genes.lihc) # 1205 3
# Identify and remove genes in lihc list missing from recount list
missing <- all.genes.lihc[-which(all.genes.lihc$EntrezID %in% gene_info$ENTREZID), ]
dim(missing) # 2 3
all.genes.lihc <- all.genes.lihc[which(all.genes.lihc$EntrezID %in% gene_info$ENTREZID), ]
dim(all.genes.lihc) # 1203 3
rm(missing)
gene_info.lihc <- unique(gene_info[which(gene_info$ENTREZID %in% all.genes.lihc$EntrezID), 1:4])
dim(gene_info.lihc) # 1205 4
# 2 duplicates
gene_info.lihc[which(gene_info.lihc$ENTREZID %in% gene_info.lihc$ENTREZID[duplicated(gene_info.lihc$ENTREZID)]), ]
# 2 genes with two Ensemble IDs. Other fields all match so assume no error and keep all.
dup.genes.lihc <- unique(pathway_genes_hepatocellular$EntrezID[duplicated(pathway_genes_hepatocellular$EntrezID)])[
  which(unique(pathway_genes_hepatocellular$EntrezID[duplicated(pathway_genes_hepatocellular$EntrezID)]) %in% gene_info.lihc$ENTREZID)
]
length(dup.genes.lihc) # 327
gene_info.lihc.dup <- gene_info.lihc[which(gene_info.lihc$ENTREZID %in% dup.genes.lihc), ]
dim(gene_info.lihc.dup) # 329 4
gene_info.lihc.dup[which(gene_info.lihc.dup$ENTREZID %in% gene_info.lihc.dup$ENTREZID[duplicated(gene_info.lihc.dup$ENTREZID)]), ]
# 2 genes with two Ensemble IDs, already considered above.
write.csv(gene_info.lihc, here(folder, "lihc_pathway_genes_recount.csv"), row.names=F)
write.csv(gene_info.lihc.dup, here(folder, "lihc_pathway_genes_recount_multiple_pathways.csv"), row.names=F)
rm(gene_info.lihc, gene_info.lihc.dup, dup.genes.lihc, all.genes.lihc, pathway_genes_hepatocellular)

# Lung squamous
gene_info.lusc <- read.csv(here(folder, "luad_pathway_genes_recount.csv"), stringsAsFactors=F)
gene_info.lusc.dup <- read.csv(here(folder, "luad_pathway_genes_recount_multiple_pathways.csv"), stringsAsFactors=F)
write.csv(gene_info.lusc, here(folder, "lusc_pathway_genes_recount.csv"), row.names=F)
write.csv(gene_info.lusc.dup, here(folder, "lusc_pathway_genes_recount_multiple_pathways.csv"), row.names=F)
rm(gene_info.lusc, gene_info.lusc.dup)

# Thyroid
pathway_genes_thyroid <- read.csv(here("Data sources/Cancer-related pathway genes", 
                                              "KEGG_pathway_genes_thyroid_cancer.csv"), 
                                         stringsAsFactors=F)
dim(pathway_genes_thyroid) # 674 3
all.genes.thca <- unique(pathway_genes_thyroid)
dim(all.genes.thca) # 607 3
# Identify and remove genes in thca list missing from recount list
missing <- all.genes.thca[-which(all.genes.thca$EntrezID %in% gene_info$ENTREZID), ]
dim(missing) # 2 3
all.genes.thca <- all.genes.thca[which(all.genes.thca$EntrezID %in% gene_info$ENTREZID), ]
dim(all.genes.thca) # 606 3
rm(missing)
gene_info.thca <- unique(gene_info[which(gene_info$ENTREZID %in% all.genes.thca$EntrezID), 1:4])
dim(gene_info.thca) # 608 4
# 2 duplicates
gene_info.thca[which(gene_info.thca$ENTREZID %in% gene_info.thca$ENTREZID[duplicated(gene_info.thca$ENTREZID)]), ]
# 2 genes with two Ensemble IDs. Other fields all match so assume no error and keep all.
dup.genes.thca <- unique(pathway_genes_thyroid$EntrezID[duplicated(pathway_genes_thyroid$EntrezID)])[
  which(unique(pathway_genes_thyroid$EntrezID[duplicated(pathway_genes_thyroid$EntrezID)]) %in% gene_info.thca$ENTREZID)
  ]
length(dup.genes.thca) # 61
gene_info.thca.dup <- gene_info.thca[which(gene_info.thca$ENTREZID %in% dup.genes.thca), ]
dim(gene_info.thca.dup) # 61 4
write.csv(gene_info.thca, here(folder, "thca_pathway_genes_recount.csv"), row.names=F)
write.csv(gene_info.thca.dup, here(folder, "thca_pathway_genes_recount_multiple_pathways.csv"), row.names=F)
rm(gene_info.thca, gene_info.thca.dup, dup.genes.thca, all.genes.thca, pathway_genes_thyroid)

# Prostate
pathway_genes_prostate <- read.csv(here("Data sources/Cancer-related pathway genes", 
                                              "KEGG_pathway_genes_prostate_cancer.csv"), 
                                         stringsAsFactors=F)
dim(pathway_genes_prostate) # 1521 3
all.genes.prad <- unique(pathway_genes_prostate)
dim(all.genes.prad) # 1178 3
# Identify and remove genes in prad list missing from recount list
missing <- all.genes.prad[-which(all.genes.prad$EntrezID %in% gene_info$ENTREZID), ]
dim(missing) # 1 3
all.genes.prad <- all.genes.prad[which(all.genes.prad$EntrezID %in% gene_info$ENTREZID), ]
dim(all.genes.prad) # 1177 3
rm(missing)
gene_info.prad <- unique(gene_info[which(gene_info$ENTREZID %in% all.genes.prad$EntrezID), 1:4])
dim(gene_info.prad) # 1180 4
# 3 duplicates
gene_info.prad[which(gene_info.prad$ENTREZID %in% gene_info.prad$ENTREZID[duplicated(gene_info.prad$ENTREZID)]), ]
# 3 genes with two Ensemble IDs. Other fields all match so assume no error and keep all.
dup.genes.prad <- unique(pathway_genes_prostate$EntrezID[duplicated(pathway_genes_prostate$EntrezID)])[
  which(unique(pathway_genes_prostate$EntrezID[duplicated(pathway_genes_prostate$EntrezID)]) %in% gene_info.prad$ENTREZID)
  ]
length(dup.genes.prad) # 253
gene_info.prad.dup <- gene_info.prad[which(gene_info.prad$ENTREZID %in% dup.genes.prad), ]
dim(gene_info.prad.dup) # 254 4
gene_info.prad.dup[which(gene_info.prad.dup$ENTREZID %in% gene_info.prad.dup$ENTREZID[duplicated(gene_info.prad.dup$ENTREZID)]), ]
# 1 genes with two Ensemble IDs, already considered above.
write.csv(gene_info.prad, here(folder, "prad_pathway_genes_recount.csv"), row.names=F)
write.csv(gene_info.prad.dup, here(folder, "prad_pathway_genes_recount_multiple_pathways.csv"), row.names=F)
rm(gene_info.prad, gene_info.prad.dup, dup.genes.prad, all.genes.prad, pathway_genes_prostate)

# Colon
pathway_genes_colon <- read.csv(here("Data sources/Cancer-related pathway genes", 
                                              "KEGG_pathway_genes_Colorectal_cancer.csv"), 
                                         stringsAsFactors=F)
dim(pathway_genes_colon) # 1474 3
all.genes.coad <- unique(pathway_genes_colon)
dim(all.genes.coad) # 1024 3
# Identify and remove genes in coad list missing from recount list
missing <- all.genes.coad[-which(all.genes.coad$EntrezID %in% gene_info$ENTREZID), ]
dim(missing) # 2 3
all.genes.coad <- all.genes.coad[which(all.genes.coad$EntrezID %in% gene_info$ENTREZID), ]
dim(all.genes.coad) # 1022 3
rm(missing)
gene_info.coad <- unique(gene_info[which(gene_info$ENTREZID %in% all.genes.coad$EntrezID), 1:4])
dim(gene_info.coad) # 1026 4
# 4 duplicates
gene_info.coad[which(gene_info.coad$ENTREZID %in% gene_info.coad$ENTREZID[duplicated(gene_info.coad$ENTREZID)]), ]
# 4 genes with two Ensemble IDs. Other fields all match so assume no error and keep all.
dup.genes.coad <- unique(pathway_genes_colon$EntrezID[duplicated(pathway_genes_colon$EntrezID)])[
  which(unique(pathway_genes_colon$EntrezID[duplicated(pathway_genes_colon$EntrezID)]) %in% gene_info.coad$ENTREZID)
  ]
length(dup.genes.coad) # 293
gene_info.coad.dup <- gene_info.coad[which(gene_info.coad$ENTREZID %in% dup.genes.coad), ]
dim(gene_info.coad.dup) # 293 4
write.csv(gene_info.coad, here(folder, "coad_pathway_genes_recount.csv"), row.names=F)
write.csv(gene_info.coad.dup, here(folder, "coad_pathway_genes_recount_multiple_pathways.csv"), row.names=F)
rm(gene_info.coad, gene_info.coad.dup, dup.genes.coad, all.genes.coad, pathway_genes_colon)

