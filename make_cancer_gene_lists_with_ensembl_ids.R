library(here)
library(recount)
library(org.Hs.eg.db)

gencode <- gsub('\\..*', '', names(recount_genes))
gene_info <- unique(select(org.Hs.eg.db, gencode, c('ENTREZID', 'SYMBOL','ENSEMBL'), 'ENSEMBL'))
folder <- "Data sources/Cancer-related genes"

# brca list
brca_sources <- readRDS(here(folder, "breast_gene_lists.rds"))
genes.brca <- c(brca_sources$cgc$Gene.Symbol, 
                brca_sources$disgenet$geneSymbol, 
                brca_sources$intogen$SYMBOL, 
                brca_sources$kegg$Symbol, 
                brca_sources$malacards$Symbol)
all.genes.brca <- unique(genes.brca)
length(all.genes.brca) # 1322
sum(all.genes.brca %in% gene_info$SYMBOL) # 1319
sum(unique(gene_info$SYMBOL) %in% all.genes.brca) # 1319
sum(gene_info$SYMBOL %in% all.genes.brca) # 1320
# 3 symbols in brca list not in recount list
# 1 symbol in brca list in recount list twice
gene_info.brca <- gene_info[which(gene_info$SYMBOL %in% all.genes.brca), ]
dim(gene_info.brca) # 1320 3
# Identify genes whose symbol in all.genes.brca isn't in gene_info symbols
# Some had only an ensemble ID
missing <- data.frame(cbind(all.genes.brca[-which(all.genes.brca %in% gene_info$SYMBOL)], 
                            c(rep("", length(all.genes.brca[-which(all.genes.brca %in% gene_info$SYMBOL)]))), 
                            c(rep("", length(all.genes.brca[-which(all.genes.brca %in% gene_info$SYMBOL)])))))
names(missing) <- names(gene_info.brca)
gene_info.brca <- rbind(gene_info.brca, missing)
dim(gene_info.brca) # 1323 3
# Identify genes whose symbol in gene_info matches more than one ensembl ID
dup <- gene_info.brca$SYMBOL[which(duplicated(gene_info.brca$SYMBOL))] # IRS4
gene_info[which(gene_info$SYMBOL %in% dup), ]
# One symbol with two ensemble IDs, partly overlapping on ensembl.org, so maybe different 
# transcripts of the same gene, but keep both.
rm(dup, missing, brca_sources)
dup.genes.brca <- unique(genes.brca[duplicated(genes.brca)])
length(dup.genes.brca) # 312
gene_info.brca.dup <- gene_info.brca[which(gene_info.brca$SYMBOL %in% dup.genes.brca), ]
dim(gene_info.brca.dup) # 312 3
saveRDS(gene_info.brca, here(folder, "brca_genes_info.rds"))
saveRDS(gene_info.brca.dup, here(folder, "brca_genes_at_least_two_dbs_info.rds"))

# kirc list
kirc_sources <- readRDS(here(folder, "renal_gene_lists.rds"))
genes.kirc <- c(kirc_sources$cgc$Gene.Symbol, 
                kirc_sources$disgenet$geneSymbol, 
                kirc_sources$intogen$SYMBOL, 
                kirc_sources$kegg$Symbol, 
                kirc_sources$malacards$Symbol)
all.genes.kirc <- unique(genes.kirc)
length(all.genes.kirc) # 297
sum(all.genes.kirc %in% gene_info$SYMBOL) # 294
sum(unique(gene_info$SYMBOL) %in% all.genes.kirc) # 294
sum(gene_info$SYMBOL %in% all.genes.kirc) # 295
# 3 symbols in kirc list not in recount list
# 1 symbol in kirc list in recount list more than once
gene_info.kirc <- gene_info[which(gene_info$SYMBOL %in% all.genes.kirc), ]
dim(gene_info.kirc) # 295 3
# Identify genes whose symbol in all.genes.kirc isn't in gene_info symbols
# Some had only an ensemble ID
missing <- data.frame(cbind(all.genes.kirc[-which(all.genes.kirc %in% gene_info$SYMBOL)], 
                            c(rep("", length(all.genes.kirc[-which(all.genes.kirc %in% gene_info$SYMBOL)]))), 
                            c(rep("", length(all.genes.kirc[-which(all.genes.kirc %in% gene_info$SYMBOL)])))))
names(missing) <- names(gene_info.kirc)
gene_info.kirc <- rbind(gene_info.kirc, missing)
dim(gene_info.kirc) # 298 3
# Identify genes whose symbol in gene_info matches more than one ensembl ID
dup <- gene_info.kirc$SYMBOL[which(duplicated(gene_info.kirc$SYMBOL))] # BUB1B-PAK6
gene_info[which(gene_info$SYMBOL %in% dup), ]
# Looks to be one gene and a read-through from that gene (PAK6) to another (BUB1B).
gene_info.kirc[grep("PAK6", gene_info.kirc$SYMBOL), ]
# Entry with symbol PAK6 has same Ensemble ID as one of the BUB1B-PAK6 entries, but a different Entez ID.
# Remove entry that seems to be PAK6 but with Entrez ID and Symbol for readthrough.
gene_info.kirc <- gene_info.kirc[-which(
  gene_info.kirc$ENSEMBL == "ENSG00000137843" & gene_info.kirc$SYMBOL == "BUB1B-PAK6"
), ]
dim(gene_info.kirc) # 297 3
rm(missing, kirc_sources)
dup.genes.kirc <- unique(genes.kirc[duplicated(genes.kirc)])
length(dup.genes.kirc) # 37
gene_info.kirc.dup <- gene_info.kirc[which(gene_info.kirc$SYMBOL %in% dup.genes.kirc), ]
dim(gene_info.kirc.dup) # 37 3
saveRDS(gene_info.kirc, here(folder, "kirc_genes_info.rds"))
saveRDS(gene_info.kirc.dup, here(folder, "kirc_genes_at_least_two_dbs_info.rds"))

# thca list
thca_sources <- readRDS(here(folder, "thyroid_gene_lists.rds"))
genes.thca <- c(thca_sources$cgc$Gene.Symbol, 
                thca_sources$disgenet$geneSymbol, 
                thca_sources$intogen$SYMBOL, 
                thca_sources$kegg$Symbol, 
                thca_sources$malacards$Symbol)
all.genes.thca <- unique(genes.thca)
length(all.genes.thca) # 265
sum(all.genes.thca %in% gene_info$SYMBOL) # 265
sum(unique(gene_info$SYMBOL) %in% all.genes.thca) # 265
sum(gene_info$SYMBOL %in% all.genes.thca) # 265
# no symbols missing or duplicated
gene_info.thca <- gene_info[which(gene_info$SYMBOL %in% all.genes.thca), ]
dim(gene_info.thca) # 265 3
rm(thca_sources)
dup.genes.thca <- unique(genes.thca[duplicated(genes.thca)])
length(dup.genes.thca) # 40
gene_info.thca.dup <- gene_info.thca[which(gene_info.thca$SYMBOL %in% dup.genes.thca), ]
dim(gene_info.thca.dup) # 40 3
saveRDS(gene_info.thca, here(folder, "thca_genes_info.rds"))
saveRDS(gene_info.thca.dup, here(folder, "thca_genes_at_least_two_dbs_info.rds"))

# luad list
luad_sources <- readRDS(here(folder, "lungad_gene_lists.rds"))
genes.luad <- c(luad_sources$cgc$Gene.Symbol, 
                luad_sources$disgenet$geneSymbol, 
                luad_sources$intogen$SYMBOL, 
                luad_sources$kegg$Symbol)
all.genes.luad <- unique(genes.luad)
length(all.genes.luad) # 306
sum(all.genes.luad %in% gene_info$SYMBOL) # 306
sum(unique(gene_info$SYMBOL) %in% all.genes.luad) # 306
sum(gene_info$SYMBOL %in% all.genes.luad) # 306
# no symbols missing or duplicated
gene_info.luad <- gene_info[which(gene_info$SYMBOL %in% all.genes.luad), ]
dim(gene_info.luad) # 306 3
dup.genes.luad <- unique(genes.luad[duplicated(genes.luad)])
length(dup.genes.luad) # 26
gene_info.luad.dup <- gene_info.luad[which(gene_info.luad$SYMBOL %in% dup.genes.luad), ]
dim(gene_info.luad.dup) # 26 3
saveRDS(gene_info.luad, here(folder, "luad_genes_info.rds"))
saveRDS(gene_info.luad.dup, here(folder, "luad_genes_at_least_two_dbs_info.rds"))

# lihc list
lihc_sources <- readRDS(here(folder, "hepatocellular_gene_lists.rds"))
genes.lihc <- c(lihc_sources$cgc$Gene.Symbol, 
                lihc_sources$disgenet$geneSymbol, 
                lihc_sources$intogen$SYMBOL, 
                lihc_sources$kegg$Symbol, 
                lihc_sources$malacards$Symbol)
all.genes.lihc <- unique(genes.lihc)
length(all.genes.lihc) # 781
sum(all.genes.lihc %in% gene_info$SYMBOL) # 771
sum(unique(gene_info$SYMBOL) %in% all.genes.lihc) # 771
sum(gene_info$SYMBOL %in% all.genes.lihc) # 773
# 10 symbols in lihc list not in recount list
# 2 duplications of symbols in lihc list in recount list
gene_info.lihc <- gene_info[which(gene_info$SYMBOL %in% all.genes.lihc), ]
dim(gene_info.lihc) # 773 3
# Identify genes whose symbol in all.genes.lihc isn't in gene_info symbols
# Some had only an ensemble ID
missing <- data.frame(cbind(all.genes.lihc[-which(all.genes.lihc %in% gene_info$SYMBOL)], 
                            c(rep("", length(all.genes.lihc[-which(all.genes.lihc %in% gene_info$SYMBOL)]))), 
                            c(rep("", length(all.genes.lihc[-which(all.genes.lihc %in% gene_info$SYMBOL)])))))
names(missing) <- names(gene_info.lihc)
gene_info.lihc <- rbind(gene_info.lihc, missing)
dim(gene_info.lihc) # 783 3
# Identify genes whose symbol in gene_info matches more than one ensembl ID
dup <- gene_info.lihc$SYMBOL[which(duplicated(gene_info.lihc$SYMBOL))] # IFNAR2, PINX1
gene_info[which(gene_info$SYMBOL %in% dup), ]
# Two symbols with two ensemble IDs, partly overlapping on ensembl.org, so maybe different 
# transcripts of the same gene, but keep both.
rm(missing, dup, lihc_sources)
dup.genes.lihc <- unique(genes.lihc[duplicated(genes.lihc)])
length(dup.genes.lihc) # 123
gene_info.lihc.dup <- gene_info.lihc[which(gene_info.lihc$SYMBOL %in% dup.genes.lihc), ]
dim(gene_info.lihc.dup) # 123 3
saveRDS(gene_info.lihc, here(folder, "lihc_genes_info.rds"))
saveRDS(gene_info.lihc.dup, here(folder, "lihc_genes_at_least_two_dbs_info.rds"))

# lusc list
lusc_sources <- readRDS(here(folder, "lungsq_gene_lists.rds"))
genes.lusc <- c(lusc_sources$cgc$Gene.Symbol, 
                lusc_sources$disgenet$geneSymbol, 
                lusc_sources$intogen$SYMBOL, 
                lusc_sources$kegg$Symbol, 
                lusc_sources$malacards$Symbol)
all.genes.lusc <- unique(genes.lusc)
length(all.genes.lusc) # 266
sum(all.genes.lusc %in% gene_info$SYMBOL) # 266
sum(unique(gene_info$SYMBOL) %in% all.genes.lusc) # 266
sum(gene_info$SYMBOL %in% all.genes.lusc) # 266
# no symbols missing or duplicated
gene_info.lusc <- gene_info[which(gene_info$SYMBOL %in% all.genes.lusc), ]
dim(gene_info.lusc) # 266 3
dup.genes.lusc <- unique(genes.lusc[duplicated(genes.lusc)])
length(dup.genes.lusc) # 46
gene_info.lusc.dup <- gene_info.lusc[which(gene_info.lusc$SYMBOL %in% dup.genes.lusc), ]
dim(gene_info.lusc.dup) # 46 3
saveRDS(gene_info.lusc, here(folder, "lusc_genes_info.rds"))
saveRDS(gene_info.lusc.dup, here(folder, "lusc_genes_at_least_two_dbs_info.rds"))

# prad list
prad_sources <- readRDS(here(folder, "prostate_gene_lists.rds"))
genes.prad <- c(prad_sources$cgc$Gene.Symbol, 
                prad_sources$disgenet$geneSymbol, 
                prad_sources$intogen$SYMBOL, 
                prad_sources$kegg$Symbol, 
                prad_sources$malacards$Symbol)
all.genes.prad <- unique(genes.prad)
length(all.genes.prad) # 1125
sum(all.genes.prad %in% gene_info$SYMBOL) # 1124
sum(unique(gene_info$SYMBOL) %in% all.genes.prad) # 1124
sum(gene_info$SYMBOL %in% all.genes.prad) # 1125
# 1 symbol in prad list not in recount list
# 1 symbol in prad list in recount list twice
gene_info.prad <- gene_info[which(gene_info$SYMBOL %in% all.genes.prad), ]
dim(gene_info.prad) # 1125 3
# Identify genes whose symbol in all.genes.prad isn't in gene_info symbols
# Some had only an ensemble ID
missing <- data.frame(cbind(all.genes.prad[-which(all.genes.prad %in% gene_info$SYMBOL)], 
                            c(rep("", length(all.genes.prad[-which(all.genes.prad %in% gene_info$SYMBOL)]))), 
                            c(rep("", length(all.genes.prad[-which(all.genes.prad %in% gene_info$SYMBOL)])))))
names(missing) <- names(gene_info.prad)
gene_info.prad <- rbind(gene_info.prad, missing)
dim(gene_info.prad) # 1126 3
# Identify genes whose symbol in gene_info matches more than one ensembl ID
dup <- gene_info.prad$SYMBOL[which(duplicated(gene_info.prad$SYMBOL))] # TBC1D3
gene_info[which(gene_info$SYMBOL %in% dup), ]
# Different ensembl IDs, same Entrez
grep(dup, prad_sources$cgc$Gene.Symbol)
grep(dup, prad_sources$disgenet$geneSymbol)
grep(dup, prad_sources$intogen$SYMBOL)
grep(dup, prad_sources$kegg)
grep(dup, prad_sources$malacards$Symbol)
# TBC1D3 only in MalaCards list so only have symbol and description
# One ensembl ID (ENSG00000274512) corresponds to a different Entrez ID and symbol on NCBI 
# (TBC1D3L), so exclude that ensembl ID.
gene_info.prad <- gene_info.prad[-which(gene_info.prad$ENSEMBL == "ENSG00000274512"), ]
dim(gene_info.prad) # 1125 3
rm(missing, dup, prad_sources)
dup.genes.prad <- unique(genes.prad[duplicated(genes.prad)])
length(dup.genes.prad) # 267
gene_info.prad.dup <- gene_info.prad[which(gene_info.prad$SYMBOL %in% dup.genes.prad), ]
dim(gene_info.prad.dup) # 267 3
saveRDS(gene_info.prad, here(folder, "prad_genes_info.rds"))
saveRDS(gene_info.prad.dup, here(folder, "prad_genes_at_least_two_dbs_info.rds"))

# hnsc not in KEGG so can't do pathway genes so didn't include for initial gene lists

# coad list
coad_sources <- readRDS(here(folder, "colon_gene_lists.rds"))
genes.coad <- c(coad_sources$cgc$Gene.Symbol, 
                coad_sources$disgenet$geneSymbol, 
                coad_sources$intogen$SYMBOL, 
                coad_sources$kegg$Symbol, 
                coad_sources$malacards$Symbol)
all.genes.coad <- unique(genes.coad)
length(all.genes.coad) # 1437
sum(all.genes.coad %in% gene_info$SYMBOL) # 1431
sum(unique(gene_info$SYMBOL) %in% all.genes.coad) # 1431
sum(gene_info$SYMBOL %in% all.genes.coad) # 1433
# 6 symbols in coad list not in recount list
# 2 duplications of symbols in coad list in recount list
gene_info.coad <- gene_info[which(gene_info$SYMBOL %in% all.genes.coad), ]
dim(gene_info.coad) # 1433 3
# Identify genes whose symbol in all.genes.kirc isn't in gene_info symbols
# Some had only an ensemble ID
missing <- data.frame(cbind(all.genes.coad[-which(all.genes.coad %in% gene_info$SYMBOL)], 
                            c(rep("", length(all.genes.coad[-which(all.genes.coad %in% gene_info$SYMBOL)]))), 
                            c(rep("", length(all.genes.coad[-which(all.genes.coad %in% gene_info$SYMBOL)])))))
names(missing) <- names(gene_info.coad)
gene_info.coad <- rbind(gene_info.coad, missing)
dim(gene_info.coad) # 1439 3
# Identify genes whose symbol in gene_info matches more than one ensembl ID
dup <- gene_info.coad$SYMBOL[which(duplicated(gene_info.coad$SYMBOL))] # IRS4, RAD54B
grep(dup[1], coad_sources$disgenet$geneSymbol)
grep(dup[2], coad_sources$disgenet$geneSymbol)
grep(dup[2], coad_sources$malacards$Symbol)
coad_sources$disgenet[which(coad_sources$disgenet$geneSymbol == dup[1]), ]
coad_sources$disgenet[which(coad_sources$disgenet$geneSymbol == dup[2]), ]
coad_sources$malacards[which(coad_sources$malacards$Symbol == dup[2]), ]
# DisGeNET only has symbol and Entrez ID, MalaCards only has symbol and description
gene_info[which(gene_info$SYMBOL %in% dup), ]
# For IRS4, ENSG00000260548 very slightly overlaps with ENSG00000133124 and has a different 
# name on ensembl.org, so exclude.
# For RAD54B, ENSG00000265817 overlaps with ENSG00000197275 but has a different name on
# ensembl.org, so exclude.
gene_info.coad <- gene_info.coad[-which(gene_info.coad$ENSEMBL %in% 
                                          c("ENSG00000260548", "ENSG00000265817")), ]
dim(gene_info.coad) # 1437 3
rm(missing, dup, coad_sources)
dup.genes.coad <- unique(genes.coad[duplicated(genes.coad)])
length(dup.genes.coad) # 280
gene_info.coad.dup <- gene_info.coad[which(gene_info.coad$SYMBOL %in% dup.genes.coad), ]
dim(gene_info.coad.dup) # 280 3
saveRDS(gene_info.coad, here(folder, "coad_genes_info.rds"))
saveRDS(gene_info.coad.dup, here(folder, "coad_genes_at_least_two_dbs_info.rds"))
