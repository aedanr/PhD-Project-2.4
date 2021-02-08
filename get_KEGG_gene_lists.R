library(here)
library(KEGGREST)

# Get list of genes for each cancer pathway ####
cancer_pathways <- c("hsa05224", "hsa05211", "hsa05223", "hsa05225", "hsa05216", "hsa05215", "hsa05210")
for (i in cancer_pathways) {
  list <- keggGet(i)[[1]]
  gene <- list$GENE
  entrez <- as.numeric(gene[seq(1, length(gene) - 1, 2)])
  symbols_names <- gene[seq(2, length(gene), 2)]
  symbols <- substr(symbols_names, 1, unlist(gregexpr(";", symbols_names)) - 1)
  names <- substr(
    symbols_names, 
    unlist(gregexpr(";", symbols_names)) + 2, 
    unlist(lapply(symbols_names, function(x) unlist(gregexpr("\\[", x))[1])) - 2
  )
  assign(
    gsub(
      " |-", "_", gsub("\\(|\\)| - Homo sapiens \\(human\\)", "", list$NAME)
    ), 
    data.frame(
      EntrezID = entrez, 
      Symbol = symbols, 
      Name = names
    )
  )
  rm(list, gene, entrez, symbols_names, symbols, names)
}
rm(i, cancer_pathways)

ls()
c(dim(Breast_cancer), dim(unique(Breast_cancer))) # 147
c(dim(Renal_cell_carcinoma), dim(unique(Renal_cell_carcinoma))) # 69
c(dim(Non_small_cell_lung_cancer), dim(unique(Non_small_cell_lung_cancer))) # 68
c(dim(Hepatocellular_carcinoma), dim(unique(Hepatocellular_carcinoma))) # 168
c(dim(Thyroid_cancer), dim(unique(Thyroid_cancer))) # 37
c(dim(Prostate_cancer), dim(unique(Prostate_cancer))) # 97
c(dim(Colorectal_cancer), dim(unique(Colorectal_cancer))) # 86

for (i in ls()) {
  write.csv(get(i), here("Data sources/Cancer-related genes/KEGG", paste0("KEGG_", i, ".csv")), row.names=F)
}
rm(list=ls())


# Get list of related pathways for each cancer pathway ####
cancer_pathways <- c("hsa05224", "hsa05211", "hsa05223", "hsa05225", "hsa05216", "hsa05215", "hsa05210")
for (i in cancer_pathways) {
  list <- keggGet(i)[[1]]
  pathway <- list$REL_PATHWAY
  assign(
    gsub(
      " |-", "_", gsub("\\(|\\)| - Homo sapiens \\(human\\)", "", list$NAME)
    ), 
    data.frame(
      KEGGID = names(pathway), 
      Pathway = unname(pathway), 
      Cancer = gsub(
        " |-", "_", gsub("\\(|\\)| - Homo sapiens \\(human\\)", "", list$NAME)
      )
    )
  )
  rm(list, pathway)
}
rm(i, cancer_pathways)

ls()
dim(Breast_cancer) # 8 3
dim(Renal_cell_carcinoma) # 6 3
dim(Non_small_cell_lung_cancer) # 7 3
dim(Hepatocellular_carcinoma) # 10 3
dim(Thyroid_cancer) # 5 3
dim(Prostate_cancer) # 8 3
dim(Colorectal_cancer) # 9 3

related_pathways <- do.call(rbind, lapply(ls(), get))
write.csv(related_pathways, here("Data sources/Cancer-related pathway genes", "KEGG cancer-related pathways.csv"), row.names=F)
rm(list=ls())


# Get list of genes for each cancer-related pathway ####
related_pathways <- read.csv(here("Data sources/Cancer-related pathway genes", "KEGG cancer-related pathways.csv"), stringsAsFactors=F)
for (i in related_pathways$KEGGID) {
  list <- keggGet(i)[[1]]
  gene <- list$GENE
  entrez <- as.numeric(gene[seq(1, length(gene) - 1, 2)])
  symbols_names <- gene[seq(2, length(gene), 2)]
  symbols <- substr(symbols_names, 1, unlist(gregexpr(";", symbols_names)) - 1)
  names <- substr(
    symbols_names, 
    unlist(gregexpr(";", symbols_names)) + 2, 
    unlist(lapply(symbols_names, function(x) unlist(gregexpr("\\[", x))[1])) - 2
  )
  assign(
    gsub(
      " |-", "_", gsub("\\(|\\)| - Homo sapiens \\(human\\)", "", list$NAME)
    ), 
    data.frame(
      EntrezID = entrez, 
      Symbol = symbols, 
      Name = names
    )
  )
  rm(list, gene, entrez, symbols_names, symbols, names)
}
rm(i, related_pathways)

for (i in ls()) {
  write.csv(get(i), here("Data sources/Cancer-related pathway genes", paste0("KEGG_", i, ".csv")), row.names=F)
}
rm(list=ls())


# Get list of genes in related pathways for each cancer ####
related_pathways <- read.csv(here("Data sources/Cancer-related pathway genes", "KEGG cancer-related pathways.csv"), stringsAsFactors=F)
for (i in related_pathways$Cancer) {
  name <- paste0("KEGG_pathway_genes_", i)
  assign(name, data.frame(EntezID=c(), Symbol=c(), Name=c()))
  pathways <- related_pathways$Pathway[which(related_pathways$Cancer == i)]
  for (j in pathways) {
    filename <- paste0("KEGG_", gsub(" |-", "_", gsub("\\(|\\)", "", j)), ".csv")
    list <- read.csv(here("Data sources/Cancer-related pathway genes", filename), stringsAsFactors=F)
    assign(name, rbind(get(name), list))
    rm(filename, list)
  }
  write.csv(get(name), here("Data sources/Cancer-related pathway genes", paste0(name, ".csv")), row.names=F)
  rm(j, name, pathways)
}
rm(list=ls())



