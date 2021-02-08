library(here)

folder <- "Data sources/Cancer-related genes"
for (i in c("brca", "kirc", "thca", "luad", "lihc", 
            "lusc", "prad", "coad")) {
  assign(paste0("known.genes.", i), 
         readRDS(here(folder, paste0(i, "_genes_info.rds"))))
}
rm(i, folder)

folder <- "Data sources/Cancer-related pathway genes"
for (i in c("brca", "kirc", "thca", "luad", "lihc", 
            "lusc", "prad", "coad")) {
  assign(paste0("pathway.genes.", i), 
         read.csv(here(folder, paste0(i, "_pathway_genes_recount.csv")), stringsAsFactors=F))
}
rm(i, folder)

for (i in c("brca", "kirc", "thca", "luad", "lihc", 
            "lusc", "prad", "coad")) {
  known <- get(paste0("known.genes.", i))
  pathway <- get(paste0("pathway.genes.", i))
  assign(paste0("pathway.minus.known.genes.", i),
         pathway[-which(pathway$ENTREZID %in% known$ENTREZID), ])
  rm(known, pathway)
}
rm(i)

ls()
c(nrow(known.genes.brca), nrow(pathway.genes.brca), nrow(pathway.minus.known.genes.brca)) # 1323 974 579
c(nrow(known.genes.kirc), nrow(pathway.genes.kirc), nrow(pathway.minus.known.genes.kirc)) # 297 624 539
c(nrow(known.genes.thca), nrow(pathway.genes.thca), nrow(pathway.minus.known.genes.thca)) # 265 608 552
c(nrow(known.genes.luad), nrow(pathway.genes.luad), nrow(pathway.minus.known.genes.luad)) # 306 924 833
c(nrow(known.genes.lihc), nrow(pathway.genes.lihc), nrow(pathway.minus.known.genes.lihc)) # 783 1205 970
c(nrow(known.genes.lusc), nrow(pathway.genes.lusc), nrow(pathway.minus.known.genes.lusc)) # 266 924 819
c(nrow(known.genes.prad), nrow(pathway.genes.prad), nrow(pathway.minus.known.genes.prad)) # 1125 1180 865
c(nrow(known.genes.coad), nrow(pathway.genes.coad), nrow(pathway.minus.known.genes.coad)) # 1437 1026 732

for (i in c("brca", "kirc", "thca", "luad", "lihc", 
            "lusc", "prad", "coad")) {
  write.csv(get(paste0("pathway.minus.known.genes.", i)), 
            here("Data sources/Cancer-related pathway genes", 
                 paste0(i, "_pathway_genes_recount_without_known_genes.csv")), 
            row.names=F)
}
rm(i)

