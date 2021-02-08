library(here)
library(dplyr)
library(recount)
library(org.Hs.eg.db)

## Recount genes ####
# Get list of gene names, symbols, IDs from recount, following 
# http://bioconductor.org/packages/release/bioc/vignettes/recount/inst/doc/recount-quickstart.html
gencode <- gsub('\\..*', '', names(recount_genes))
gene_info <- select(org.Hs.eg.db, 
                    gencode, 
                    c('ENTREZID', 'GENENAME', 'SYMBOL','ENSEMBL', 'ALIAS'),
                    # c('ENTREZID', 'GENENAME', 'SYMBOL','ENSEMBL'),
                    'ENSEMBL')
# dim(gene_info) # 115974 5
# dim(gene_info) # 58236 4


## Cancer Gene Census ####
cgc_all <- read.csv2(here("Data sources/Cancer-related genes/Cancer Gene Census", 
                          "Census_allTue Apr 21 04_16_09 2020.csv"), 
                     sep=",", stringsAsFactors=F)
# names(cgc_all)
# dim(cgc_all) # 723 20
# dim(unique(cgc_all)) # 723 20
# length(unique(cgc_all$Gene.Symbol)) # 723

# Try to match any genes whose symbols aren't in Recount
# sum(!(cgc_all$Gene.Symbol %in% gene_info$SYMBOL)) # 11
missing <- cgc_all %>% 
  filter(!(cgc_all$Gene.Symbol %in% gene_info$SYMBOL)) %>% 
  dplyr::select(c(Gene.Symbol, Name, Entrez.GeneId, Synonyms))
# dim(missing) # 11 4
# dim(unique(missing)) # 11 4
found <- gene_info[match(missing$Entrez.GeneId, gene_info$ENTREZID, nomatch=0), ] %>% 
  filter(!is.na(ENTREZID))
# dim(found) # 3 5
cgc_all$Gene.Symbol[match(found$ENTREZID, cgc_all$Entrez.GeneId)] <- found$SYMBOL
# sum(!(cgc_all$Gene.Symbol %in% gene_info$SYMBOL)) # 8
missing <- cgc_all %>% 
  filter(!(cgc_all$Gene.Symbol %in% gene_info$SYMBOL)) %>% 
  dplyr::select(c(Gene.Symbol, Name, Entrez.GeneId, Synonyms))
# dim(missing) # 8 4
found <- gene_info[match(toupper(missing$Name), toupper(gene_info$GENENAME), nomatch=0), ] %>% 
  filter(!is.na(GENENAME))
# dim(found) # 0 5
found <- gene_info[match(missing$Gene.Symbol, gene_info$ALIAS, nomatch=0), ] %>% 
  filter(!is.na(ALIAS))
# dim(found) # 1 5
cgc_all$Gene.Symbol[match(found$ALIAS, cgc_all$Gene.Symbol)] <- found$SYMBOL
# sum(!(cgc_all$Gene.Symbol %in% gene_info$SYMBOL)) # 7
# No names or Entrez IDs match anything in Recount and none have synonyms, so assume they're 
# not in Recount and remove from CGC list.
cgc_all <- cgc_all %>% filter(cgc_all$Gene.Symbol %in% gene_info$SYMBOL)
# dim(cgc_all) # 716 20
cgc_all$Tumour.Types.Combined <- paste0(cgc_all$Tumour.Types.Somatic., ", ", cgc_all$Tumour.Types.Germline.)
# dim(cgc_all) # 716 21

# Breast
cgc_breast <- cgc_all %>% 
  filter(grepl("breast", cgc_all$Tumour.Types.Combined, ignore.case=T))
# dim(cgc_breast) # 43 21
# table(cgc_breast$Tumour.Types.Combined)
cgc_breast <- cgc_breast %>% 
  filter(!(grepl("secretory breast", cgc_breast$Tumour.Types.Combined, ignore.case=T) | 
             grepl("lobular breast", cgc_breast$Tumour.Types.Combined, ignore.case=T) | 
             grepl("luminal A breast", cgc_breast$Tumour.Types.Combined, ignore.case=T) | 
             grepl("phyllodes tumour of the breast", cgc_breast$Tumour.Types.Combined, ignore.case=T)))
# dim(cgc_breast) # 38 21
# table(cgc_breast$Tumour.Types.Combined)
cgc_genes_breast <- dplyr::select(cgc_breast, c(Gene.Symbol, Name, Entrez.GeneId, Synonyms))
# dim(cgc_genes_breast) # 38 4

# Renal
cgc_renal <- cgc_all %>% 
  filter(grepl("clear cell renal", cgc_all$Tumour.Types.Combined, ignore.case=T) | 
           grepl("CCRCC", cgc_all$Tumour.Types.Combined, ignore.case=T))
# dim(cgc_renal) # 9 21
# table(cgc_renal$Tumour.Types.Combined)
cgc_genes_renal <- dplyr::select(cgc_renal, c(Gene.Symbol, Name, Entrez.GeneId, Synonyms))
# dim(cgc_genes_renal) # 9 4

# Thyroid
cgc_thyroid <- cgc_all %>% 
  filter(grepl("papillary thyroid", cgc_all$Tumour.Types.Combined, ignore.case=T))
# dim(cgc_thyroid) # 18 21
# table(cgc_thyroid$Tumour.Types.Combined)
cgc_genes_thyroid <- dplyr::select(cgc_thyroid, c(Gene.Symbol, Name, Entrez.GeneId, Synonyms))
# dim(cgc_genes_thyroid) # 18 4

# Lung adenocarcinoma
cgc_lungad <- cgc_all %>% 
  filter(grepl("lung", cgc_all$Tumour.Types.Combined, ignore.case=T))
# dim(cgc_lungad) # 29 21
# table(cgc_lungad$Tumour.Types.Combined)
cgc_lungad <- cgc_lungad %>% 
  filter(!(grepl("small cell lung", cgc_lungad$Tumour.Types.Combined, ignore.case=T) | 
             grepl("lung SCC", cgc_lungad$Tumour.Types.Combined, ignore.case=T)))
# dim(cgc_lungad) # 25 21
# table(cgc_lungad$Tumour.Types.Combined)
cgc_genes_lungad <- dplyr::select(cgc_lungad, c(Gene.Symbol, Name, Entrez.GeneId, Synonyms))
# dim(cgc_genes_lungad) # 25 4

# Hepatocellular
cgc_hepatocellular <- cgc_all %>% 
  filter(grepl("hepatocellular", cgc_all$Tumour.Types.Combined, ignore.case=T))
# dim(cgc_hepatocellular) # 13 21
# table(cgc_hepatocellular$Tumour.Types.Combined)
cgc_hepatocellular <- cgc_hepatocellular %>% 
  filter(!grepl("fibrolamellar", cgc_hepatocellular$Tumour.Types.Combined))
# dim(cgc_hepatocellular) # 11 21
# table(cgc_hepatocellular$Tumour.Types.Combined)
cgc_genes_hepatocellular <- dplyr::select(cgc_hepatocellular, c(Gene.Symbol, Name, Entrez.GeneId, Synonyms))
# dim(cgc_genes_hepatocellular) # 11 4

# Lung squamous
cgc_lungsq <- cgc_all %>% 
  filter(grepl("lung scc", cgc_all$Tumour.Types.Combined, ignore.case=T))
# dim(cgc_lungsq) # 2 21
# table(cgc_lungsq$Tumour.Types.Combined)
cgc_genes_lungsq <- dplyr::select(cgc_lungsq, c(Gene.Symbol, Name, Entrez.GeneId, Synonyms))
# dim(cgc_genes_lungsq) # 2 4

# Prostate
cgc_prostate <- cgc_all %>% 
  filter(grepl("prostate", cgc_all$Tumour.Types.Combined, ignore.case=T))
# dim(cgc_prostate) # 27 21
# table(cgc_prostate$Tumour.Types.Combined)
cgc_genes_prostate <- dplyr::select(cgc_prostate, c(Gene.Symbol, Name, Entrez.GeneId, Synonyms))
# dim(cgc_genes_prostate) # 27 4

# Colon
cgc_colon <- cgc_all %>% 
  filter(grepl("colon", cgc_all$Tumour.Types.Combined, ignore.case=T) | 
           grepl("colorectal", cgc_all$Tumour.Types.Combined, ignore.case=T) | 
           grepl(" crc", cgc_all$Tumour.Types.Combined, ignore.case=T) | 
           grepl("crc,", cgc_all$Tumour.Types.Combined, ignore.case=T)) # space before or comma after crc excludes ccrcc
# dim(cgc_colon) # 77 21
# table(cgc_colon$Tumour.Types.Combined)
cgc_genes_colon <- dplyr::select(cgc_colon, c(Gene.Symbol, Name, Entrez.GeneId, Synonyms))
# dim(cgc_genes_colon) # 77 4


## DisGeNET ####
disgenet_all <- read.delim(here("Data sources/Cancer-related genes/DisGeNET", 
                                "curated_gene_disease_associations (downloaded 2020-04-22).tsv"), 
                           stringsAsFactors=F)
# names(disgenet_all)
# dim(disgenet_all) # 81746 16
# dim(unique(disgenet_all)) # 81746 16
# length(unique(disgenet_all$diseaseId)) # 10370
# length(unique(disgenet_all$geneSymbol)) # 9411

# Try to match any genes whose symbols aren't in Recount
# sum(!(disgenet_all$geneSymbol %in% gene_info$SYMBOL)) # 985
# sum(!(unique(disgenet_all$geneSymbol) %in% gene_info$SYMBOL)) # 228
missing <- disgenet_all %>% 
  filter(!(disgenet_all$geneSymbol %in% gene_info$SYMBOL)) %>% 
  dplyr::select(c(geneId, geneSymbol))
# dim(missing) # 985 2
found <- gene_info[match(missing$geneId, gene_info$ENTREZID, nomatch=0), ] %>% 
  filter(!is.na(ENTREZID))
# dim(found) # 163 5
# dim(unique(found)) # 30 5
disgenet_all$geneSymbol[which(disgenet_all$geneId %in% found$ENTREZID)] <- found$SYMBOL
# sum(!(disgenet_all$geneSymbol %in% gene_info$SYMBOL)) # 822
# sum(!(unique(disgenet_all$geneSymbol) %in% gene_info$SYMBOL)) # 198
missing <- disgenet_all %>% 
  filter(!(disgenet_all$geneSymbol %in% gene_info$SYMBOL)) %>% 
  dplyr::select(c(geneId, geneSymbol))
# dim(missing) # 822 2
found <- gene_info[match(missing$geneSymbol, gene_info$ALIAS, nomatch=0), ] %>% 
  filter(!is.na(ALIAS))
# dim(found) # 64 5
# dim(unique(found)) # 15 5
disgenet_all$geneSymbol[which(disgenet_all$geneSymbol %in% found$ALIAS)] <- found$SYMBOL
# sum(!(disgenet_all$geneSymbol %in% gene_info$SYMBOL)) # 758
# sum(!(unique(disgenet_all$geneSymbol) %in% gene_info$SYMBOL)) # 183

# DisGeNET only has Entrez ID and symbol, so nothing to match genes whose IDs and symbols 
# aren't in Recount by, so remove.
disgenet_all <- disgenet_all %>% filter(disgenet_all$geneSymbol %in% gene_info$SYMBOL)
# dim(disgenet_all) # 80988 16
# length(unique(disgenet_all$geneSymbol)) # 9214

# Breast
# sort(table(droplevels(disgenet_all %>%
#                         dplyr::select(diseaseName) %>%
#                         filter(grepl("breast", disgenet_all$diseaseName, ignore.case=T)))))
disgenet_breast <- filter(disgenet_all, diseaseName %in% c("Breast Carcinoma", 
                                                           "Breast adenocarcinoma"))
# dim(disgenet_breast) # 524 16
# length(unique(disgenet_breast$geneSymbol)) # 481
# length(unique(paste0(disgenet_breast$geneSymbol, disgenet_breast$diseaseName))) # 524
# Overlap is only because of same gene coming under carcinoma and adenocarcinoma.
# Surprisingly little overlap between them.
# length(unique(disgenet_breast$geneSymbol[which(disgenet_breast$diseaseName == "Breast Carcinoma")])) # 481
# length(unique(disgenet_breast$geneSymbol[which(disgenet_breast$diseaseName == "Breast adenocarcinoma")])) # 43
# Only 43 genes for adenocarcinoma and they are all also in carcinoma, so only need carcinoma
disgenet_breast <- filter(disgenet_all, diseaseName == "Breast Carcinoma")
# dim(disgenet_breast) # 481 16
# length(unique(disgenet_breast$geneSymbol)) # 481
disgenet_genes_breast <- unique(dplyr::select(disgenet_breast, c(geneId, geneSymbol)))
# dim(disgenet_genes_breast) # 481 2

# Renal
# sort(table(droplevels(disgenet_all %>%
#                         dplyr::select(diseaseName) %>%
#                         filter(grepl("renal", disgenet_all$diseaseName, ignore.case=T)))))
disgenet_renal <- filter(disgenet_all, diseaseName == "Conventional (Clear Cell) Renal Cell Carcinoma")
# dim(disgenet_renal) # 143 16
disgenet_genes_renal <- unique(dplyr::select(disgenet_renal, c(geneId, geneSymbol)))
# dim(disgenet_genes_renal) # 143 2

# Thyroid
# sort(table(droplevels(disgenet_all %>%
#                         dplyr::select(diseaseName) %>%
#                         filter(grepl("thyroid", disgenet_all$diseaseName, ignore.case=T)))))
disgenet_thyroid <- filter(disgenet_all, diseaseName == "Papillary thyroid carcinoma")
# dim(disgenet_thyroid) # 83 16
disgenet_genes_thyroid <- unique(dplyr::select(disgenet_thyroid, c(geneId, geneSymbol)))
# dim(disgenet_genes_thyroid) # 83 2

# Lung adenocarcinoma
# sort(table(droplevels(disgenet_all %>%
#                         dplyr::select(diseaseName) %>%
#                         filter(grepl("lung", disgenet_all$diseaseName, ignore.case=T)))))
disgenet_lungad <- filter(disgenet_all, diseaseName == "Adenocarcinoma of lung (disorder)")
# dim(disgenet_lungad) # 210 16
disgenet_genes_lungad <- unique(dplyr::select(disgenet_lungad, c(geneId, geneSymbol)))
# dim(disgenet_genes_lungad) # 210 2

# Hepatocellular
# sort(table(droplevels(disgenet_all %>%
#                         dplyr::select(diseaseName) %>%
#                         filter(grepl("hepatocellular", disgenet_all$diseaseName, ignore.case=T)))))
disgenet_hepatocellular <- filter(disgenet_all, diseaseName == "Adult Hepatocellular Carcinoma")
# dim(disgenet_hepatocellular) # 9 16
disgenet_genes_hepatocellular <- unique(dplyr::select(disgenet_hepatocellular, c(geneId, geneSymbol)))
# dim(disgenet_genes_hepatocellular) # 9 2

# Lung squamous
# sort(table(droplevels(disgenet_all %>%
#                         dplyr::select(diseaseName) %>%
#                         filter(grepl("lung", disgenet_all$diseaseName, ignore.case=T) & 
#                                  (grepl("squamous", disgenet_all$diseaseName, ignore.case=T) | 
#                                  grepl("scc", disgenet_all$diseaseName, ignore.case=T))))))
disgenet_lungsq <- filter(disgenet_all, diseaseName == "Squamous cell carcinoma of lung")
# dim(disgenet_lungsq) # 32 16
disgenet_genes_lungsq <- unique(dplyr::select(disgenet_lungsq, c(geneId, geneSymbol)))
# dim(disgenet_genes_lungsq) # 32 2

# Prostate
# sort(table(droplevels(disgenet_all %>%
#                         dplyr::select(diseaseName) %>%
#                         filter(grepl("prostate", disgenet_all$diseaseName, ignore.case=T)))))
disgenet_prostate <- filter(disgenet_all, diseaseName %in% c("Malignant neoplasm of prostate", 
                                                             "Adenocarcinoma of prostate"))
# dim(disgenet_prostate) # 619 16
# length(unique(disgenet_prostate$geneSymbol)) # 602
disgenet_genes_prostate <- unique(dplyr::select(disgenet_prostate, c(geneId, geneSymbol)))
# dim(disgenet_genes_prostate) # 602 2

# Colon
# sort(table(droplevels(disgenet_all %>%
#                         dplyr::select(diseaseName) %>%
#                         filter(grepl("colorectal", disgenet_all$diseaseName, ignore.case=T) |
#                                  grepl("colon", disgenet_all$diseaseName, ignore.case=T) |
#                                  grepl("crc", disgenet_all$diseaseName, ignore.case=T)))))
disgenet_colon <- filter(disgenet_all, diseaseName %in% c("Colorectal Carcinoma", 
                                                          "Malignant tumor of colon", 
                                                          "Colorectal Cancer", 
                                                          "Colon Carcinoma", 
                                                          "Adenocarcinoma of colon"))
# dim(disgenet_colon) # 1093 16
# length(unique(disgenet_colon$geneSymbol)) # 788
disgenet_genes_colon <- unique(dplyr::select(disgenet_colon, c(geneId, geneSymbol)))
# dim(disgenet_genes_colon) # 788 2


## IntOGen ####
intogen_filtered <- read.delim(here("Data sources/Cancer-related genes/IntOGen/2020-02-02_IntOGen-Drivers-20200213", 
                                    "Compendium_Cancer_Genes.tsv"), 
                               stringsAsFactors=F)
# names(intogen_filtered)
# dim(intogen_filtered) # 3333 18
# dim(unique(intogen_filtered)) # 3333 18
# length(unique(intogen_filtered$SYMBOL)) # 568
intogen_filtered$SYMBOL_COHORT <- paste0(intogen_filtered$SYMBOL, "_", intogen_filtered$COHORT)
# table(table(intogen_filtered$SYMBOL_COHORT))
# Entries are uniquely defined by SYMBOL_COHORT, i.e. gene and cohort
# dim(intogen_filtered) # 3333 19

# Try to match any genes whose symbols aren't in Recount
# Can only try to match by alias
# sum(!(intogen_filtered$SYMBOL %in% gene_info$SYMBOL)) # 4
# sum(!(unique(intogen_filtered$SYMBOL) %in% gene_info$SYMBOL)) # 2
missing <- intogen_filtered %>% 
  filter(!(intogen_filtered$SYMBOL %in% gene_info$SYMBOL)) %>% 
  dplyr::select(c(TRANSCRIPT, SYMBOL))
# dim(missing) # 4 2
found <- gene_info[match(missing$SYMBOL, gene_info$ALIAS, nomatch=0), ] %>% 
  filter(!is.na(ALIAS))
# dim(found) # 4 5
# dim(unique(found)) # 2 5
intogen_filtered$SYMBOL[which(intogen_filtered$SYMBOL %in% found$ALIAS)] <- found$SYMBOL
# sum(!(intogen_filtered$geneSymbol %in% gene_info$SYMBOL)) # 0

intogen_raw <- read.delim(here("Data sources/Cancer-related genes/IntOGen/2020-02-02_IntOGen-Drivers-20200213", 
                               "Unfiltered_driver_results_05.tsv"), 
                          stringsAsFactors=F)
# names(intogen_raw)
# dim(intogen_raw) # 4164 21
# dim(unique(intogen_raw)) # 4164 21
# length(unique(intogen_raw$SYMBOL)) # 976
intogen_raw$SYMBOL_COHORT <- paste0(intogen_raw$SYMBOL, "_", intogen_raw$COHORT)
# table(table(intogen_raw$SYMBOL_COHORT))
# Entries are uniquely defined by SYMBOL_COHORT, i.e. gene and cohort
# dim(intogen_raw) # 4164 22

# Try to match any genes whose symbols aren't in Recount
# Can only try to match by alias
# sum(!(intogen_raw$SYMBOL %in% gene_info$SYMBOL)) # 7
# sum(!(unique(intogen_raw$SYMBOL) %in% gene_info$SYMBOL)) # 4
missing <- intogen_raw %>% 
  filter(!(intogen_raw$SYMBOL %in% gene_info$SYMBOL)) %>% 
  dplyr::select(c(SYMBOL))
# dim(missing) # 7 1
found <- gene_info[match(missing$SYMBOL, gene_info$ALIAS, nomatch=0), ] %>% 
  filter(!is.na(ALIAS))
# dim(found) # 7 5
# dim(unique(found)) # 4 5
intogen_raw$SYMBOL[which(intogen_raw$SYMBOL %in% found$ALIAS)] <- found$SYMBOL
# sum(!(intogen_raw$geneSymbol %in% gene_info$SYMBOL)) # 0
# dim(intogen_raw) # 4164 22
intogen_raw_filtered <- filter(intogen_raw, FILTER=="PASS")
# dim(intogen_raw_filtered) # 3252 22
# length(unique(intogen_raw_filtered$SYMBOL)) # 503

# There are some entries in the filtered list that don't have "PASS" under "FILTER" in the raw list. Not sure if 
# this is an error or if there is a reason they don't correspond, but should make sure I don't select any genes 
# that don't have "PASS" under "FILTER" in the raw list. To be safe, I won't keep any genes that have "PASS" 
# under "FILTER" in the raw list but aren't in the filtered list, i.e. only select genes that are in filtered 
# list and have "PASS" under "FILTER" in the raw list. But each gene may appear multiple times, and only needs 
# to pass once to be ok.

# Breast
intogen_raw_breast <- filter(intogen_raw, CANCER_TYPE == "BRCA")
# dim(intogen_raw_breast) # 250 22
intogen_raw_filtered_breast <- filter(intogen_raw_filtered, CANCER_TYPE=="BRCA")
# dim(intogen_raw_filtered_breast) # 222 22
intogen_filtered_breast <- filter(intogen_filtered, CANCER_TYPE == "BRCA")
# dim(intogen_filtered_breast) # 224 19
intogen_filtered_breast <- filter(intogen_filtered_breast, SYMBOL %in% intogen_raw_filtered_breast$SYMBOL)
# dim(intogen_filtered_breast) # 218 19
intogen_genes_breast <- intogen_filtered_breast[match(
  unique(intogen_filtered_breast$SYMBOL), intogen_filtered_breast$SYMBOL), c("SYMBOL", "TRANSCRIPT")]
# dim(intogen_genes_breast) # 95 2

# Renal
intogen_raw_renal <- filter(intogen_raw, CANCER_TYPE == "RCCC")
# dim(intogen_raw_renal) # 51 22
intogen_raw_filtered_renal <- filter(intogen_raw_filtered, CANCER_TYPE=="RCCC")
# dim(intogen_raw_filtered_renal) # 47 22
intogen_filtered_renal <- filter(intogen_filtered, CANCER_TYPE == "RCCC")
# dim(intogen_filtered_renal) # 47 19
intogen_filtered_renal <- filter(intogen_filtered_renal, SYMBOL %in% intogen_raw_filtered_renal$SYMBOL)
# dim(intogen_filtered_renal) # 47 19
intogen_genes_renal <- intogen_filtered_renal[match(
  unique(intogen_filtered_renal$SYMBOL), intogen_filtered_renal$SYMBOL), c("SYMBOL", "TRANSCRIPT")]
# dim(intogen_genes_renal) # 27 2

# Thyroid
intogen_raw_thyroid <- filter(intogen_raw, CANCER_TYPE == "THCA")
# dim(intogen_raw_thyroid) # 92 22
intogen_raw_filtered_thyroid <- filter(intogen_raw_filtered, CANCER_TYPE=="THCA")
# dim(intogen_raw_filtered_thyroid) # 48 22
intogen_filtered_thyroid <- filter(intogen_filtered, CANCER_TYPE == "THCA")
# dim(intogen_filtered_thyroid) # 50 19
intogen_filtered_thyroid <- filter(intogen_filtered_thyroid, SYMBOL %in% intogen_raw_filtered_thyroid$SYMBOL)
# dim(intogen_filtered_thyroid) # 42 19
intogen_genes_thyroid <- intogen_filtered_thyroid[match(
  unique(intogen_filtered_thyroid$SYMBOL), intogen_filtered_thyroid$SYMBOL), c("SYMBOL", "TRANSCRIPT")]
# dim(intogen_genes_thyroid) # 32 2

# Lung adenocarcinoma
intogen_raw_lungad <- filter(intogen_raw, CANCER_TYPE == "LUAD")
# dim(intogen_raw_lungad) # 71 22
intogen_raw_filtered_lungad <- filter(intogen_raw_filtered, CANCER_TYPE=="LUAD")
# dim(intogen_raw_filtered_lungad) # 57 22
intogen_filtered_lungad <- filter(intogen_filtered, CANCER_TYPE == "LUAD")
# dim(intogen_filtered_lungad) # 60 19
intogen_filtered_lungad <- filter(intogen_filtered_lungad, SYMBOL %in% intogen_raw_filtered_lungad$SYMBOL)
# dim(intogen_filtered_lungad) # 56 19
intogen_genes_lungad <- intogen_filtered_lungad[match(
  unique(intogen_filtered_lungad$SYMBOL), intogen_filtered_lungad$SYMBOL), c("SYMBOL", "TRANSCRIPT")]
# dim(intogen_genes_lungad) # 38 2

# Hepatocellular
# Listed as hepatic cancer, but all 9 cohorts look to be hepatocellular carcinoma, although some specific 
# variants for ICGC cohorts, e.g. virus-associated, secondary to alcohol and adiposity. 
intogen_raw_hepatocellular <- filter(intogen_raw, CANCER_TYPE == "HC")
# dim(intogen_raw_hepatocellular) # 189 22
intogen_raw_filtered_hepatocellular <- filter(intogen_raw_filtered, CANCER_TYPE=="HC")
# dim(intogen_raw_filtered_hepatocellular) # 157 22
intogen_filtered_hepatocellular <- filter(intogen_filtered, CANCER_TYPE == "HC")
# dim(intogen_filtered_hepatocellular) # 157 19
intogen_filtered_hepatocellular <- filter(intogen_filtered_hepatocellular, SYMBOL %in% intogen_raw_filtered_hepatocellular$SYMBOL)
# dim(intogen_filtered_hepatocellular) # 155 19
intogen_genes_hepatocellular <- intogen_filtered_hepatocellular[match(
  unique(intogen_filtered_hepatocellular$SYMBOL), intogen_filtered_hepatocellular$SYMBOL), c("SYMBOL", "TRANSCRIPT")]
# dim(intogen_genes_hepatocellular) # 74 2

# Lung squamous
intogen_raw_lungsq <- filter(intogen_raw, CANCER_TYPE == "LUSC")
# dim(intogen_raw_lungsq) # 104 22
intogen_raw_filtered_lungsq <- filter(intogen_raw_filtered, CANCER_TYPE=="LUSC")
# dim(intogen_raw_filtered_lungsq) # 65 22
intogen_filtered_lungsq <- filter(intogen_filtered, CANCER_TYPE == "LUSC")
# dim(intogen_filtered_lungsq) # 75 19
intogen_filtered_lungsq <- filter(intogen_filtered_lungsq, SYMBOL %in% intogen_raw_filtered_lungsq$SYMBOL)
# dim(intogen_filtered_lungsq) # 64 19
intogen_genes_lungsq <- intogen_filtered_lungsq[match(
  unique(intogen_filtered_lungsq$SYMBOL), intogen_filtered_lungsq$SYMBOL), c("SYMBOL", "TRANSCRIPT")]
# dim(intogen_genes_lungsq) # 50 2

# Prostate
intogen_raw_prostate <- filter(intogen_raw, CANCER_TYPE == "PRAD")
# dim(intogen_raw_prostate) # 278 22
intogen_raw_filtered_prostate <- filter(intogen_raw_filtered, CANCER_TYPE=="PRAD")
# dim(intogen_raw_filtered_prostate) # 212 22
intogen_filtered_prostate <- filter(intogen_filtered, CANCER_TYPE == "PRAD")
# dim(intogen_filtered_prostate) # 215 19
intogen_filtered_prostate <- filter(intogen_filtered_prostate, SYMBOL %in% intogen_raw_filtered_prostate$SYMBOL)
# dim(intogen_filtered_prostate) # 208 19
intogen_genes_prostate <- intogen_filtered_prostate[match(
  unique(intogen_filtered_prostate$SYMBOL), intogen_filtered_prostate$SYMBOL), c("SYMBOL", "TRANSCRIPT")]
# dim(intogen_genes_prostate) # 76 2

# Colon
intogen_raw_colon <- filter(intogen_raw, CANCER_TYPE == "COREAD")
# dim(intogen_raw_colon) # 198 22
intogen_raw_filtered_colon <- filter(intogen_raw_filtered, CANCER_TYPE=="COREAD")
# dim(intogen_raw_filtered_colon) # 146 22
intogen_filtered_colon <- filter(intogen_filtered, CANCER_TYPE == "COREAD")
# dim(intogen_filtered_colon) # 147 19
intogen_filtered_colon <- filter(intogen_filtered_colon, SYMBOL %in% intogen_raw_filtered_colon$SYMBOL)
# dim(intogen_filtered_colon) # 143 19
intogen_genes_colon <- intogen_filtered_colon[match(
  unique(intogen_filtered_colon$SYMBOL), intogen_filtered_colon$SYMBOL), c("SYMBOL", "TRANSCRIPT")]
# dim(intogen_genes_colon) # 69 2


## KEGG ####
# Only have symbols so just have to remove any that don't match Recount

# Breast
kegg_pathway_genes_breast <- read.csv(here("Data sources/Cancer-related genes/KEGG", 
                                           "KEGG_Breast_cancer.csv"), 
                                      header=T, stringsAsFactors=F)
# dim(kegg_pathway_genes_breast) # 147 3
kegg_pathway_genes_breast <- kegg_pathway_genes_breast[which(kegg_pathway_genes_breast$Symbol %in% gene_info$SYMBOL), ]
# dim(kegg_pathway_genes_breast) # 147 3

# Renal
kegg_pathway_genes_renal <- read.csv(here("Data sources/Cancer-related genes/KEGG", 
                                          "KEGG_Renal_cell_carcinoma.csv"), 
                                     header=T, stringsAsFactors=F)
# dim(kegg_pathway_genes_renal) # 69 3
kegg_pathway_genes_renal <- kegg_pathway_genes_renal[which(kegg_pathway_genes_renal$Symbol %in% gene_info$SYMBOL), ]
# dim(kegg_pathway_genes_renal) # 69 3

# Lung adenocarcinoma
kegg_pathway_genes_lungad <- read.csv(here("Data sources/Cancer-related genes/KEGG", 
                                             "KEGG_Non_small_cell_lung_cancer.csv"), 
                                        header=T, stringsAsFactors=F)
# dim(kegg_pathway_genes_lungad) # 68 3
kegg_pathway_genes_lungad <- kegg_pathway_genes_lungad[which(kegg_pathway_genes_lungad$Symbol %in% gene_info$SYMBOL), ]
# dim(kegg_pathway_genes_lungad) # 68 3

# Hepatocellular
kegg_pathway_genes_hepatocellular <- read.csv(here("Data sources/Cancer-related genes/KEGG", 
                                                   "KEGG_Hepatocellular_carcinoma.csv"), 
                                              header=T, stringsAsFactors=F)
# dim(kegg_pathway_genes_hepatocellular) # 168 3
kegg_pathway_genes_hepatocellular[
  -which(kegg_pathway_genes_hepatocellular$Symbol %in% gene_info$SYMBOL), ]
gene_info[grep("2952", gene_info$ENTREZID), ]
gene_info[grep("ENSG00000277656", gene_info$ENSEMBL), ]
gene_info[grep("GSTT1", gene_info$ALIAS), ]
gene_info[grep("glutathione", gene_info$GENENAME), ]
# Entrez ID, Ensemble ID, symbol and name don't appear anywhere in recount genes, so remove
kegg_pathway_genes_hepatocellular <- kegg_pathway_genes_hepatocellular[
  which(kegg_pathway_genes_hepatocellular$Symbol %in% gene_info$SYMBOL), ]
# dim(kegg_pathway_genes_hepatocellular) # 167 3

# Lung squamous
kegg_pathway_genes_lungsq <- kegg_pathway_genes_lungad
# dim(kegg_pathway_genes_lungsq) # 68 3

# Thyroid
kegg_pathway_genes_thyroid <- read.csv(here("Data sources/Cancer-related genes/KEGG", 
                                            "KEGG_Thyroid_cancer.csv"), 
                                       header=T, stringsAsFactors=F)
# dim(kegg_pathway_genes_thyroid) # 37 3
kegg_pathway_genes_thyroid <- kegg_pathway_genes_thyroid[which(kegg_pathway_genes_thyroid$Symbol %in% gene_info$SYMBOL), ]
# dim(kegg_pathway_genes_thyroid) # 37 3

# Prostate
kegg_pathway_genes_prostate <- read.csv(here("Data sources/Cancer-related genes/KEGG", 
                                             "KEGG_Prostate_cancer.csv"), 
                                        header=T, stringsAsFactors=F)
# dim(kegg_pathway_genes_prostate) # 97 3
kegg_pathway_genes_prostate <- kegg_pathway_genes_prostate[which(kegg_pathway_genes_prostate$Symbol %in% gene_info$SYMBOL), ]
# dim(kegg_pathway_genes_prostate) # 97 3

# Colon
kegg_pathway_genes_colon <- read.csv(here("Data sources/Cancer-related genes/KEGG", 
                                          "KEGG_Colorectal_cancer.csv"), 
                                     header=T, stringsAsFactors=F)
# dim(kegg_pathway_genes_colon) # 86 3
kegg_pathway_genes_colon <- kegg_pathway_genes_colon[which(kegg_pathway_genes_colon$Symbol %in% gene_info$SYMBOL), ]
# dim(kegg_pathway_genes_colon) # 86 3


## MalaCards ####
# Have gene names in Description so can try to match, assuming description is same as gene name, and 
# some have ensembl ID as symbol. Some that have ensembl ID as symbol don't have a symbol in gene_info, 
# so just keep those as the ensembl ID.

# Breast
malacards_genes_breast_adeno <- read.csv(here("Data sources/Cancer-related genes/MalaCards", 
                                              "MalaCards_Breast_Adenocarcinoma_genes.csv"), 
                                         stringsAsFactors=F)
malacards_genes_breast_can <- read.csv(here("Data sources/Cancer-related genes/MalaCards", 
                                              "MalaCards_Breast_Cancer_genes.csv"), 
                                       stringsAsFactors=F)
# dim(malacards_genes_breast_adeno) # 108 2
# dim(malacards_genes_breast_can) # 990 2
malacards_genes_breast <- unique(rbind(malacards_genes_breast_adeno, 
                                       malacards_genes_breast_can))
# dim(malacards_genes_breast) # 1022 2
# sum(!(malacards_genes_breast$Symbol %in% gene_info$SYMBOL)) # 50
missing <- malacards_genes_breast %>% 
  filter(!(malacards_genes_breast$Symbol %in% gene_info$SYMBOL))
# dim(missing) # 50 2
found <- gene_info[match(toupper(missing$Symbol), toupper(gene_info$GENENAME), nomatch=0), ] %>% 
  filter(!is.na(GENENAME))
# dim(found) # 0 5
found <- gene_info[match(missing$Symbol, gene_info$ALIAS, nomatch=0), ] %>% 
  filter(!is.na(ALIAS))
# dim(found) # 1 5
malacards_genes_breast$Symbol[match(found$ALIAS, malacards_genes_breast$Symbol)] <- found$SYMBOL
# sum(!(malacards_genes_breast$Symbol %in% gene_info$SYMBOL)) # 49
missing <- malacards_genes_breast %>% 
  filter(!(malacards_genes_breast$Symbol %in% gene_info$SYMBOL))
# dim(missing) # 49 2
found <- gene_info[match(missing$Symbol, gene_info$ENSEMBL, nomatch=0), ] %>% 
  filter(!is.na(SYMBOL))
# dim(found) # 1 5
malacards_genes_breast$Symbol[match(found$ENSEMBL, malacards_genes_breast$Symbol)] <- found$SYMBOL
# sum(!(malacards_genes_breast$Symbol %in% gene_info$SYMBOL)) # 48
missing <- malacards_genes_breast %>% 
  filter(!(malacards_genes_breast$Symbol %in% gene_info$SYMBOL))
# dim(missing) # 48 2
malacards_genes_breast <- malacards_genes_breast %>% 
  filter(malacards_genes_breast$Symbol %in% gene_info$SYMBOL | 
           malacards_genes_breast$Symbol %in% gene_info$ENSEMBL)
# dim(malacards_genes_breast) # 977 2

# Renal
malacards_genes_renal <- read.csv(here("Data sources/Cancer-related genes/MalaCards", 
                                       "MalaCards_CCRCC_genes.csv"), 
                                  stringsAsFactors=F)
# dim(malacards_genes_renal) # 113 2
# sum(!(malacards_genes_renal$Symbol %in% gene_info$SYMBOL)) # 9
missing <- malacards_genes_renal %>% 
  filter(!(malacards_genes_renal$Symbol %in% gene_info$SYMBOL))
# dim(missing) # 9 2
found <- gene_info[match(toupper(missing$Symbol), toupper(gene_info$GENENAME), nomatch=0), ] %>% 
  filter(!is.na(GENENAME))
# dim(found) # 0 5
found <- gene_info[match(missing$Symbol, gene_info$ALIAS, nomatch=0), ] %>% 
  filter(!is.na(ALIAS))
# dim(found) # 0 5
found <- gene_info[match(missing$Symbol, gene_info$ENSEMBL, nomatch=0), ] %>% 
  filter(!is.na(SYMBOL))
# dim(found) # 0 5
malacards_genes_renal <- malacards_genes_renal %>% 
  filter(malacards_genes_renal$Symbol %in% gene_info$SYMBOL | 
           malacards_genes_renal$Symbol %in% gene_info$ENSEMBL)
# dim(malacards_genes_renal) # 107 2

# Thyroid
malacards_genes_thyroid_glca <- read.csv(here("Data sources/Cancer-related genes/MalaCards", 
                                              "MalaCards_Thyroid_Gland_Cancer_genes.csv"), 
                                         stringsAsFactors=F)
malacards_genes_thyroid_car <- read.csv(here("Data sources/Cancer-related genes/MalaCards", 
                                             "MalaCards_Thyroid_Carcinoma_genes.csv"), 
                                        stringsAsFactors=F)
# dim(malacards_genes_thyroid_glca) # 166 2
# dim(malacards_genes_thyroid_car) # 31 2
malacards_genes_thyroid <- unique(rbind(malacards_genes_thyroid_glca, 
                                        malacards_genes_thyroid_car), 
                                  stringsAsFactors=F)
# dim(malacards_genes_thyroid) # 177 2
# sum(!(malacards_genes_thyroid$Symbol %in% gene_info$SYMBOL)) # 16
missing <- malacards_genes_thyroid %>% 
  filter(!(malacards_genes_thyroid$Symbol %in% gene_info$SYMBOL))
# dim(missing) # 16 2
found <- gene_info[match(toupper(missing$Symbol), toupper(gene_info$GENENAME), nomatch=0), ] %>% 
  filter(!is.na(GENENAME))
# dim(found) # 0 5
found <- gene_info[match(missing$Symbol, gene_info$ALIAS, nomatch=0), ] %>% 
  filter(!is.na(ALIAS))
# dim(found) # 0 5
found <- gene_info[match(missing$Symbol, gene_info$ENSEMBL, nomatch=0), ] %>% 
  filter(!is.na(SYMBOL))
# dim(found) # 0 5
malacards_genes_thyroid <- malacards_genes_thyroid %>% 
  filter(malacards_genes_thyroid$Symbol %in% gene_info$SYMBOL | 
           malacards_genes_thyroid$Symbol %in% gene_info$ENSEMBL)
# dim(malacards_genes_thyroid) # 161 2

# Lung adenocarcinoma
# No categories matching appropriately; there is "non-squamous non-small cell", but the entry 
# has no genes.

# Hepatocellular
malacards_genes_hepatocellular <- read.csv(here("Data sources/Cancer-related genes/MalaCards", 
                                                "MalaCards_Hepatocellular_Carcinoma_genes.csv"), 
                                           stringsAsFactors=F)
# dim(malacards_genes_hepatocellular) # 745 2
# sum(!(malacards_genes_hepatocellular$Symbol %in% gene_info$SYMBOL)) # 77
missing <- malacards_genes_hepatocellular %>% 
  filter(!(malacards_genes_hepatocellular$Symbol %in% gene_info$SYMBOL))
# dim(missing) # 77 2
found <- gene_info[match(toupper(missing$Symbol), toupper(gene_info$GENENAME), nomatch=0), ] %>% 
  filter(!is.na(GENENAME))
# dim(found) # 0 5
found <- gene_info[match(missing$Symbol, gene_info$ALIAS, nomatch=0), ] %>% 
  filter(!is.na(ALIAS))
# dim(found) # 0 5
found <- gene_info[match(missing$Symbol, gene_info$ENSEMBL, nomatch=0), ] %>% 
  filter(!is.na(SYMBOL))
# dim(found) # 2 5
malacards_genes_hepatocellular$Symbol[match(found$ENSEMBL, malacards_genes_hepatocellular$Symbol)] <- found$SYMBOL
# sum(!(malacards_genes_hepatocellular$Symbol %in% gene_info$SYMBOL)) # 75
malacards_genes_hepatocellular <- malacards_genes_hepatocellular %>% 
  filter(malacards_genes_hepatocellular$Symbol %in% gene_info$SYMBOL | 
           malacards_genes_hepatocellular$Symbol %in% gene_info$ENSEMBL)
# dim(malacards_genes_hepatocellular) # 680 2

# Lung squamous
malacards_genes_lungsq <- read.csv(here("Data sources/Cancer-related genes/MalaCards", 
                                        "MalaCards_Lung_Squamous_Cell_Carcinoma_genes.csv"), 
                                   stringsAsFactors=F)
# dim(malacards_genes_lungsq) # 175 2
# sum(!(malacards_genes_lungsq$Symbol %in% gene_info$SYMBOL)) # 7
missing <- malacards_genes_lungsq %>% 
  filter(!(malacards_genes_lungsq$Symbol %in% gene_info$SYMBOL))
# dim(missing) # 7 2
found <- gene_info[match(toupper(missing$Symbol), toupper(gene_info$GENENAME), nomatch=0), ] %>% 
  filter(!is.na(GENENAME))
# dim(found) # 0 5
found <- gene_info[match(missing$Symbol, gene_info$ALIAS, nomatch=0), ] %>% 
  filter(!is.na(ALIAS))
# dim(found) # 1 5
malacards_genes_lungsq$Symbol[match(found$ALIAS, malacards_genes_lungsq$Symbol)] <- found$SYMBOL
missing <- malacards_genes_lungsq %>% 
  filter(!(malacards_genes_lungsq$Symbol %in% gene_info$SYMBOL))
# dim(missing) # 6 2
found <- gene_info[match(missing$Symbol, gene_info$ENSEMBL, nomatch=0), ] %>% 
  filter(!is.na(SYMBOL))
# dim(found) # 0 5
malacards_genes_lungsq <- malacards_genes_lungsq %>% 
  filter(malacards_genes_lungsq$Symbol %in% gene_info$SYMBOL | 
           malacards_genes_lungsq$Symbol %in% gene_info$ENSEMBL)
# dim(malacards_genes_lungsq) # 169 2

# Prostate
malacards_genes_prostate <- read.csv(here("Data sources/Cancer-related genes/MalaCards", 
                                          "MalaCards_Prostate_Cancer_genes.csv"), 
                                     stringsAsFactors=F)
# dim(malacards_genes_prostate) # 726 2
# sum(!(malacards_genes_prostate$Symbol %in% gene_info$SYMBOL)) # 58
missing <- malacards_genes_prostate %>% 
  filter(!(malacards_genes_prostate$Symbol %in% gene_info$SYMBOL))
# dim(missing) # 58 2
found <- gene_info[match(toupper(missing$Symbol), toupper(gene_info$GENENAME), nomatch=0), ] %>% 
  filter(!is.na(GENENAME))
# dim(found) # 0 5
found <- gene_info[match(missing$Symbol, gene_info$ALIAS, nomatch=0), ] %>% 
  filter(!is.na(ALIAS))
# dim(found) # 3 5
malacards_genes_prostate$Symbol[match(found$ALIAS, malacards_genes_prostate$Symbol)] <- found$SYMBOL
missing <- malacards_genes_prostate %>% 
  filter(!(malacards_genes_prostate$Symbol %in% gene_info$SYMBOL))
# dim(missing) # 55 2
found <- gene_info[match(missing$Symbol, gene_info$ENSEMBL, nomatch=0), ] %>% 
  filter(!is.na(SYMBOL))
# dim(found) # 1 5
malacards_genes_prostate$Symbol[match(found$ENSEMBL, malacards_genes_prostate$Symbol)] <- found$SYMBOL
missing <- malacards_genes_prostate %>% 
  filter(!(malacards_genes_prostate$Symbol %in% gene_info$SYMBOL))
# dim(missing) # 54 2
malacards_genes_prostate <- malacards_genes_prostate %>% 
  filter(malacards_genes_prostate$Symbol %in% gene_info$SYMBOL | 
           malacards_genes_prostate$Symbol %in% gene_info$ENSEMBL)
# dim(malacards_genes_prostate) # 673 2

# Colon
malacards_genes_colon_ca <- read.csv(here("Data sources/Cancer-related genes/MalaCards", 
                                          "MalaCards_Colon_Cancer_genes.csv"), 
                                     stringsAsFactors=F)
malacards_genes_colon_ad <- read.csv(here("Data sources/Cancer-related genes/MalaCards", 
                                          "MalaCards_Colon_Adenocarcinoma_genes.csv"), 
                                     stringsAsFactors=F)
malacards_genes_colorectal <- read.csv(here("Data sources/Cancer-related genes/MalaCards", 
                                            "MalaCards_Colorectal_Adenocarcinoma_genes.csv"), 
                                       stringsAsFactors=F)
# dim(malacards_genes_colon_ca) # 838 2
# dim(malacards_genes_colon_ad) # 96 2
# dim(malacards_genes_colorectal) # 49 2
malacards_genes_colon <- unique(rbind(malacards_genes_colon_ca, 
                                       malacards_genes_colon_ad, 
                                       malacards_genes_colorectal))
# dim(malacards_genes_colon) # 866 2
# sum(!(malacards_genes_colon$Symbol %in% gene_info$SYMBOL)) # 82
missing <- malacards_genes_colon %>% 
  filter(!(malacards_genes_colon$Symbol %in% gene_info$SYMBOL))
# dim(missing) # 82 2
found <- gene_info[match(toupper(missing$Symbol), toupper(gene_info$GENENAME), nomatch=0), ] %>% 
  filter(!is.na(GENENAME))
# dim(found) # 0 5
found <- gene_info[match(missing$Symbol, gene_info$ALIAS, nomatch=0), ] %>% 
  filter(!is.na(ALIAS))
# dim(found) # 1 5
malacards_genes_colon$Symbol[match(found$ALIAS, malacards_genes_colon$Symbol)] <- found$SYMBOL
missing <- malacards_genes_colon %>% 
  filter(!(malacards_genes_colon$Symbol %in% gene_info$SYMBOL))
# dim(missing) # 81 2
found <- gene_info[match(missing$Symbol, gene_info$ENSEMBL, nomatch=0), ] %>% 
  filter(!is.na(SYMBOL))
# dim(found) # 0 5
malacards_genes_colon <- malacards_genes_colon %>% 
  filter(malacards_genes_colon$Symbol %in% gene_info$SYMBOL | 
           malacards_genes_colon$Symbol %in% gene_info$ENSEMBL)
# dim(malacards_genes_colon) # 791 2


## Combined ####

# Breast
breast_genes_combined <- c(cgc_genes_breast$Gene.Symbol, 
                           disgenet_genes_breast$geneSymbol, 
                           intogen_genes_breast$SYMBOL, 
                           kegg_pathway_genes_breast$Symbol, 
                           malacards_genes_breast$Symbol)
breast_genes <- list(cgc = cgc_genes_breast, 
                     disgenet = disgenet_genes_breast,
                     intogen = intogen_genes_breast, 
                     kegg = kegg_pathway_genes_breast, 
                     malacards = malacards_genes_breast)
# lapply(breast_genes, nrow)
# CGC 38 Breast cancer, breast carcinoma
# DisGeNET 481 Breast carcinoma, breast adenocarcinoma
# IntOGen 95 Breast adenocarcinoma
# KEGG 147 Breast cancer
# MalaCards 977 Breast cancer, breast adenocarcinoma
# length(breast_genes_combined) # 1738
# length(unique(breast_genes_combined)) # 1322
# table(table(breast_genes_combined)) # 1010, 240, 48, 16, 8; 312 in at least 2, 72 in at least 3


# Renal
renal_genes_combined <- c(cgc_genes_renal$Gene.Symbol, 
                          disgenet_genes_renal$geneSymbol, 
                          intogen_genes_renal$SYMBOL, 
                          kegg_pathway_genes_renal$Symbol, 
                          malacards_genes_renal$Symbol)
renal_genes <- list(cgc = cgc_genes_renal, 
                    disgenet = disgenet_genes_renal, 
                    intogen = intogen_genes_renal, 
                    kegg = kegg_pathway_genes_renal, 
                    malacards = malacards_genes_renal)
# lapply(renal_genes, nrow)
# CGC 9 Clear cell renal carcinoma, CCRCC
# DisGeNET 143 Conventional (clear cell) renal call carcinoma
# IntOGen 27 Renal clear cell carcinoma
# KEGG 69 Renal cell carcinoma
# MalaCards 107 Clear cell renal cell carcinoma
# length(renal_genes_combined) # 355
# length(unique(renal_genes_combined)) # 297
# table(table(renal_genes_combined)) # 260, 22, 9, 6; 37 in at least 2, 15 in at least 3

# Thyroid
thyroid_genes_combined <- c(cgc_genes_thyroid$Gene.Symbol, 
                            disgenet_genes_thyroid$geneSymbol, 
                            intogen_genes_thyroid$SYMBOL, 
                            kegg_pathway_genes_thyroid$Symbol, 
                            malacards_genes_thyroid$Symbol)
thyroid_genes <- list(cgc = cgc_genes_thyroid, 
                      disgenet = disgenet_genes_thyroid, 
                      intogen = intogen_genes_thyroid, 
                      kegg = kegg_pathway_genes_thyroid, 
                      malacards = malacards_genes_thyroid)
# lapply(thyroid_genes, nrow)
# CGC 18 Papillary thyroid
# DisGeNET 83 Papillary thyroid carcinoma
# IntOGen 32 Thyroid adenocarcinoma
# KEGG 37 Thyroid cancer
# MalaCards 161 Thyroid gland cancer, thyroid carcinoma
# length(thyroid_genes_combined) # 331
# length(unique(thyroid_genes_combined)) # 265
# table(table(thyroid_genes_combined)) # 225, 25, 6, 7, 2; 40 in at least 2, 15 in at least 3

# Lung adenocarcinoma
lungad_genes_combined <- c(cgc_genes_lungad$Gene.Symbol, 
                           disgenet_genes_lungad$geneSymbol, 
                           intogen_genes_lungad$SYMBOL, 
                           kegg_pathway_genes_lungad$Symbol)
lungad_genes <- list(cgc = cgc_genes_lungad, 
                     disgenet = disgenet_genes_lungad, 
                     intogen = intogen_genes_lungad, 
                     kegg = kegg_pathway_genes_lungad)
# lapply(lungad_genes, nrow)
# CGC 25 Lung cancer, lung adenocarcinoma
# DisGeNet 210 Adenocarcinoma of lung
# IntOGen 38 Lung adenocarcinoma
# KEGG 68 Non-small cell lung cancer
# MalaCards 0
# length(lungad_genes_combined) # 341
# length(unique(lungad_genes_combined)) # 306
# table(table(lungad_genes_combined)) # 280, 18, 7, 1; 26 in at least 2, 8 in at least 3

# Hepatocellular
hepatocellular_genes_combined <- c(cgc_genes_hepatocellular$Gene.Symbol, 
                                   disgenet_genes_hepatocellular$geneSymbol, 
                                   intogen_genes_hepatocellular$SYMBOL, 
                                   kegg_pathway_genes_hepatocellular$Symbol, 
                                   malacards_genes_hepatocellular$Symbol)
hepatocellular_genes <- list(cgc = cgc_genes_hepatocellular, 
                             disgene = disgenet_genes_hepatocellular, 
                             intogen = intogen_genes_hepatocellular, 
                             kegg = kegg_pathway_genes_hepatocellular, 
                             malacards = malacards_genes_hepatocellular)
# lapply(hepatocellular_genes, nrow)
# CGC 11 Hepatocellular carcinoma
# DisGeNET 9 Adult hepatocellular carcinoma
# IntOGen 74 Hepatic cancer
# KEGG 167 Hepatocellular carcinoma
# MalaCards 680 Hepatocellular carcinoma
# length(hepatocellular_genes_combined) # 941
# length(unique(hepatocellular_genes_combined)) # 781
# table(table(hepatocellular_genes_combined)) # 656, 98, 20, 6, 1; 118 in at least 2, 27 in at least 3

# Lung squamous
lungsq_genes_combined <- c(cgc_genes_lungsq$Gene.Symbol, 
                           disgenet_genes_lungsq$geneSymbol, 
                           intogen_genes_lungsq$SYMBOL, 
                           kegg_pathway_genes_lungsq$Symbol, 
                           malacards_genes_lungsq$Symbol)
lungsq_genes <- list(cgc = cgc_genes_lungsq, 
                     disgenet = disgenet_genes_lungsq, 
                     intogen = intogen_genes_lungsq, 
                     kegg = kegg_pathway_genes_lungsq, 
                     malacards = malacards_genes_lungsq)
# lapply(lungsq_genes, nrow)
# CGC 2 Lung SCC
# DisGeNET 32 Squamous cell carcinoma of lung
# IntOGen 50 Lung squamous cell carcinoma
# KEGG 68 Non-small cell lung cancer
# MalaCards 169 Lung squamous cell carcinoma
# length(lungsq_genes_combined) # 321
# length(unique(lungsq_genes_combined)) # 266
# table(table(lungsq_genes_combined)) # 220, 38, 7, 1; 46 in at least 2, 8 in at least 3

# Prostate
prostate_genes_combined <- c(cgc_genes_prostate$Gene.Symbol, 
                             disgenet_genes_prostate$geneSymbol, 
                             intogen_genes_prostate$SYMBOL, 
                             kegg_pathway_genes_prostate$Symbol, 
                             malacards_genes_prostate$Symbol)
prostate_genes <- list(cgc = cgc_genes_prostate, 
                       disgenet = disgenet_genes_prostate, 
                       intogen = intogen_genes_prostate, 
                       kegg = kegg_pathway_genes_prostate, 
                       malacards = malacards_genes_prostate)
# lapply(prostate_genes, nrow)
# CGC 27 Prostate cancer, prostate carcinoma
# DisGeNET 602 Adenocarcinoma of prostate, malignant neoplasm of prostate
# IntOGen 76 Prostate adenocarcinoma
# KEGG 97 Prostate cancer
# MalaCards 673 Prostate cancer
# length(prostate_genes_combined) # 1475
# length(unique(prostate_genes_combined)) # 1125
# table(table(prostate_genes_combined)) # 858, 206, 43, 14, 4; 267 in at least 2, 61 in at least 3

# Colon
colon_genes_combined <- c(cgc_genes_colon$Gene.Symbol, 
                          disgenet_genes_colon$geneSymbol, 
                          intogen_genes_colon$SYMBOL, 
                          kegg_pathway_genes_colon$Symbol, 
                          malacards_genes_colon$Symbol)
colon_genes <- list(cgc = cgc_genes_colon, 
                    disgenet = disgenet_genes_colon, 
                    intogen = intogen_genes_colon, 
                    kegg = kegg_pathway_genes_colon, 
                    malacards = malacards_genes_colon)
# lapply(colon_genes, nrow)
# CGC 77 CRC, colon cancer/carcinoma, colorectal cancer/carcinoma
# DisGeNET 788 Colorectal cancer/carcinoma, colorectal cancer, adenocarcinoma of colon, malignant tumor of colon
# IntOGen 69 Colorectal adenocarcinoma
# KEGG 86 Colorectal cancer
# MalaCards 791 Colon cancer, colon adenocarcinoma, coloerctal adenocarcinoma
# length(colon_genes_combined) # 1811
# length(unique(colon_genes_combined)) # 1437
# table(table(colon_genes_combined)) # 1157, 222, 33, 14, 11; 280 in at least 2, 58 in at least 3


## Save lists ####
folder <- "Data sources/Cancer-related genes"
saveRDS(breast_genes_combined, here(folder, "all_breast_genes.rds"))
saveRDS(renal_genes_combined, here(folder, "all_hepatocellular_genes.rds"))
saveRDS(thyroid_genes_combined, here(folder, "all_thyroid_genes.rds"))
saveRDS(lungad_genes_combined, here(folder, "all_lungad_genes.rds"))
saveRDS(hepatocellular_genes_combined, here(folder, "all_hepatocellular_genes.rds"))
saveRDS(lungsq_genes_combined, here(folder, "all_lungsq_genes.rds"))
saveRDS(prostate_genes_combined, here(folder, "all_prostate_genes.rds"))
saveRDS(colon_genes_combined, here(folder, "all_colon_genes.rds"))
saveRDS(breast_genes, here(folder, "breast_gene_lists.rds"))
saveRDS(renal_genes, here(folder, "renal_gene_lists.rds"))
saveRDS(thyroid_genes, here(folder, "thyroid_gene_lists.rds"))
saveRDS(lungad_genes, here(folder, "lungad_gene_lists.rds"))
saveRDS(hepatocellular_genes, here(folder, "hepatocellular_gene_lists.rds"))
saveRDS(lungsq_genes, here(folder, "lungsq_gene_lists.rds"))
saveRDS(prostate_genes, here(folder, "prostate_gene_lists.rds"))
saveRDS(colon_genes, here(folder, "colon_gene_lists.rds"))




