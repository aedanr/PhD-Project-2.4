# Select GO terms in levels 2 to 6 which are also in org.Hs.eg.db and in similarity matrices computed by GOSemSim () (i.e. those 
# for which information content can be calculated).
library(GO.db)
library(org.Hs.eg.db)
# org.Hs.eg.db has two GO ID keys - "GO" and "GOALL". GOALL has more terms so using that.

# Get all GO IDs in org.Hs.eg.db
org.Hs_IDs <- keys(org.Hs.eg.db, keytype="GOALL")
# length(org.Hs_IDs) # 22746

# Confirm that all GO IDs in org.Hs.eg.db are in GO.db (since GO.db will be used to identify levels and ontologies)
# sum(org.Hs_IDs %in% GOID(GOTERM)) # 22746

# Create subsets of terms for each ontology in levels 2 to 6 (defined by highest level term appears in, where level 1 is level 
# below root, i.e. "BP", "MF", "CC" are level 0), using GO.db
getAllChildren <- function(goids, ontology)
{
  ans <- unique(unlist(mget(goids, get(paste0("GO", ontology, "CHILDREN"))), use.names=FALSE))
  ans <- ans[!is.na(ans)]
}

level1BP <- getAllChildren("GO:0008150", "BP") # GO:0008150 is level 0 ID for BP
level2BP <- getAllChildren(level1BP, "BP")
level3BP <- getAllChildren(level2BP, "BP")
level4BP <- getAllChildren(level3BP, "BP")
level5BP <- getAllChildren(level4BP, "BP")
level6BP <- getAllChildren(level5BP, "BP")
GO.db_IDs_BP <- unique(c(level2BP, level3BP, level4BP, level5BP, level6BP))[-which(
  unique(c(level2BP, level3BP, level4BP, level5BP, level6BP)) %in% c("GO:0008150", level1BP)
)]
# length(GO.db_IDs_BP) # 25554
rm(level1BP, level2BP, level3BP, level4BP, level5BP, level6BP)

level1MF <- getAllChildren("GO:0003674", "MF") # GO:0003674 is level 0 ID for MF
level2MF <- getAllChildren(level1MF, "MF")
level3MF <- getAllChildren(level2MF, "MF")
level4MF <- getAllChildren(level3MF, "MF")
level5MF <- getAllChildren(level4MF, "MF")
level6MF <- getAllChildren(level5MF, "MF")
GO.db_IDs_MF <- unique(c(level2MF, level3MF, level4MF, level5MF, level6MF))[-which(
  unique(c(level2MF, level3MF, level4MF, level5MF, level6MF)) %in% c("GO:0003674", level1MF)
)]
# length(GO.db_IDs_MF) # 0
# No values returned because there is no overlap between levels 0-1 and the rest, so which() returns nothing
GO.db_IDs_MF <- unique(c(level2MF, level3MF, level4MF, level5MF, level6MF))
# length(GO.db_IDs_MF) # 10067
rm(level1MF, level2MF, level3MF, level4MF, level5MF, level6MF)

level1CC <- getAllChildren("GO:0005575", "CC") # GO:0005575 is level 0 ID for CC
level2CC <- getAllChildren(level1CC, "CC")
level3CC <- getAllChildren(level2CC, "CC")
level4CC <- getAllChildren(level3CC, "CC")
level5CC <- getAllChildren(level4CC, "CC")
level6CC <- getAllChildren(level5CC, "CC")
GO.db_IDs_CC <- unique(c(level2CC, level3CC, level4CC, level5CC, level6CC))[-which(
  unique(c(level2CC, level3CC, level4CC, level5CC, level6CC)) %in% c("GO:0005575", level1CC)
)]
# length(GO.db_IDs_CC) # 4130
rm(level1CC, level2CC, level3CC, level4CC, level5CC, level6CC)

# Subset level 2 to 6 terms to only include those in org.Hs.eg.db
IDs_to_match_BP <- intersect(GO.db_IDs_BP, org.Hs_IDs)
IDs_to_match_MF <- intersect(GO.db_IDs_MF, org.Hs_IDs)
IDs_to_match_CC <- intersect(GO.db_IDs_CC, org.Hs_IDs)
# length(IDs_to_match_BP) # 14202
# length(IDs_to_match_MF) # 4092
# length(IDs_to_match_CC) # 1960
rm(org.Hs_IDs, GO.db_IDs_BP, GO.db_IDs_MF, GO.db_IDs_CC)

# Get terms for IDs to be used from GO.db using Term(GOTERM), which has IDs as names
terms_to_match_BP <- Term(GOTERM)[which(names(Term(GOTERM)) %in% IDs_to_match_BP)]
terms_to_match_MF <- Term(GOTERM)[which(names(Term(GOTERM)) %in% IDs_to_match_MF)]
terms_to_match_CC <- Term(GOTERM)[which(names(Term(GOTERM)) %in% IDs_to_match_CC)]
rm(IDs_to_match_BP, IDs_to_match_MF, IDs_to_match_CC)

# Strip formatting of terms - make all upper case; remove brackets, commas, points, spaces, primes, hyphens, forward slashes, 
# colons, >, =; change + to PLUS
# For MF, first change "NAD(P)" to "NADBRPBR" so that it doesn't end up the same as "NADP".
terms_to_match_formatted_BP <- toupper(terms_to_match_BP)
terms_to_match_formatted_BP <- gsub("\\(", "", gsub("\\)", "", gsub("\\[", "", gsub("\\]", "", 
                                                                                    terms_to_match_formatted_BP))))
terms_to_match_formatted_BP <- gsub("\\.", "", gsub(",", "", terms_to_match_formatted_BP))
terms_to_match_formatted_BP <- gsub(" ", "", gsub("'", "", terms_to_match_formatted_BP))
terms_to_match_formatted_BP <- gsub("-", "", gsub("/", "", terms_to_match_formatted_BP))
terms_to_match_formatted_BP <- gsub("\\+", "PLUS", terms_to_match_formatted_BP)
terms_to_match_formatted_BP <- gsub(":", "", terms_to_match_formatted_BP)
terms_to_match_formatted_BP <- gsub(">", "", terms_to_match_formatted_BP)
terms_to_match_formatted_BP <- gsub("=", "", terms_to_match_formatted_BP)

terms_to_match_formatted_MF <- toupper(terms_to_match_MF)
terms_to_match_formatted_MF <- gsub("NAD\\(P\\)", "NADBRPBR", terms_to_match_formatted_MF)
terms_to_match_formatted_MF <- gsub("\\(", "", gsub("\\)", "", gsub("\\[", "", gsub("\\]", "", 
                                                                                    terms_to_match_formatted_MF))))
terms_to_match_formatted_MF <- gsub("\\.", "", gsub(",", "", terms_to_match_formatted_MF))
terms_to_match_formatted_MF <- gsub(" ", "", gsub("'", "", terms_to_match_formatted_MF))
terms_to_match_formatted_MF <- gsub("-", "", gsub("/", "", terms_to_match_formatted_MF))
terms_to_match_formatted_MF <- gsub("\\+", "PLUS", terms_to_match_formatted_MF)
terms_to_match_formatted_MF <- gsub(":", "", terms_to_match_formatted_MF)
terms_to_match_formatted_MF <- gsub(">", "", terms_to_match_formatted_MF)
terms_to_match_formatted_MF <- gsub("=", "", terms_to_match_formatted_MF)

terms_to_match_formatted_CC <- toupper(terms_to_match_CC)
terms_to_match_formatted_CC <- gsub("\\(", "", gsub("\\)", "", gsub("\\[", "", gsub("\\]", "", 
                                                                                    terms_to_match_formatted_CC))))
terms_to_match_formatted_CC <- gsub("\\.", "", gsub(",", "", terms_to_match_formatted_CC))
terms_to_match_formatted_CC <- gsub(" ", "", gsub("'", "", terms_to_match_formatted_CC))
terms_to_match_formatted_CC <- gsub("-", "", gsub("/", "", terms_to_match_formatted_CC))
terms_to_match_formatted_CC <- gsub("\\+", "PLUS", terms_to_match_formatted_CC)
terms_to_match_formatted_CC <- gsub(":", "", terms_to_match_formatted_CC)
terms_to_match_formatted_CC <- gsub(">", "", terms_to_match_formatted_CC)
terms_to_match_formatted_CC <- gsub("=", "", terms_to_match_formatted_CC)

# Import MSigDB terms from gmt files using GSEABase function getGmt()
library(GSEABase)
GSEA_data_BP <- getGmt("GSEA data/c5.bp.v7.1.symbols.gmt", geneIdType=SymbolIdentifier(), collectionType=GOCollection())
msigdb_terms_BP <- names(GSEA_data_BP)
GSEA_data_MF <- getGmt("GSEA data/c5.mf.v7.1.symbols.gmt", geneIdType=SymbolIdentifier(), collectionType=GOCollection())
msigdb_terms_MF <- names(GSEA_data_MF)
GSEA_data_CC <- getGmt("GSEA data/c5.cc.v7.1.symbols.gmt", geneIdType=SymbolIdentifier(), collectionType=GOCollection())
msigdb_terms_CC <- names(GSEA_data_CC)

# Strip formatting of terms - remove "GO" from start and remove underscores
# For MF, first change "NAD_P_" to "NADBRPBR" to match formatting of GO.db/org.Hs.eg.db terms.
msigdb_terms_formatted_BP <- gsub("_", "", msigdb_terms_BP)
msigdb_terms_formatted_BP <- gsub("^GO", "", msigdb_terms_formatted_BP)
# length(msigdb_terms_formatted_BP) # 7530
msigdb_terms_formatted_MF <- gsub("NAD_P_", "NADBRPBR", msigdb_terms_MF)
msigdb_terms_formatted_MF <- gsub("_", "", msigdb_terms_formatted_MF)
msigdb_terms_formatted_MF <- gsub("^GO", "", msigdb_terms_formatted_MF)
# length(msigdb_terms_formatted_MF) # 1663
msigdb_terms_formatted_CC <- gsub("_", "", msigdb_terms_CC)
msigdb_terms_formatted_CC <- gsub("^GO", "", msigdb_terms_formatted_CC)
# length(msigdb_terms_formatted_CC) # 999

# Create data frame matching between IDs, GO.db terms, MSigDB term and stripped terms for each ontology
# length(terms_to_match_BP) # 14202, IDs as names
# length(terms_to_match_formatted_BP) # 14202, IDs as names
# length(msigdb_terms_BP) # 7530
# length(msigdb_terms_formatted_BP) # 7530
terms_to_include_BP <- data.frame(
  ID = names(terms_to_match_BP), 
  term = unname(terms_to_match_BP), 
  stripped_term = terms_to_match_formatted_BP
)
terms_to_include_MF <- data.frame(
  ID = names(terms_to_match_MF), 
  term = unname(terms_to_match_MF), 
  stripped_term = terms_to_match_formatted_MF
)
terms_to_include_CC <- data.frame(
  ID = names(terms_to_match_CC), 
  term = unname(terms_to_match_CC), 
  stripped_term = terms_to_match_formatted_CC
)
msigdb_terms_BP <- data.frame(
  term = msigdb_terms_BP, 
  stripped_term = msigdb_terms_formatted_BP
)
msigdb_terms_MF <- data.frame(
  term = msigdb_terms_MF, 
  stripped_term = msigdb_terms_formatted_MF
)
msigdb_terms_CC <- data.frame(
  term = msigdb_terms_CC, 
  stripped_term = msigdb_terms_formatted_CC
)
GOdb_terms_to_include_in_matrix_BP <- terms_to_include_BP[
  which(terms_to_include_BP$stripped_term %in% msigdb_terms_BP$stripped_term), 
]
GOdb_terms_to_include_in_matrix_MF <- terms_to_include_MF[
  which(terms_to_include_MF$stripped_term %in% msigdb_terms_MF$stripped_term), 
]
GOdb_terms_to_include_in_matrix_CC <- terms_to_include_CC[
  which(terms_to_include_CC$stripped_term %in% msigdb_terms_CC$stripped_term), 
]
msigdb_terms_to_include_in_matrix_BP <- msigdb_terms_BP[
  which(msigdb_terms_BP$stripped_term %in% GOdb_terms_to_include_in_matrix_BP$stripped_term), 
]
msigdb_terms_to_include_in_matrix_MF <- msigdb_terms_MF[
  which(msigdb_terms_MF$stripped_term %in% GOdb_terms_to_include_in_matrix_MF$stripped_term), 
]
msigdb_terms_to_include_in_matrix_CC <- msigdb_terms_CC[
  which(msigdb_terms_CC$stripped_term %in% GOdb_terms_to_include_in_matrix_CC$stripped_term), 
]
# c(nrow(GOdb_terms_to_include_in_matrix_BP), nrow(msigdb_terms_to_include_in_matrix_BP)) # 6783 6783
# c(nrow(GOdb_terms_to_include_in_matrix_MF), nrow(msigdb_terms_to_include_in_matrix_MF)) # 1499 1499
# c(nrow(GOdb_terms_to_include_in_matrix_CC), nrow(msigdb_terms_to_include_in_matrix_CC)) # 986 986
rm(terms_to_match_BP, terms_to_match_MF, terms_to_match_CC, 
   terms_to_match_formatted_BP, terms_to_match_formatted_MF, terms_to_match_formatted_CC, 
   msigdb_terms_formatted_BP, msigdb_terms_formatted_MF, msigdb_terms_formatted_CC, 
   terms_to_include_BP, terms_to_include_MF, terms_to_include_CC, 
   msigdb_terms_BP, msigdb_terms_MF, msigdb_terms_CC)

# Create similarity matrices
library(GOSemSim)
library(rrvgo)
BP_semdata <- godata(OrgDb="org.Hs.eg.db", ont="BP")
MF_semdata <- godata(OrgDb="org.Hs.eg.db", ont="MF")
CC_semdata <- godata(OrgDb="org.Hs.eg.db", ont="CC")
BP_matrix <- calculateSimMatrix(x=GOdb_terms_to_include_in_matrix_BP$ID, 
                                orgdb=org.Hs.eg.db, 
                                ont="BP", method="Rel", 
                                semdata=BP_semdata)
MF_matrix <- calculateSimMatrix(x=GOdb_terms_to_include_in_matrix_MF$ID, 
                                orgdb=org.Hs.eg.db, 
                                ont="MF", method="Rel", 
                                semdata=MF_semdata)
CC_matrix <- calculateSimMatrix(x=GOdb_terms_to_include_in_matrix_CC$ID, 
                                orgdb=org.Hs.eg.db, 
                                ont="CC", method="Rel", 
                                semdata=CC_semdata)
# dim(BP_matrix) # 6567 6567; IC couldn't be computed for 216
# dim(MF_matrix) # 1450 1450; IC couldn't be computed for 49
# dim(CC_matrix) # 929 929; IC couldn't be computed for 57
rm(BP_semdata, MF_semdata, CC_semdata)

# Create data frames with final lists of terms to be included in GSEA analysis
GOdb_terms_to_include_in_gmt_BP <- GOdb_terms_to_include_in_matrix_BP[
  which(GOdb_terms_to_include_in_matrix_BP$ID %in% rownames(BP_matrix)), 
]
GOdb_terms_to_include_in_gmt_MF <- GOdb_terms_to_include_in_matrix_MF[
  which(GOdb_terms_to_include_in_matrix_MF$ID %in% rownames(MF_matrix)), 
]
GOdb_terms_to_include_in_gmt_CC <- GOdb_terms_to_include_in_matrix_CC[
  which(GOdb_terms_to_include_in_matrix_CC$ID %in% rownames(CC_matrix)), 
]
msigdb_terms_to_include_in_gmt_BP <- msigdb_terms_to_include_in_matrix_BP[
  which(msigdb_terms_to_include_in_matrix_BP$stripped_term %in% GOdb_terms_to_include_in_gmt_BP$stripped_term), 
]
msigdb_terms_to_include_in_gmt_MF <- msigdb_terms_to_include_in_matrix_MF[
  which(msigdb_terms_to_include_in_matrix_MF$stripped_term %in% GOdb_terms_to_include_in_gmt_MF$stripped_term), 
]
msigdb_terms_to_include_in_gmt_CC <- msigdb_terms_to_include_in_matrix_CC[
  which(msigdb_terms_to_include_in_matrix_CC$stripped_term %in% GOdb_terms_to_include_in_gmt_CC$stripped_term), 
]
# c(nrow(GOdb_terms_to_include_in_gmt_BP), nrow(msigdb_terms_to_include_in_gmt_BP)) # 6567 6567
# c(nrow(GOdb_terms_to_include_in_gmt_MF), nrow(msigdb_terms_to_include_in_gmt_MF)) # 1450 1450
# c(nrow(GOdb_terms_to_include_in_gmt_CC), nrow(msigdb_terms_to_include_in_gmt_CC)) # 929 929
rm(GOdb_terms_to_include_in_matrix_BP, GOdb_terms_to_include_in_matrix_MF, GOdb_terms_to_include_in_matrix_CC, 
   msigdb_terms_to_include_in_matrix_BP, msigdb_terms_to_include_in_matrix_MF, msigdb_terms_to_include_in_matrix_CC)

# Create gmt files and save
GSEA_data_matched_BP <- GSEA_data_BP[which(names(GSEA_data_BP) %in% msigdb_terms_to_include_in_gmt_BP$term)]
toGmt(GSEA_data_matched_BP, "GSEA data/GSEA_data_BP_2to6_in_sim_mat.gmt")
GSEA_data_matched_MF <- GSEA_data_MF[which(names(GSEA_data_MF) %in% msigdb_terms_to_include_in_gmt_MF$term)]
toGmt(GSEA_data_matched_MF, "GSEA data/GSEA_data_MF_2to6_in_sim_mat.gmt")
GSEA_data_matched_CC <- GSEA_data_CC[which(names(GSEA_data_CC) %in% msigdb_terms_to_include_in_gmt_CC$term)]
toGmt(GSEA_data_matched_CC, "GSEA data/GSEA_data_CC_2to6_in_sim_mat.gmt")
rm(GSEA_data_BP, GSEA_data_MF, GSEA_data_CC, 
   GSEA_data_matched_BP, GSEA_data_matched_MF, GSEA_data_matched_CC)

# Import GSEA results
folders <- c(
  "BP_levels2to6_in_sim_mat_lnHMdisp.brca.GseaPreranked.1595201934330", 
  "BP_levels2to6_in_sim_mat_lnHMdisp.coad.GseaPreranked.1595206728483", 
  "BP_levels2to6_in_sim_mat_lnHMdisp.kirc.GseaPreranked.1595202835930", 
  "BP_levels2to6_in_sim_mat_lnHMdisp.lihc.GseaPreranked.1595203284354", 
  "BP_levels2to6_in_sim_mat_lnHMdisp.luad.GseaPreranked.1595207734594", 
  "BP_levels2to6_in_sim_mat_lnHMdisp.lusc.GseaPreranked.1595208769887",
  "BP_levels2to6_in_sim_mat_lnHMdisp.prad.GseaPreranked.1595209253858", 
  "BP_levels2to6_in_sim_mat_lnHMdisp.thca.GseaPreranked.1595209905298", 
  "BP_levels2to6_in_sim_mat_voom.brca.GseaPreranked.1595203995426", 
  "BP_levels2to6_in_sim_mat_voom.coad.GseaPreranked.1595210756317", 
  "BP_levels2to6_in_sim_mat_voom.kirc.GseaPreranked.1595205183010", 
  "BP_levels2to6_in_sim_mat_voom.lihc.GseaPreranked.1595205671251",
  "BP_levels2to6_in_sim_mat_voom.luad.GseaPreranked.1595212166982", 
  "BP_levels2to6_in_sim_mat_voom.lusc.GseaPreranked.1595213302808", 
  "BP_levels2to6_in_sim_mat_voom.prad.GseaPreranked.1595214193989", 
  "BP_levels2to6_in_sim_mat_voom.thca.GseaPreranked.1595214604890", 
  "CC_levels2to6_in_sim_mat_lnHMdisp.brca.GseaPreranked.1595215046508", 
  "CC_levels2to6_in_sim_mat_lnHMdisp.coad.GseaPreranked.1595215738416",
  "CC_levels2to6_in_sim_mat_lnHMdisp.kirc.GseaPreranked.1595215129137", 
  "CC_levels2to6_in_sim_mat_lnHMdisp.lihc.GseaPreranked.1595215241436", 
  "CC_levels2to6_in_sim_mat_lnHMdisp.luad.GseaPreranked.1595215833114", 
  "CC_levels2to6_in_sim_mat_lnHMdisp.lusc.GseaPreranked.1595215935372", 
  "CC_levels2to6_in_sim_mat_lnHMdisp.prad.GseaPreranked.1595216014039", 
  "CC_levels2to6_in_sim_mat_lnHMdisp.thca.GseaPreranked.1595216095416",
  "CC_levels2to6_in_sim_mat_voom.brca.GseaPreranked.1595215330317", 
  "CC_levels2to6_in_sim_mat_voom.coad.GseaPreranked.1595216230413", 
  "CC_levels2to6_in_sim_mat_voom.kirc.GseaPreranked.1595215399806", 
  "CC_levels2to6_in_sim_mat_voom.lihc.GseaPreranked.1595215647362", 
  "CC_levels2to6_in_sim_mat_voom.luad.GseaPreranked.1595216299294", 
  "CC_levels2to6_in_sim_mat_voom.lusc.GseaPreranked.1595216383598",
  "CC_levels2to6_in_sim_mat_voom.prad.GseaPreranked.1595216472362", 
  "CC_levels2to6_in_sim_mat_voom.thca.GseaPreranked.1595216540455", 
  "MF_levels2to6_in_sim_mat_lnHMdisp.brca.GseaPreranked.1595217549397", 
  "MF_levels2to6_in_sim_mat_lnHMdisp.coad.GseaPreranked.1595218097755", 
  "MF_levels2to6_in_sim_mat_lnHMdisp.kirc.GseaPreranked.1595217633976", 
  "MF_levels2to6_in_sim_mat_lnHMdisp.lihc.GseaPreranked.1595217711530",
  "MF_levels2to6_in_sim_mat_lnHMdisp.luad.GseaPreranked.1595218200991", 
  "MF_levels2to6_in_sim_mat_lnHMdisp.lusc.GseaPreranked.1595218309641", 
  "MF_levels2to6_in_sim_mat_lnHMdisp.prad.GseaPreranked.1595218414721", 
  "MF_levels2to6_in_sim_mat_lnHMdisp.thca.GseaPreranked.1595218503196", 
  "MF_levels2to6_in_sim_mat_voom.brca.GseaPreranked.1595217782492", 
  "MF_levels2to6_in_sim_mat_voom.coad.GseaPreranked.1595218613350",
  "MF_levels2to6_in_sim_mat_voom.kirc.GseaPreranked.1595217860651", 
  "MF_levels2to6_in_sim_mat_voom.lihc.GseaPreranked.1595217955363", 
  "MF_levels2to6_in_sim_mat_voom.luad.GseaPreranked.1595218694458", 
  "MF_levels2to6_in_sim_mat_voom.lusc.GseaPreranked.1595218787179", 
  "MF_levels2to6_in_sim_mat_voom.prad.GseaPreranked.1595218980388", 
  "MF_levels2to6_in_sim_mat_voom.thca.GseaPreranked.1595219108692"
)
library(stringr)
files <- character(48)
for (i in 1:48) {
  files[i] <- paste0(
    "gsea_report_for_na_pos_", 
    str_extract_all(folders[i], "\\(?[0-9]+\\)?")[[1]][3], 
    ".xls"
    )
}
analyses <- character(0)
for (i in c("BP", "CC", "MF")) {
  for (j in c("lnHMdisp", "voom")) {
    for (k in c("brca", "coad", "kirc", "lihc", "luad", "lusc", "prad", "thca")) {
      analyses <- c(analyses, paste0(i, "_", j, "_", k))
    }
  }
}
for (i in 1:48) {
  assign(
    analyses[i], 
    read.delim(
      paste0(
        "Results/GSEA results June 2020/jul20/", 
        folders[i], 
        "/", 
        files[i]
      )
    )
  )
  assign(
    analyses[i], 
    get(analyses[i])[c("NAME", "SIZE", "NES", "NOM.p.val", "FDR.q.val")]
  )
}
rm(folders, files, i)

# Add GO IDs to results and sort by increasing FDR
all_msigdb_terms_in_GSEA_analysis <- rbind(
  msigdb_terms_to_include_in_gmt_BP, 
  msigdb_terms_to_include_in_gmt_MF, 
  msigdb_terms_to_include_in_gmt_CC
)
all_GOdb_terms_in_GSEA_analysis <- rbind(
  GOdb_terms_to_include_in_gmt_BP, 
  GOdb_terms_to_include_in_gmt_MF, 
  GOdb_terms_to_include_in_gmt_CC
)
for (i in analyses) {
  temp <- get(i)
  temp$NAME_stripped <- all_msigdb_terms_in_GSEA_analysis$stripped_term[
    match(temp$NAME, all_msigdb_terms_in_GSEA_analysis$term)
  ]
  temp$GO_ID <- all_GOdb_terms_in_GSEA_analysis$ID[
    match(temp$NAME_stripped, all_GOdb_terms_in_GSEA_analysis$stripped_term)
  ]
  temp <- temp[order(temp$FDR.q.val), ]
  assign(i, temp)
  rm(temp)
}

# Create score vectors from nominal p-values from GSEA results
for (i in analyses) {
  assign(
    paste0("scores_", i), 
    setNames(
      -log10(get(i)$NOM.p.val), get(i)$GO_ID
    )
  )
}

# Run rrvgo
for (i in c("BP", "MF", "CC")) {
  mat <- get(paste0(i, "_matrix"))
  for (j in c("lnHMdisp", "voom")) {
    for (k in c("brca", "coad", "kirc", "lihc", "luad", "lusc", "prad", "thca")) {
      assign(
        paste0("rrvgo_", i, "_", j, "_", k), 
        reduceSimMatrix(
          simMatrix=mat[
            rownames(mat) %in% names(get(paste0("scores_", i, "_", j, "_", k))), 
            colnames(mat) %in% names(get(paste0("scores_", i, "_", j, "_", k)))
          ], 
          scores=get(paste0("scores_", i, "_", j, "_", k)), 
          threshold=0.5, 
          orgdb=org.Hs.eg.db
        )
      )
      terms_to_keep <- get(paste0("rrvgo_", i, "_", j, "_", k))$go[
        which(
          get(paste0("rrvgo_", i, "_", j, "_", k))$go == get(paste0("rrvgo_", i, "_", j, "_", k))$parent
        )
      ]
      assign(
        paste0(i, "_", j, "_", k, "_reduced"), 
        get(paste0(i, "_", j, "_", k))[
          which(get(paste0(i, "_", j, "_", k))$GO_ID %in% terms_to_keep), 
        ]
      )
    }
  }
}
for (i in c("BP", "MF", "CC")) {
  for (j in c("lnHMdisp", "voom")) {
    for (k in c("brca", "coad", "kirc", "lihc", "luad", "lusc", "prad", "thca")) {
      rm(list=paste0("scores_", i, "_", j, "_", k))
    }
  }
}
rm(i,j,k,mat,terms_to_keep)

# Get number of significant terms at 0.05 level for each rrvgo analysis
rrvgo_analyses <- character(0)
for (i in c("BP", "MF", "CC")) {
  for (j in c("voom", "lnHMdisp")) {
    for (k in c("brca", "coad", "kirc", "lihc", "luad", "lusc", "prad", "thca")) {
      rrvgo_analyses <- c(rrvgo_analyses, paste0(i, "_", j, "_", k, "_reduced"))
    }
  }
}
for (i in 1:48) {
  assign(
    paste0("sig_0.05_", rrvgo_analyses[i]), 
    sum(get(rrvgo_analyses[i])$FDR.q.val < 0.05)
  )
}
rm(i,j,k)
sig_0.05 <- data.frame(
  rbind(
    c(sig_0.05_BP_voom_brca_reduced, sig_0.05_MF_voom_brca_reduced, sig_0.05_CC_voom_brca_reduced, 
      sig_0.05_BP_lnHMdisp_brca_reduced, sig_0.05_MF_lnHMdisp_brca_reduced, sig_0.05_CC_lnHMdisp_brca_reduced), 
    c(sig_0.05_BP_voom_coad_reduced, sig_0.05_MF_voom_coad_reduced, sig_0.05_CC_voom_coad_reduced, 
      sig_0.05_BP_lnHMdisp_coad_reduced, sig_0.05_MF_lnHMdisp_coad_reduced, sig_0.05_CC_lnHMdisp_coad_reduced), 
    c(sig_0.05_BP_voom_kirc_reduced, sig_0.05_MF_voom_kirc_reduced, sig_0.05_CC_voom_kirc_reduced, 
      sig_0.05_BP_lnHMdisp_kirc_reduced, sig_0.05_MF_lnHMdisp_kirc_reduced, sig_0.05_CC_lnHMdisp_kirc_reduced), 
    c(sig_0.05_BP_voom_lihc_reduced, sig_0.05_MF_voom_lihc_reduced, sig_0.05_CC_voom_lihc_reduced, 
      sig_0.05_BP_lnHMdisp_lihc_reduced, sig_0.05_MF_lnHMdisp_lihc_reduced, sig_0.05_CC_lnHMdisp_lihc_reduced), 
    c(sig_0.05_BP_voom_luad_reduced, sig_0.05_MF_voom_luad_reduced, sig_0.05_CC_voom_luad_reduced, 
      sig_0.05_BP_lnHMdisp_luad_reduced, sig_0.05_MF_lnHMdisp_luad_reduced, sig_0.05_CC_lnHMdisp_luad_reduced), 
    c(sig_0.05_BP_voom_lusc_reduced, sig_0.05_MF_voom_lusc_reduced, sig_0.05_CC_voom_lusc_reduced, 
      sig_0.05_BP_lnHMdisp_lusc_reduced, sig_0.05_MF_lnHMdisp_lusc_reduced, sig_0.05_CC_lnHMdisp_lusc_reduced), 
    c(sig_0.05_BP_voom_prad_reduced, sig_0.05_MF_voom_prad_reduced, sig_0.05_CC_voom_prad_reduced, 
      sig_0.05_BP_lnHMdisp_prad_reduced, sig_0.05_MF_lnHMdisp_prad_reduced, sig_0.05_CC_lnHMdisp_prad_reduced), 
    c(sig_0.05_BP_voom_thca_reduced, sig_0.05_MF_voom_thca_reduced, sig_0.05_CC_voom_thca_reduced, 
      sig_0.05_BP_lnHMdisp_thca_reduced, sig_0.05_MF_lnHMdisp_thca_reduced, sig_0.05_CC_lnHMdisp_thca_reduced)
  )
)
for (i in c("BP", "MF", "CC")) {
  for (j in c("lnHMdisp", "voom")) {
    for (k in c("brca", "coad", "kirc", "lihc", "luad", "lusc", "prad", "thca")) {
      rm(list=paste0("sig_0.05_", i, "_", j, "_", k, "_reduced"))
    }
  }
}
sig_0.05 <- cbind(
  c("brca", "coad", "kirc", "lihc", "luad", "lusc", "prad", "thca"), sig_0.05
)
names(sig_0.05) <- c("cancer", "voom BP", "voom MF", "voom CC", "lnHMdisp BP", "lnHMdisp MF", "lnHMdisp CC")
sig_0.05
#   cancer voom BP voom MF voom CC lnHMdisp BP lnHMdisp MF lnHMdisp CC
# 1   brca     175      35      38          15           2           5
# 2   coad     122      44      31          31           7          23
# 3   kirc     187      53      29          45           7          11
# 4   lihc     158      41      17          70          35          18
# 5   luad     164      45      30          31          22          30
# 6   lusc     169      37      34          70          17          33
# 7   prad     166      58      28           0           0           0
# 8   thca     182      49      41         112          45          28

# Save results
for (i in c("BP", "MF", "CC")) {
  for (j in c("lnHMdisp", "voom")) {
    for (k in c("brca", "coad", "kirc", "lihc", "luad", "lusc", "prad", "thca")) {
      full_res <- get(paste0(i, "_", j, "_", k))
      rrvgo_res <- get(paste0("rrvgo_", i, "_", j, "_", k))
      reduced_res <- get(paste0(i, "_", j, "_", k, "_reduced"))
      write.csv(
        full_res, 
        paste0(
          "Results/GSEA results June 2020/full_results_", i, "_", j, "_", k, ".csv"
          ), 
        row.names=F
      )
      write.csv(
        rrvgo_res, 
        paste0(
          "Results/GSEA results June 2020/rrvgo_redundancy_results_", i, "_", j, "_", k, ".csv"
        ), 
        row.names=F
      )
      write.csv(
        reduced_res, 
        paste0(
          "Results/GSEA results June 2020/redundant_terms_removed_", i, "_", j, "_", k, ".csv"
        ), 
        row.names=F
      )
    }
  }
}

