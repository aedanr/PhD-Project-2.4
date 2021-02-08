library(topGO)
bp <- annFUN.org("BP", mapping = "org.Hs.eg.db", ID = "symbol")
class(bp)
length(bp) # 11964
head(names(bp))
head(bp$`GO:0000002`)
mf <- annFUN.org("MF", mapping = "org.Hs.eg.db", ID = "symbol")
length(mf) # 3903
head(names(mf))
cc <- annFUN.org("CC", mapping = "org.Hs.eg.db", ID = "symbol")
length(cc) # 1629
head(names(cc))

library(GSEABase)
generic <- getOBOCollection("GSEA data/goslim_generic.obo")
slim_terms <- generic@ids
length(slim_terms) # 145
head(slim_terms)
sum(slim_terms %in% names(bp)) # 56
sum(slim_terms %in% names(mf)) # 36
sum(slim_terms %in% names(cc)) # 27
mean(unlist(lapply(bp[which(names(bp) %in% slim_terms)], length))) # 110.5
mean(unlist(lapply(mf[which(names(mf) %in% slim_terms)], length))) # 173.75
mean(unlist(lapply(cc[which(names(cc) %in% slim_terms)], length))) # 1129.037


go_list <- mapIds(org.Hs.eg.db, keys(org.Hs.eg.db, "GO"),
                  "SYMBOL", "GO", multiVals = "list")
length(go_list)
head(names(go_list))
head(go_list$`GO:0002576`)
sum(slim_terms %in% names(go_list)) # 119
slim_terms_bp <- slim_terms[which(slim_terms %in% names(bp))]
slim_terms_mf <- slim_terms[slim_terms %in% names(mf)]
slim_terms_cc <- slim_terms[slim_terms %in% names(cc)]
length(slim_terms_bp) # 56
length(slim_terms_mf) # 36
length(slim_terms_cc) # 27
slim_list_bp <- go_list[which(names(go_list) %in% slim_terms_bp)]
slim_list_mf <- go_list[which(names(go_list) %in% slim_terms_mf)]
slim_list_cc <- go_list[which(names(go_list) %in% slim_terms_cc)]

for (i in 1:length(slim_list_bp)) {
  assign(gsub("GO:", "BP", names(slim_list_bp)[i]), GeneSet(unique(unlist(slim_list_bp[i])), setName=names(slim_list_bp)[i]))
}
bp_sets_list <- list()
for (i in 1:length(slim_list_bp)) {
  bp_sets_list[[i]] <- get(gsub("GO:", "BP", names(slim_list_bp)[i]))
}
bp_collection <- GeneSetCollection(bp_sets_list)
toGmt(bp_collection, "GSEA data/bp_org.Hs.eg.db.gmt")

for (i in 1:length(slim_list_mf)) {
  assign(gsub("GO:", "MF", names(slim_list_mf)[i]), GeneSet(unique(unlist(slim_list_mf[i])), setName=names(slim_list_mf)[i]))
}
mf_sets_list <- list()
for (i in 1:length(slim_list_mf)) {
  mf_sets_list[[i]] <- get(gsub("GO:", "MF", names(slim_list_mf)[i]))
}
mf_collection <- GeneSetCollection(mf_sets_list)
toGmt(mf_collection, "GSEA data/mf_org.Hs.eg.db.gmt")

for (i in 1:length(slim_list_cc)) {
  assign(gsub("GO:", "CC", names(slim_list_cc)[i]), GeneSet(unique(unlist(slim_list_cc[i])), setName=names(slim_list_cc)[i]))
}
cc_sets_list <- list()
for (i in 1:length(slim_list_cc)) {
  cc_sets_list[[i]] <- get(gsub("GO:", "CC", names(slim_list_cc)[i]))
}
cc_collection <- GeneSetCollection(cc_sets_list)
toGmt(cc_collection, "GSEA data/cc_org.Hs.eg.db.gmt")











goslim <- read.delim("GSEA data/goslim_generic.obo", stringsAsFactors=F)
goslim_names <- goslim[grep("name:", goslim[, 1]), ]
length(goslim_names) # 155
goslim_names <- goslim_names[1:145]
goslim_ontology <- goslim[grep("namespace:", goslim[, 1]), ]
goslim_ontology <- goslim_ontology[1:145]
table(goslim_ontology)
goslim_names <- tolower(gsub("\\(", "", gsub("\\)", "", gsub(",", "", gsub("-", " ", gsub("name: ", "", goslim_names))))))
goslim_names_bp <- goslim_names[which(goslim_ontology == "namespace: biological_process")]
goslim_names_mf <- goslim_names[which(goslim_ontology == "namespace: molecular_function")]
goslim_names_cc <- goslim_names[which(goslim_ontology == "namespace: cellular_component")]
bp_names <- tolower(gsub("_", " ", gsub("GO_", "", names(bp))))
sum(goslim_names_bp %in% bp_names) # 44; 27 missing
goslim_names_bp[-which(goslim_names_bp %in% bp_names)]
grep("immune system process", bp_names, value=T)


