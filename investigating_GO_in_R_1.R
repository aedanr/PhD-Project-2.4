BiocManager::install(c("GSEABase", "topGO", "GOSemSim", "clusterProfiler", "GOSim", "GO.db"))
library(GSEABase)
library(topGO)
library(GOSemSim)
library(clusterProfiler) # namespace ‘vctrs’ 0.2.4 is already loaded, but >= 0.3.0 is required
library(GOSim)
library(GO.db)

## GSEABase
# toGmt() - export to gmt file
# getGmt()
# goSlim() - map from a given source of GO terms to a given set of slim terms
# getOBOCollection() - parse a file or url encoded using GO OBO specification
bp <- getGmt("GSEA data/c5.bp.v7.1.symbols.gmt", geneIdType=SymbolIdentifier(), collectionType=GOCollection())
bp
toGmt(bp, "GSEA data/bp.gmt")
generic <- getOBOCollection("GSEA data/goslim_generic.obo")
generic # 145 GO terms
bp_slim <- goSlim(bp, generic, "BP")
class(bp) # GeneSetCollection
class(generic) # OBOCollection
subsets(generic)
# want to go from OBO collection to geneset collection and save that as gmt
collectionType(generic) <- GOCollection()
# unable to find an inherited method for function ‘collectionType<-’ for signature ‘"OBOCollection", "GOCollection"’
generic@ids # list of GO terms in generic, which is an OBOCollection

# might be able to use goSlim if I can get list of GO IDs from the gmt file, which only has names
?GOTERM
?Term
terms <- Term(GOTERM)
length(terms)
bp
head(names(bp))
head(terms)
gmt_names <- tolower(gsub("_", " ", gsub("GO_", "", names(bp))))
terms <- tolower(gsub("\\.", " ", gsub("\\)", "", gsub("\\(", "", gsub("/", " ", gsub("'", "", gsub(",", "", gsub("-", " ", terms))))))))
terms <- gsub("\\[", "", gsub("\\]", "", terms))
gmt_names <- gsub("v j ", "vj ", gsub("v d j ", "vdj ", gmt_names))
gmt_names <- gsub("2 6 ", "26 ", gmt_names)
gmt_names <- gsub("nad p h", "nadph", gmt_names)
length(gmt_names) # 7530
sum(gmt_names %in% terms) # 7465
gmt_ids <- names(terms)[match(gmt_names, terms)]
head(gmt_ids)
head(gmt_names)
names(bp) <- gmt_ids
bp
head(names(bp))
missing <- which(is.na(gmt_ids))
gmt_ids_trim <- gmt_ids[-missing]
bp_slim <- goSlim(GOCollection(gmt_ids_trim), generic, "BP")
class(bp_slim)
dim(bp_slim)
# don't see how this is of any use. just gives count for each goslim term.


generic <- getOBOCollection("GSEA data/goslim_generic.obo")
terms <- Term(GOTERM)

bp <- getGmt("GSEA data/c5.bp.v7.1.symbols.gmt", geneIdType=SymbolIdentifier(), collectionType=GOCollection())
gmt_names <- tolower(gsub("_", " ", gsub("GO_", "", names(bp))))
terms <- tolower(gsub("\\.", " ", gsub("\\)", "", gsub("\\(", "", gsub("/", " ", gsub("'", "", gsub(",", "", gsub("-", " ", terms))))))))
terms <- gsub("\\[", "", gsub("\\]", "", terms))
gmt_names <- gsub("v j ", "vj ", gsub("v d j ", "vdj ", gmt_names))
gmt_names <- gsub("2 6 ", "26 ", gmt_names)
gmt_names <- gsub("nad p h", "nadph", gmt_names)
length(gmt_names) # 7530
sum(gmt_names %in% terms) # 7465
gmt_ids <- names(terms)[match(gmt_names, terms)]
# head(gmt_ids)
# head(gmt_names)
sum(generic@ids %in% gmt_ids) # 44
bp_cut <- GeneSetCollection(bp[which(gmt_ids %in% generic@ids)])
length(bp_cut)
toGmt(bp_cut, "GSEA data/bp.gmt")
mean(unlist(lapply(geneIds(bp_cut), length))) # 795.6364

cc <- getGmt("GSEA data/c5.cc.v7.1.symbols.gmt", geneIdType=SymbolIdentifier(), collectionType=GOCollection())
gmt_names <- tolower(gsub("_", " ", gsub("GO_", "", names(cc))))
terms <- tolower(gsub("\\.", " ", gsub("\\)", "", gsub("\\(", "", gsub("/", " ", gsub("'", "", gsub(",", "", gsub("-", " ", terms))))))))
terms <- gsub("\\[", "", gsub("\\]", "", terms))
gmt_names <- gsub("v j ", "vj ", gsub("v d j ", "vdj ", gmt_names))
gmt_names <- gsub("2 6 ", "26 ", gmt_names)
gmt_names <- gsub("nad p h", "nadph", gmt_names)
length(gmt_names) # 999
sum(gmt_names %in% terms) # 980
gmt_ids <- names(terms)[match(gmt_names, terms)]
# head(gmt_ids)
# head(gmt_names)
sum(generic@ids %in% gmt_ids) # 15
cc_cut <- GeneSetCollection(cc[which(gmt_ids %in% generic@ids)])
length(cc_cut)
toGmt(cc_cut, "GSEA data/cc.gmt")
mean(unlist(lapply(geneIds(cc_cut), length))) # 856

mf <- getGmt("GSEA data/c5.mf.v7.1.symbols.gmt", geneIdType=SymbolIdentifier(), collectionType=GOCollection())
gmt_names <- tolower(gsub("_", " ", gsub("GO_", "", names(mf))))
terms <- tolower(gsub("\\.", " ", gsub("\\)", "", gsub("\\(", "", gsub("/", " ", gsub("'", "", gsub(",", "", gsub("-", " ", terms))))))))
terms <- gsub("\\[", "", gsub("\\]", "", terms))
gmt_names <- gsub("v j ", "vj ", gsub("v d j ", "vdj ", gmt_names))
gmt_names <- gsub("2 6 ", "26 ", gmt_names)
gmt_names <- gsub("nad p h", "nadph", gmt_names)
gmt_names <- gsub("nadplus", "nad\\+", gmt_names)
length(gmt_names) # 1663
sum(gmt_names %in% terms) # 1554
gmt_ids <- names(terms)[match(gmt_names, terms)]
# head(gmt_ids)
# head(gmt_names)
sum(generic@ids %in% gmt_ids) # 31
mf_cut <- GeneSetCollection(mf[which(gmt_ids %in% generic@ids)])
length(mf_cut)
toGmt(mf_cut, "GSEA data/mf.gmt")
mean(unlist(lapply(geneIds(mf_cut), length))) # 437.1613

goslim <- read.delim("GSEA data/goslim_generic.obo", stringsAsFactors=F)
goslim_names <- goslim[grep("name:", goslim[, 1]), ]
length(goslim_names)
goslim_names <- tolower(gsub("\\(", "", gsub("\\)", "", gsub(",", "", gsub("-", " ", gsub("name: ", "", goslim_names))))))
sum(goslim_names %in% gmt_names) # 31
goslim_names[which(goslim_names %in% gmt_names)]
generic@ids[which(generic@ids %in% gmt_ids)]

