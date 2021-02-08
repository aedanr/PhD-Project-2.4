## Map GO terms from GSEA gmt files to identifiers, using information from GO.db or directly from GO ####
go <- read.delim("GSEA data/go.obo", stringsAsFactors=F)
go_names <- go[grep("name:", go[, 1]), ]
go_ids <- go[grep("^id:", go[, 1]), ]
length(go_names) # 47368
length(go_ids) # 47368
go_names <- gsub("name: ", "", go_names[1:47358])
go_ids <- gsub("id: ", "", go_ids[1:47358])
go_table <- data.frame(go_ids, go_names)
dim(go_table) # 47358 2
head(go_table)
go_table <- go_table[-grep("^obsolete", go_table$go_names), ]
dim(go_table) # 44411 2
names(go_table) <- c("id", "term")

library(GO.db)
go.db.terms <- Term(GOTERM)
length(go.db.terms) # 45050
head(go.db.terms)

sum(go.db.terms %in% go_table$term) # 43337
sum(names(go.db.terms) %in% go_table$id) # 43883

head(go.db.terms[-which(go.db.terms %in% go_table$term)], 20)
go_table[which(go_table$id %in% names(head(go.db.terms[-which(go.db.terms %in% go_table$term)], 20))), ]
# All of first 20 either obsolete or changes in terminology - no differences in punctuation, etc.

# Try matching both sources to GSEA GO terms
# Make all upper case:
go.db.terms <- toupper(go.db.terms)
go_table$term <- toupper(go_table$term)
# Remove brackets:
go.db.terms <- gsub("\\(", "", gsub("\\)", "", gsub("\\[", "", gsub("\\]", "", go.db.terms))))
go_table$term <- gsub("\\(", "", gsub("\\)", "", gsub("\\[", "", gsub("\\]", "", go_table$term))))
# Remove commas and points:
go.db.terms <- gsub("\\.", "", gsub(",", "", go.db.terms))
go_table$term <- gsub("\\.", "", gsub(",", "", go_table$term))
# Remove spaces and primes:
go.db.terms <- gsub(" ", "", gsub("'", "", go.db.terms))
go_table$term <- gsub(" ", "", gsub("'", "", go_table$term))
# Remove hyphens and forward slashes:
go.db.terms <- gsub("-", "", gsub("/", "", go.db.terms))
go_table$term <- gsub("-", "", gsub("/", "", go_table$term))
# Change + to "PLUS":
go.db.terms <- gsub("\\+", "PLUS", go.db.terms)
go_table$term <- gsub("\\+", "PLUS", go_table$term)
# Remove colons:
go.db.terms <- gsub(":", "", go.db.terms)
go_table$term <- gsub(":", "", go_table$term)
# Remove ">":
go.db.terms <- gsub(">", "", go.db.terms)
go_table$term <- gsub(">", "", go_table$term)
# Remove "=":
go.db.terms <- gsub("=", "", go.db.terms)
go_table$term <- gsub("=", "", go_table$term)

library(GSEABase)
go.gsea <- names(getGmt("GSEA data/c5.all.v7.1.symbols.gmt", geneIdType=SymbolIdentifier(), collectionType=GOCollection()))
# Remove "GO" and underscores from GSEA list
go.gsea <- gsub("_", "", go.gsea)
go.gsea <- gsub("^GO", "", go.gsea)

# Match
length(go.gsea) # 10192
sum(go.gsea %in% go.db.terms) # 10087 (105 missing)
sum(go.gsea %in% go_table$term) # 10106 (86 missing)
sum(go.gsea %in% unique(c(go.db.terms, go_table$term))) # 10189 (3 missing)
# Looks like best way to match terms used by GSEA is to combine the two sources
# (i.e. GO.db and .obo file directly downloaded from GO). Still means I'm using 
# a list that clearly isn't the most up to date, but it should be reasonable since 
# GSEA was updated in Jan 2020, and definitely defensible since a lot of people 
# use GSEA.
# But whether combining the lists is the most useful thing to do depends on how I 
# narrow down my list of GO terms. e.g. if I use GO.db to identify parent/child 
# terms then there's no point in using terms from the .obo file that aren't in 
# GO.db.
go.gsea[-which(go.gsea %in% unique(c(go.db.terms, go_table$term)))]


## Does it make sense to narrow down GO terms by level using GO.db? ####
library(GO.db)
# Example from https://www.biostars.org/p/367404/ for moving down levels
getAllBPChildren <- function(goids)
{
  ans <- unique(unlist(mget(goids, GOBPCHILDREN), use.names=FALSE))
  ans <- ans[!is.na(ans)]
}
level3_terms <- getAllBPChildren("GO:0002376")
level4_terms <- getAllBPChildren(level3_terms)
length(intersect(level3_terms, level4_terms)) # how many terms are in both lists
# Example from ?GOBPCHILDREN
## Objects in this package can be accessed using the select() interface
## from the AnnotationDbi package. See ?select for details.
xx <- as.list(GOBPCHILDREN)
# Remove GO IDs that do not have any children
xx <- xx[!is.na(xx)]
if(length(xx) > 0){
  # Get the parent GO IDs for the first elents of xx
  goids <- xx[[1]]
  # Find out the GO terms for the first parent goid
  GOID(GOTERM[[goids[1]]])
  Term(GOTERM[[goids[1]]])
  Synonym(GOTERM[[goids[1]]])
  Secondary(GOTERM[[goids[1]]])
  Definition(GOTERM[[goids[1]]])
  Ontology(GOTERM[[goids[1]]])
}
# 
# First want to identify terms at top, so those that don't have any parents (or ancestors)
goids <- GOID(GOTERM)
ans <- unique(unlist(mget(goids, GOBPPARENTS), use.names=F))
length(ans)
ans
goids
# 
AllBPParents <- unique(unlist(mget(GOID(GOTERM), GOBPPARENTS), use.names=F))
# Error because terms in MF, CC don't appear in BP objects. Need to select BP objects only. 
# Need to find out how to do that. Think should be able to use select() but don't know how.
# toTable(GOTERM) converts GOTERM to a table, with one line for each unique combination of 
# ID, term, ontology, definition, synonym, secondary. Could use this to get a list of ids 
# for BP. Probably a better way, but this should work.
BPids <- unique(toTable(GOTERM)$go_id[which(toTable(GOTERM)$Ontology == "BP")])
length(BPids) # 29699
MFids <- unique(toTable(GOTERM)$go_id[which(toTable(GOTERM)$Ontology == "MF")])
length(MFids) # 11148
CCids <- unique(toTable(GOTERM)$go_id[which(toTable(GOTERM)$Ontology == "CC")])
length(CCids) # 4202
AllBPParents <- unlist(mget(GOID(BPids), GOBPPARENTS))
length(AllBPParents) # 71746
head(AllBPParents)
sum(is.na(AllBPParents))
# Manually find BP ID and see what its parents are
mget(GOID("GO:0008150"), GOBPPARENTS) # is_a "all" - so base ontologies still have parents
head(as.list(GOTERM))
# Can also use this format:
GOBPPARENTS$"GO:0008150"
GOBPCHILDREN$"GO:0008150"
# Easiest might be to follow format of first example above and just rely on finding base 
# ontology IDs manually.
getAllChildren <- function(goids, ontology)
{
  ans <- unique(unlist(mget(goids, get(paste0("GO", ontology, "CHILDREN"))), use.names=FALSE))
  ans <- ans[!is.na(ans)]
}
# If I want to avoid having very general terms, I should include terms at each level that 
# don't appear at a lower level, i.e. each term should be only in its lowest level.

## Investigate BP levels ####
level1BP <- "GO:0008150"
level2BP <- getAllChildren(level1BP, "BP")
level1_2BP <- unique(c(level1BP, level2BP))
level3BPall <- getAllChildren(level2BP, "BP")
level3BPmin <- level3BPall[-which(level3BPall %in% level1_2BP)]
level1_3BP <- unique(c(level1_2BP, level3BPall))
level4BPall <- getAllChildren(level3BPmin, "BP")
level4BPmin <- level4BPall[-which(level4BPall %in% level1_3BP)]
level1_4BP <- unique(c(level1_3BP, level4BPall))
level5BPall <- getAllChildren(level4BPmin, "BP")
level5BPmin <- level5BPall[-which(level5BPall %in% level1_4BP)]
level1_5BP <- unique(c(level1_4BP, level5BPall))
level6BPall <- getAllChildren(level5BPmin, "BP")
level6BPmin <- level6BPall[-which(level6BPall %in% level1_5BP)]
level1_6BP <- unique(c(level1_5BP, level6BPall))
level7BPall <- getAllChildren(level6BPmin, "BP")
level7BPmin <- level7BPall[-which(level7BPall %in% level1_6BP)]
level1_7BP <- unique(c(level1_6BP, level7BPall))
level8BPall <- getAllChildren(level7BPmin, "BP")
level8BPmin <- level8BPall[-which(level8BPall %in% level1_7BP)]
level1_8BP <- unique(c(level1_7BP, level8BPall))
level9BPall <- getAllChildren(level8BPmin, "BP")
level9BPmin <- level9BPall[-which(level9BPall %in% level1_8BP)]
level1_9BP <- unique(c(level1_8BP, level9BPall))
level10BPall <- getAllChildren(level9BPmin, "BP")
level10BPmin <- level10BPall[-which(level10BPall %in% level1_9BP)]
level1_10BP <- unique(c(level1_9BP, level10BPall))
level1_9BP <- unique(c(level1_8BP, level9BPall))
level11BPall <- getAllChildren(level10BPmin, "BP")
level11BPmin <- level11BPall[-which(level11BPall %in% level1_10BP)]
level1_11BP <- unique(c(level1_10BP, level11BPall))
level12BPall <- getAllChildren(level11BPmin, "BP")
level12BPmin <- level12BPall[-which(level12BPall %in% level1_11BP)]
level1_12BP <- unique(c(level1_11BP, level12BPall))
level13BPall <- getAllChildren(level12BPmin, "BP")
level13BPmin <- level13BPall[-which(level13BPall %in% level1_12BP)]
level1_13BP <- unique(c(level1_12BP, level13BPall))

length(level2BP) # 32
length(level1_2BP) # 33
length(level3BPall) # 558
length(level3BPmin) # 554
length(level1_3BP) # 587
length(level4BPall) # 3759
length(level4BPmin) # 3591
length(level1_4BP) # 4178
length(level5BPall) # 9993
length(level5BPmin) # 8171
length(level1_5BP) # 12349
length(level6BPall) # 14396
length(level6BPmin) # 8618
length(level1_6BP) # 20967
length(level7BPall) # 11979
length(level7BPmin) # 5020
length(level1_7BP) # 25987
length(level8BPall) # 6793
length(level8BPmin) # 2479
length(level1_8BP) # 28466
length(level9BPall) # 3071
length(level9BPmin) # 989
length(level1_9BP) # 29455
length(level10BPall) # 893
length(level10BPmin) # 208
length(level1_10BP) # 29663
length(level11BPall) # 175
length(level11BPmin) # 34
length(level1_11BP) # 29697
length(level12BPall) # 25
length(level12BPmin) # 2
length(level1_12BP) # 29699
length(level13BPall) # 0
length(level13BPmin) # 0
length(level1_13BP) # 29699
length(unique(toTable(GOTERM)$go_id[which(toTable(GOTERM)$Ontology == "BP")])) # 29699

# Look at some examples of terms in each level
for (i in sample(length(level2BP), 10)) {
  print(Term(GOTERM[[level2BP[i]]]))
} # 32
# rhythmic process; regulation of biological process; localization; carbon utilization; 
# multi-organism process; locomotion; cellular process; immune system process; growth; signaling
grep("RNA", unname(Term(GOTERM)[which(GOID(GOTERM) %in% level2BP)]), value=T)
# none

for (i in sample(length(level3BPmin), 10)) {
  print(Term(GOTERM[[level3BPmin[i]]]))
} # 554
# fertilization, exchange of chromosomal proteins; positive regulation of fibrinolysis; 
# septum digestion after cytokinesis; leukocyte proliferation; immunological memory process; 
# multicellular organism growth; positive regulation of protein localization to cell-cell adherens junction; 
# leukocyte activation; embryo implantation; negative regulation of DNA-binding transcription factor activity
grep("RNA", unname(Term(GOTERM)[which(GOID(GOTERM) %in% level3BPmin)]), value=T)
# Some terms that are quite specific, e.g. "regulation of transcription from RNA polymerase II promoter 
# involved in spermatogenesis", which doesn't have any children, so excluding this level would risk missing 
# some terms completely.

for (i in sample(length(level4BPmin), 10)) {
  print(Term(GOTERM[[level4BPmin[i]]]))
} # 3591
# head segmentation; primary spermatocyte growth; epithelial cell proliferation involved in salivary gland morphogenesis; 
# positive regulation of cellular response to hepatocyte growth factor stimulus; 
# negative regulation of c-di-cGMP signaling; rhombomere development; cell wall integrity MAPK cascade; 
# microtubule bundle formation involved in horsetail-astral microtubule organization; 
# regulation of entry into reproductive diapause; response to herbicide
grep("RNA", unname(Term(GOTERM)[which(GOID(GOTERM) %in% level4BPmin)]), value=T)
# Mix of general and specific terms, e.g. "tRNA folding", "negative regulation of mitochondrial unfolded protein response 
# by negative regulation of transcription from RNA polymeras II promoter", so shouldn't exclude, unless can be 
# confident that child terms of these terms will be informative, which surely will be the case.

for (i in sample(length(level5BPmin), 10)) {
  print(Term(GOTERM[[level5BPmin[i]]]))
} # 8171
# negative regulation of extrachromosomal rDNA circle accumulation involved in replicative cell aging; 
# ureter morphogenesis; response to muscle stretch; positive regulation of trypanothione biosynthetic process; 
# negative regulation of histamine uptake; embryonic crystal cell differentiation; 
# positive regulation of VCP-NPL4-UFD1 AAA ATPase complex assembly; positive regulation of glycolytic process; 
# border follicle cell migration; positive regulatiokn of phosphatidylserine exposure on apoptotic cell surface
grep("RNA", unname(Term(GOTERM)[which(GOID(GOTERM) %in% level5BPmin)]), value=T)
# Again a mix of general and specific terms.

for (i in sample(length(level6BPmin), 10)) {
  print(Term(GOTERM[[level6BPmin[i]]]))
} # 8618
# L-aspartate transmembrane transport; lateral line ganglion development; Golgi vesicle budding; 
# regulation of sensory neuron axon guidance; 
# smoothened signaling pathway involved in growh plate cartilage chondrocyte development; 
# positive regulation of sterol import; photosynthetic phosphorylation; RNA (guanin-N7)-methylation; 
# regulation of myeloid leukocyte mediated immunity; metabolism by symbiont of host xylan
grep("RNA", unname(Term(GOTERM)[which(GOID(GOTERM) %in% level6BPmin)]), value=T)
# Again a mix of general and specific terms.

for (i in sample(length(level7BPmin), 10)) {
  print(Term(GOTERM[[level7BPmin[i]]]))
} # 5020
# 2-deoxyribose 1-phosphate catabolic process; helper T cell chemotaxis; 
# positive regulation of sodium ion export across plasma membrane; filtration diaphragm assembly; 
# mini excitatory postsynaptic potential; N-terminal peptidyl-arginine acetulation; regulation of ERK5 cascade; 
# positive regulation of symbiont of host adenylate cyclase activity; modulation by symbiont of host immune reponse; 
# activation of microtubule activation
grep("RNA", unname(Term(GOTERM)[which(GOID(GOTERM) %in% level7BPmin)]), value=T)
# Again a mix of general and specific terms, but a lot look like they would probably be adequately covered by parent terms.

for (i in sample(length(level8BPmin), 10)) {
  print(Term(GOTERM[[level8BPmin[i]]]))
} # 2479
# regulation of dense core cranule exocytosis; positive regulation of striated muscle cell apoptotic process; 
# dCDP phosphorylation; cellular response to diacul bacterial lipopeptide; D-leucine catabolic process; 
# enzyme active site formation via O-sulfo-L-threonine; tRNA wobble base cytosine methylation; 
# neuromast hair cell differentiation involved in neuromast regeneration; CMP phosphorylation; O-glycan processing, core 2
grep("RNA", unname(Term(GOTERM)[which(GOID(GOTERM) %in% level8BPmin)]), value=T) # none
# Most look like they would probably be adequately covered by parent terms, but still some fairly general, e.g. "mRNA 
# cleavage", whose parent terms aren't very informative: "mRNA metabolic process", "RNA phosphodiester bond hydrolysis".

for (i in sample(length(level9BPmin), 10)) {
  print(Term(GOTERM[[level9BPmin[i]]]))
} # 989
# C-terminal peptidyl-arginine amidation; group A colicin transport; propanediol catabolic process; 
# regulation of cAMP-dependent protein kinase activity; GDP-fucose import into Golgi lumen; hypoxanthine oxidation; 
# peptidyl-aspartagine methylation; DNM1L-mediated stimulation of mitophagy in response to mitochondrial depolarization; 
# laminarabiose transport; hexuronate transmembrane transport
grep("RNA", unname(Term(GOTERM)[which(GOID(GOTERM) %in% level9BPmin)]), value=T)
# Generally look like would have informative parent terms.

for (i in sample(length(level10BPmin), 10)) {
  print(Term(GOTERM[[level10BPmin[i]]]))
} # 208
# passive induction of host innate immune response by virus; ferrous iron transmembrane transport; 
# negative regulation of cytosolic calcium ion concentration; store-operated calcium entry; 
# activation of protein kinase A activity; epinephrine-mediated vasodilation; 
# nuclear tRNA 3'-trailer cleavage, endonucleolytic; positive regulation of cytosolic calcium ion concentration; 
# smooth endoplasmic reticulum calcium ion homeostasis; protein monoubiquitination
grep("RNA", unname(Term(GOTERM)[which(GOID(GOTERM) %in% level10BPmin)]), value=T)
# All quite specific.

for (i in sample(length(level11BPmin), 10)) {
  print(Term(GOTERM[[level11BPmin[i]]]))
} # 34
# alpha-tubulin acetylation; T follicular helper cell differentiation; protein localization to microtubule plus-end; 
# iron assimilation by capture and transport; protein K69-linked ufmylation; peptidyl-lysine N6-acetylation; 
# capsanthin catabolic process; D-glucuronate transmembrane transport; calcium ion export across plasma membrane; 
# negative regulatiokn of smooth endoplasmic reticulum calcium ion concentration
grep("RNA", unname(Term(GOTERM)[which(GOID(GOTERM) %in% level11BPmin)]), value=T)
# none

for (i in 1:length(level12BPmin)) {
  print(Term(GOTERM[[level12BPmin[i]]]))
} # 2
# asymmetric protein localization to old mitotic spindle pole body; 
# asymmetric protein localization to new mitotic spindle pole body


## Investigate MF levels ####
level1MF <- "GO:0003674"
level2MF <- getAllChildren(level1MF, "MF")
level1_2MF <- unique(c(level1MF, level2MF))
level3MF <- getAllChildren(level2MF, "MF") # No level 3 terms in level 2
level1_3MF <- unique(c(level1_2MF, level3MF))
level4MFall <- getAllChildren(level3MF, "MF")
level4MFmin <- level4MFall[-which(level4MFall %in% level1_3MF)]
level1_4MF <- unique(c(level1_3MF, level4MFall))
level5MFall <- getAllChildren(level4MFmin, "MF")
level5MFmin <- level5MFall[-which(level5MFall %in% level1_4MF)]
level1_5MF <- unique(c(level1_4MF, level5MFall))
level6MFall <- getAllChildren(level5MFmin, "MF")
level6MFmin <- level6MFall[-which(level6MFall %in% level1_5MF)]
level1_6MF <- unique(c(level1_5MF, level6MFall))
level7MFall <- getAllChildren(level6MFmin, "MF")
level7MFmin <- level7MFall[-which(level7MFall %in% level1_6MF)]
level1_7MF <- unique(c(level1_6MF, level7MFall))
level8MFall <- getAllChildren(level7MFmin, "MF")
level8MFmin <- level8MFall[-which(level8MFall %in% level1_7MF)]
level1_8MF <- unique(c(level1_7MF, level8MFall))
level9MFall <- getAllChildren(level8MFmin, "MF")
level9MFmin <- level9MFall[-which(level9MFall %in% level1_8MF)]
level1_9MF <- unique(c(level1_8MF, level9MFall))
level10MFall <- getAllChildren(level9MFmin, "MF")
level10MFmin <- level10MFall[-which(level10MFall %in% level1_9MF)]
level1_10MF <- unique(c(level1_9MF, level10MFall))
level1_9MF <- unique(c(level1_8MF, level9MFall))
level11MFall <- getAllChildren(level10MFmin, "MF")
level11MFmin <- level11MFall[-which(level11MFall %in% level1_10MF)]
level1_11MF <- unique(c(level1_10MF, level11MFall))
level12MFall <- getAllChildren(level11MFmin, "MF")
level12MFmin <- level12MFall[-which(level12MFall %in% level1_11MF)]
level1_12MF <- unique(c(level1_11MF, level12MFall))
level13MFall <- getAllChildren(level12MFmin, "MF")
level13MFmin <- level13MFall[-which(level13MFall %in% level1_12MF)]
level1_13MF <- unique(c(level1_12MF, level13MFall))

length(level2MF) # 15
length(level1_2MF) # 16
length(level3MF) # 146
length(level1_3MF) # 162
length(level4MFall) # 871
length(level4MFmin) # 869
length(level1_4MF) # 1031
length(level5MFall) # 2233
length(level5MFmin) # 2098
length(level1_5MF) # 3129
length(level6MFall) # 5388
length(level6MFmin) # 5054
length(level1_6MF) # 8183
length(level7MFall) # 2332
length(level7MFmin) # 1926
length(level1_7MF) # 10109
length(level8MFall) # 997
length(level8MFmin) # 727
length(level1_8MF) # 10836
length(level9MFall) # 271
length(level9MFmin) # 201
length(level1_9MF) # 11037
length(level10MFall) # 95
length(level10MFmin) # 79
length(level1_10MF) # 11116
length(level11MFall) # 21
length(level11MFmin) # 13
length(level1_11MF) # 11129
length(level12MFall) # 22
length(level12MFmin) # 19
length(level1_12MF) # 11148
length(level13MFall) # 0
length(level13MFmin) # 0
length(level1_13MF) # 11148
length(unique(toTable(GOTERM)$go_id[which(toTable(GOTERM)$Ontology == "MF")])) # 11148

# Look at some examples of terms in each level
for (i in sample(length(level2MF), 10)) {
  print(Term(GOTERM[[level2MF[i]]]))
} # 15
# translation regulator activity; molecular carrier activity; transporter activity; structural molecule activity; 
# toxin activity; cargo receptor activity; hijacked molecular function; antioxidant activity; 
# molecular function regulator; molecular transducer activity
grep("RNA", unname(Term(GOTERM)[which(GOID(GOTERM) %in% level2MF)]), value=T)
# none

for (i in sample(length(level3MF), 10)) {
  print(Term(GOTERM[[level3MF[i]]]))
} # 146
# metal cluster binding; amide binding; microfibril binding; neurotransmitter transporter activity; 
# transforming growth factor beta receptor, cytoplasmic mediator activity; general transcription initiation factor activity; 
# general transcription initiation factor activity; dinitrosyl-iron complex binding; structural constituent of cytoskeleton; 
# structural constituent of eye lens; apolipoprotein receptor activity
grep("RNA", unname(Term(GOTERM)[which(GOID(GOTERM) %in% level3MF)]), value=T)
# Two out of 3 terms don't have any children, so really can't exclude this level.

for (i in sample(length(level4MFmin), 10)) {
  print(Term(GOTERM[[level4MFmin[i]]]))
} # 869
# sn-glycerol 3-phosphate binding; tyrosine binding; nitric-oxide synthase regulator activity; clathrin adaptor activity; 
# carbon phosphorus lysase activity; TFIIH-class transcription factor complex binding; riboflavin binding; 
# violaxanthin de-epoxidase activity; phosphopantetheine binding; 4-alpha-glucanotransferase activity
grep("RNA", unname(Term(GOTERM)[which(GOID(GOTERM) %in% level4MFmin)]), value=T)
# Mix of general and specific terms.

for (i in sample(length(level5MFmin), 10)) {
  print(Term(GOTERM[[level5MFmin[i]]]))
} # 2098
# methylation-dependent protein binding; talin binding; sepiapterin deaminase activity; trichloro-p-hydroquinone reductive 
# dehalogenase activity; ion channel activity; chalcone isomerase activity; ATP-dependent RNA helicase activity; 
# diacyl lipopeptide binding; oxidoreductase activity, acting on single donors with incorporation of molecular oxygen, 
# incorporation of two atoms of oxygen;dehydroascorbic acid transmembrane transporter activity
grep("RNA", unname(Term(GOTERM)[which(GOID(GOTERM) %in% level5MFmin)]), value=T)
# Some very specific terms, but still some fairly general.

for (i in sample(length(level6MFmin), 10)) {
  print(Term(GOTERM[[level6MFmin[i]]]))
} # 5054
# DNA-dependent protein kinase activity; peptidase activator activity involved in apoptotic process; single-stranded DNA 
# binding; cyclomaltodextrin glucanotransferase activity; uracil:catio symporter activity; 6C-naringenin dibenzoylmethane 
# tautomer glucoosyltransferase activity; ferredoxin-NAD(P) reductase activity; 
# 4alpha-formyl-5alpha-cholesta-7,24-dien-3beta-ol-4alpha-methyl oxidase activity; glucose-6-phosphate 3-dehydrogenase 
# activity
grep("RNA", unname(Term(GOTERM)[which(GOID(GOTERM) %in% level6MFmin)]), value=T)
# Again a mix of general and specific terms.

for (i in sample(length(level7MFmin), 10)) {
  print(Term(GOTERM[[level7MFmin[i]]]))
} # 1926
# miRNA binding; glucuronoarabinoxylan endo-1,4-beta-xylanase activity; diisopropyl-fluorophosphatase activity; 
# karmpferol 3-0-galactosyltransferase activity; eoxin E4 synthase activity; very-long-chain-(S)-2-hydroxy-acid oxidase 
# activity; undecaprenyl-phosphate mannosyltransferase activity; plastid single-subunit type RNA polymerase biding; 
# ecdysteroid-phosphate phosphatase activity
grep("RNA", unname(Term(GOTERM)[which(GOID(GOTERM) %in% level7MFmin)]), value=T)
# Again a mix of general and specific terms.

for (i in sample(length(level8MFmin), 10)) {
  print(Term(GOTERM[[level8MFmin[i]]]))
} # 727
# N-hydroxyarylamine O-acetyltransferase activity; store-operated calcium channel activity; autotransporter activity; 
# palmitoyl-CoA hydrolase activity; N6-hydroxylysine O-acetyltransferase activity; long-chain acyl-CoA hydrolase activity; 
# ADP reductase activity; acetyl CoA:(Z)-3-hexen-1-ol acetyltransferase activity; translation release factor activity; 
# 5-aminolevulinate synthase activity
grep("RNA", unname(Term(GOTERM)[which(GOID(GOTERM) %in% level8MFmin)]), value=T)
# Mostly very specific terms that surely would be adequately covered by parent terms.

for (i in sample(length(level9MFmin), 10)) {
  print(Term(GOTERM[[level9MFmin[i]]]))
} # 201
# inositol 1,4,5-trisphosphate-sensitive calcium-release channel activity; RNA polymerase II proximal promoter 
# sequence-specific DNA binding; insulin-responsive glucose:proton symporter activity; translation release factor 
# activity, codon nonspecific; voltage-gated sodium channel activity involved in cardiac muscle cell action potential; 
# G-protein activated inward rectifier potassium channel activity; dinucleotide insertion or deletion binding; 
# inositol phosphorylceramide phospholipase activity; phosphatidylinositol-3,4-bisphosphate 4-phosphatase activity; 
# inositol-2,4,5-triphosphate 5-phosphatase activity
grep("RNA", unname(Term(GOTERM)[which(GOID(GOTERM) %in% level9MFmin)]), value=T)
# All look like very specific terms.

for (i in sample(length(level10MFmin), 10)) {
  print(Term(GOTERM[[level10MFmin[i]]]))
} # 79
# carbohydrate response element binding"
# single-stranded DNA-dependent CTPase activity; tubulin N-acetyltransferase activity; mannosyl-inositol 
# phosphorylceramide phospholipase activity; tubulin-dependent ATPase activity; 
# ryanodine-sensitive calcium-release channel activity involved in regulation of postsynaptic cytosolic calcium levels; 
# voltage-gated calcium channel activity involved SA node cell action potential; ATP:3'-cytidine-cytidine-tRNA 
# adenylyltransferase activity; GTP-dependent helicase activity; peptide-glutamate-N-acetyltransferase activity
grep("RNA", unname(Term(GOTERM)[which(GOID(GOTERM) %in% level10MFmin)]), value=T)
# Only 2 terms, probably have informative parents.

for (i in sample(length(level11MFmin), 10)) {
  print(Term(GOTERM[[level11MFmin[i]]]))
} # 13
# cohesin ATPase activity; 
# ATP-dependent microtubule motor activity, plus-end-directed; RNA translocase activity; H3 histone acetyltransferase 
# activity; ATP-dependent microtubule motor activity, minus-end-directed; single-stranded DNA-dependent ATPase activity; 
# protein-DNA loading ATPase activity; H2A histone acetyltransferase activity; H2B histone acetyltransferase activity; 
# H4 histone acetyltransferase activity
grep("RNA", unname(Term(GOTERM)[which(GOID(GOTERM) %in% level11MFmin)]), value=T)
# Only 1 term, "RNA translocase activity", which has parent "RNA-dependent ATPase activity", which has parent "ATPase 
# activity, coupled", which has parent "ATPase activity", which has parent "nucleoside-triphosphatase activity", which has 
# parent "pyrophosphatase activity", which has parent "hydrolase activity, acting on acid anhydrides, in 
# phosphorus-containing anhydrides", which has parent "hydrolase activity, acting on acid anhydrides", which has parent 
# "hydrolase activity", which has parent "catalytic activity", which has parent "molecular function". It's really 
# impossible to say which level would be the most useful. It's difficult to justify excluding even level 11 terms based on 
# this.

for (i in 1:length(level12MFmin)) {
  print(Term(GOTERM[[level12MFmin[i]]]))
} # 19
# histone acetyltransferase activity (H3-K56 specific); histone acetyltransferase activity (H3-K14 specific); etc.
grep("RNA", unname(Term(GOTERM)[which(GOID(GOTERM) %in% level12MFmin)]), value=T)
# "RNA translocase activity involved in viral RNA genome packaging" definitely very specific.

## Investigate CC levels ####
level1CC <- "GO:0005575"
level2CC <- getAllChildren(level1CC, "CC")
level1_2CC <- unique(c(level1CC, level2CC))
level3CCall <- getAllChildren(level2CC, "CC")
level3CCmin <- level3CCall[-which(level3CCall %in% level1_2CC)]
level1_3CC <- unique(c(level1_2CC, level3CCall))
level4CCall <- getAllChildren(level3CCmin, "CC")
level4CCmin <- level4CCall[-which(level4CCall %in% level1_3CC)]
level1_4CC <- unique(c(level1_3CC, level4CCall))
level5CCall <- getAllChildren(level4CCmin, "CC")
level5CCmin <- level5CCall[-which(level5CCall %in% level1_4CC)]
level1_5CC <- unique(c(level1_4CC, level5CCall))
level6CCall <- getAllChildren(level5CCmin, "CC")
level6CCmin <- level6CCall[-which(level6CCall %in% level1_5CC)]
level1_6CC <- unique(c(level1_5CC, level6CCall))
level7CCall <- getAllChildren(level6CCmin, "CC")
level7CCmin <- level7CCall[-which(level7CCall %in% level1_6CC)]
level1_7CC <- unique(c(level1_6CC, level7CCall))
level8CCall <- getAllChildren(level7CCmin, "CC")
level8CCmin <- level8CCall[-which(level8CCall %in% level1_7CC)]
level1_8CC <- unique(c(level1_7CC, level8CCall))
level9CCall <- getAllChildren(level8CCmin, "CC")
level9CCmin <- level9CCall[-which(level9CCall %in% level1_8CC)]
level1_9CC <- unique(c(level1_8CC, level9CCall))
level10CCall <- getAllChildren(level9CCmin, "CC")
level10CCmin <- level10CCall[-which(level10CCall %in% level1_9CC)]
level1_10CC <- unique(c(level1_9CC, level10CCall))
level1_9CC <- unique(c(level1_8CC, level9CCall))
level11CCall <- getAllChildren(level10CCmin, "CC")
level11CCmin <- level11CCall[-which(level11CCall %in% level1_10CC)]
level1_11CC <- unique(c(level1_10CC, level11CCall))
level12CCall <- getAllChildren(level11CCmin, "CC")
level12CCmin <- level12CCall[-which(level12CCall %in% level1_11CC)]
level1_12CC <- unique(c(level1_11CC, level12CCall))
level13CCall <- getAllChildren(level12CCmin, "CC")
level13CCmin <- level13CCall[-which(level13CCall %in% level1_12CC)]
level1_13CC <- unique(c(level1_12CC, level13CCall))

length(level2CC) # 21
length(level1_2CC) # 22
length(level3CCall) # 755
length(level3CCmin) # 748
length(level1_3CC) # 770
length(level4CCall) # 1381
length(level4CCmin) # 1118
length(level1_4CC) # 1888
length(level5CCall) # 2136
length(level5CCmin) # 1383
length(level1_5CC) # 3271
length(level6CCall) # 1403
length(level6CCmin) # 690
length(level1_6CC) # 3961
length(level7CCall) # 555
length(level7CCmin) # 195
length(level1_7CC) # 4156
length(level8CCall) # 126
length(level8CCmin) # 41
length(level1_8CC) # 4197
length(level9CCall) # 18
length(level9CCmin) # 4
length(level1_9CC) # 4201
length(level10CCall) # 6
length(level10CCmin) # 1
length(level1_10CC) # 4202
length(level11CCall) # 0
length(level11CCmin) # 0
length(level1_11CC) # 4202
length(unique(toTable(GOTERM)$go_id[which(toTable(GOTERM)$Ontology == "CC")])) # 4202

# Look at some examples of terms in each level
for (i in sample(length(level2CC), 10)) {
  print(Term(GOTERM[[level2CC[i]]]))
} # 21
# supramolecular complex; symplast; membrane part; cell junction; organelle part; protein-containing complex; 
# mitochondrion-associated adherens complex; virion; membrane; synapse
grep("RNA", unname(Term(GOTERM)[which(GOID(GOTERM) %in% level2CC)]), value=T)
# none

for (i in sample(length(level3CCmin), 10)) {
  print(Term(GOTERM[[level3CCmin[i]]]))
} # 748
# enterobactin synthetase complex; mediator complex; intrinsic component of synaptic membrane; sulfite reductase complex 
# (NADPH); BID-BCL-cl complex; retromer, tubulation complex; eukaryotic translation initiation factor 4F complex; 
# ventral disc lateral crest; postsynaptic specialization membrane; cell pole
grep("RNA", unname(Term(GOTERM)[which(GOID(GOTERM) %in% level3CCmin)]), value=T)
# Some very general but some very specific terms, all involving "complex"

for (i in sample(length(level4CCmin), 10)) {
  print(Term(GOTERM[[level4CCmin[i]]]))
} # 1118
# pollen wall; coated vesicle membrane; cell training edge membrane; nucleotide-excision repair complex; mitochondrial 
# part; bacterial-type flagellum hook-filament junction; gas vesicle shell; mannosyl-oligosaccharide 1,2-alpha-mannosidase 
# complex; Mitochondria-associated ER Membrane; spine mat
grep("RNA", unname(Term(GOTERM)[which(GOID(GOTERM) %in% level4CCmin)]), value=T) # none
# Level of detail that wouldn't want to exclude.

for (i in sample(length(level5CCmin), 10)) {
  print(Term(GOTERM[[level5CCmin[i]]]))
} # 1383
# cell cortex part; Cdc48p-Npl4p-Vms1p AAA ATPase complex; Fused-Smurf ubiquitin ligase complex; myosin III complex; 
# transmembrane actin-associated (TAN) line; anchored component of postsynaptic density membrane; lamin filament; 
# fascia adherens; alphaPDGFR-SHP-2-complex; intrinsic component of mycolate outer membrane
grep("RNA", unname(Term(GOTERM)[which(GOID(GOTERM) %in% level5CCmin)]), value=T) # none
# Generally quite specific terms.

for (i in sample(length(level6CCmin), 10)) {
  print(Term(GOTERM[[level6CCmin[i]]]))
} # 690
# septin cytoskeleton; CatSper complex; Rpd3S complex; chloroplast thylakoid lumen; thiazole synthase complex; 
# extrinsic component of stromal side of plastic inner membrane; nucleolar exosome (RNase complex); chromatoid body; 
# protein kinase complex; RNA polymerase I upstream activating factor complex
grep("RNA", unname(Term(GOTERM)[which(GOID(GOTERM) %in% level6CCmin)]), value=T) # none
# Generally very specific terms. Would expect all to have informative parent terms.

for (i in sample(length(level7CCmin), 10)) {
  print(Term(GOTERM[[level7CCmin[i]]]))
} # 195
# terminal cisterna lumen; SCF-Ufo1/Pof10 ubiquitin ligase complex; host cell rough endoplasmic reticulum; Z chromosome; 
# early recombination nodule; hypolemmal cisterna; MLL1 complex; Golgi cis cisterna membrane; iridosome; 
# U4/U6 x US tri-snRNP complex
grep("RNA", unname(Term(GOTERM)[which(GOID(GOTERM) %in% level7CCmin)]), value=T)
# none

for (i in sample(length(level8CCmin), 10)) {
  print(Term(GOTERM[[level8CCmin[i]]]))
} # 41
# glial limiting end-foot; amylin receptor complex 3; condensed chromatin of inactivated sex chromosome; 
# NuA3b histon acetyltransferase complex; chromatin of active sex chromosome; heterochromatin domain; intercalary 
# heterochromatin; Mst2 hisone acetyltransferase complex; bacteroid-containing symbiosome; alpha-heterochromatin
grep("RNA", unname(Term(GOTERM)[which(GOID(GOTERM) %in% level8CCmin)]), value=T)
# none

for (i in sample(length(level9CCmin), 10)) {
  print(Term(GOTERM[[level9CCmin[i]]]))
} # 4
# C-terminal peptidyl-arginine amidation; group A colicin transport; propanediol catabolic process; 
# regulation of cAMP-dependent protein kinase activity; GDP-fucose import into Golgi lumen; hypoxanthine oxidation; 
# peptidyl-aspartagine methylation; DNM1L-mediated stimulation of mitophagy in response to mitochondrial depolarization; 
# laminarabiose transport; hexuronate transmembrane transport
grep("RNA", unname(Term(GOTERM)[which(GOID(GOTERM) %in% level9CCmin)]), value=T)
# none

for (i in 1:length(level10CCmin)) {
  print(Term(GOTERM[[level10CCmin[i]]]))
} # SAS acetyltransferase complex


## Create subsets based on various level thresholds ####
# May be best to keep all terms, but test a couple of thresholds: discarding levels 1-2 and 9 upwards, which should really 
# only exclude terms that are better represented by parents or children but won't reduce the number of terms by much, and 
# more stringently discard levels 1-3 and 7 upwards, which risks excluding informative terms that don't have informative 
# parents or children and/or don't have children.
# Note that what I'm calling level 1 (i.e. BP, MF, CC) seems to be called level 2 by at least some people (which makes 
# sense really, making level 1 the overall root). It will probably be best to use that terminology if I talk about levels 
# in paper or thesis, but not really any need to change it in all the code now.

length(unique(toTable(GOTERM)$go_id[which(toTable(GOTERM)$Ontology == "BP")])) # 29699
length(unique(toTable(GOTERM)$go_id[which(toTable(GOTERM)$Ontology == "MF")])) # 11148
length(unique(toTable(GOTERM)$go_id[which(toTable(GOTERM)$Ontology == "CC")])) # 4202
length(Term(GOTERM)) # 45050 = 29699 + 11148 + 4202 + 1 (1 entry with id and term "all" has Ontology "universal")
length(unique(toTable(GOTERM)$go_id)) # 45050
dim(unique(toTable(GOTERM)[, c("go_id", "Term")])) # 45050 2
identical(GOID(GOTERM), unique(toTable(GOTERM)$go_id)) # FALSE
head(GOID(GOTERM))
head(unique(toTable(GOTERM)$go_id))
identical(sort(unname(GOID(GOTERM))), sort(unique(toTable(GOTERM)$go_id))) # TRUE
identical(sort(unname(Term(GOTERM))), sort(unique(toTable(GOTERM)$Term))) # TRUE
# Different order for each but same information. Doesn't matter which I use as long as it's consistent and allows me 
# to map between terms and IDs.
all_GO.db_terms <- unique(toTable(GOTERM)[c("go_id", "Term", "Ontology")])
dim(all_GO.db_terms) # 45050 3

# Add column with manipulated terms for matching to MSigDB terms
GO.db_stripped <- toupper(all_GO.db_terms$Term)
# Remove brackets:
GO.db_stripped <- gsub("\\(", "", gsub("\\)", "", gsub("\\[", "", gsub("\\]", "", GO.db_stripped))))
# Remove commas and points:
GO.db_stripped <- gsub("\\.", "", gsub(",", "", GO.db_stripped))
# Remove spaces and primes:
GO.db_stripped <- gsub(" ", "", gsub("'", "", GO.db_stripped))
# Remove hyphens and forward slashes:
GO.db_stripped <- gsub("-", "", gsub("/", "", GO.db_stripped))
# Change + to "PLUS":
GO.db_stripped <- gsub("\\+", "PLUS", GO.db_stripped)
# Remove colons:
GO.db_stripped <- gsub(":", "", GO.db_stripped)
# Remove ">":
GO.db_stripped <- gsub(">", "", GO.db_stripped)
# Remove "=":
GO.db_stripped <- gsub("=", "", GO.db_stripped)

all_GO.db_terms <- cbind(all_GO.db_terms, GO.db_stripped)

all_GO.db_BP <- all_GO.db_terms[which(all_GO.db_terms$Ontology == "BP"), ]
all_GO.db_MF <- all_GO.db_terms[which(all_GO.db_terms$Ontology == "MF"), ]
all_GO.db_CC <- all_GO.db_terms[which(all_GO.db_terms$Ontology == "CC"), ]
level3to8_BP <- all_GO.db_BP[which((all_GO.db_BP$go_id %in% level1_8BP) & !(all_GO.db_BP$go_id %in% level1_2BP)), ]
level4to6_BP <- all_GO.db_BP[which((all_GO.db_BP$go_id %in% level1_6BP) & !(all_GO.db_BP$go_id %in% level1_3BP)), ]
level3to8_MF <- all_GO.db_MF[which((all_GO.db_MF$go_id %in% level1_8MF) & !(all_GO.db_MF$go_id %in% level1_2MF)), ]
level4to6_MF <- all_GO.db_MF[which((all_GO.db_MF$go_id %in% level1_6MF) & !(all_GO.db_MF$go_id %in% level1_3MF)), ]
level3to8_CC <- all_GO.db_CC[which((all_GO.db_CC$go_id %in% level1_8CC) & !(all_GO.db_CC$go_id %in% level1_2CC)), ]
level4to6_CC <- all_GO.db_CC[which((all_GO.db_CC$go_id %in% level1_6CC) & !(all_GO.db_CC$go_id %in% level1_3CC)), ]
level3to8_combined <- rbind(level3to8_BP, level3to8_MF, level3to8_CC)
level4to6_combined <- rbind(level4to6_BP, level4to6_MF, level4to6_CC)
rbind(c(nrow(all_GO.db_terms), nrow(all_GO.db_BP), nrow(all_GO.db_MF), nrow(all_GO.db_CC)), 
      c(nrow(level3to8_combined), nrow(level3to8_BP), nrow(level3to8_MF), nrow(level3to8_CC)), 
      c(nrow(level4to6_combined), nrow(level4to6_BP), nrow(level4to6_MF), nrow(level4to6_CC)))
#       [,1]  [,2]  [,3] [,4]
# [1,] 45050 29699 11148 4202
# [2,] 43428 28433 10820 4175
# [3,] 31592 20380  8021 3191
# Excluding levels 1-2 and 9+ removes only 3.6% of terms overall, and excluding levels 1-3 and 7+ removes 30%.

# Want to run a few examples to see which version gives the most informative categories. Run each with and without 
# applying some sort of reduction of redundancy. Methods like GOSemSim identify pairwise semantic similarity (i.e. 
# similarity of meaning rather than just words) between terms, but don't give a way of choosing which to keep. 
# REVIGO takes a list of GO IDs and p-values and chooses terms to keep (I think based on keeping the term with the 
# lowest p-value), and seems to do more than just pairwise comparisons, although it's probably just combining 
# pairwise comparisons based on the given similarity threshold). REVIGO is browser based, but there is a Bioconductor 
# package, rrvgo, based on REVIGO. There is also a function, simplify(), in clusterProfiler, which uses GOSemSim to 
# identify redundancy and remove redundant terms. It looks like it chooses which term to keep based on p-value, so 
# should be similar to REVIGO/rrvgo. For both, you need to choose a similarity threshold, so should also try a few 
# different values to see which gives the most useful terms.
# Will end up with probably 9 different variants to test - 3 sets of levels and 3 thresholds - for each of the 3 
# ontologies, but will use the same set of levels and same threshold for each ontology. Actually will probably look 
# at each ontology separately and all together - or possibly just all together initially at least.

## Map full and reduced GO term lists to MSigDB list and export ####
GO.db_stripped_combined_full <- as.character(all_GO.db_terms$GO.db_stripped)
GO.db_stripped_combined_3to8 <- as.character(level3to8_combined$GO.db_stripped)
GO.db_stripped_combined_4to6 <- as.character(level4to6_combined$GO.db_stripped)
GO.db_stripped_BP_full <- as.character(all_GO.db_BP$GO.db_stripped)
GO.db_stripped_BP_3to8 <- as.character(level3to8_BP$GO.db_stripped)
GO.db_stripped_BP_4to6 <- as.character(level4to6_BP$GO.db_stripped)
GO.db_stripped_MF_full <- as.character(all_GO.db_MF$GO.db_stripped)
GO.db_stripped_MF_3to8 <- as.character(level3to8_MF$GO.db_stripped)
GO.db_stripped_MF_4to6 <- as.character(level4to6_MF$GO.db_stripped)
GO.db_stripped_CC_full <- as.character(all_GO.db_CC$GO.db_stripped)
GO.db_stripped_CC_3to8 <- as.character(level3to8_CC$GO.db_stripped)
GO.db_stripped_CC_4to6 <- as.character(level4to6_CC$GO.db_stripped)

library(GSEABase)
GSEA_data <- getGmt("GSEA data/c5.all.v7.1.symbols.gmt", geneIdType=SymbolIdentifier(), collectionType=GOCollection())
go.gsea <- names(GSEA_data)
# Remove "GO" and underscores from GSEA list
go.gsea_stripped <- gsub("_", "", go.gsea)
go.gsea_stripped <- gsub("^GO", "", go.gsea_stripped)
GO.GSEA_stripped_combined_full <- go.gsea[which(go.gsea_stripped %in% GO.db_stripped_combined_full)]
GO.GSEA_stripped_combined_3to8 <- go.gsea[which(go.gsea_stripped %in% GO.db_stripped_combined_3to8)]
GO.GSEA_stripped_combined_4to6 <- go.gsea[which(go.gsea_stripped %in% GO.db_stripped_combined_4to6)]
GO.GSEA_stripped_BP_full <- go.gsea[which(go.gsea_stripped %in% GO.db_stripped_BP_full)]
GO.GSEA_stripped_BP_3to8 <- go.gsea[which(go.gsea_stripped %in% GO.db_stripped_BP_3to8)]
GO.GSEA_stripped_BP_4to6 <- go.gsea[which(go.gsea_stripped %in% GO.db_stripped_BP_4to6)]
GO.GSEA_stripped_MF_full <- go.gsea[which(go.gsea_stripped %in% GO.db_stripped_MF_full)]
GO.GSEA_stripped_MF_3to8 <- go.gsea[which(go.gsea_stripped %in% GO.db_stripped_MF_3to8)]
GO.GSEA_stripped_MF_4to6 <- go.gsea[which(go.gsea_stripped %in% GO.db_stripped_MF_4to6)]
GO.GSEA_stripped_CC_full <- go.gsea[which(go.gsea_stripped %in% GO.db_stripped_CC_full)]
GO.GSEA_stripped_CC_3to8 <- go.gsea[which(go.gsea_stripped %in% GO.db_stripped_CC_3to8)]
GO.GSEA_stripped_CC_4to6 <- go.gsea[which(go.gsea_stripped %in% GO.db_stripped_CC_4to6)]
rbind(c(length(GO.GSEA_stripped_combined_full), length(GO.GSEA_stripped_BP_full), 
        length(GO.GSEA_stripped_MF_full), length(GO.GSEA_stripped_CC_full)), 
      c(length(GO.GSEA_stripped_combined_3to8), length(GO.GSEA_stripped_BP_3to8), 
        length(GO.GSEA_stripped_MF_3to8), length(GO.GSEA_stripped_CC_3to8)), 
      c(length(GO.GSEA_stripped_combined_4to6), length(GO.GSEA_stripped_BP_4to6), 
        length(GO.GSEA_stripped_MF_4to6), length(GO.GSEA_stripped_CC_4to6)))
#       [,1] [,2] [,3] [,4]
# [1,] 10087 7472 1622  993
# [2,]  9827 7273 1567  987
# [3,]  7402 5529 1158  715

GSEA_data_GO.db_combined_full <- GSEA_data[which(names(GSEA_data) %in% GO.GSEA_stripped_combined_full)]
GSEA_data_GO.db_combined_3to8 <- GSEA_data[which(names(GSEA_data) %in% GO.GSEA_stripped_combined_3to8)]
GSEA_data_GO.db_combined_4to6 <- GSEA_data[which(names(GSEA_data) %in% GO.GSEA_stripped_combined_4to6)]
GSEA_data_GO.db_BP_full <- GSEA_data[which(names(GSEA_data) %in% GO.GSEA_stripped_BP_full)]
GSEA_data_GO.db_BP_3to8 <- GSEA_data[which(names(GSEA_data) %in% GO.GSEA_stripped_BP_3to8)]
GSEA_data_GO.db_BP_4to6 <- GSEA_data[which(names(GSEA_data) %in% GO.GSEA_stripped_BP_4to6)]
GSEA_data_GO.db_MF_full <- GSEA_data[which(names(GSEA_data) %in% GO.GSEA_stripped_MF_full)]
GSEA_data_GO.db_MF_3to8 <- GSEA_data[which(names(GSEA_data) %in% GO.GSEA_stripped_MF_3to8)]
GSEA_data_GO.db_MF_4to6 <- GSEA_data[which(names(GSEA_data) %in% GO.GSEA_stripped_MF_4to6)]
GSEA_data_GO.db_CC_full <- GSEA_data[which(names(GSEA_data) %in% GO.GSEA_stripped_CC_full)]
GSEA_data_GO.db_CC_3to8 <- GSEA_data[which(names(GSEA_data) %in% GO.GSEA_stripped_CC_3to8)]
GSEA_data_GO.db_CC_4to6 <- GSEA_data[which(names(GSEA_data) %in% GO.GSEA_stripped_CC_4to6)]
rbind(c(length(GSEA_data_GO.db_combined_full), length(GSEA_data_GO.db_BP_full), 
        length(GSEA_data_GO.db_MF_full), length(GSEA_data_GO.db_CC_full)), 
      c(length(GSEA_data_GO.db_combined_3to8), length(GSEA_data_GO.db_BP_3to8), 
        length(GSEA_data_GO.db_MF_3to8), length(GSEA_data_GO.db_CC_3to8)), 
      c(length(GSEA_data_GO.db_combined_4to6), length(GSEA_data_GO.db_BP_4to6), 
        length(GSEA_data_GO.db_MF_4to6), length(GSEA_data_GO.db_CC_4to6)))
# Match expected lengths.
toGmt(GSEA_data_GO.db_combined_full, "GSEA data/GSEA_data_GO.db_combined_full.gmt")
toGmt(GSEA_data_GO.db_combined_3to8, "GSEA data/GSEA_data_GO.db_combined_3to8.gmt")
toGmt(GSEA_data_GO.db_combined_4to6, "GSEA data/GSEA_data_GO.db_combined_4to6.gmt")
toGmt(GSEA_data_GO.db_BP_full, "GSEA data/GSEA_data_GO.db_BP_full.gmt")
toGmt(GSEA_data_GO.db_BP_3to8, "GSEA data/GSEA_data_GO.db_BP_3to8.gmt")
toGmt(GSEA_data_GO.db_BP_4to6, "GSEA data/GSEA_data_GO.db_BP_4to6.gmt")
toGmt(GSEA_data_GO.db_MF_full, "GSEA data/GSEA_data_GO.db_MF_full.gmt")
toGmt(GSEA_data_GO.db_MF_3to8, "GSEA data/GSEA_data_GO.db_MF_3to8.gmt")
toGmt(GSEA_data_GO.db_MF_4to6, "GSEA data/GSEA_data_GO.db_MF_4to6.gmt")
toGmt(GSEA_data_GO.db_CC_full, "GSEA data/GSEA_data_GO.db_CC_full.gmt")
toGmt(GSEA_data_GO.db_CC_3to8, "GSEA data/GSEA_data_GO.db_CC_3to8.gmt")
toGmt(GSEA_data_GO.db_CC_4to6, "GSEA data/GSEA_data_GO.db_CC_4to6.gmt")

## Ran each of these on voom and lnHMdisp brca datasets (so 24 runs altogether).
## Import results, look at stats and test out rrvgo ####
folders <- c(
  "all_full_lnHMdisp.brca.GseaPreranked.1594381710427", 
  "all_full_voom.brca.GseaPreranked.1594382319270", 
  "all_levels_3to8_lnHMdisp.brca.GseaPreranked.1594378947490", 
  "all_levels_3to8_voom.brca.GseaPreranked.1594379546460", 
  "all_levels_4to6_lnHMdisp.brca.GseaPreranked.1594381108970", 
  "all_levels_4to6_voom.brca.GseaPreranked.1594380594371", 
  "BP_full_lnHMdisp.brca.GseaPreranked.1594367212385", 
  "BP_full_voom.brca.GseaPreranked.1594368037535", 
  "BP_levels_3to8_lnHMdisp.brca.GseaPreranked.1594364121920", 
  "BP_levels_3to8_voom.brca.GseaPreranked.1594364995851", 
  "BP_levels_4to6_lnHMdisp.brca.GseaPreranked.1594366506792", 
  "BP_levels_4to6_voom.brca.GseaPreranked.1594365869348", 
  "CC_full_lnHMdisp.brca.GseaPreranked.1594377336137", 
  "CC_full_voom.brca.GseaPreranked.1594377984307", 
  "CC_levels_3to8_lnHMdisp.brca.GseaPreranked.1594376754526", 
  "CC_levels_3to8_voom.brca.GseaPreranked.1594376104051", 
  "CC_levels_4to6_lnHMdisp.brca.GseaPreranked.1594376860767", 
  "CC_levels_4to6_voom.brca.GseaPreranked.1594376929826", 
  "MF_full_lnHMdisp.brca.GseaPreranked.1594378816723", 
  "MF_full_voom.brca.GseaPreranked.1594378730650", 
  "MF_levels_3to8_lnHMdisp.brca.GseaPreranked.1594378455073", 
  "MF_levels_3to8_voom.brca.GseaPreranked.1594378303515", 
  "MF_levels_4to6_lnHMdisp.brca.GseaPreranked.1594378560064", 
  "MF_levels_4to6_voom.brca.GseaPreranked.1594378628283"
)
files <- c(
  "gsea_report_for_na_pos_1594381710427.xls", 
  "gsea_report_for_na_pos_1594382319270.xls", 
  "gsea_report_for_na_pos_1594378947490.xls", 
  "gsea_report_for_na_pos_1594379546460.xls", 
  "gsea_report_for_na_pos_1594381108970.xls", 
  "gsea_report_for_na_pos_1594380594371.xls", 
  "gsea_report_for_na_pos_1594367212385.xls", 
  "gsea_report_for_na_pos_1594368037535.xls", 
  "gsea_report_for_na_pos_1594364121920.xls", 
  "gsea_report_for_na_pos_1594364995851.xls", 
  "gsea_report_for_na_pos_1594366506792.xls", 
  "gsea_report_for_na_pos_1594365869348.xls", 
  "gsea_report_for_na_pos_1594377336137.xls", 
  "gsea_report_for_na_pos_1594377984307.xls", 
  "gsea_report_for_na_pos_1594376754526.xls", 
  "gsea_report_for_na_pos_1594376104051.xls", 
  "gsea_report_for_na_pos_1594376860767.xls", 
  "gsea_report_for_na_pos_1594376929826.xls", 
  "gsea_report_for_na_pos_1594378816723.xls", 
  "gsea_report_for_na_pos_1594378730650.xls", 
  "gsea_report_for_na_pos_1594378455073.xls", 
  "gsea_report_for_na_pos_1594378303515.xls", 
  "gsea_report_for_na_pos_1594378560064.xls", 
  "gsea_report_for_na_pos_1594378628283.xls"
)
analyses <- c(
  "all_full_lnHMdisp", 
  "all_full_voom", 
  "all_3to8_lnHMdisp", 
  "all_3to8_voom", 
  "all_4to6_lnHMdisp", 
  "all_4to6_voom", 
  "BP_full_lnHMdisp", 
  "BP_full_voom", 
  "BP_3to8_lnHMdisp", 
  "BP_3to8_voom", 
  "BP_4to6_lnHMdisp", 
  "BP_4to6_voom", 
  "CC_full_lnHMdisp", 
  "CC_full_voom", 
  "CC_3to8_lnHMdisp", 
  "CC_3to8_voom", 
  "CC_4to6_lnHMdisp", 
  "CC_4to6_voom", 
  "MF_full_lnHMdisp", 
  "MF_full_voom", 
  "MF_3to8_lnHMdisp", 
  "MF_3to8_voom", 
  "MF_4to6_lnHMdisp", 
  "MF_4to6_voom"
)
for (i in 1:24) {
  assign(
    analyses[i], 
    read.delim(
      paste0(
        "Results/GSEA results June 2020/jul10/", 
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

for (i in 1:24) {
  assign(
    paste0("number_q_0.05_", analyses[i]), 
    sum(get(analyses[i])$FDR.q.val < 0.05)
  )
  # assign(
  #   paste0("number_q_0.01_", analyses[i]), 
  #   sum(get(analyses[i])$FDR.q.val < 0.01)
  # )
}
c(number_q_0.05_all_full_lnHMdisp, number_q_0.05_all_full_voom) # 44 746
c(number_q_0.05_all_3to8_lnHMdisp, number_q_0.05_all_3to8_voom) # 38 742
c(number_q_0.05_all_4to6_lnHMdisp, number_q_0.05_all_4to6_voom) # 28 599
c(number_q_0.01_all_full_lnHMdisp, number_q_0.01_all_full_voom) # 5 401
c(number_q_0.01_all_3to8_lnHMdisp, number_q_0.01_all_3to8_voom) # 5 394
c(number_q_0.01_all_4to6_lnHMdisp, number_q_0.01_all_4to6_voom) # 3 313

library(rrvgo)
library(clusterProfiler)
library(GOSemSim)
# simplify() in clusterProfiler only takes output of clusterProfiler enrichment functions 
# as input, so easiest to use rrvgo (and no real reason to choose one or the other anyway).
# First need to create similarity matrix using calculateSimMatrix(), which takes a list of 
# GO terms, an OrgDb object for an organism, the ontology of interest, and which method to 
# use to calculate similarity scores. Then use reduceSimMatrix() with similarity matrix, 
# scores and threshold to group terms. Scores should be transformed so that higher is 
# better (authors suggest -log, but I don't think it matters as long as the direction is 
# right and the order is maintained).
# GOSemSim::godata() us used by default to obtain semdata argument of calculateSimMatrix(), 
# an "object with prepared GO DATA for measuring semantic similarity" - i.e. a GOSemSimDATA 
# object. Similarity measures are "Resnik", "Lin", "Rel", "Jiang" or "Wang", as implemented 
# in GOSemSim.
BP_semdata <- GOSemSim::godata(OrgDb="org.Hs.eg.db", ont="BP")
MF_semdata <- GOSemSim::godata(OrgDb="org.Hs.eg.db", ont="MF")
CC_semdata <- GOSemSim::godata(OrgDb="org.Hs.eg.db", ont="CC")
# Need to map GO terms from GSEA back to ... whatever list of terms GOSemSim/rrvgo use. 
# That is presumably whatever comes from the OrgDb object, in this case org.Hs.eg.db, which 
# has "GOSOURCEURL: http://current.geneontology.org/ontology/go-basic.obo".
head(keys(org.Hs.eg.db, keytype="GO"))
head(keys(org.Hs.eg.db, keytype="GOALL"))
# Not sure what the difference is between these, but GOALL contains all entries in GO. 
# Both are GO IDs rather than names, so looks like I'll need to map MSigDB GO terms back to 
# IDs from GO itself.
# Get GO terms and IDs directly downloaded from GO:
go <- read.delim("GSEA data/go.obo", stringsAsFactors=F)
go_names <- go[grep("name:", go[, 1]), ]
go_names <- gsub("name: ", "", go_names[1:47358])
go_ids <- go[grep("^id:", go[, 1]), ]
go_ids <- gsub("id: ", "", go_ids[1:47358])
go_table <- data.frame(go_ids, go_names)
go_table <- go_table[-grep("^obsolete", go_table$go_names), ]
names(go_table) <- c("id", "term")
# This version is from 2020-06-01. org.Hs.eg.db has GOSOURCEDATE 2020-05-02. Doesn't seem 
# to be an obvious source to download old versions from, so I'll just see how well the two 
# versions match for how and hope that's ok.
# Manipulate terms to allow matching with MSigDB terms:
# Add column with manipulated terms for matching to MSigDB terms
go_term_stripped <- toupper(go_table$term)
# Remove brackets:
go_term_stripped <- gsub("\\(", "", gsub("\\)", "", gsub("\\[", "", gsub("\\]", "", go_term_stripped))))
# Remove commas and points:
go_term_stripped <- gsub("\\.", "", gsub(",", "", go_term_stripped))
# Remove spaces and primes:
go_term_stripped <- gsub(" ", "", gsub("'", "", go_term_stripped))
# Remove hyphens and forward slashes:
go_term_stripped <- gsub("-", "", gsub("/", "", go_term_stripped))
# Change + to "PLUS":
go_term_stripped <- gsub("\\+", "PLUS", go_term_stripped)
# Remove colons:
go_term_stripped <- gsub(":", "", go_term_stripped)
# Remove ">":
go_term_stripped <- gsub(">", "", go_term_stripped)
# Remove "=":
go_term_stripped <- gsub("=", "", go_term_stripped)
go_table <- cbind(go_table, go_term_stripped)
# Do same to MSigDB terms:
library(GSEABase)
GSEA_data <- getGmt("GSEA data/c5.all.v7.1.symbols.gmt", geneIdType=SymbolIdentifier(), collectionType=GOCollection())
go.gsea <- names(GSEA_data)
go.gsea_stripped <- gsub("_", "", go.gsea)
go.gsea_stripped <- gsub("^GO", "", go.gsea_stripped)
gsea_matched_terms <- data.frame(go.gsea, go.gsea_stripped)
length(org.db_stripped) # 44411
length(go.gsea_stripped) # 10192
sum(go.gsea_stripped %in% org.db_stripped) # 10106
# Now have stripped terms to link terms in GSEA results to IDs in GO, which should match IDs 
# in org.Hs.eg.db and hence in GOSemSim/rrvgo. Need to add a column to GSEA results with 
# those IDs (and then to keep things neater, create new results objects that just have IDs 
# and -log10 p-values).
for (i in analyses) {
  temp <- get(i)
  temp$NAME_stripped <- gsea_matched_terms$go.gsea_stripped[
      match(temp$NAME, gsea_matched_terms$go.gsea)
    ]
  temp$GO_ID <- go_table$id[
    match(temp$NAME_stripped, go_table$go_term_stripped)
  ]
  assign(i, temp)
  rm(temp)
}
# Now have results with columns needed for GOSemSim/rrvgo, but there are some NAs in IDs 
# because not all terms in MSigDB could be matched to current GO ontology. I'll see whether 
# I can still run rrvgo like this, in which case this there shouldn't be any need to re-run 
# GSEA with those terms excluded since I'm using raw p-values and they won't change given 
# a different list of terms. It's a bit messy though, so might be better to re-run GSEA 
# only including terms that can be matched anyway. This also hasn't accounted for any 
# mismatch between the current GO ontology (from 2020-06-01, downloaded directly), and the 
# version in org.Hs.eg.db, from 2020-05-02.

# Now have GO terms, orgdb and GOSemSim object, so can run calculateSimMatrix() to get a 
# similarity matrix for each ontology. Instead of getting a similarity matrix for each set 
# of results, will just create overall similarity matrices for each entire ontology, which 
# can then be applied to all analyses for that ontology (I think). (This also seems to mean 
# that I can only look at each ontology separately, but I don't think that's really an issue). 
# Test run using method="Resnik", but need to decide which to use.
BP_matrix <- calculateSimMatrix(x=keys(org.Hs.eg.db, keytype="GOALL"), 
                                orgdb="org.Hs.eg.db", 
                                ont="BP", 
                                method="Resnik")
# Warning: Removed 7483 terms that were not found in orgdb for BP
dim(BP_matrix) # 15262 15262
15262+7483 # 22745
MF_matrix <- calculateSimMatrix(x=keys(org.Hs.eg.db, keytype="GOALL"), 
                                orgdb="org.Hs.eg.db", 
                                ont="MF", 
                                method="Resnik")
# Warning: Removed 18368 terms that were not found in orgdb for MF
dim(MF_matrix) # 4377 4377
4377+18368 # 22745
CC_matrix <- calculateSimMatrix(x=keys(org.Hs.eg.db, keytype="GOALL"), 
                                orgdb="org.Hs.eg.db", 
                                ont="CC", 
                                method="Resnik")
# Warning: Removed 20914 terms that were not found in orgdb for CC
dim(CC_matrix) # 1831 1831
1831+20914 # 22745
15262+4377+1831 # 21470
length(keys(org.Hs.eg.db, keytype="GOALL")) # 22746
22746-21470 # 1276 missing
sum(keys(org.Hs.eg.db, keytype="GOALL") %in% go_table$id) # 22715

# Now test applying rrvgo with these simulation matrices
?reduceSimMatrix
# Will need to create named vectors of scores for each analysis. For now just make one so 
# can run test analysis. Default threshold is 0.7 (which authors call medium).
scores_BP_3to8_lnHMdisp <- -log10(BP_3to8_lnHMdisp$NOM.p.val)
names(scores_BP_3to8_lnHMdisp) <- BP_3to8_lnHMdisp$GO_ID
# Obvious problem of infinite -logp for p=0. Let's see if it allows infinite scores.
reduce_BP_3to8_lnHMdisp <- reduceSimMatrix(simMatrix=BP_matrix, 
                                           scores=scores_BP_3to8_lnHMdisp, 
                                           threshold=0.7, 
                                           orgdb=org.Hs.eg.db)
# Error: scores vector does not contain all terms in the similarity matrix - so need to 
# either run calculateSimMatrix for each analysis, or subset the general matrices.
reduce_BP_3to8_lnHMdisp <- reduceSimMatrix(simMatrix=BP_matrix[
  which(rownames(BP_matrix) %in% names(scores_BP_3to8_lnHMdisp)), 
  which(colnames(BP_matrix) %in% names(scores_BP_3to8_lnHMdisp))], 
  scores=scores_BP_3to8_lnHMdisp, 
  threshold=0.7, 
  orgdb=org.Hs.eg.db)
# Seems to be ok with infinite scores.
class(reduce_BP_3to8_lnHMdisp) # data.frame
dim(reduce_BP_3to8_lnHMdisp) # 1097 8
dim(BP_3to8_lnHMdisp) # 1127 7
# 30 out of 1127 terms removed.
head(reduce_BP_3to8_lnHMdisp)
names(reduce_BP_3to8_lnHMdisp)
sum(reduce_BP_3to8_lnHMdisp$go %in% BP_3to8_lnHMdisp$GO_ID) # 1097
identical(reduce_BP_3to8_lnHMdisp$score, 
          -log10(BP_3to8_lnHMdisp[which(
            BP_3to8_lnHMdisp$GO_ID %in% reduce_BP_3to8_lnHMdisp$go), ]$NOM.p.val)) # FALSE
identical(sort(reduce_BP_3to8_lnHMdisp$score), 
          sort(-log10(BP_3to8_lnHMdisp[which(
            BP_3to8_lnHMdisp$GO_ID %in% reduce_BP_3to8_lnHMdisp$go), ]$NOM.p.val))) # TRUE
# Order is changed, but original (transformed) p-values are retained.
# Looking at help, results still contain all terms, so reducing from 239 to 214 just reflects 
# the number of GO IDs in the results that could be matched to org.Hs.eg.db.
# Output doesn't really seem to match what the help says it should be though: "a data.frame 
# with all terms and its 'reducer' (NA if the term was not reduced)". "reducer" must be the 
# "parent" column, but it's never NA - it looks like if a term isn't reduced, this column has 
# that term's ID.
dim(reduce_BP_3to8_lnHMdisp[which(reduce_BP_3to8_lnHMdisp$go == reduce_BP_3to8_lnHMdisp$parent), ]) # 175 8
# So a lot removed - from 1097 (matching) terms to 175.
terms_to_keep <- reduce_BP_3to8_lnHMdisp$go[which(reduce_BP_3to8_lnHMdisp$go == reduce_BP_3to8_lnHMdisp$parent)]
# Match back to MSigDB terms:
BP_3to8_lnHMdisp_reduced <- BP_3to8_lnHMdisp[which(BP_3to8_lnHMdisp$GO_ID %in% terms_to_keep), ]
head(BP_3to8_lnHMdisp_reduced)
head(BP_3to8_lnHMdisp)
head(reduce_BP_3to8_lnHMdisp)
sum(BP_3to8_lnHMdisp$FDR.q.val < 0.05) # 25
sum(BP_3to8_lnHMdisp_reduced$FDR.q.val < 0.05) # 16
# Number of terms significant at 0.05 not changed by much.
scores_BP_3to8_voom <- -log10(BP_3to8_voom$NOM.p.val)
names(scores_BP_3to8_voom) <- BP_3to8_voom$GO_ID
reduce_BP_3to8_voom <- reduceSimMatrix(simMatrix=BP_matrix[
  which(rownames(BP_matrix) %in% names(scores_BP_3to8_voom)), 
  which(colnames(BP_matrix) %in% names(scores_BP_3to8_voom))], 
  scores=scores_BP_3to8_voom, 
  threshold=0.7, 
  orgdb=org.Hs.eg.db)
dim(reduce_BP_3to8_voom) # 2491 8
dim(BP_3to8_voom) # 2569 7
# 78 out of 2569 terms not matched.
dim(reduce_BP_3to8_voom[which(reduce_BP_3to8_voom$go == reduce_BP_3to8_voom$parent), ]) # 278 8
# 2491 matching terms reduced to 278.
terms_to_keep <- reduce_BP_3to8_voom$go[which(reduce_BP_3to8_voom$go == reduce_BP_3to8_voom$parent)]
# Match back to MSigDB terms:
BP_3to8_voom_reduced <- BP_3to8_voom[which(BP_3to8_voom$GO_ID %in% terms_to_keep), ]
sum(BP_3to8_voom$FDR.q.val < 0.05) # 571
sum(BP_3to8_voom_reduced$FDR.q.val < 0.05) # 138
# Number of significant terms at 0.05 level reduced by quite a lot, but still a  lot more for voom than lnHMdisp.
BP_3to8_voom_reduced_IDs <- BP_3to8_voom_reduced$GO_ID[which(BP_3to8_voom_reduced$FDR.q.val < 0.05)]
BP_3to8_lnHMdisp_reduced_IDs <- BP_3to8_lnHMdisp_reduced$GO_ID[which(BP_3to8_lnHMdisp_reduced$FDR.q.val < 0.05)]
length(BP_3to8_voom_reduced_IDs)
length(BP_3to8_lnHMdisp_reduced_IDs)
sum(BP_3to8_lnHMdisp_reduced_IDs %in% BP_3to8_voom_reduced_IDs)
# No overlap in reduced lists at FDR 0.05
sum(BP_3to8_lnHMdisp$FDR.q.val < 0.01) # 9
sum(BP_3to8_lnHMdisp_reduced$FDR.q.val < 0.01) # 5
sum(BP_3to8_voom$FDR.q.val < 0.01) # 298
sum(BP_3to8_voom_reduced$FDR.q.val < 0.01) # 85

# Try more stringent threshold of 0.5
reduce_BP_3to8_lnHMdisp <- reduceSimMatrix(simMatrix=BP_matrix[
  which(rownames(BP_matrix) %in% names(scores_BP_3to8_lnHMdisp)), 
  which(colnames(BP_matrix) %in% names(scores_BP_3to8_lnHMdisp))], 
  scores=scores_BP_3to8_lnHMdisp, 
  threshold=0.5, 
  orgdb=org.Hs.eg.db)
dim(reduce_BP_3to8_lnHMdisp[which(reduce_BP_3to8_lnHMdisp$go == reduce_BP_3to8_lnHMdisp$parent), ]) # 265 8
# 562 terms remaining compared to 175 with threshold=0.7, which doesn't make sense, because a lower threshold should 
# mean that more terms are removed.
terms_to_keep <- reduce_BP_3to8_lnHMdisp$go[which(reduce_BP_3to8_lnHMdisp$go == reduce_BP_3to8_lnHMdisp$parent)]
BP_3to8_lnHMdisp_reduced <- BP_3to8_lnHMdisp[which(BP_3to8_lnHMdisp$GO_ID %in% terms_to_keep), ]
reduce_BP_3to8_voom <- reduceSimMatrix(simMatrix=BP_matrix[
  which(rownames(BP_matrix) %in% names(scores_BP_3to8_voom)), 
  which(colnames(BP_matrix) %in% names(scores_BP_3to8_voom))], 
  scores=scores_BP_3to8_voom, 
  threshold=0.5, 
  orgdb=org.Hs.eg.db)
dim(reduce_BP_3to8_voom[which(reduce_BP_3to8_voom$go == reduce_BP_3to8_voom$parent), ]) # 1053 8
# Again far fewer terms removed with lower threshold.
# Raised an issue on github about this and author replied to say that the algorithm has been changed but the 
# documentation hasn't been updated yet, so the threshold now means the point at which a hierarchical clustering 
# tree is cut, so a higher threshold means fewer clusters and therefore more terms considered redundant. Author 
# also confirmed that the redundant terms are those for which reducedTerms$go != reducedTerms$parent.

# Try a more stringent threshold of 0.9
reduce_BP_3to8_lnHMdisp <- reduceSimMatrix(simMatrix=BP_matrix[
  which(rownames(BP_matrix) %in% names(scores_BP_3to8_lnHMdisp)), 
  which(colnames(BP_matrix) %in% names(scores_BP_3to8_lnHMdisp))], 
  scores=scores_BP_3to8_lnHMdisp, 
  threshold=0.9, 
  orgdb=org.Hs.eg.db)
dim(reduce_BP_3to8_lnHMdisp[which(reduce_BP_3to8_lnHMdisp$go == reduce_BP_3to8_lnHMdisp$parent), ]) # 35 8
# 35 terms remaining compared to 175 with threshold=0.7
terms_to_keep <- reduce_BP_3to8_lnHMdisp$go[which(reduce_BP_3to8_lnHMdisp$go == reduce_BP_3to8_lnHMdisp$parent)]
BP_3to8_lnHMdisp_reduced <- BP_3to8_lnHMdisp[which(BP_3to8_lnHMdisp$GO_ID %in% terms_to_keep), ]
reduce_BP_3to8_voom <- reduceSimMatrix(simMatrix=BP_matrix[
  which(rownames(BP_matrix) %in% names(scores_BP_3to8_voom)), 
  which(colnames(BP_matrix) %in% names(scores_BP_3to8_voom))], 
  scores=scores_BP_3to8_voom, 
  threshold=0.9, 
  orgdb=org.Hs.eg.db)
dim(reduce_BP_3to8_voom[which(reduce_BP_3to8_voom$go == reduce_BP_3to8_voom$parent), ]) # 48 8
# 48 terms remaining compared to 278 with threshold=0.7
terms_to_keep <- reduce_BP_3to8_voom$go[which(reduce_BP_3to8_voom$go == reduce_BP_3to8_voom$parent)]
BP_3to8_voom_reduced <- BP_3to8_voom[which(BP_3to8_voom$GO_ID %in% terms_to_keep), ]
sum(BP_3to8_lnHMdisp$FDR.q.val < 0.05) # 25
sum(BP_3to8_lnHMdisp_reduced$FDR.q.val < 0.05) # 7
sum(BP_3to8_voom$FDR.q.val < 0.05) # 571
sum(BP_3to8_voom_reduced$FDR.q.val < 0.05) # 30
BP_3to8_lnHMdisp_reduced_IDs <- BP_3to8_lnHMdisp_reduced$GO_ID[which(BP_3to8_lnHMdisp_reduced$FDR.q.val < 0.05)]
BP_3to8_voom_reduced_IDs <- BP_3to8_voom_reduced$GO_ID[which(BP_3to8_voom_reduced$FDR.q.val < 0.05)]
length(BP_3to8_lnHMdisp_reduced_IDs) # 7
length(BP_3to8_voom_reduced_IDs) # 30
sum(BP_3to8_lnHMdisp_reduced_IDs %in% BP_3to8_voom_reduced_IDs)
# No overlap in reduced lists at FDR 0.05


## Test different thresholds for similarity with different level terms excluded ####
# Have, for each ontology, full, 3to8 and 4to6, for brca, for voom and lnHMdisp. 
# Really don't think there's much chance I'll end up using 4to6, but will include here anyway. Want to end up with 
# some number of significant terms that is manageable, in that I can look at each significant term manually, but 
# they also need to be informative, and I need to be sure that I'm not discarding terms that could be more 
# informative.
# For now, just apply rrvgo to GSEA analyses and look at remaining terms with FDRs given by GSEA, but think I will 
# want to run again with the redundant terms excluded. I think it will really be better to not use GSEA at all and 
# figure out how to run an enrichment analysis using say fgsea or something else in R, and integrate the reduction 
# process properly, without having to keep matching between databases and probably losing information (since GSEA 
# terms don't always match GO.db and GO sources, and also have to match gene symbols - if can find a way that just 
# uses ENSEMBL terms, I can include everything in my results and let the functions worry about matching).
# Also still need to decide which similarity measure to use.

# Look at BP initially
rbind(c(number_q_0.05_BP_full_voom, number_q_0.05_BP_3to8_voom, number_q_0.05_BP_4to6_voom), 
      c(number_q_0.05_BP_full_lnHMdisp, number_q_0.05_BP_3to8_lnHMdisp, number_q_0.05_BP_4to6_lnHMdisp))
# Restricting levels doesn't change number of terms much for voom, but does for lnHMdisp. A lot more terms identified 
# for voom, so might be difficult to find a set of criteria that works well for both. Number of significant terms 
# ranges from 495 to 585 for voom, and 14 to 27 for lnHMdisp.

library(rrvgo)
library(GOSemSim)
library(GSEABase)
BP_semdata <- GOSemSim::godata(OrgDb="org.Hs.eg.db", ont="BP")
go <- read.delim("GSEA data/go.obo", stringsAsFactors=F)
go_names <- go[grep("name:", go[, 1]), ]
go_names <- gsub("name: ", "", go_names[1:47358])
go_ids <- go[grep("^id:", go[, 1]), ]
go_ids <- gsub("id: ", "", go_ids[1:47358])
go_table <- data.frame(go_ids, go_names)
go_table <- go_table[-grep("^obsolete", go_table$go_names), ]
names(go_table) <- c("id", "term")
go_term_stripped <- toupper(go_table$term)
go_term_stripped <- gsub("\\(", "", gsub("\\)", "", gsub("\\[", "", gsub("\\]", "", go_term_stripped))))
go_term_stripped <- gsub("\\.", "", gsub(",", "", go_term_stripped))
go_term_stripped <- gsub(" ", "", gsub("'", "", go_term_stripped))
go_term_stripped <- gsub("-", "", gsub("/", "", go_term_stripped))
go_term_stripped <- gsub("\\+", "PLUS", go_term_stripped)
go_term_stripped <- gsub(":", "", go_term_stripped)
go_term_stripped <- gsub(">", "", go_term_stripped)
go_term_stripped <- gsub("=", "", go_term_stripped)
go_table <- cbind(go_table, go_term_stripped)
GSEA_data <- getGmt("GSEA data/c5.all.v7.1.symbols.gmt", geneIdType=SymbolIdentifier(), collectionType=GOCollection())
go.gsea <- names(GSEA_data)
go.gsea_stripped <- gsub("_", "", go.gsea)
go.gsea_stripped <- gsub("^GO", "", go.gsea_stripped)
gsea_matched_terms <- data.frame(go.gsea, go.gsea_stripped)
for (i in analyses) {
  temp <- get(i)
  temp$NAME_stripped <- gsea_matched_terms$go.gsea_stripped[
    match(temp$NAME, gsea_matched_terms$go.gsea)
  ]
  temp$GO_ID <- go_table$id[
    match(temp$NAME_stripped, go_table$go_term_stripped)
  ]
  assign(i, temp)
  rm(temp)
}
BP_matrix <- calculateSimMatrix(x=keys(org.Hs.eg.db, keytype="GOALL"), 
                                orgdb="org.Hs.eg.db", 
                                ont="BP", 
                                method="Resnik")
for (i in c("full", "3to8", "4to6")) {
  for (j in c("voom", "lnHMdisp")) {
    assign(paste0("scores_BP_", i, "_", j), 
           setNames(-log10(get(paste0("BP_", i, "_", j))$NOM.p.val), get(paste0("BP_", i, "_", j))$GO_ID))
  }
}

# Try rrvgo with thresholds of 0.5, 0.7 and 0.9
for (i in c("full", "3to8", "4to6")) {
  for (j in c("voom", "lnHMdisp")) {
    for (t in c(0.5,0.7,0.9)) {
      assign(paste0("rrvgo_BP_", i, "_", j, "_", t),
             reduceSimMatrix(
               simMatrix=BP_matrix[which(rownames(BP_matrix) %in% names(get(paste0("scores_BP_", i, "_", j)))),
                                   which(colnames(BP_matrix) %in% names(get(paste0("scores_BP_", i, "_", j))))],
               scores=get(paste0("scores_BP_", i, "_", j)),
               threshold=t,
               orgdb=org.Hs.eg.db
             ))
      terms_to_keep <- get(paste0("rrvgo_BP_", i, "_", j, "_", t))$go[
        which(get(paste0("rrvgo_BP_", i, "_", j, "_", t))$go == get(paste0("rrvgo_BP_", i, "_", j, "_", t))$parent)
      ]
      assign(paste0("BP_", i, "_", j, "_reduced_", t), 
             get(paste0("BP_", i, "_", j))[which(get(paste0("BP_", i, "_", j))$GO_ID %in% terms_to_keep), ])
    }
  }
}
rm(i,j,t,terms_to_keep)

rrvgo_analyses <- character(0)
for (i in c("full", "3to8", "4to6")) {
  for (j in c("voom", "lnHMdisp")) {
    for (t in c(0.5,0.7,0.9)) {
      rrvgo_analyses <- c(rrvgo_analyses, paste0("BP_", i, "_", j, "_reduced_", t))
    }
  }
}
for (i in 1:18) {
  assign(
    paste0("number_q_0.05_", rrvgo_analyses[i]), 
    sum(get(rrvgo_analyses[i])$FDR.q.val < 0.05)
  )
}
rm(i,j,t)

results <- as.data.frame(
  rbind(c(number_q_0.05_BP_full_voom, number_q_0.05_BP_3to8_voom, number_q_0.05_BP_4to6_voom, 
          number_q_0.05_BP_full_lnHMdisp, number_q_0.05_BP_3to8_lnHMdisp, number_q_0.05_BP_4to6_lnHMdisp), 
        c(number_q_0.05_BP_full_voom_reduced_0.5, number_q_0.05_BP_3to8_voom_reduced_0.5, 
          number_q_0.05_BP_4to6_voom_reduced_0.5, 
          number_q_0.05_BP_full_lnHMdisp_reduced_0.5, number_q_0.05_BP_3to8_lnHMdisp_reduced_0.5, 
          number_q_0.05_BP_4to6_lnHMdisp_reduced_0.5), 
        c(number_q_0.05_BP_full_voom_reduced_0.7, number_q_0.05_BP_3to8_voom_reduced_0.7, 
          number_q_0.05_BP_4to6_voom_reduced_0.7, 
          number_q_0.05_BP_full_lnHMdisp_reduced_0.7, number_q_0.05_BP_3to8_lnHMdisp_reduced_0.7, 
          number_q_0.05_BP_4to6_lnHMdisp_reduced_0.7), 
        c(number_q_0.05_BP_full_voom_reduced_0.9, number_q_0.05_BP_3to8_voom_reduced_0.9, 
          number_q_0.05_BP_4to6_voom_reduced_0.9, 
          number_q_0.05_BP_full_lnHMdisp_reduced_0.9, number_q_0.05_BP_3to8_lnHMdisp_reduced_0.9, 
          number_q_0.05_BP_4to6_lnHMdisp_reduced_0.9))
)
names(results) <- c("voom full", "voom 3to8", "voom 4to6", "lnHMdisp full", "lnHMdisp 3to8", "lnHMdisp 4to6")
results$redundancy_threshold <- c(0, 0.5, 0.7, 0.9)
results
#   voom full voom 3to8 voom 4to6 lnHMdisp full lnHMdisp 3to8 lnHMdisp 4to6 redundancy_threshold
# 1       585       571       495            27            25            14                  0.0
# 2       373       362       324            23            22            12                  0.5
# 3       137       138       124            18            16            10                  0.7
# 4        29        30        25             7             7             7                  0.9
# Difference that restricting levels makes decreases as redundancy threshold increases, and difference between voom and 
# lnHMdisp decreases as well. In terms of manageable numbers of significant terms, really only the highest threshold 
# tested, 0.9, is any good. Still need to see whether the terms that are left are informative and whether there are 
# informative terms that have been discarded. If find that there aren't useful terms with threshold 0.9, could use a 
# lower threshold and a stricter FDR cutoff.

# Look at top 10 terms for each threshold (for voom initially)
BP_full_voom$NAME[1:10]
# regulation of vasculature development, extracellular structure organisation, circulatory system process, 
# regulation of system process, cell cell adhesion via plasma membrane adhesion molecules, homophilic cell cell adhesion 
# via plasma membrane adhesion molecules, regulation of blood pressure, connective tissue development, regulation of ion 
# transmembrane transport, DNA packaging
BP_full_voom_reduced_0.5$NAME[1:10]
# regulation of vasculature development, extracellular structure organisation, circulatory system process, 
# regulation of system process, cell cell adhesion via plasma membrane adhesion molecules, regulation of ion transmembrane 
# transport, DNA packaging, regulation of transmembrane transport, chromatin assembly, ameboidal type cell migration
BP_full_voom_reduced_0.7$NAME[1:10]
# regulation of ion transmembrane transport, regulation of membrane potential, regulation of hormone levels, sensory organ 
# development, G protein coupled receptor signalling pathway couple to cyclic nucleotide second messenger, multicellular 
# organism signalling, DNA conformation change, adenylate cyclase activating G protein coupled receptor signalling pathway, 
# striated muscle cell differentiation, antibacterial humoral response
BP_full_voom_reduced_0.9$NAME[1:10]
# regulation of hormone levels, DNA conformation change, adenylate cyclase activating G protein coupled receptor signalling 
# pathway, antibacterial humoral response, sister chromatid segregation, regulation of transporter activity, behaviour, 
# rhythmic process, microtubule cytoskeleton organisation involved in mitosis, negative regulation of mitotic cell cycle
# Categories seem to get broader as threshold increases. e.g. "behaviour", "rhythmic process" are completely useless. 
# Maybe excluding levels 1 and 2 will help with this, although they're both pretty small, so might need to go down to 
# excluding level 3.
# Same thing but for levels 3 to 8 only:
BP_3to8_voom$NAME[1:10]
# regulation of vasculature development, circulatory system process, extracellular structure organisation, regulation of 
# system process, homophilic cell adhesion via plasma membrane adhesion molecules, regulation of blood pressure, cell cell 
# adhesion via plasma membrane adhesion molecules, regulation of ion transmembrane transport, connective tissue development, 
# DNA packaging
BP_3to8_voom_reduced_0.5$NAME[1:10]
# regulation of vasculature development, circulatory system process, extracellular structure organisation, regulation of 
# system process, homophilic adhesion via plasma membrane adhesion molecules, regulation of blood pressure, regulation of 
# ion transmembrane transport, connective tissue development, regulation of ion transmembrane transport, connective tissue 
# development, DNA packaging, antimicrobial humoral response
BP_3to8_voom_reduced_0.7$NAME[1:10]
# regulation of ion transmembrane transport, connective tissue development, G protein coupled receptor signalling pathway 
# coupled to cyclic nucleotide second messenger, sensory organ development, response to corticosteroid, skeletal system 
# development, multicellular organismal signalling, second messenger mediated signalling, DNA conformational change, 
# nuclear chromosome segregation
BP_3to8_voom_reduced_0.9$NAME[1:10]
# nuclear chromosome segregation, regulation of body fluid levels, humoral immune response, regulation of cation channel 
# activity, epithelial cell proliferation, ovulation cycle, negative regulation of mitotic cell cycle, cellular response to 
# peptide hormone stimulus, sex differentiation, regulation of cell division
# Still some very broad terms, but a lot that look clearly relevant, especially with threshold 0.9. 
# Same thing but for levels 4 to 6 only:
BP_4to6_voom$NAME[1:10]
# regulation of vasculature development, extracellular structure organisation, circulatory system process, regulation of 
# system process, regulation of ion transmembrane transport, cell cell adhesion via plasma membrane adhesion molecules, 
# regulation of blood pressure, homophilic cell adhesion via plasma membrane adhesion molecules, connective tissue 
# development, cyclic nucleotide mediated signalling
BP_4to6_voom_reduced_0.5$NAME[1:10]
# regulation of vasculature development, extracellular structure organisation, circulatory system process, regulation of 
# system process, regulation of ion transmembrane transport, regulation of blood pressure, homophilic cell adhesion via 
# plasma membrane adhesion molecules, cyclin nucleotide mediated signalling, regulation of transmembrane transport, 
# chromatin assembly
BP_4to6_voom_reduced_0.7$NAME[1:10]
# cyclic nucleotide mediated signalling, chromatin assembly, G protein coupled receptor signalling pathway coupled to cyclic 
# nucleotide second messenger, sensory organ development, drug transport, mitotic sister chromatid segregation, regulation 
# of hormone levels, striated muscle cell development, cartilage development
BP_4to6_voom_reduced_0.9$NAME[1:10]
# G protein coupled receptor signalling pathway coupled to cyclic nucleotide second messenger, mitotic sister chromatid
# segregation, response to oxygen levels, neuron projection guidance, cell substrate adhesion, actomyosin structure 
# organisation, negative regulation of hydrolase activity, sex differentiation, regulation of cell division, 
# gland development
# Seems to be less really broad terms when only including levels 4 to 6, but also a lot of really specific ones.
# Comparing different excluded levels for each threshold, there seems to be very little difference for 0.5 except in the 
# order in which terms appear, and I'm not sure why that should be - rrvgo obviously changes the order, but it's not clear 
# how it sorts. Also not much difference in terms for 0.7, and some very general and some very specific terms remain with 
# levels 4 to 6 only. For 0.9, some very general terms removed going from full to levels 3 to 8 - behaviour, rhythmic 
# process (verified that they're not there using grep).

# See how lnHMdisp compares
BP_full_lnHMdisp$NAME[1:10]
# homophilic cell adhesion via plasma membrane adhesion molecules, gene silencing by RNA, double strand break repair, ncRNA 
# metabolic process, nucleic acid phosphodiester bond hydrolysis, RNA 3' end processing, macromolecule methylation, gene 
# silencing, cell cell adhesion via plasma membrane adhesion molecules, DNA recombination
BP_full_lnHMdisp_reduced_0.5$NAME[1:10]
# homophilic cell adhesion via plasma membrane adhesion molecules, gene silencing by RNA, double strand break repair, ncRNA 
# metabolic process, nucleic acid phosphodiester bond hydrolysis, RNA 3' end processing, macromolecule methylation, gene 
# silencing, DNA recombination, tRNA metabolic process
BP_full_lnHMdisp_reduced_0.7$NAME[1:10]
# homophilic cell adhesion via plasma membrane adhesion molecules, gene silencing by RNA, double strand break repair, ncRNA 
# metabolic process, nucleic acid phosphodiester bond hydrolysis, RNA 3' end processing, macromolecule methylation, DNA 
# recombination, RNA localisation telomere organisation
BP_full_lnHMdisp_reduced_0.9$NAME[1:10]
# homophilic cell adhesion via plasma membrane adhesion molecules, nucleic acid phosphodiester bond hydrolysis, 
# macromolecule methylation, RNA localisation, immunoglobulin production, smoothened signalling pathway, muscle filament 
# sliding, protein DNA complex subunit organisation, skeletal muscle adaptation dorsal ventral neural tube patterning
# Not much difference until threshold 0.9, and most of the terms that are removed then look useful; some clearly redundant 
# terms removed with threshold 0.5 and 0.7.
# Same thing but for levels 3 to 8 only:
BP_3to8_lnHMdisp$NAME[1:10]
# homophilic cell adhesion via plasma membrane adhesion molecules, gene silencing by RNA, double strand break repair, ncRNA 
# metabolic process, macromolecule methylation, RNA 3' end processing, gene silencing, cell cell adhesion via plasma 
# membrane adhesion molecules, nucleic acid phosphodiester bond hydrolysis, DNA recombination
BP_3to8_lnHMdisp_reduced_0.5$NAME[1:10]
# homophilic cell adhesion via plasma membrane adhesion molecules, gene silencing by RNA, double strand break repair, ncRNA 
# metabolic process, macromolecule methylation, RNA 3' end processing, gene silencing, nucleic acid phosphodiester bond 
# hydrolysis, DNA recombination, tRNA metabolic process
BP_3to8_lnHMdisp_reduced_0.7$NAME[1:10]
# homophilic cell adhesion via plasma membrane adhesion molecules, gene silencing by RNA, double strand break repair, 
# RNA 3' end processing, gene silencing, nucleic acid phosphodiester bond hydrolysis, DNA recombination, tRNA metabolic 
# process, mRNA processing, RNA localisation, RNA methylation
BP_3to8_lnHMdisp_reduced_0.9$NAME[1:10]
# homophilic cell adhesion via plasma membrane adhesion molecules, double strand break repair, gene silencing, DNA 
# recombination, protein localisation to cilium, methylation, smoothened signalling pathway, muscle filament sliding, 
# immunoglobulin production, protein DNA complex subunit organisation
# Again not much difference until threshold 0.9, and then a few useful looking terms removed, some of which have broader 
# parent terms that may be less informative, and some of which don't have obvious parent terms, and again, the few terms 
# that are removed with threshold 0.5 and 0.7 are generally clearly redundant.
# Same thing but for levels 4 to 6 only:
BP_4to6_lnHMdisp$NAME[1:10]
# homophilic cell adhesion via plasma membrane adhesion molecules, gene silencing by RNA, nucleic acid phosphodiester bond 
# hydrolysis, macromolecule methylation, cell cell adhesion via plasma membrane adhesion molecules, telomere organisation, 
# DNA recombination, RNA localisation, RNA methylation, RNA modification
BP_4to6_lnHMdisp_reduced_0.5$NAME[1:10]
# gene silencing by RNA, nucleic acid phosphodiester bond hydrolysis, macromolecule methylation, cell cell adhesion via 
# plasma membrane adhesion molecules, telomere organisation, DNA recombination, RNA localisation, RNA methylation, DNA 
# strand elongation, regulation of gene expression epigenetic
BP_4to6_lnHMdisp_reduced_0.7$NAME[1:10]
# gene silencing by RNA, nucleic acid phosphodiester bond hydrolysis, macromolecule methylation, cell cell adhesion via 
# plasma membrane adhesion molecules, telomere organisation, DNA recombination, RNA methylation, regulation of gene 
# expression epigenetic, smoothened signalling pathway, regulation of response to DNA damage stimulus
BP_4to6_lnHMdisp_reduced_0.9$NAME[1:10]
# gene silencing by RNA, nucleic acid phosphodiester bond hydrolysis, cell cell adhesion via plasma membrane adhesion 
# molecules, telomere organisation, regulation of gene expression epigenetic, smoothened signalling pathway, regulation of 
# response to DNA damage stimulus, muscle filament sliding, immunoglobulin production, skeletal muscle adaptation
# Again terms that are removed with threshold 0.5 and 0.7 generally look redundant, but still some clear redundancy with 
# 0.7 - macromolecule methylation and RNA methylation.
# Comparing different excluded levels for each threshold, no difference between full and levels 3 to 8 for 0.5, but some 
# quite general terms removed with levels 4 to 6. Not much difference for 0.7, but possibly some general terms that aren't 
# very informative removed with levels 3 to 8. Also not much difference for 0.9.

# Seems clear that some removal of redundant terms helps, but that 0.9 is probably too high a threshold. Either 0.5 or 0.7 
# should be ok. No obviously useful terms removed between 0.5 and 0.7 that I can see, so think should use 0.7 to reduce 
# sizes of lists as much as reasonably possible.
# Exclusion of levels is less clear. Probably want to be conservative with this, although there isn't any clear evidence 
# that removing levels 1 and 2 and below 8 removes useful terms. Level 3 is reasonably small, so can probably inspect it 
# manually to see if it's ok to remove.
# Also need to remember this is all just for BP, and it could be different for MF and CC. That's another reason to be 
# conservative, but I think it's safe to assume that redundancy removal will be similar, so I'm comfortable with using 0.7.
c(length(level3BPmin), length(level3MF), length(level3CCmin)) # 523 151 978
# Actually a lot of terms to look through manually, and very different distribution between ontologies:
c(length(level3BPmin) / length(level1_13BP), 
  length(level3MF) / length(level1_13MF), 
  length(level3CCmin) / length(level1_13CC))
# 1-2% of terms in BP and MF are in level 3 as their highest level, but 23% of terms in CC.
# Safest to keep level 3 terms.
# Didn't seem to be much difference between excluding levels 7 and 8 or not, but excluding them means only a very small 
# proportion of terms are excluded. Looking at samples of terms from each level above, seems safe to exclude level 8, and 
# probably level 7, but to be conservative and still remove a reasonable number of mostly irrelevant terms, will exclude 
# level 8, so keep only levels 3 to 7.

# So will use rrvgo with threshold 0.7, and exclude terms that only appear in levels 1 or 2, or below 7.
# Want to avoid going back and forth between sources removing non-matching terms, so should start with only terms that I 
# can match to the list of terms that will be used by GSEA. For removing levels, I'm using GO.db. For GOSemSim and rrvgo, 
# I'm using org.Hs.eg.db - specifically (presumably) keys(org.Hs.eg.db, keytype="GOALL"). So need to get a list of terms 
# common to GO.db and keys(org.Hs.eg.db, keytype="GOALL"), match them to the terms used by GSEA, and exclude any from GSEA 
# analysis that aren't in both GO.db and keys(org.Hs.eg.db, keytype="GOALL"). keys(org.Hs.eg.db, keytype="GOALL") has IDs 
# only, so the only matching of actual terms will be between GO.db and MSigDB.
library(GO.db)
library(org.Hs.eg.db)
org.Hs_IDs <- keys(org.Hs.eg.db, keytype="GOALL")
length(org.Hs_IDs) # 22746 GOALL, 18107 GO
GO.db_IDs <- GOID(GOTERM)
length(GO.db_IDs) # 44509
IDs_to_match <- intersect(GO.db_IDs, org.Hs_IDs)
length(IDs_to_match) # 22746 GOALL, 18107 GO
# All IDs in org.Hs.eg.db are in GO.db

# Separate into ontologies and keep only terms in levels 3 to 7
# But change terminology to make level below root "level 1" - so I'm only getting rid of level 2 and 
# below level 6.
getAllChildren <- function(goids, ontology)
{
  ans <- unique(unlist(mget(goids, get(paste0("GO", ontology, "CHILDREN"))), use.names=FALSE))
  ans <- ans[!is.na(ans)]
}
level1BP <- getAllChildren("GO:0008150", "BP")
level2BP <- getAllChildren(level1BP, "BP")
level3BP <- getAllChildren(level2BP, "BP")
level4BP <- getAllChildren(level3BP, "BP")
level5BP <- getAllChildren(level4BP, "BP")
level6BP <- getAllChildren(level5BP, "BP")
GO.db_IDs_BP <- unique(c(level2BP, level3BP, level4BP, level5BP, level6BP))[-which(
  unique(c(level2BP, level3BP, level4BP, level5BP, level6BP)) %in% c("GO:0008150", level1BP)
)]
length(GO.db_IDs_BP) # 25554
length(unique(toTable(GOTERM)$go_id[which(toTable(GOTERM)$Ontology == "BP")])) # 29211
# 1271 terms removed (13%)

level1MF <- getAllChildren("GO:0003674", "MF")
level2MF <- getAllChildren(level1MF, "MF")
level3MF <- getAllChildren(level2MF, "MF")
level4MF <- getAllChildren(level3MF, "MF")
level5MF <- getAllChildren(level4MF, "MF")
level6MF <- getAllChildren(level5MF, "MF")
GO.db_IDs_MF <- unique(c(level2MF, level3MF, level4MF, level5MF, level6MF))[-which(
  unique(c(level2MF, level3MF, level4MF, level5MF, level6MF)) %in% c("GO:0003674", level1MF)
)]
# no overlap between levels 0 and 1 and the rest
GO.db_IDs_MF <- unique(c(level2MF, level3MF, level4MF, level5MF, level6MF))
length(GO.db_IDs_MF) # 10067
length(unique(toTable(GOTERM)$go_id[which(toTable(GOTERM)$Ontology == "MF")])) # 11113
# 1046 terms removed (9%)

level1CC <- getAllChildren("GO:0005575", "CC")
level2CC <- getAllChildren(level1CC, "CC")
level3CC <- getAllChildren(level2CC, "CC")
level4CC <- getAllChildren(level3CC, "CC")
level5CC <- getAllChildren(level4CC, "CC")
level6CC <- getAllChildren(level5CC, "CC")
GO.db_IDs_CC <- unique(c(level2CC, level3CC, level4CC, level5CC, level6CC))[-which(
  unique(c(level2CC, level3CC, level4CC, level5CC, level6CC)) %in% c("GO:0005575", level1CC)
)]
length(GO.db_IDs_CC) # 4130
length(unique(toTable(GOTERM)$go_id[which(toTable(GOTERM)$Ontology == "CC")])) # 4184
# 54 terms removed (1%)

IDs_to_match_BP <- intersect(GO.db_IDs_BP, org.Hs_IDs)
IDs_to_match_MF <- intersect(GO.db_IDs_MF, org.Hs_IDs)
IDs_to_match_CC <- intersect(GO.db_IDs_CC, org.Hs_IDs)
length(IDs_to_match_BP) # 14202 GOALL, 10681 GO
length(IDs_to_match_MF) # 4092 GOALL, 3640 GO
length(IDs_to_match_CC) # 1960 GOALL, 1717 GO
# Total 20254, compared to 22746 without splitting and refining for GOALL.
# 12% removed, so consistent with proportions removed before matching to org.Hs.eg.db.
# For GO, total 16038 compared to 18107, so 11% removed.

# Now identify terms from GO.db that have IDs in refined org.Hs.eg.db lists
terms_to_match_BP <- Term(GOTERM)[which(names(Term(GOTERM)) %in% IDs_to_match_BP)]
terms_to_match_MF <- Term(GOTERM)[which(names(Term(GOTERM)) %in% IDs_to_match_MF)]
terms_to_match_CC <- Term(GOTERM)[which(names(Term(GOTERM)) %in% IDs_to_match_CC)]
identical(sort(IDs_to_match_BP), names(terms_to_match_BP)) # TRUE
identical(sort(IDs_to_match_MF), names(terms_to_match_MF)) # TRUE
identical(sort(IDs_to_match_CC), names(terms_to_match_CC)) # TRUE

# Remove formatting from terms to match from GO.db (upper case; remove brackets, commas, points, spaces, primes, 
# hyphens, forward slashes, colons, >, =; change + to PLUS):
terms_to_match_formatted_BP <- toupper(terms_to_match_BP)
terms_to_match_formatted_BP <- toupper(terms_to_match_formatted_BP)
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
terms_to_match_formatted_MF <- toupper(terms_to_match_formatted_MF)
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
terms_to_match_formatted_CC <- toupper(terms_to_match_formatted_CC)
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

# Remove "GO" and underscores to match formatting of terms to match
msigdb_terms_formatted_BP <- gsub("_", "", msigdb_terms_BP)
msigdb_terms_formatted_BP <- gsub("^GO", "", msigdb_terms_formatted_BP)
length(msigdb_terms_formatted_BP) # 7530
sum(msigdb_terms_formatted_BP %in% terms_to_match_formatted_BP) # 6783 GOALL, 5463 GO
# 747 (10%) terms in MSigDB either not in org.Hs.eg.db or in excluded levels for GOALL, 2067 (27%) for GO
msigdb_terms_formatted_MF <- gsub("_", "", msigdb_terms_MF)
msigdb_terms_formatted_MF <- gsub("^GO", "", msigdb_terms_formatted_MF)
length(msigdb_terms_formatted_MF) # 1663
sum(msigdb_terms_formatted_MF %in% terms_to_match_formatted_MF) # 1499 GOALL, 1294 GO
# 164 (10%) terms in MSigDB either not in org.Hs.eg.db or in excluded levels for GOALL, 205 (14%) for GO
msigdb_terms_formatted_CC <- gsub("_", "", msigdb_terms_CC)
msigdb_terms_formatted_CC <- gsub("^GO", "", msigdb_terms_formatted_CC)
length(msigdb_terms_formatted_CC) # 999
sum(msigdb_terms_formatted_CC %in% terms_to_match_formatted_CC) # 986 GOALL, 868 GO
# 13 (1%) terms in MSigDB either not in org.Hs.eg.db or in excluded levels for GOALL, 131 (13%) for GO

# Create similarity matrices
library(GOSemSim)
library(rrvgo)
BP_semdata <- godata(OrgDb="org.Hs.eg.db", ont="BP")
MF_semdata <- godata(OrgDb="org.Hs.eg.db", ont="MF")
CC_semdata <- godata(OrgDb="org.Hs.eg.db", ont="CC")
BP_matrix <- calculateSimMatrix(x=IDs_to_match_BP, orgdb=org.Hs.eg.db, ont="BP", method="Wang", semdata=BP_semdata)
dim(BP_matrix) # 13431 GOALL, 10246 GO
MF_matrix <- calculateSimMatrix(x=IDs_to_match_MF, orgdb=org.Hs.eg.db, ont="MF", method="Resnik", semdata=MF_semdata)
dim(MF_matrix) # 3872 GOALL, 3465 GO
CC_matrix <- calculateSimMatrix(x=IDs_to_match_CC, orgdb=org.Hs.eg.db, ont="CC", method="Resnik", semdata=CC_semdata)
dim(CC_matrix) # 1817 GOALL, 1615 GO
# Not sure why some terms (771 BP, 220 MF, 143 CC) weren't found in orgdb since I've already matched to it.
# Maybe GOSemSim uses GO not GOALL from org.Hs.eg.db.
# Doesn't look like it, because (435, 175, 102) weren't found when I used GO.
sum(IDs_to_match_BP %in% keys(org.Hs.eg.db, keytype="GOALL")) # 14202 (= 13431 + 771)
sum(IDs_to_match_MF %in% keys(org.Hs.eg.db, keytype="GOALL")) # 4092 (= 3872 + 220)
sum(IDs_to_match_CC %in% keys(org.Hs.eg.db, keytype="GOALL")) # 1960 (= 1817 + 143)
sum(IDs_to_match_BP %in% keys(org.Hs.eg.db, keytype="GO")) # 10681 (= 10246 + 435)
sum(IDs_to_match_MF %in% keys(org.Hs.eg.db, keytype="GO")) # 3640 (= 3465 + 175)
sum(IDs_to_match_CC %in% keys(org.Hs.eg.db, keytype="GO")) # 1717 (= 1615 + 102)
# Maybe issue comes from godata()
length(unique(BP_semdata@geneAnno$GO)) + 
  length(unique(MF_semdata@geneAnno$GO)) + 
  length(unique(CC_semdata@geneAnno$GO)) # 18107
sum(unique(c(BP_semdata@geneAnno$GO, MF_semdata@geneAnno$GO, CC_semdata@geneAnno$GO)) %in% 
      keys(org.Hs.eg.db, keytype="GO")) # 18107
sum(unique(c(BP_semdata@geneAnno$GO, MF_semdata@geneAnno$GO, CC_semdata@geneAnno$GO)) %in% 
      keys(org.Hs.eg.db, keytype="GOALL")) # 18107
# godata() seems to get terms from GO entry of org.db (has a GO entry but not a GOALL one, and the GO entry has all  terms 
# in org.HS.eg.db GO slot), but that doesn't make sense because calculateSimMatrix() returns more terms when using GOALL 
# than when using GO.
length(c(rownames(BP_matrix), rownames(MF_matrix), rownames(CC_matrix))) # 19120
sum(c(rownames(BP_matrix), rownames(MF_matrix), rownames(CC_matrix)) %in% 
      c(BP_semdata@geneAnno$GO, MF_semdata@geneAnno$GO, CC_semdata@geneAnno$GO)) # 15326
# So there are 3794 terms in the similarity matrices that aren't in the godata object.
sum(c(IDs_to_match_BP, IDs_to_match_MF, IDs_to_match_CC) %in% keys(org.Hs.eg.db, keytype="GO")) # 16038
missing_IDs <- c(IDs_to_match_BP, IDs_to_match_MF, IDs_to_match_CC)[
  -which(c(IDs_to_match_BP, IDs_to_match_MF, IDs_to_match_CC) %in% 
           c(rownames(BP_matrix), rownames(MF_matrix), rownames(CC_matrix)))
]
length(missing_IDs) # 1134 (= 771 + 220 + 143)
sum(missing_IDs %in% keys(org.Hs.eg.db, keytype="GOALL")) # 1134
# All terms that calculateSimMatrix() says are missing from orgdb are definitely there.
BP_semdata@metadata # GOSOURCEDATE 2020-05-02, same as in org.Hs.eg.db; DBSCHEMAVERSION also same
# looks like godata() takes these directly from the orgdb object supplied.

# So there are two things happening: calculateSimMatrix() is returning results for terms that aren't in the godata object 
# supplied for measuring similarity, and it isn't finding some terms that are in the orgdb object.

length(IDs_to_match_BP) # 14202
length(unique(IDs_to_match_BP)) # 14202
sum(IDs_to_match_BP %in% names(BP_semdata@IC)) # 13431
# So this is where the missing terms are coming from. For some reason, some terms don't appear in the IC (information 
# content) vector produced by godata().

# Found code for calculateSimMatrix(): https://rdrr.io/bioc/rrvgo/src/R/rrvgo.R. Confirms that IC is used to filter out 
# "terms that were not found in orgdb" (along with whether or not terms have ancestors, which all of mine should):
x <- unique(x)
found <- x %in% names(semdata@IC)
hasAncestor <- !is.na(sapply(x, function(x) tryCatch(GOSemSim:::getAncestors(ont)[x], error=function(e) NA)))
if(all(!found)) {
  warning("No terms were found in orgdb for ", ont,
          "\nCheck that the organism and ontology match the ones provided by orgdb")
  return(NA)
} else if(!all(found)) {
  warning("Removed ", length(x) - sum(found), " terms that were not found in orgdb for ", ont)
}
x <- x[found & hasAncestor]
# Similarity matrix is then created using GOSemSim::goSim:
m <- matrix(GOSemSim::goSim(x, x, semData=semdata, measure=method),
            ncol=length(x), dimnames=list(x, x))

# removing terms which the similarity couldn't be calculated
out <- apply(m, 2, function(x) all(is.na(x)))
m[!out, !out]

# I should be able to create the similarity matrix directly using goSim and avoid throwing out terms that are in 
# orgdb but not in the IC vector
simMatrix <- function(x, orgdb, semdata, ont=c("BP", "MF", "CC"), method=c("Resnik", "Lin", "Rel", "Jiang", "Wang")) {
  ont <- match.arg(ont)
  method <- match.arg(method)
  m <- matrix(GOSemSim::goSim(x, x, semData=semdata, measure=method),
              ncol=length(x), dimnames=list(x, x))
  out <- apply(m, 2, function(x) all(is.na(x)))
  m[!out, !out]
}
simMat_BP <- simMatrix(x=IDs_to_match_BP, orgdb=org.Hs.eg.db, ont="BP", method="Resnik", semdata=BP_semdata)
dim(simMat_BP) # 13431 (14202 if don't remove NAs in simMatrix)
simMat_MF <- simMatrix(x=IDs_to_match_MF, orgdb=org.Hs.eg.db, ont="MF", method="Resnik", semdata=MF_semdata)
dim(simMat_MF) # 3872 (4092 if don't remove NAs in simMatrix)
simMat_CC <- simMatrix(x=IDs_to_match_CC, orgdb=org.Hs.eg.db, ont="CC", method="Resnik", semdata=CC_semdata)
dim(simMat_CC) # 1817 (1960 if don't remove NAs in simMatrix)
# So looks there is a reason to remove terms that don't appear in IC, but they end up with similarity NA anyway so 
# are removed in that step. But the similarity measure Wang, which is the default in goSim, is the only one that 
# doesn't use the information content, and running again with Wang shows that it obtains a similarity value for all 
# pairs. It's also the most recently published and most highly cited of the methods. Schlicker is the only other 
# method that was designed specifically for GO terms - the others were all originally developed for natural 
# language taxonomies rather than semantic similarity of GO terms. Found a small number of papers that compared 
# some methods, and apart from proposing new ones that aren't options in goSim, there doesn't seem to be a clear 
# choice between them, but certainly no good reason to not use Wang since that also lets me keep more terms.

# Tried running calculateSimMatrix with method="Wang", and it still discards terms with no IC even though Wang 
# doesn't use IC, so I'll continue with applying goSim directly using Wang. Also don't need to compute IC in 
# godata(), which might speed it up a bit (doesn't matter while I've got the data in the workspace, but will be 
# faster if I need to do it again; might also run to see if it makes any difference to output of simMatrix).
BP_semdata <- godata(OrgDb="org.Hs.eg.db", ont="BP", computeIC=F)
MF_semdata <- godata(OrgDb="org.Hs.eg.db", ont="MF", computeIC=F)
CC_semdata <- godata(OrgDb="org.Hs.eg.db", ont="CC", computeIC=F)
simMat_BP <- simMatrix(x=IDs_to_match_BP, orgdb=org.Hs.eg.db, ont="BP", method="Wang", semdata=BP_semdata)
simMat_MF <- simMatrix(x=IDs_to_match_MF, orgdb=org.Hs.eg.db, ont="MF", method="Wang", semdata=MF_semdata)
simMat_CC <- simMatrix(x=IDs_to_match_CC, orgdb=org.Hs.eg.db, ont="CC", method="Wang", semdata=CC_semdata)
dim(simMat_BP) # 13431 (14202 if don't remove NAs in simMatrix)
dim(simMat_MF) # 3872 (4092 if don't remove NAs in simMatrix)
dim(simMat_CC) # 1817 (1960 if don't remove NAs in simMatrix)

# Takes a long time to run using Wang, and seems to use a lot of memory. Trying both via rrgvo and directly using 
# simMatrix(). My way doesn't need IC, whereas rrvgo does, even though it's not used for similarity calculation, 
# because calculateSimMatrix removes all terms without an IC value. If I end up using rrvgo, the main advantage of 
# using Wang is gone since I still won't be able to get similarities for terms that don't have an IC, so in that 
# case I'd be better to switch to Schlicker.
# Using calculateSimMatrix with Wang for BP, even with fresh R session and minimal workspace, runs for a couple of 
# hours and then exits with a "memory exhausted" error. Same with simMatrix, so looks like I just can't use Wang 
# unless I break the list of terms down into subsets and combine the similarity matrices, but I don't think I've 
# got enough reason to be so committed to Wang as to do that, so I'll try Schlicker (which is called "Rel") and 
# revert to calculateSimMatrix().
BP_semdata <- godata(OrgDb="org.Hs.eg.db", ont="BP")
MF_semdata <- godata(OrgDb="org.Hs.eg.db", ont="MF")
CC_semdata <- godata(OrgDb="org.Hs.eg.db", ont="CC")
BP_matrix <- calculateSimMatrix(x=IDs_to_match_BP, orgdb=org.Hs.eg.db, ont="BP", method="Rel", semdata=BP_semdata)
MF_matrix <- calculateSimMatrix(x=IDs_to_match_MF, orgdb=org.Hs.eg.db, ont="MF", method="Rel", semdata=MF_semdata)
CC_matrix <- calculateSimMatrix(x=IDs_to_match_CC, orgdb=org.Hs.eg.db, ont="CC", method="Rel", semdata=CC_semdata)
dim(BP_matrix) # 13431
dim(MF_matrix) # 3872
dim(CC_matrix) # 1817

# Next step is to re-check specificity of terms included/excluded using whichever similarity measure I end up using, 
# since the threshold of 0.7 was chosen using one method (Resnik), and there's no reason to assume that others will 
# give similar results (although one paper I saw showed correlations between similarity measures of around 0.8, so 
# do expect them to be at least roughly similar, but can't assume they'll be similar enough to rely on inferences 
# from a different measure). But since I'm also using a different set of levels than before, I might as well run 
# GSEA with the new set of terms that I've decided on first, and look at different similarity thresholds for those 
# results. It will also give me a better insight if I only include the terms that will be kept after removing 
# redundancy (i.e. there's no point in including terms in GSEA that will be removed anyway because the IC can't be 
# calculated for them), so need to match remaining terms back to MSigDB terms and export .gmt files to run in GSEA.
# Just do each ontology separately, don't bother combining for overall analysis.

length(msigdb_terms_BP) # 7530 actual BP terms in msigdb
length(msigdb_terms_formatted_BP) # 7530 stripped BP terms in msigdb
length(rownames(BP_matrix)) # 13431 BP IDs to match with msigdb to keep in gmt file (levels 2 to 6 and have IC)
length(terms_to_match_formatted_BP) # 14202 links between IDs and stripped terms through names
sum(terms_to_match_formatted_BP %in% msigdb_terms_formatted_BP) # 6783/7530 msigdb terms in levels 2 to 6 and in orgdb
sum(names(terms_to_match_formatted_BP) %in% rownames(BP_matrix)) # 13431/13431 IDs in similarity matrix can be matched

# Need to select MSigDB terms that match IDs in similarity matrix. Do for BP first.
# First match IDs in similarity matrix to stripped terms
stripped_terms_to_keep_BP <- terms_to_match_formatted_BP[which(
  names(terms_to_match_formatted_BP) %in% rownames(BP_matrix)
)]
length(stripped_terms_to_keep_BP) # 13431 terms, as expected - all terms in BP similarity matrix
head(stripped_terms_to_keep_BP) # stripped terms with IDs as names
# Next select stripped terms in MSigDB to keep in gmt files by matching to stripped terms in similarity matrix
stripped_terms_to_keep_in_MSigDB_BP <- msigdb_terms_formatted_BP[which(
  msigdb_terms_formatted_BP %in% stripped_terms_to_keep_BP
)]
length(stripped_terms_to_keep_in_MSigDB_BP) # 6567/6783 BP MSigDB terms in levels 2 to 6 and in orgdb are in sim matrix
head(stripped_terms_to_keep_in_MSigDB_BP) # stripped terms, as expected
# Now match stripped terms to keep in gmt files to MSigDB terms
terms_to_keep_in_gmt_BP <- msigdb_terms_BP[which(msigdb_terms_formatted_BP %in% stripped_terms_to_keep_in_MSigDB_BP)]
length(terms_to_keep_in_gmt_BP) # 6567, as expected
cbind(head(terms_to_keep_in_gmt_BP), head(stripped_terms_to_keep_in_MSigDB_BP)) # First few terms match
# Subset GSEABase file to reduce to terms want to keep
GSEA_data_matched_BP <- GSEA_data_BP[which(names(GSEA_data_BP) %in% terms_to_keep_in_gmt_BP)]
c(length(GSEA_data_BP), length(GSEA_data_matched_BP)) # 7530 6567, as expected
# Finally, export as gmt file
toGmt(GSEA_data_matched_BP, "GSEA data/GSEA_data_BP_2to6_in_sim_mat.gmt")

# Same for MF
# Match IDs in similarity matrix to stripped terms
stripped_terms_to_keep_MF <- terms_to_match_formatted_MF[which(
  names(terms_to_match_formatted_MF) %in% rownames(MF_matrix)
)]
length(stripped_terms_to_keep_MF) # 3872 terms, as expected - all terms in MF similarity matrix
head(stripped_terms_to_keep_MF) # stripped terms with IDs as names
# Select stripped terms in MSigDB to keep in gmt files by matching to stripped terms in similarity matrix
stripped_terms_to_keep_in_MSigDB_MF <- msigdb_terms_formatted_MF[which(
  msigdb_terms_formatted_MF %in% stripped_terms_to_keep_MF
)]
length(stripped_terms_to_keep_in_MSigDB_MF) # 1450/1500 MF MSigDB terms in levels 2 to 6 and in orgdb are in sim matrix
# But they're not unique terms. Not sure why - would have to go back and investigate - but shouldn't make a difference 
# to the data ultimately used.
head(stripped_terms_to_keep_in_MSigDB_MF) # stripped terms, as expected
# Match stripped terms to keep in gmt files to MSigDB terms
terms_to_keep_in_gmt_MF <- msigdb_terms_MF[which(msigdb_terms_formatted_MF %in% stripped_terms_to_keep_in_MSigDB_MF)]
length(terms_to_keep_in_gmt_MF) # 1450, as expected
cbind(head(terms_to_keep_in_gmt_MF), head(stripped_terms_to_keep_in_MSigDB_MF)) # First few terms match
# Subset GSEABase file to reduce to terms want to keep
GSEA_data_matched_MF <- GSEA_data_MF[which(names(GSEA_data_MF) %in% terms_to_keep_in_gmt_MF)]
c(length(GSEA_data_MF), length(GSEA_data_matched_MF)) # 1663 1450, as expected
# Export as gmt file
toGmt(GSEA_data_matched_MF, "GSEA data/GSEA_data_MF_2to6_in_sim_mat.gmt")

# Same for CC
# Match IDs in similarity matrix to stripped terms
stripped_terms_to_keep_CC <- terms_to_match_formatted_CC[which(
  names(terms_to_match_formatted_CC) %in% rownames(CC_matrix)
)]
length(stripped_terms_to_keep_CC) # 1817 terms, as expected - all terms in CC similarity matrix
head(stripped_terms_to_keep_CC) # stripped terms with IDs as names
# Select stripped terms in MSigDB to keep in gmt files by matching to stripped terms in similarity matrix
stripped_terms_to_keep_in_MSigDB_CC <- msigdb_terms_formatted_CC[which(
  msigdb_terms_formatted_CC %in% stripped_terms_to_keep_CC
)]
length(stripped_terms_to_keep_in_MSigDB_CC) # 929/986 CC MSigDB terms in levels 2 to 6 and in orgdb are in sim matrix
head(stripped_terms_to_keep_in_MSigDB_CC) # stripped terms, as expected
# Match stripped terms to keep in gmt files to MSigDB terms
terms_to_keep_in_gmt_CC <- msigdb_terms_CC[which(msigdb_terms_formatted_CC %in% stripped_terms_to_keep_in_MSigDB_CC)]
length(terms_to_keep_in_gmt_CC) # 929, as expected
cbind(head(terms_to_keep_in_gmt_CC), head(stripped_terms_to_keep_in_MSigDB_CC)) # First few terms match
# Subset GSEABase file to reduce to terms want to keep
GSEA_data_matched_CC <- GSEA_data_CC[which(names(GSEA_data_CC) %in% terms_to_keep_in_gmt_CC)]
c(length(GSEA_data_CC), length(GSEA_data_matched_CC)) # 999 929, as expected
# Export as gmt file
toGmt(GSEA_data_matched_CC, "GSEA data/GSEA_data_CC_2to6_in_sim_mat.gmt")

# Import GSEA results
folders <- c(
  "BP_levels2to6_in_sim_mat_lnHMdisp.brca.GseaPreranked.1594982623706", 
  "BP_levels2to6_in_sim_mat_voom.brca.GseaPreranked.1594983137435", 
  "CC_levels2to6_in_sim_mat_lnHMdisp.brca.GseaPreranked.1594983681808", 
  "CC_levels2to6_in_sim_mat_voom.brca.GseaPreranked.1594983599033", 
  "MF_levels2to6_in_sim_mat_lnHMdisp.brca.GseaPreranked.1594983813669", 
  "MF_levels2to6_in_sim_mat_voom.brca.GseaPreranked.1594983893361"
)
files <- c(
  "gsea_report_for_na_pos_1594982623706.xls", 
  "gsea_report_for_na_pos_1594983137435.xls", 
  "gsea_report_for_na_pos_1594983681808.xls", 
  "gsea_report_for_na_pos_1594983599033.xls", 
  "gsea_report_for_na_pos_1594983813669.xls", 
  "gsea_report_for_na_pos_1594983893361.xls"
)
analyses <- c(
  "BP_2to6_lnHMdisp", 
  "BP_2to6_voom", 
  "CC_2to6_lnHMdisp", 
  "CC_2to6_voom", 
  "MF_2to6_lnHMdisp", 
  "MF_2to6_voom"
)
for (i in 1:6) {
  assign(
    analyses[i], 
    read.delim(
      paste0(
        "Results/GSEA results June 2020/jul17/", 
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

# Apply redundancy removal with thresholds of 0.5, 0.7 and 0.9. Can use similarity matrices already produced.
# Add columns with stripped terms and IDs to GSEA results
terms_to_keep_in_gmt_all <- c(terms_to_keep_in_gmt_BP, terms_to_keep_in_gmt_MF, terms_to_keep_in_gmt_CC)
stripped_terms_to_keep_in_MSigDB_all <- c(
  stripped_terms_to_keep_in_MSigDB_BP, 
  stripped_terms_to_keep_in_MSigDB_MF, 
  stripped_terms_to_keep_in_MSigDB_CC
)
stripped_terms_with_IDs <- c(stripped_terms_to_keep_BP, stripped_terms_to_keep_MF, stripped_terms_to_keep_CC)
for (i in analyses) {
  temp <- get(i)
  temp$NAME_stripped <- stripped_terms_to_keep_in_MSigDB_all[
    match(temp$NAME, terms_to_keep_in_gmt_all)
  ]
  temp$GO_ID <- names(stripped_terms_with_IDs)[
    match(temp$NAME_stripped, stripped_terms_with_IDs)
  ]
  assign(i, temp)
  rm(temp)
}

# Create score vectors from nominal p-values from GSEA results
for (i in c("BP", "MF", "CC")) {
  for (j in c("voom", "lnHMdisp")) {
    assign(paste0("scores_", i, "_2to6_", j), 
           setNames(-log10(get(paste0(i, "_2to6_", j))$NOM.p.val), get(paste0(i, "_2to6_", j))$GO_ID))
  }
}

# Run rrvgo for each ontology, each threshold, voom and lnHMdisp
for (i in c("BP", "MF", "CC")) {
  mat <- get(paste0(i, "_matrix"))
  for (j in c("voom", "lnHMdisp")) {
    for (t in c(0.5,0.7,0.9)) {
      assign(paste0("rrvgo_", i, "_2to6_", j, "_", t),
             reduceSimMatrix(
               simMatrix=mat[which(rownames(mat) %in% names(get(paste0("scores_", i, "_2to6_", j)))),
                             which(colnames(mat) %in% names(get(paste0("scores_", i, "_2to6_", j))))],
               scores=get(paste0("scores_", i, "_2to6_", j)),
               threshold=t,
               orgdb=org.Hs.eg.db
             ))
      terms_to_keep <- get(paste0("rrvgo_", i, "_2to6_", j, "_", t))$go[
        which(get(paste0("rrvgo_", i, "_2to6_", j, "_", t))$go == get(paste0("rrvgo_", i, "_2to6_", j, "_", t))$parent)
      ]
      assign(paste0(i, "_2to6_", j, "_reduced_", t), 
             get(paste0(i, "_2to6_", j))[which(get(paste0(i, "_2to6_", j))$GO_ID %in% terms_to_keep), ])
    }
  }
}
rm(i,j,t,terms_to_keep)

# Get number of significant terms at 0.05 level for each rrvgo analysis
rrvgo_analyses <- character(0)
for (i in c("BP", "MF", "CC")) {
  for (j in c("voom", "lnHMdisp")) {
    for (t in c(0.5,0.7,0.9)) {
      rrvgo_analyses <- c(rrvgo_analyses, paste0(i, "_2to6_", j, "_reduced_", t))
    }
  }
}
for (i in 1:18) {
  assign(
    paste0("number_q_0.05_", rrvgo_analyses[i]), 
    sum(get(rrvgo_analyses[i])$FDR.q.val < 0.05)
  )
}
rm(i,j,t)
significant <- as.data.frame(
  rbind(c(number_q_0.05_BP_2to6_voom_reduced_0.5, number_q_0.05_MF_2to6_voom_reduced_0.5, 
          number_q_0.05_CC_2to6_voom_reduced_0.5, 
          number_q_0.05_BP_2to6_lnHMdisp_reduced_0.5, number_q_0.05_MF_2to6_lnHMdisp_reduced_0.5, 
          number_q_0.05_CC_2to6_lnHMdisp_reduced_0.5), 
        c(number_q_0.05_BP_2to6_voom_reduced_0.7, number_q_0.05_MF_2to6_voom_reduced_0.7, 
          number_q_0.05_CC_2to6_voom_reduced_0.7, 
          number_q_0.05_BP_2to6_lnHMdisp_reduced_0.7, number_q_0.05_MF_2to6_lnHMdisp_reduced_0.7, 
          number_q_0.05_CC_2to6_lnHMdisp_reduced_0.7), 
        c(number_q_0.05_BP_2to6_voom_reduced_0.9, number_q_0.05_MF_2to6_voom_reduced_0.9, 
          number_q_0.05_CC_2to6_voom_reduced_0.9, 
          number_q_0.05_BP_2to6_lnHMdisp_reduced_0.9, number_q_0.05_MF_2to6_lnHMdisp_reduced_0.9, 
          number_q_0.05_CC_2to6_lnHMdisp_reduced_0.9))
)
names(significant) <- c("voom BP", "voom MF", "voom CC", "lnHMdisp BP", "lnHMdisp MF", "lnHMdisp CC")
significant$redundancy_threshold <- c(0.5, 0.7, 0.9)
significant
#   voom BP voom MF voom CC lnHMdisp BP lnHMdisp MF lnHMdisp CC redundancy_threshold
# 1     178      34      39          17           2           6                  0.5
# 2      93      22      30          13           2           5                  0.7
# 3      30      18      15           6           1           4                  0.9
# Again only threshold 0.9 really reduces down to a manageable number of terms for voom, and even then it's still quite a lot 
# for BP. May need to either restrict to a fixed number of top ranked terms or use a stricter FDR threshold.

# Look at top 10 terms for each threshold for BP
BP_2to6_voom_reduced_0.5$NAME[1:10]
# homophilic cell adhesion via plasma membrane adhesion molecules; chromatin assembly; ameboidal type cell migration; 
# cyclic nucleotide mediated signalling; regulation of hormone levels; mitotic sister chromatid segregation; 
# sensory organ development; muscle cell differentiation; G protein coupled receptor signalling pathway coupled to cyclic 
# nucleotide second messenger; cognition
BP_2to6_voom_reduced_0.7$NAME[1:10]
# homophilic cell adhesion via plasma membrane adhesion molecules; chromatin assembly; cyclic nucleotide mediated signalling; 
# regulation of hormone levels; mitotic sister chromatid segregation; sensory organ development; muscle cell differentiation; 
# G protein couple receptor signalling pathway coupled to cyclic nucleotide second messenger; 
# multicellular organismal signalling; negative regulation of vasculature development
# MISSING: ameboidal type cell migration; cognition
rrvgo_BP_2to6_voom_0.7[
  which(
    rrvgo_BP_2to6_voom_0.7$go %in% BP_2to6_voom_reduced_0.5$GO_ID[1:10][
      -which(BP_2to6_voom_reduced_0.5$GO_ID[1:10] %in% BP_2to6_voom_reduced_0.7$GO_ID[1:10])
      ]
    ), c("term", "parentTerm")
  ]
# REPLACED BY: regulation of locomotion; regulation of excretion
# But cognition and regulation of excretion only have similarity 0.475
BP_2to6_voom_reduced_0.9$NAME[1:10]
# homophilic cell adhesion via plasma membrane adhesion molecules; regulation of hormone levels; 
# mitotic sister chromatid segregation; regulation of transporter activity; antibacterial humoral response; ovulation cycle; 
# positive regulation of ion transport; actoymyosin structure organisation; mitotic DNA replication; 
# positive regulation of fatty acid metabolic process
# MISSING: chromatin assembly; cyclic nucleotide mediated signalling; sensory organ development; muscle cell differentiation; 
# G protein couple receptor signalling pathway coupled to cyclic nucleotide second messenger; multicellular organismal 
# signalling; negative regulation of vasculature development
rrvgo_BP_2to6_voom_0.9[
  which(
    rrvgo_BP_2to6_voom_0.9$go %in% BP_2to6_voom_reduced_0.7$GO_ID[1:10][
      -which(BP_2to6_voom_reduced_0.7$GO_ID[1:10] %in% BP_2to6_voom_reduced_0.9$GO_ID[1:10])
    ]
  ), c("term", "parentTerm")
]
# REPLACED BY: actin filament organisation; SMAD protein signal transduction; negative regulation of morphogenesis involved in 
# differentiation; negative regulation of morphogenesis involved in differentiation; SMAD protein signal transduction; 
# SMAD protein signal transduction; negative regulation of cell morphogenesis involved in differentiation
BP_2to6_lnHMdisp_reduced_0.5$NAME[1:10]
# homophilic cell adhesion via plasma membrane adhesion molecules; gene silencing by RNA; RNA 3' end processing; 
# nucleic acid phosphodiester bond hydrolysis; DNA recombination; telomere organisation; RNA methylation; mRNA processing; 
# DNA strand elongation; protein localisation to cilium
BP_2to6_lnHMdisp_reduced_0.7$NAME[1:10]
# homophilic cell adhesion via plasma membrane adhesion molecules; gene silencing by RNA; RNA 3' end processing; 
# nucleic acid phosphodiester bond hydrolysis; telomere organisation; RNA methylation; 
# regulation of gene expression epigenetic; mRNA transport
# MISSING: DNA recombination; mRNA processing; protein localisation to cilium
rrvgo_BP_2to6_lnHMdisp_0.7[
  which(
    rrvgo_BP_2to6_lnHMdisp_0.7$go %in% BP_2to6_lnHMdisp_reduced_0.5$GO_ID[1:10][
      -which(BP_2to6_lnHMdisp_reduced_0.5$GO_ID[1:10] %in% BP_2to6_lnHMdisp_reduced_0.7$GO_ID[1:10])
    ]
  ), c("term", "parentTerm")
]
# REPLACED BY: DNA strand elongation; RNA 3' end processing; mRNA transport
# protein localisation to cilium and mRNA transport only have similarity 0.447
BP_2to6_lnHMdisp_reduced_0.9$NAME[1:10]
# homophilic cell adhesion via plasma membrane adhesion molecules; nucleic acid phosphodiester bond hydrolysis; 
# RNA methylation; mRNA transport; smoothened signalling pathway; protein DNA complex subunit organisation; 
# muscle filament sliding; skeletal muscle adaptation; regulation of viral genome replication; protein folding
# MISSING: gene silencing by RNA; RNA 3' end processing; telomere organisation; DNA strand elongation; methylation; 
# regulation of gene expression epigenetic
rrvgo_BP_2to6_lnHMdisp_0.9[
  which(
    rrvgo_BP_2to6_lnHMdisp_0.9$go %in% BP_2to6_lnHMdisp_reduced_0.7$GO_ID[1:10][
      -which(BP_2to6_lnHMdisp_reduced_0.7$GO_ID[1:10] %in% BP_2to6_lnHMdisp_reduced_0.9$GO_ID[1:10])
    ]
  ), c("term", "parentTerm")
]
# REPLACED BY: RNA methylation; nucleic acid phosphodiester bond hydrolysis; protein-DNA complex subunit organisation; 
# nucleic acid phosphodiester bond hydrolysis; nucleic acid phosphodiester bond hydrolysis; RNA methylation

# Look at top 10 terms for each threshold for MF
MF_2to6_voom_reduced_0.5$NAME[1:10]
# passive transmembrane transporter activity; heparin binding; metal ion transmembrane transporter activity; 
# G protein coupled receptor activity; hormone binding; integrin binding; extracellular matrix structural constituent 
# conferring compression resistance; transmembrane receptor protein tyrosine kinase activity; sulphur compound binding; 
# chloride transmembrane transporter activity
MF_2to6_voom_reduced_0.7$NAME[1:10]
# passive transmembrane transporter activity; heparin binding; hormone binding; integrin binding; extracellular matrix 
# structural constituent conferring compression resistance; transmembrane receptor protein tyrosine kinase activity; 
# sulphur compound binding; G protein coupled receptor binding; peptidase regulator activity; 
# protein heterodimerisation activity
# MISSING: metal ion transmembrane transporter activity; G protein coupled receptor activity; 
# chloride transmembrane transporter activity
rrvgo_MF_2to6_voom_0.7[
  which(
    rrvgo_MF_2to6_voom_0.7$go %in% MF_2to6_voom_reduced_0.5$GO_ID[1:10][
      -which(MF_2to6_voom_reduced_0.5$GO_ID[1:10] %in% MF_2to6_voom_reduced_0.7$GO_ID[1:10])
    ]
  ), c("term", "parentTerm")
]
# REPLACED BY: passive transmembrane transporter activity; peptide receptor activity; passive transmembrane transporter 
# activity
MF_2to6_voom_reduced_0.9$NAME[1:10]
# passive transmembrane transporter activity; heparin binding; hormone binding; integrin binding; extracellular matrix 
# structural constituent conferring compression resistance; sulphur compound binding; G protein coupled receptor binding; 
# peptidase regulator activity; retinol dehydrogenase activity; peptide receptor activity
# MISSING: transmembrane receptor protein tyrosine kinase activity; protein heterodimerisation activity
rrvgo_MF_2to6_voom_0.9[
  which(
    rrvgo_MF_2to6_voom_0.9$go %in% MF_2to6_voom_reduced_0.7$GO_ID[1:10][
      -which(MF_2to6_voom_reduced_0.7$GO_ID[1:10] %in% MF_2to6_voom_reduced_0.9$GO_ID[1:10])
    ]
  ), c("term", "parentTerm")
]
# REPLACED BY: retinol dehydrogenase activity; integrin binding
MF_2to6_lnHMdisp_reduced_0.5$NAME[1:10]
# exonuclease activity active with either ribo or deoxyribonucleic acids and producing 5' phosphomonoesters; 
# histone binding; RNA methyltransferase activity; syntaxin 1 binding; single stranded DNA binding; chromatin binding; 
# telomeric DNA binding; ubiquitin like protein ligase activity; basal transcription machinery binding; 
# RNA polymerase core enzyme binding
MF_2to6_lnHMdisp_reduced_0.7$NAME[1:10]
# exonuclease activity active with either ribo or deoxyribonucleic acids and producing 5' phosphomonoesters; 
# histone binding; RNA methyltransferase activity; syntaxin 1 binding; chromatin binding; telomeric DNA binding; basal 
# transcription machinery binding; RNA polymerase core enzyme binding; structural constituent of nuclear pore; 
# peptide antigen binding
# MISSING: single stranded DNA binding; ubiquitin like protein ligase activity
rrvgo_MF_2to6_lnHMdisp_0.7[
  which(
    rrvgo_MF_2to6_lnHMdisp_0.7$go %in% MF_2to6_lnHMdisp_reduced_0.5$GO_ID[1:10][
      -which(MF_2to6_lnHMdisp_reduced_0.5$GO_ID[1:10] %in% MF_2to6_lnHMdisp_reduced_0.7$GO_ID[1:10])
    ]
  ), c("term", "parentTerm")
]
# REPLACED BY: telomeric DNA binding; RNA methyltransferase activity
MF_2to6_lnHMdisp_reduced_0.9$NAME[1:10]
# exonuclease activity active with either ribo or deoxyribonucleic acids and producing 5' phosphomonoesters; 
# histone binding; syntaxin 1 binding; chromatin binding; telomeric DNA binding; basal transcription machinery binding; 
# RNA polymerase core enzyme binding; structural component of nuclear pore; peptide antigen binding; 
# neurotransmitter receptor activity involved in regulation of postsynaptic membrane potential
# MISSING: RNA methyltransferase activity
rrvgo_MF_2to6_lnHMdisp_0.9[
  which(
    rrvgo_MF_2to6_lnHMdisp_0.9$go %in% MF_2to6_lnHMdisp_reduced_0.7$GO_ID[1:10][
      -which(MF_2to6_lnHMdisp_reduced_0.7$GO_ID[1:10] %in% MF_2to6_lnHMdisp_reduced_0.9$GO_ID[1:10])
    ]
  ), c("term", "parentTerm")
]
# REPLACED BY: exonuclease activity, active with either ribo or deoxyribonucleic acids and producing 5' phosphomonoesters

# Look at top 10 terms for each threshold for CC
CC_2to6_voom_reduced_0.5$NAME[1:10]
# extracellular matrix; DNA packaging complex; transporter complex; immunoglobulin complex; endoplasmic reticulum lumen; 
# sarcolemma; protein DNA complex; contractile fibre; postsynaptic membrane; cell cell junction
CC_2to6_voom_reduced_0.7$NAME[1:10]
# DNA packaging complex; transporter complex; immunoglobulin complex; endoplasmic reticulum lumen; sarcolemma; 
# contractile fibre; postsynaptic membrane; blood microparticle; apical part of cell; collagen trimer
# MISSING: extracellular matrix; protein DNA complex; cell cell junction
rrvgo_CC_2to6_voom_0.7[
  which(
    rrvgo_CC_2to6_voom_0.7$go %in% CC_2to6_voom_reduced_0.5$GO_ID[1:10][
      -which(CC_2to6_voom_reduced_0.5$GO_ID[1:10] %in% CC_2to6_voom_reduced_0.7$GO_ID[1:10])
    ]
  ), c("term", "parentTerm")
]
# REPLACED BY: blood microparticle; DNA packaging complex; presynapse
CC_2to6_voom_reduced_0.9$NAME[1:10]
# DNA packaging complex; endoplasmic reticulum lumen; contractile fibre; postsynaptic membrane; blood microparticle; 
# apical part of cell; cell body membrane; caveola; presynapse; chromosome centromeric region
# MISSING: transporter complex; immunoglobulin complex; sarcolemma; collagen trimer
rrvgo_CC_2to6_voom_0.9[
  which(
    rrvgo_CC_2to6_voom_0.9$go %in% CC_2to6_voom_reduced_0.7$GO_ID[1:10][
      -which(CC_2to6_voom_reduced_0.7$GO_ID[1:10] %in% CC_2to6_voom_reduced_0.9$GO_ID[1:10])
    ]
  ), c("term", "parentTerm")
]
# REPLACED BY: DNA packaging complex; DNA packaging complex; postsynaptic membrane; DNA packaging complex
CC_2to6_lnHMdisp_reduced_0.5$NAME[1:10]
# chromosomal region; endoplasmic reticulum quality control compartment; cilium; nuclear pore; replication fork; 
# transferase complex transferring phosphorus containing groups; integrator complex; intrinsic component of golgi membrane; 
# transport vesicle membrane; mediator complex
CC_2to6_lnHMdisp_reduced_0.7$NAME[1:10]
# chromosomal region; endoplasmic reticulum quality control compartment; cilium; nuclear pore; transferase complex 
# transferring phosphorus containing groups; integrator complex; intrinsic component of golgi membrane; 
# transport vesicle membrane; mediator complex; intraciliary transport particle
# MISSING: replication fork
rrvgo_CC_2to6_lnHMdisp_0.7[
  which(
    rrvgo_CC_2to6_lnHMdisp_0.7$go %in% CC_2to6_lnHMdisp_reduced_0.5$GO_ID[1:10][
      -which(CC_2to6_lnHMdisp_reduced_0.5$GO_ID[1:10] %in% CC_2to6_lnHMdisp_reduced_0.7$GO_ID[1:10])
    ]
  ), c("term", "parentTerm")
]
# REPLACED BY: chromosomal region
CC_2to6_lnHMdisp_reduced_0.9$NAME[1:10]
# chromosomal region; endoplasmic reticulum quality control compartment; cilium; transferase complex transferring 
# phosphorus containing groups; intrinsic component of golgi membrane; transport vesicle membrane; myofilament; 
# trans golgi network membrane; host cellular component; phagophore assembly site membrane
# nuclear pore; integrator complex; mediator complex; intraciliary transport particle
rrvgo_CC_2to6_lnHMdisp_0.9[
  which(
    rrvgo_CC_2to6_lnHMdisp_0.9$go %in% CC_2to6_lnHMdisp_reduced_0.7$GO_ID[1:10][
      -which(CC_2to6_lnHMdisp_reduced_0.7$GO_ID[1:10] %in% CC_2to6_lnHMdisp_reduced_0.9$GO_ID[1:10])
    ]
  ), c("term", "parentTerm")
]
# REPLACED BY: endoplasmic reticulum quality control compartment; transferase complex, transferring phosphorus-containing 
# groups; endoplasmic reticulum quality control compartment; transferase complex, transferring phosphorus-containing groups

# Definitely shouldn't go as high as 0.9 as a lot of terms are replaced by terms that don't seem related, but that also 
# happens going from 0.5 to 0.7, so maybe shouldn't even trust 0.7. But then how do I know I can even trust 0.5?

# Is it viable to use a stricter FDR threshold and a looser similarity threshold?
for (i in 1:18) {
  assign(
    paste0("number_q_0.01_", rrvgo_analyses[i]), 
    sum(get(rrvgo_analyses[i])$FDR.q.val < 0.01)
  )
}
sig_0.01 <- as.data.frame(
  rbind(c(number_q_0.01_BP_2to6_voom_reduced_0.5, number_q_0.01_MF_2to6_voom_reduced_0.5, 
          number_q_0.01_CC_2to6_voom_reduced_0.5, 
          number_q_0.01_BP_2to6_lnHMdisp_reduced_0.5, number_q_0.01_MF_2to6_lnHMdisp_reduced_0.5, 
          number_q_0.01_CC_2to6_lnHMdisp_reduced_0.5), 
        c(number_q_0.01_BP_2to6_voom_reduced_0.7, number_q_0.01_MF_2to6_voom_reduced_0.7, 
          number_q_0.01_CC_2to6_voom_reduced_0.7, 
          number_q_0.01_BP_2to6_lnHMdisp_reduced_0.7, number_q_0.01_MF_2to6_lnHMdisp_reduced_0.7, 
          number_q_0.01_CC_2to6_lnHMdisp_reduced_0.7), 
        c(number_q_0.01_BP_2to6_voom_reduced_0.9, number_q_0.01_MF_2to6_voom_reduced_0.9, 
          number_q_0.01_CC_2to6_voom_reduced_0.9, 
          number_q_0.01_BP_2to6_lnHMdisp_reduced_0.9, number_q_0.01_MF_2to6_lnHMdisp_reduced_0.9, 
          number_q_0.01_CC_2to6_lnHMdisp_reduced_0.9))
)
names(sig_0.01) <- c("voom BP", "voom MF", "voom CC", "lnHMdisp BP", "lnHMdisp MF", "lnHMdisp CC")
sig_0.01$redundancy_threshold <- c(0.5, 0.7, 0.9)
sig_0.01
#   voom BP voom MF voom CC lnHMdisp BP lnHMdisp MF lnHMdisp CC redundancy_threshold
# 1     111      23      25           5           0           1                  0.5
# 2      65      16      21           4           0           1                  0.7
# 3      24      12      11           2           0           1                  0.9
significant
#   voom BP voom MF voom CC lnHMdisp BP lnHMdisp MF lnHMdisp CC redundancy_threshold
# 1     178      34      39          17           2           6                  0.5
# 2      93      22      30          13           2           5                  0.7
# 3      30      18      15           6           1           4                  0.9
# Using an FDR threshold of 0.01 instead of 0.05 doesn't reduce the number of significant terms by much for voom, but 
# does for lnHMdisp.

# Really doubting whether I can trust the similarity values given some of the terms with relatively high similarity 
# that were identified with "Rel" (Schlicker). But this is also the default in REVIGO, and they looked at some examples 
# that suggested that it is reliable. Then can I be sure that I can trust the clustering in rrvgo? With a threshold of 
# 0.7, there were terms with similarities around 0.5 that were considered redundant. Maybe I'll just go with 0.7 as the 
# default in reduceSimMatrix() and assume that there won't be many cases where terms that don't seem obviously related 
# are clustered together.
mean(BP_matrix > 0.5) # 0.014
mean(MF_matrix > 0.5) # 0.012
mean(CC_matrix > 0.5) # 0.018
mean(BP_matrix > 0.7) # 0.002
mean(MF_matrix > 0.7) # 0.002
mean(CC_matrix > 0.7) # 0.004
# Very few pairs with high similarity, so really should be safe, if the similarity measure is reliable, to remove terms 
# with similarity higher than say 0.5. But the other problem is that rrvgo doesn't work as simply as that, and now uses 
# the threshold in the opposite way, so that setting a higher threshold removes more terms. It's really not clear then 
# how the threshold set in rrvgo relates to the actual similarity between terms that will be considered redundant. But 
# reading the author's reply to my issue on github (https://github.com/ssayols/rrvgo/issues/1), he says it uses 
# (1 - simMatrix), so I should be able to use that to set the threshold in an informed way. He also says that the 
# threshold set should be taken as "the expected similarity of terms within a group", but that "this is not entirely 
# correct, and you'll see similarities below this threshold being put in the same group" - which is what looks to have 
# happened above when pairs with similarity less than 0.5 were clustered together with the threshold set at 0.7. Looking 
# at code for reduceSimMatrix(), as the author said, it really just applies hierarchical clustering using (1 - similarity) 
# as the distance, and uses the threshold as the point at which to cut the tree, using cutree(), which takes a tree and 
# either a number of groups or a cut height as input.
# So should roughly consider threshold meanings as 1 - what is described in the help (which looks to be taken directly 
# from REVIGO): so "large" 0.1 (similarity above 0.9 considered redundant), "medium" 0.3 (similarity above 0.7 considered 
# redundant), "small" 0.5 (similarity above 0.5 considered redundant), "tiny" 0.6 (similarity above 0.4 considered 
# redundant). Looking at it this way, it's not surprising that there are terms considered redundant with threshold 0.7, 
# even though the actual similarity is below 0.5, but this implies that a threshold of 0.7 would be extremely strict, and 
# probably the highest that should be considered is 0.6, and even then it should be done with caution (REVIGO site has a 
# warning for using threshold 0.4).
# Overall then, most sensible way to go is to use threshold of 0.5. I'll just need to only look at say up to the top 10 
# terms.

# Start new script to apply the settings I've decided on. Only need to do anything new from the point of running GSEA 
# onwards, since I already have gmt files for each ontology restricted to levels 2 to 6 and including only terms that are 
# in similarity matrix. Include the code for this in new script though so that I have it all together, and not hidden among 
# other code figuring out how to do it. Then start fresh from importing GSEA results for each cancer onwards, and add the 
# final step of re-running GSEA with gmt files that only contain the terms included after removing redundancy (which I expect 
# to be different for each cancer, since the p-values for each term will be different - so even though the similarities used 
# and the clustering results will be the same, the term that is chosen to keep for each cluster won't be the same as it's the 
# term with the lowest p-value (highest score in terms of rrvgo input) that is kept for each cluster).

# But maybe I should try selecting redundant terms to remove without reference to GSEA p-values. That would mean I would be 
# doing the same thing for all cancers, so would only have to run rrvgo once, and would mean there wouldn't be a possible 
# criticism about selecting terms to include based on p-value biasing FDRs if I re-run GSEA analyses after removing redundant 
# terms, which is what I planned to do. The default for rrvgo if scores aren't supplied is to keep the term with the most 
# genes, which would generally mean the broadest term among each cluster. This might end up being good since I think there are 
# probably a lot more too-specific than too-general terms.

# So repeat above analysis but with reduceSimMatrix applied just once for each threshold (instead of once for each threshold 
# for both voom and lnHMdisp), with no scores supplied.

# Run rrvgo for each ontology and each threshold, using full similarity matrices
# There is a problem running reduceSimMatrix without scores - it results in an error in hclust() because of NA/NaN/Inf 
# values in as.dist(1 - simMatrix), which must happen when the similarity matrix is reordered, possibly because it's ordered 
# based on scores before checking if scores is NULL and setting the sizes of the terms as the scores. Looks like I can get 
# around the problem by taking the function that gets the set sizes outside of rrvgo and manually setting the scores as the 
# set sizes.
getGoSize <- function(terms, orgdb) {
  if(all(is(orgdb) != "OrgDb")) {
    orgdb <- loadOrgdb(orgdb)
  }
  
  # get all GO terms with genes associated
  go <- suppressMessages(
    AnnotationDbi::select(orgdb,
                          keytype="ENTREZID",
                          columns=c("GO", "ONTOLOGY"),
                          keys=AnnotationDbi::keys(orgdb, keytype="ENTREZID")))
  go <- go[!is.na(go$GO), ]
  go <- go[go$GO %in% terms, ]
  
  # count
  counts   <- table(go$GO)
  go <- go[go$GO %in% terms, ]
  empty    <- terms[!(terms %in% names(counts))]
  nocounts <- setNames(rep(0, length(empty)), empty)
  
  c(counts, nocounts)
}

# Need to only include terms in similarity matrices that are also in MSigDB to be able to compare properly with previous 
# results using scores to decide which terms to remove, since they only considered terms in MSigDB.
formatted_terms_to_keep_in_BP_matrix <- msigdb_terms_formatted_BP[
  which(msigdb_terms_BP %in% terms_to_keep_in_gmt_BP)
]
IDs_to_keep_in_BP_matrix <- names(terms_to_match_formatted_BP)[
  which(terms_to_match_formatted_BP %in% formatted_terms_to_keep_in_BP_matrix)
]
formatted_terms_to_keep_in_MF_matrix <- msigdb_terms_formatted_MF[
  which(msigdb_terms_MF %in% terms_to_keep_in_gmt_MF)
]
IDs_to_keep_in_MF_matrix <- names(terms_to_match_formatted_MF)[
  which(terms_to_match_formatted_MF %in% formatted_terms_to_keep_in_MF_matrix)
]
formatted_terms_to_keep_in_CC_matrix <- msigdb_terms_formatted_CC[
  which(msigdb_terms_CC %in% terms_to_keep_in_gmt_CC)
]
IDs_to_keep_in_CC_matrix <- names(terms_to_match_formatted_CC)[
  which(terms_to_match_formatted_CC %in% formatted_terms_to_keep_in_CC_matrix)
]

for (i in c("BP", "MF", "CC")) {
    for (t in c(0.5)) {
      mat <- get(paste0(i, "_matrix"))[
        rownames(get(paste0(i, "_matrix"))) %in% get(paste0("IDs_to_keep_in_", i, "_matrix")), 
        colnames(get(paste0(i, "_matrix"))) %in% get(paste0("IDs_to_keep_in_", i, "_matrix"))
      ]
      assign(paste0("rrvgo_", i, "_2to6_", t),
             reduceSimMatrix(
               simMatrix=mat, 
               scores=getGoSize(rownames(mat), org.Hs.eg.db), 
               threshold=t, 
               orgdb=org.Hs.eg.db
             ))
      terms_to_keep <- get(paste0("rrvgo_", i, "_2to6_", t))$go[
        which(get(paste0("rrvgo_", i, "_2to6_", t))$go == get(paste0("rrvgo_", i, "_2to6_", t))$parent)
      ]
      for (j in c("voom", "lnHMdisp")) { 
      assign(paste0(i, "_2to6_", j, "_reduced_", t), 
             get(paste0(i, "_2to6_", j))[which(get(paste0(i, "_2to6_", j))$GO_ID %in% terms_to_keep), ])
      }
    }
}
rm(i,t,mat,terms_to_keep,j)

# Get number of significant terms at 0.05 level for each rrvgo analysis
rrvgo_analyses <- character(0)
for (i in c("BP", "MF", "CC")) {
  for (j in c("voom", "lnHMdisp")) {
    for (t in c(0.5)) {
      rrvgo_analyses <- c(rrvgo_analyses, paste0(i, "_2to6_", j, "_reduced_", t))
    }
  }
}
for (i in 1:6) {
  assign(
    paste0("number_q_0.05_", rrvgo_analyses[i]), 
    sum(get(rrvgo_analyses[i])$FDR.q.val < 0.05)
  )
}
rm(i,j,t)
significant <- as.data.frame(
  rbind(c(number_q_0.05_BP_2to6_voom_reduced_0.5, number_q_0.05_MF_2to6_voom_reduced_0.5, 
          number_q_0.05_CC_2to6_voom_reduced_0.5, 
          number_q_0.05_BP_2to6_lnHMdisp_reduced_0.5, number_q_0.05_MF_2to6_lnHMdisp_reduced_0.5, 
          number_q_0.05_CC_2to6_lnHMdisp_reduced_0.5))
)
names(significant) <- c("voom BP", "voom MF", "voom CC", "lnHMdisp BP", "lnHMdisp MF", "lnHMdisp CC")
significant$redundancy_threshold <- c(0.5)
significant
#   voom BP voom MF voom CC lnHMdisp BP lnHMdisp MF lnHMdisp CC redundancy_threshold
# 1      75      19      24           9           1           5                  0.5
# Compare to results from using scores to decide which terms to keep:
#   voom BP voom MF voom CC lnHMdisp BP lnHMdisp MF lnHMdisp CC redundancy_threshold
# 1     178      34      39          17           2           6                  0.5
# Far fewer significant terms when sizes of sets are used to decide which terms to keep. Unless something else has changed 
# that I haven't realised, that means that there were previously a lot of significant terms in clusters which are now not 
# being selected as the representatives of those clusters. The question is whether that means that by doing it this way, I'm 
# risking missing a lot of significantly enriched terms, or removing likely false positives (since this would mean that 
# closely related terms are not significant), or if it's possible that it doesn't really make much difference at all since I 
# would be re-running GSEA with fewer terms, so that terms that weren't significant at FDR 0.05 previously might become 
# significant.

for (i in 1:6) {
  assign(
    paste0("number_q_0.01_", rrvgo_analyses[i]), 
    sum(get(rrvgo_analyses[i])$FDR.q.val < 0.01)
  )
}
sig_0.01 <- as.data.frame(
  rbind(c(number_q_0.01_BP_2to6_voom_reduced_0.5, number_q_0.01_MF_2to6_voom_reduced_0.5, 
          number_q_0.01_CC_2to6_voom_reduced_0.5, 
          number_q_0.01_BP_2to6_lnHMdisp_reduced_0.5, number_q_0.01_MF_2to6_lnHMdisp_reduced_0.5, 
          number_q_0.01_CC_2to6_lnHMdisp_reduced_0.5))
)
names(sig_0.01) <- c("voom BP", "voom MF", "voom CC", "lnHMdisp BP", "lnHMdisp MF", "lnHMdisp CC")
sig_0.01$redundancy_threshold <- c(0.5)
sig_0.01
#   voom BP voom MF voom CC lnHMdisp BP lnHMdisp MF lnHMdisp CC redundancy_threshold
# 1      36      12      16           4           0           1                  0.5
# Compare to version using score to decide which terms to remove:
#   voom BP voom MF voom CC lnHMdisp BP lnHMdisp MF lnHMdisp CC redundancy_threshold
# 1     111      23      25           5           0           1                  0.5
# Again far fewer significant terms when removing without reference to scores. Now close to a manageable number of significant 
# terms at 0.01 level for voom, but still very few or none for lnHMdisp.

# Look at terms that are left using scores versus using set sizes
rrvgo_BP_2to6_0.5_terms_in_lnHMdisp_results <- rrvgo_BP_2to6_0.5[
  which(rrvgo_BP_2to6_0.5$go %in% rrvgo_BP_2to6_lnHMdisp_0.5$go), 
]
rrvgo_BP_2to6_0.5_terms_in_voom_results <- rrvgo_BP_2to6_0.5[
  which(rrvgo_BP_2to6_0.5$go %in% rrvgo_BP_2to6_voom_0.5$go), 
]
length(unique(rrvgo_BP_2to6_lnHMdisp_0.5$parent)) # 234
length(unique(rrvgo_BP_2to6_0.5_terms_in_lnHMdisp_results$parent)) # 342
length(unique(rrvgo_BP_2to6_voom_0.5$parent)) # 395
length(unique(rrvgo_BP_2to6_0.5_terms_in_voom_results$parent)) # 533
unique(rrvgo_BP_2to6_lnHMdisp_0.5$parentTerm)[
  -which(unique(rrvgo_BP_2to6_lnHMdisp_0.5$parentTerm) %in% 
           unique(rrvgo_BP_2to6_0.5_terms_in_lnHMdisp_results$parentTerm))
] # Kept by using scores, removed by using set sizes
unique(rrvgo_BP_2to6_0.5_terms_in_lnHMdisp_results$parentTerm)[
  -which(unique(rrvgo_BP_2to6_0.5_terms_in_lnHMdisp_results$parentTerm) %in% 
           unique(rrvgo_BP_2to6_lnHMdisp_0.5$parentTerm))
] # Kept by using set sizes, removed by using scores
# Too many terms to really look at properly. 
unique(rrvgo_BP_2to6_voom_0.5$parentTerm)[
  -which(unique(rrvgo_BP_2to6_voom_0.5$parentTerm) %in% 
           unique(rrvgo_BP_2to6_0.5_terms_in_voom_results$parentTerm))
] # Kept by using scores, removed by using set sizes
unique(rrvgo_BP_2to6_0.5_terms_in_voom_results$parentTerm)[
  -which(unique(rrvgo_BP_2to6_0.5_terms_in_voom_results$parentTerm) %in% 
           unique(rrvgo_BP_2to6_voom_0.5$parentTerm))
] # Kept by using set sizes, removed by using scores
# Too many terms to really look at properly. 

rrvgo_MF_2to6_0.5_terms_in_lnHMdisp_results <- rrvgo_MF_2to6_0.5[
  which(rrvgo_MF_2to6_0.5$go %in% rrvgo_MF_2to6_lnHMdisp_0.5$go), 
]
rrvgo_MF_2to6_0.5_terms_in_voom_results <- rrvgo_MF_2to6_0.5[
  which(rrvgo_MF_2to6_0.5$go %in% rrvgo_MF_2to6_voom_0.5$go), 
]
length(unique(rrvgo_MF_2to6_lnHMdisp_0.5$parent)) # 114
length(unique(rrvgo_MF_2to6_0.5_terms_in_lnHMdisp_results$parent)) # 129
length(unique(rrvgo_MF_2to6_voom_0.5$parent)) # 156
length(unique(rrvgo_MF_2to6_0.5_terms_in_voom_results$parent)) # 177
unique(rrvgo_MF_2to6_lnHMdisp_0.5$parentTerm)[
  -which(unique(rrvgo_MF_2to6_lnHMdisp_0.5$parentTerm) %in% 
           unique(rrvgo_MF_2to6_0.5_terms_in_lnHMdisp_results$parentTerm))
] # Kept by using scores, removed by using set sizes
unique(rrvgo_MF_2to6_0.5_terms_in_lnHMdisp_results$parentTerm)[
  -which(unique(rrvgo_MF_2to6_0.5_terms_in_lnHMdisp_results$parentTerm) %in% 
           unique(rrvgo_MF_2to6_lnHMdisp_0.5$parentTerm))
] # Kept by using set sizes, removed by using scores
# Too many terms to really look at properly. 
unique(rrvgo_MF_2to6_voom_0.5$parentTerm)[
  -which(unique(rrvgo_MF_2to6_voom_0.5$parentTerm) %in% 
           unique(rrvgo_MF_2to6_0.5_terms_in_voom_results$parentTerm))
] # Kept by using scores, removed by using set sizes
unique(rrvgo_MF_2to6_0.5_terms_in_voom_results$parentTerm)[
  -which(unique(rrvgo_MF_2to6_0.5_terms_in_voom_results$parentTerm) %in% 
           unique(rrvgo_MF_2to6_voom_0.5$parentTerm))
] # Kept by using set sizes, removed by using scores
# Too many terms to really look at properly. 

rrvgo_CC_2to6_0.5_terms_in_lnHMdisp_results <- rrvgo_CC_2to6_0.5[
  which(rrvgo_CC_2to6_0.5$go %in% rrvgo_CC_2to6_lnHMdisp_0.5$go), 
]
rrvgo_CC_2to6_0.5_terms_in_voom_results <- rrvgo_CC_2to6_0.5[
  which(rrvgo_CC_2to6_0.5$go %in% rrvgo_CC_2to6_voom_0.5$go), 
]
length(unique(rrvgo_CC_2to6_lnHMdisp_0.5$parent)) # 77
length(unique(rrvgo_CC_2to6_0.5_terms_in_lnHMdisp_results$parent)) # 93
length(unique(rrvgo_CC_2to6_voom_0.5$parent)) # 78
length(unique(rrvgo_CC_2to6_0.5_terms_in_voom_results$parent)) # 88
unique(rrvgo_CC_2to6_lnHMdisp_0.5$parentTerm)[
  -which(unique(rrvgo_CC_2to6_lnHMdisp_0.5$parentTerm) %in% 
           unique(rrvgo_CC_2to6_0.5_terms_in_lnHMdisp_results$parentTerm))
] # Kept by using scores, removed by using set sizes
unique(rrvgo_CC_2to6_0.5_terms_in_lnHMdisp_results$parentTerm)[
  -which(unique(rrvgo_CC_2to6_0.5_terms_in_lnHMdisp_results$parentTerm) %in% 
           unique(rrvgo_CC_2to6_lnHMdisp_0.5$parentTerm))
] # Kept by using set sizes, removed by using scores
# Can't really say one set of terms is more useful than the other.
unique(rrvgo_CC_2to6_voom_0.5$parentTerm)[
  -which(unique(rrvgo_CC_2to6_voom_0.5$parentTerm) %in% 
           unique(rrvgo_CC_2to6_0.5_terms_in_voom_results$parentTerm))
] # Kept by using scores, removed by using set sizes
unique(rrvgo_CC_2to6_0.5_terms_in_voom_results$parentTerm)[
  -which(unique(rrvgo_CC_2to6_0.5_terms_in_voom_results$parentTerm) %in% 
           unique(rrvgo_CC_2to6_voom_0.5$parentTerm))
] # Kept by using set sizes, removed by using scores
# Can't really say one set of terms is more useful than the other.

# On reflection, there's not really anything to gain by trying to make sure I'm using the same set of non-redundant terms to 
# start with for each analysis. It's completely valid, and more likely to give useful results, to run each analysis and then 
# remove redundant terms, keeping only the term for each cluster with the lowest p-value. There's no need to run the analyses 
# again with the redundant terms removed to get a new FDR estimate.







