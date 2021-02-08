## Import data and break down by project and specific diagnoses ####
library(here)
library(recount)
data <- readRDS(here("recount data/TCGA", "tcga_matched_pairs.rds"))
data <- scale_counts(data)

# Check sample sizes by ICD10 diagnosis
sort(table(data$gdc_cases.diagnoses.primary_diagnosis), dec=T)
identical(toupper(data$gdc_cases.diagnoses.primary_diagnosis), data$xml_icd_10)
# Both agree
# C50.9 212 breast, unspecified
# C64.9 194 kidney, except renal pelvis (C64)
# C34.1 128 bronchus and lung: upper lobe, bronchus or lung
# C73   118 thyroid gland
# C22.0 100 liver and intrahepatic bile ducts: liver cell carcinoma
# C61    92 prostate
# C34.3  64 bronchus and lung: lower lobe, bronchus or lung
# C64.1  64 kidney, excep renal pelvis (C64)
# C54.1  44 corpus uteri: endometrium
# C18.9  28 colon, unspecified
# C02.9  26 tongue, unspecified
# C14.8  26 overlapping lesio of lip, oral cavity and pharynx
# C16.0  26 stomach: cardia
# C16.2  26 stomach: body of stomach
# C18.2  26 colon: ascending colon
# C15.5  22 oesophagus: lower third of aoesophagus
# C32.9  22 larynx, unspecified
table(data$xml_icd_10[which(data$gdc_cases.project.project_id == "TCGA-BRCA")])
# All C50.x, different areas of breast
table(data$xml_icd_10[which(data$gdc_cases.project.project_id == "TCGA-KIRC")])
# All C64.9 (not in ICD-10, only has C64, kidney except renal pelvis)
table(data$xml_icd_10[which(data$gdc_cases.project.project_id == "TCGA-THCA")])
# All C73, thyroid gland
table(data$xml_icd_10[which(data$gdc_cases.project.project_id == "TCGA-LUAD")])
# All C34.x, different areas of bronchus or lung
table(data$xml_icd_10[which(data$gdc_cases.project.project_id == "TCGA-LIHC")])
# All C22.0, liver cell carcinoma
table(data$xml_icd_10[which(data$gdc_cases.project.project_id == "TCGA-LUSC")])
# All C34.x, different areas of bronchus or lung
table(data$xml_icd_10[which(data$gdc_cases.project.project_id == "TCGA-PRAD")])
# All C61, prostate
# ICD-10 codes don't help - they only specify organ/location.

brca <- data[, data$gdc_cases.project.project_id == "TCGA-BRCA"]
kirc <- data[, data$gdc_cases.project.project_id == "TCGA-KIRC"]
thca <- data[, data$gdc_cases.project.project_id == "TCGA-THCA"]
luad <- data[, data$gdc_cases.project.project_id == "TCGA-LUAD"]
lihc <- data[, data$gdc_cases.project.project_id == "TCGA-LIHC"]
lusc <- data[, data$gdc_cases.project.project_id == "TCGA-LUSC"]
prad <- data[, data$gdc_cases.project.project_id == "TCGA-PRAD"]
hnsc <- data[, data$gdc_cases.project.project_id == "TCGA-HNSC"]
coad <- data[, data$gdc_cases.project.project_id == "TCGA-COAD"]

# Remove NA metadata fields
for (j in c("brca", "kirc", "thca", "luad", "lihc", "lusc", "prad", "hnsc", "coad")) {
  remove <- numeric(0)
  temp <- get(j)
  for (i in 1:ncol(colData(data))) {if (mean(is.na(temp[[i]])) == 1) {remove <- c(remove, i)}}
  colData(temp) <- colData(temp)[, -remove]
  assign(j, temp)
}
rm(remove, temp)
diag.fields <- grep("diag", unique(c(names(colData(brca)), names(colData(kirc)), names(colData(thca)), 
                                     names(colData(luad)), names(colData(lihc)), names(colData(lusc)), 
                                     names(colData(prad)), names(colData(hnsc)), names(colData(coad)))), 
                    value=T)
# gdc_cases.diagnoses.tumor_stage:
# Not reported, Stage I, Ia, Ib, II, IIa, IIb, III, IIIa, IIIb, IIIc, IV, IVa, IVb, IVc
# gdc_cases.diagnoses.morphology:
# Format 8xxx/3; ICD-O codes. /3 means malignant, primary site
# BRCA 85xx, KIRC all 8310, THCA 82xx or 83xx, LUAD 81xx-84xx, LIHC 817x, LUSC 807x or 808x, 
# PRAD 81xx or 85xx, HNSC 807x, COAD 81xx or 84xx
table(brca$gdc_cases.diagnoses.morphology) # 178/212 8500/3 infiltrating duct carcinoma NOS
table(kirc$gdc_cases.diagnoses.morphology) # 144/144 8310/3 clear cell adenocarcinoma NOS
table(thca$gdc_cases.diagnoses.morphology) # 98/118 8260/3 papillary adenocarcinoma NOS
table(luad$gdc_cases.diagnoses.morphology) # 84/116 8140/3 adenocarcinoma NOS
table(lihc$gdc_cases.diagnoses.morphology) # 98/100 8170/3 hepatocellular carcinoma
table(lusc$gdc_cases.diagnoses.morphology) # 88/98 8070/3 squamous cell carcinoma NOS
table(prad$gdc_cases.diagnoses.morphology) # 86/92 8550/3 acinar cell carcinoma
table(hnsc$gdc_cases.diagnoses.morphology) # 80/86 8070/3 squamous cell carcinoma NOS
table(coad$gdc_cases.diagnoses.morphology) # 78/82 8140/3 adenocarcinoma NOS
# cgc_case_histological_diagnosis
table(brca$cgc_case_histological_diagnosis) # 182/224 infiltrating ductal carcinoma
# Others inflitrating lobular, medullary, mucinous, mixed histology, other
table(kirc$cgc_case_histological_diagnosis) # 144/144 kidney clear cell renal carcinoma
table(thca$cgc_case_histological_diagnosis) # 98/118 classical/usual
# Others follicular (>- 99% follicular patterned), tall cell (>= 50% tall cell features)
table(luad$cgc_case_histological_diagnosis) # 86/116 NOS
# Others mixed subtype, briocioalveolar nonmucinous, clear cell, mucinous, papillary, mucinous (colloid)
table(lihc$cgc_case_histological_diagnosis) # 98/100 hepatocellular carcinoma
# Others hepatocholangiocarcinoma mixed
table(lusc$cgc_case_histological_diagnosis) # 88/98 lung squamous cell carcinoma NOS
# Others basaloid, small cell
table(prad$cgc_case_histological_diagnosis) # 90/92 acinar type
# Others other subtype
table(hnsc$cgc_case_histological_diagnosis) # 86/86 head and neck squamous cell carcinoma
table(coad$cgc_case_histological_diagnosis) # 78/82 colon adenocarcinoma
# Others mucinous
# cgc_case_other_histological_diagnosis:
table(brca$cgc_case_other_histological_diagnosis)
# 24 samples with entries, all infiltrating carcinoma; ductal and/or lobular, some apocrine and/or clear cell
# xml_diagnosis:
table(luad$xml_diagnosis) # 116/116 lung adenocarcinoma
table(lusc$xml_diagnosis) # 98/98 lung squamous cell carcinoma

grep("icd", names(colData(data)), value=T)
# cgc_case_icd_o3site - ICD codes
# cgc_case_icd_o3histology - ICD-O codes
# cgc_case_icd10 - ICD codes
# xml_icd_o_3_site - ICD codes
# xml_icd_o_3_histology - ICD-O codes
# xml_icd_10 - ICD codes
grep("morph", names(colData(data)), value=T)
# gdc_cases.diagnoses.morphology - ICD-O codes
# xml_primary_pathology_tumor_morphology_list - percentages of epithelioid cells, no data
# xml_leukemia_french_american_british_morphology_code - no data
grep("histol", names(colData(data)), value=T)
# cgc_case_histological_diagnosis - cancer type descriptions
# cgc_case_other_histological_diagnosis - descriptions for infiltrating carcinoma, 20 entries
# xml_primary_pathology_histological_type - descriptions for various cancer types, 50 entries
# xml_primary_pathology_leiomyosarcoma_histologic_subtype - leiomyosarcoma types, no data
# xml_histological_type - cancer type descriptions
# xml_histological_type_other - cancer type descriptions, 24 entries
# xml_neoplasm_histologic_grade - high/low grade, G1/2/3/4/X/B
# xml_histologic_grading_tier_category - three or four tier, 6 entries
# xml_primary_pathology_neoplasm_histologic_grade - G1/2/3/4/X, 44 entries
# xml_primary_pathology_bone_marrow_sample_histology - concordant/discordant, no data
# xml_primary_pathology_synchronous_tumor_histology_list - no data
# xml_primary_pathology_histology_list - seminoma types, no data
# xml_primary_pathology_histological_type_list - thymoma types, 4 entries

## ICD - site
# gdc_cases.diagnoses.primary_diagnosis
# cgc_case_icd_o3site
# cgc_case_icd10
# xml_icd_o_3_site
# xml_icd_10
sum(toupper(data$gdc_cases.diagnoses.primary_diagnosis) != data$cgc_case_icd_o3site, na.rm=T)
# 308 differences
sum(toupper(data$gdc_cases.diagnoses.primary_diagnosis) != data$cgc_case_icd10, na.rm=T)
# No differences
sum(toupper(data$gdc_cases.diagnoses.primary_diagnosis) != data$xml_icd_o_3_site, na.rm=T)
# 308 differences
sum(toupper(data$gdc_cases.diagnoses.primary_diagnosis) != data$xml_icd_10, na.rm=T)
# No differences
sum(is.na(data$gdc_cases.diagnoses.primary_diagnosis)) # 0
sum(is.na(data$cgc_case_icd10)) # 4
sum(is.na(data$xml_icd_10)) # 0
# So gdc_cases.diagnoses.primary_diagnosis, cgc_case_icd10 and xml_icd_10 are identical except for case and 
# 4 NA in cgc_case_icd10. Use xml_icd_10 since no NA, upper case like most others, and short name.
sum(data$cgc_case_icd_o3site != data$xml_icd_o_3_site, na.rm=T)
# No differences
sum(is.na(data$cgc_case_icd_o3site)) # 4
sum(is.na(data$xml_icd_o_3_site)) # 0
# cgc_case_icd_o3site and xml_icd_o_3_site identical except for 4 NA in cgc_case_icd_o3site, so keep 
# xml_icd_o_3_site.
table(data$xml_icd_o_3_site[which(data$xml_icd_o_3_site != data$xml_icd_10)], 
      data$xml_icd_10[which(data$xml_icd_o_3_site != data$xml_icd_10)])
# Slight differences that probably don't mean anything. Mostly codes ending in .9 in xml_icd_o_3_site which 
# don't have the decimal in xml_icd_10, but also (C34.3, C34.30), (C34.9, C34.1), (C64.9, C64.1). May be slight 
# differences in coding between ICD-O and ICD, or different versions used.
# ICD-O                         ICD    
# C01.9 Base of tongue NOS      C01    Base of tongue
# C19.9 Rectosigmoid junction   C19    Rectosigmoid junction
# C20.9 Rectum NOS              C20    Rectum
# C34.3 Lower lobe, lung        C34.30 (may be error; no C34.30; C34.3 Lower lobe, bronchus or lung)
# C34.9 Lung NOS                C34.1  Bronchus or lung, unspecified
# C37.9 Thymus                  C37    Thymus
# C61.9 Prostate                C61    Prostate
# C64.9 Kidney NOS              C64.1  (may be error; no C64.1; C64 Kidney, except renal pelvis)
# C73.9 Thyroid                 C73    Thyroid
# Essentially no differences. Will use ICD 10 for site, so xml_icd_10.

## ICD-O - histology
# gdc_cases.diagnoses.morphology
# cgc_case_icd_o3histology
# xml_icd_o_3_histology
sum(data$gdc_cases.diagnoses.morphology != data$cgc_case_icd_o3histology, na.rm=T)
# 2 differences
sum(data$gdc_cases.diagnoses.morphology != data$xml_icd_o_3_histology, na.rm=T)
# No differences
identical(data$gdc_cases.diagnoses.morphology, data$xml_icd_o_3_histology)
# TRUE
# gdc_cases.diagnoses.morphology and xml_icd_o_3_histology identical. Use xml_icd_o_3_histology.
sum(data$xml_icd_o_3_histology != data$cgc_case_icd_o3histology, na.rm=T)
# 2 differences
sum(is.na(data$xml_icd_o_3_histology)) # 0
sum(is.na(data$cgc_case_icd_o3histology)) # 4
table(data$cgc_case_icd_o3histology[which(data$xml_icd_o_3_histology != data$cgc_case_icd_o3histology)], 
      data$xml_icd_o_3_histology[which(data$xml_icd_o_3_histology != data$cgc_case_icd_o3histology)])
# Two cases listed as 8500/3 in xml_icd_o_3_histology and 8490/3 in cgc_case_icd_o3histology.
# 8500/3 Infiltrating duct carcinoma NOS
# 8490/3 Signet ring cell carcinoma
data$cgc_case_histological_diagnosis[which(data$xml_icd_o_3_histology != data$cgc_case_icd_o3histology)]
data$xml_histological_type[which(data$xml_icd_o_3_histology != data$cgc_case_icd_o3histology)]
# Both listed as infiltrating ductal carcinoma in both histological description fields, as assume that is right 
# and so use xml_icd_o_3_histology.
table(data$xml_icd_o_3_histology[which(is.na(data$cgc_case_icd_o3histology))])
# Four cases listed as 8070/3 in xml_icd_o_3_histology NA in cgc_case_icd_o3histology.

## Descriptions
# cgc_case_histological_diagnosis
# xml_histological_type
sum(toupper(data$cgc_case_histological_diagnosis) != toupper(data$xml_histological_type), na.rm=T)
# 162 differences, not due do differences in case
unique(paste(data$cgc_case_histological_diagnosis, data$xml_histological_type)[
  which(data$cgc_case_histological_diagnosis != data$xml_histological_type)])
# 10 combinations, 9 of which only differ in punctuation and spacing.
sum(data$cgc_case_histological_diagnosis == "Other  specify" & 
      data$xml_histological_type == "Mixed Histology (please specify)", na.rm=T)
# Two cases.
data$cgc_case_other_histological_diagnosis[which(
  data$cgc_case_histological_diagnosis == "Other  specify" & 
    data$xml_histological_type == "Mixed Histology (please specify)")]
# Both "infiltrating carcinoma with ductal and lobular features", so other and mixed both look right.
sum(is.na(data$cgc_case_histological_diagnosis))
# 60
sum(is.na(data$xml_histological_type))
# 56
# Four cases with entries for xml_histological_type and NA for cgc_case_histological_diagnosis
# Use xml_histological_type since fewer NA, also shorter name and more info, i.e. mixed rather than other.


## Create combined diagnosis field to ensure consistent sample selection, save separate files ####
brca$combined_diag <- paste(brca$xml_icd_o_3_site, brca$xml_icd_o_3_histology, brca$xml_histological_type)
kirc$combined_diag <- paste(kirc$xml_icd_o_3_site, kirc$xml_icd_o_3_histology, kirc$xml_histological_type)
thca$combined_diag <- paste(thca$xml_icd_o_3_site, thca$xml_icd_o_3_histology, thca$xml_histological_type)
luad$combined_diag <- paste(luad$xml_icd_o_3_site, luad$xml_icd_o_3_histology, luad$xml_histological_type)
lihc$combined_diag <- paste(lihc$xml_icd_o_3_site, lihc$xml_icd_o_3_histology, lihc$xml_histological_type)
lusc$combined_diag <- paste(lusc$xml_icd_o_3_site, lusc$xml_icd_o_3_histology, lusc$xml_histological_type)
prad$combined_diag <- paste(prad$xml_icd_o_3_site, prad$xml_icd_o_3_histology, prad$xml_histological_type)
hnsc$combined_diag <- paste(hnsc$xml_icd_o_3_site, hnsc$xml_icd_o_3_histology, hnsc$xml_histological_type)
coad$combined_diag <- paste(coad$xml_icd_o_3_site, coad$xml_icd_o_3_histology, coad$xml_histological_type)

saveRDS(brca, file=here("recount data/TCGA", "brca_matched_pairs.rds"))
saveRDS(kirc, file=here("recount data/TCGA", "kirc_matched_pairs.rds"))
saveRDS(thca, file=here("recount data/TCGA", "thca_matched_pairs.rds"))
saveRDS(luad, file=here("recount data/TCGA", "luad_matched_pairs.rds"))
saveRDS(lihc, file=here("recount data/TCGA", "lihc_matched_pairs.rds"))
saveRDS(lusc, file=here("recount data/TCGA", "lusc_matched_pairs.rds"))
saveRDS(prad, file=here("recount data/TCGA", "prad_matched_pairs.rds"))
saveRDS(hnsc, file=here("recount data/TCGA", "hnsc_matched_pairs.rds"))
saveRDS(coad, file=here("recount data/TCGA", "coad_matched_pairs.rds"))

## Load data and assess which diagnosis categories to use ####
library(here)
library(recount)
brca <- readRDS(here("recount data/TCGA", "brca_matched_pairs.rds"))
kirc <- readRDS(here("recount data/TCGA", "kirc_matched_pairs.rds"))
thca <- readRDS(here("recount data/TCGA", "thca_matched_pairs.rds"))
luad <- readRDS(here("recount data/TCGA", "luad_matched_pairs.rds"))
lihc <- readRDS(here("recount data/TCGA", "lihc_matched_pairs.rds"))
lusc <- readRDS(here("recount data/TCGA", "lusc_matched_pairs.rds"))
prad <- readRDS(here("recount data/TCGA", "prad_matched_pairs.rds"))
hnsc <- readRDS(here("recount data/TCGA", "hnsc_matched_pairs.rds"))
coad <- readRDS(here("recount data/TCGA", "coad_matched_pairs.rds"))

table(brca$combined_diag)
# 168 C50.9 8500/3 Infiltrating Ductal Carcinoma
#  18 C50.9 8524/3 Mixed Histology (please specify)
#  12 C50.9 8522/3 C50.9 8520/3 Infiltrating Lobular Carcinoma
#   4 C50.3 8500/3 Infiltrating Ductal Carcinoma
#   4 C50.4 8500/3 Infiltrating Ductal Carcinoma
#   4 C50.9 8510/3 Medullary Carcinoma
#   2 C50.9 8524/3 Infiltrating Ductal Carcinoma
#   2 C50.3 8523/3 Mucinous Carcinoma
#   2 C50.8 8500/3 Infiltrating Ductal Carcinoma
#   2 C50.9 8541/3 Infiltrating Ductal Carcinoma
#   2 C50.9 8520/3 Mixed Histology (please specify)
#   2 C50.9 8523/3 Other, specify
#   2 C50.9 8575/3 Other, specify
table(kirc$combined_diag)
# 144 C64.9 8310/3 Kidney Clear Cell Renal Carcinoma
table(thca$combined_diag)
# 98 C73.9 8260/3 Thyroid Papillary Carcinoma - Classical/usual
# 14 C73.9 8340/3 Thyroid Papillary Carcinoma - Follicular (>= 99% follicular patterned)
# 6 C73.9 8344/3 Thyroid Papillary Carcinoma - Tall Cell (>= 50% tall cell features)
table(luad$combined_diag)
# 58 C34.1 8140/3 Lung Adenocarcinoma- Not Otherwise Specified (NOS)
# 18 C34.3 8140/3 Lung Adenocarcinoma- Not Otherwise Specified (NOS)
#  6 C34.2 8140/3 Lung Adenocarcinoma- Not Otherwise Specified (NOS)
#  6 C34.3 8252/3 Lung Bronchioloalveolar Carcinoma Nonmucinous
#  4 C34.1 8255/3 Lung Adenocarcinoma Mixed Subtype
#  4 C34.1 8260/3 Lung Papillary Adenocarcinoma
#  4 C34.3 8255/3 Lung Adenocarcinoma Mixed Subtype
#  4 C34.3 8480/3 Lung Mucinous Adenocarcinoma
#  4 C34.3 8480/3 Mucinous (Colloid) Carcinoma
#  2 C34.1 8260/3 Lung Adenocarcinoma- Not Otherwise Specified (NOS)
#  2 C34.1 8310/3 Lung Clear Cell Adenocarcinoma
#  2 C34.1 8480/3 Mucinous (Colloid) Carcinoma
#  2 C34.9 8140/3 Lung Adenocarcinoma- Not Otherwise Specified (NOS)
table(lihc$combined_diag)
# 98 C22.0 8170/3 Hepatocellular Carcinoma
#  2 C22.0 8171/3 Hepatocholangiocarcinoma (Mixed)
table(lusc$combined_diag)
# 52 C34.1 8070/3 Lung Squamous Cell Carcinoma- Not Otherwise Specified (NOS)
# 32 C34.3 8070/3 Lung Squamous Cell Carcinoma- Not Otherwise Specified (NOS)
#  2 C34.1 8083/3 Lung Basaloid Squamous Cell Carcinoma
#  2 C34.2 8070/3 Lung Squamous Cell Carcinoma- Not Otherwise Specified (NOS)
#  2 C34.2 8073/3 Lung Small Cell Squamous Cell Carcinoma
#  2 C34.3 8071/3 Lung Squamous Cell Carcinoma- Not Otherwise Specified (NOS)
#  2 C34.3 8083/3 Lung Basaloid Squamous Cell Carcinoma
#  2 C34.9 8070/3 Lung Squamous Cell Carcinoma- Not Otherwise Specified (NOS)
#  2 C34.9 8071/3 Lung Squamous Cell Carcinoma- Not Otherwise Specified (NOS)
table(prad$combined_diag)
# 86 C61.9 8550/3 Prostate Adenocarcinoma Acinar Type
#  4 C61.9 8140/3 Prostate Adenocarcinoma Acinar Type
#  2 C61.9 8500/3 Prostate Adenocarcinoma, Other Subtype
table(hnsc$combined_diag)
# 26 C14.8 8070/3 Head & Neck Squamous Cell Carcinoma
# 24 C02.9 8070/3 Head & Neck Squamous Cell Carcinoma
# 22 C32.9 8070/3 Head & Neck Squamous Cell Carcinoma
#  4 C01.9 8070/3 Head & Neck Squamous Cell Carcinoma
#  4 C04.9 8070/3 Head & Neck Squamous Cell Carcinoma
#  2 C02.9 8071/3 Head & Neck Squamous Cell Carcinoma
#  2 C04.9 8071/3 Head & Neck Squamous Cell Carcinoma
#  2 C06.9 8071/3 Head & Neck Squamous Cell Carcinoma
table(coad$combined_diag)
# 26 C18.2 8140/3 Colon Adenocarcinoma
# 24 C18.9 8140/3 Colon Adenocarcinoma
# 14 C18.7 8140/3 Colon Adenocarcinoma
#  6 C18.0 8140/3 Colon Adenocarcinoma
#  2 C18.3 8140/3 Colon Adenocarcinoma
#  2 C18.5 8140/3 Colon Adenocarcinoma
#  2 C18.4 8140/3 Colon Adenocarcinoma
#  2 C18.6 8140/3 Colon Adenocarcinoma
#  2 C18.7 8480/3 Colon Mucinous Adenocarcinoma
#  2 C18.9 8480/3 Colon Mucinous Adenocarcinoma

