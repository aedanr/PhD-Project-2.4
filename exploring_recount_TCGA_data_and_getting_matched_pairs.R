library(here)
library(recount)

## Explore TCGA data through recount ####
load(here("recount data/TCGA", "rse_gene.Rdata"))
dim(rse_gene)
# 58037 11284
dim(colData(rse_gene))
# 11284 864
head(names(colData(rse_gene)), 21) # Standard 21 entries of colData for recount files
# Looks like "project" is all "TCGA", "paired_end" is all "TRUE", and all other fields except 
# "reads_downloaded", "mapped_read_count", "auc" and "bigwig_file" are NA.

# Other metadata columns:
# Columns 22-158 begin gdc_ and include info on file names and sizes, data types, platform, 
# centres, submitter, projects, demographic info, tissue sources, diagnoses, outcomes, 
# treatment, other clinical info.
# Columns 159-259 begin cgc_ and look to contain similar info.
# Columns 260-864 begin xml_ and look to contain far more detailed clinical info.
# From f1000research.com/articles/6-1558/v1 - recount workflow:
# "GTEx and TCGA rse objects include additional metadata as available from the raw sources. 
# In particular, we compiled metadata for GTEx using the v6 phenotype information available 
# at gtexportal.org, and we put together a large table of TCGA case and sample information by 
# combining information accumulated across Seven Bridges' Cancer Genomics Cloud and 
# TCGAbiolinks."
# From recount R package manual:
# "For TCGA we acquired metaata information from 3 different sources:
# - GDC: via a json query
# - CGC: via json queries and a custom script to merge the tables
# - TCGAbiolinks: we used to parse GDC's XML files.
# For more information, check 
# https://github.com/leekgroup/recount-website/tree/master/metadata/tcga_prep"
# GDC is NIH Genomic Data Commons: https://portal.gdc.cancer.gov/
# More information for Seven Bridges (CGC) at: 
# https://docs.cancergenomicscloud.org/v1.0/docs/tcga-metadata
info <- colData(rse_gene)
grep("category", names(info), value=T)
# gdc_cases.samples.annotations.category - some reasons to possibly exclude samples
grep("sample", names(info), value=T)
# gdc_cases.samples.oct_embedded
# gdc_cases.samples.is_ffpe / cgc_sample_is_ffpe
# gdc_cases.samples.sample_type / cgc_sample_sample_type
# gdc_cases.samples.portions.annotations.notes - data for one sample only which was originally 
# wrongly identified as primary but is metatstatic. Has "Barcode incorrect" in 
# gdc_cases.samples.annotations.category.
# gdc_cases.samples.annotations.notes - extensive notes on a small number of samples with 
# possible issues.
grep("project", names(info), value=T)
# gdc_cases.project.name / gdc_cases.tissue_source_site.project - more or less specific cancer 
# types, identical to each other except for capitalisation and spacing
# gdc_cases.project.primary_site - site of primary tumour
# gdc_cases.project.project_id - 3/4 letter TCGA project codes
grep("type", names(info), value=T)
# cgc_file_disease_type - identical to gdc_cases.project.name except for one NA
# xml_histological_type
# xml_histological_type_other
# xml_tumor_type - only for a few samples; Primary, Type 1 or Type 2
# xml_diagnosis_subtype - only for some samples; Non-Papillary or Papillary
grep("cancer", names(info), value=T)
# xml_person_neoplasm_cancer_status
# xml_colorectal_cancer
grep("tumor", names(info), value=T)
# gdc_cases.diagnoses.tumor_stage
# cgc_case_tumor_status
# xml_primary_pathology_tumor_tissue_sites
# xml_primary_pathology_tumor_tissue_site_other
# xml_tumor_tissue_site
# xml_primary_pathology_tumor_tissue_site
# xml_tumor_location - for brain tumours only
grep("id", names(info), value=T)
# gdc_cases.case_id
# gdc_cases.demographic.demographic_id
# gdc_cases.project.project_id
# gdc_cases.project.program.program_id
# gdc_cases.tissue_source_site.tissue_source_site_id
# gdc_cases.samples.sample_type_id
# gdc_cases.samples.sample_id
# cgc_case_id
# xml_patient_id
names(info)[-unique(c(
  grep("category", names(info)), 
  grep("sample", names(info)), 
  grep("project", names(info)), 
  grep("type", names(info)), 
  grep("cancer", names(info)), 
  grep("tumor", names(info)), 
  grep("id", names(info))
))]
# gdc_platform
# gdc_cases.diagnoses.primary_diagnosis / gdc_cases.diagnoses.tissue_or_organ_of_origin / 
# cgc_case_icd_o3site / cgc_case_icd10 / xml_icd_o_3_site
# cgc_file_investigation - 3/4 letter TCGA codes
# cgc_case_primary_site - mostly matches gdc_cases.project.primary_site; check
# cgc_case_histological_diagnosis - mostly matches xml_histological_type; check
# cgc_case_clinical_stage
# cgc_case_pathologic_stage - mostly matches gdc_cases.diagnoses.tumor_stage
# cgc_case_icd_o3histology / xml_icd_o_3_histology
# cgc_case_other_histological_diagnosis
# xml_other_dx
# xml_distant_metastasis_present_ind2
# xml_diagnosis - only for lung cancer, Squamous Cell Carcinoma or Adenocarcinoma
# xml_metastatic_site - mostly empty, must only be for a small subset
# xml_other_metastatic_site - mostly empty, must only be for a small subset
# xml_metastatic_site_list - mostly empty, must only be for a small subset

case_project_info <- info[, unique(c(
  grep("case", names(info)), 
  grep("project", names(info)), 
  grep("program", names(info)), 
  grep("sample", names(info)), 
  grep("type", names(info)), 
  grep("patient", names(info))
))]
dim(case_project_info)
# 211 columns
# Relevant for finding and verifying cases with multiple samples:
case_project_info[1, c('gdc_cases.case_id', # long alphanumeric identifier
                       'gdc_cases.submitter_id', # TCGA-xx-xxxx
                       'gdc_cases.project.name', # cancer type
                       'gdc_cases.project.project_id', # TCGA-xxxx
                       'gdc_cases.samples.sample_type_id', # two digit number
                       'gdc_cases.samples.sample_type', # sample type description
                       'cgc_file_case', # long alphanumeric identifier
                       'cgc_case_sample', # long alphanumeric identifier
                       'cgc_case_id', # TCGA-xx-xxxx
                       'cgc_case_program', # TCGA
                       'cgc_file_disease_type', # cancer type
                       'cgc_sample_sample_type_code', # integer
                       'cgc_sample_sample_type', # sample type description
                       'xml_bcr_patient_barcode', # TCGA-xx-xxxx
                       'xml_patient_id')] # xxxx
# gdc_cases.case_id and gdc_cases.submitter_id have a one-to-one correspondence
# gdc_cases.project.name and gdc_cases.project.project_id have a one-to-one correspondence
# gdc_cases.samples.sample_type_id and gdc_cases.samples.sample_type have a one-to-one correspondence
# cgc_file_case, cgc_case_sample and cgc_case_id have a one-to-one correspondence and one NA
# cgc_case_program has only one value
# cgc_sample_sample_type_code and cgc_sample_sample_type have a one-to-one correspondence and one NA
# xml_bcr_patient_barcode and xml_patient_id have a one-to-one correspondence and 39 NA
case_info <- info[, c('gdc_cases.submitter_id', # Patient: TCGA-xx-xxxx
                      'gdc_cases.project.name', # Disease: cancer type
                      'gdc_cases.project.project_id', # Disease: TCGA-xxxx
                      'gdc_cases.samples.sample_type', # Sample type: sample type description
                      'cgc_case_id', # Patient: TCGA-xx-xxxx
                      'cgc_file_disease_type', # Disease: cancer type
                      'cgc_sample_sample_type', # Sample type: sample type description
                      'xml_bcr_patient_barcode')] # Patient: TCGA-xx-xxxx
# Three separate sources for patient, two for cancer/project, two for sample type

sum(case_info$gdc_cases.submitter_id != case_info$cgc_case_id, na.rm=T)
sum(case_info$gdc_cases.submitter_id != case_info$xml_bcr_patient_barcode, na.rm=T)
sum(case_info$cgc_case_id != case_info$cgc_case_id, na.rm=T)
# Only differences in patient IDs are NAs
table(table(case_info$gdc_cases.submitter_id))
table(table(case_info$cgc_case_id))
table(table(case_info$xml_bcr_patient_barcode))
table(table(paste(case_info$gdc_cases.submitter_id, 
                  case_info$cgc_case_id)))
# One patient ID in gdc, which has one sample, is in NA in cgc
table(table(paste(case_info$gdc_cases.submitter_id, 
                  case_info$xml_bcr_patient_barcode)))
# 37 patients with one sample and one patient with two samples in gdc are NA in xml

sum(case_info$gdc_cases.project.name != case_info$cgc_file_disease_type, na.rm=T)
# Only differences in cancer/project names are NAs
table(table(case_info$gdc_cases.project.name))
table(table(case_info$cgc_file_disease_type))
table(table(paste(case_info$gdc_cases.project.name, 
                  case_info$cgc_file_disease_type)))
case_info$gdc_cases.project.name[which(is.na(case_info$cgc_file_disease_type))]
# Stomach Adenocarcinoma has one entry in gdc which is NA in cgc

sum(case_info$gdc_cases.samples.sample_type != case_info$cgc_sample_sample_type, na.rm=T)
# Only differences in sample types are NAs
table(table(case_info$gdc_cases.samples.sample_type))
table(table(case_info$cgc_sample_sample_type))
table(table(paste(case_info$gdc_cases.samples.sample_type, 
                  case_info$cgc_sample_sample_type)))
case_info$gdc_cases.samples.sample_type[which(is.na(case_info$cgc_sample_sample_type))]
# One primary tumour in gdc has NA as sample type in cgc
which(is.na(case_info$cgc_file_disease_type))
which(is.na(case_info$cgc_sample_sample_type))
# Same entry has NA for cancer type and sample type in cgc; a stomach adenocarcinoma 
# primary tumour sample from a patient with one sample

c(sum(table(table(case_info$gdc_cases.submitter_id))[2:6]), 
  sum(table(table(case_info$cgc_case_id))[2:6]), 
  sum(table(table(case_info$xml_bcr_patient_barcode))[2:6]))
# 853 cases with more than one sample, one of which is NA in xml


## Remove samples with possible issues ####
which(!is.na(rse_gene$gdc_cases.samples.portions.annotations.notes))
which(!is.na(rse_gene$gdc_cases.samples.annotations.category))
which(!is.na(rse_gene$gdc_cases.samples.annotations.notes))
# Total 15 samples with possible issues. Don't need to consider 
# gdc_cases.samples.portions.annotations.notes separately as only one sample, which is 
# included in gdc_cases.samples.annotations.category and notes.
# gdc_cases.samples.annotations.category and notes go together.
paste(rse_gene$gdc_cases.samples.annotations.category, 
      rse_gene$gdc_cases.samples.annotations.notes)[
        -which(is.na(rse_gene$gdc_cases.samples.annotations.notes))]
# 1. "Barcode incorrect": Was classed as primary then reclassified as metastatic, but correct 
# now except for barcode. Keep.
# 2. "BCR Notification": Either possible normal/tumour swaps or tumour clusters with normals, 
# indicating that tumour purity is very low. Exclude.
# 3. "Tumor class but appears normal": "Genomic data resembles normal data". Exclude.
# 4. "Item does not meet study protocol": "Normal sample should not have shipped for 
# methylation, please use exome and SNP only". Not sure if relevant. Exclude to be safe.
# 5. "Item is noncanonical": Seems to just be a barcode naming anomaly. Keep.
# 6. "Pathology outside specification": Normal tissue contains tumour. Exclude.
table(rse_gene$gdc_cases.samples.annotations.category)
# 13 samples to exclude, so should be left with 11271
rse_gene <- rse_gene[, -which(
  rse_gene$gdc_cases.samples.annotations.category %in% c("BCR Notification", 
                                                         "Tumor class but appears normal", 
                                                         "Item does not meet study protocol", 
                                                         "Pathology outside specification")
)]
dim(rse_gene) # 58037 11271
info <- colData(rse_gene)


## Cases with more than one sample, primary tumour and normal only ####
table(rse_gene$gdc_cases.samples.sample_type)
# Probably would be ok to kep "Additional - New Primary" but only 11 of them so to be really safe 
# will only keep "Primary Tumor" for tumour, and "Solid Tissue Normal" and 
# "Primary Blood Derived Cancer - Peripheral Blood" for normal.
tcga_normal_tumour <- rse_gene[, which(info$gdc_cases.samples.sample_type %in% 
                                         c("Primary Blood Derived Cancer - Peripheral Blood", 
                                           "Primary Tumor", 
                                           "Solid Tissue Normal"))]
dim(tcga_normal_tumour) # 58037 10815

# Keep only cases that have multiple samples
mult_samples <- unique(tcga_normal_tumour$gdc_cases.submitter_id[
  duplicated(tcga_normal_tumour$gdc_cases.submitter_id)])
tcga_nt_mult_samples <- tcga_normal_tumour[, which(
  tcga_normal_tumour$gdc_cases.submitter_id %in% mult_samples)]
dim(tcga_nt_mult_samples) # 58037 1622
table(tcga_nt_mult_samples$gdc_cases.samples.sample_type)
# 917 primary tumour, 705 solid tissue normal, 0 peripheral blood
length(unique(tcga_nt_mult_samples$gdc_cases.submitter_id))
# 771 patients
# Must be some with multiple primary tumour samples and no normal
table(table(tcga_nt_mult_samples$gdc_cases.submitter_id))
# 718 with 2 samples, 33 with 3, 14 with 4, 5 with 5, 1 with 6
table(table(tcga_nt_mult_samples$gdc_cases.submitter_id[which(
  tcga_nt_mult_samples$gdc_cases.samples.sample_type == "Primary Tumor")]))
# 687 patients with 1 tumour sample, 34 with 2, 43 with 3, 2 with 4, 5 with 5
table(table(tcga_nt_mult_samples$gdc_cases.submitter_id[which(
  tcga_nt_mult_samples$gdc_cases.samples.sample_type == "Solid Tissue Normal")]))
# 705 patients with 1 solid tissue normal, none with >1

# Keep only cases that have a normal sample
cases_with_normal <- tcga_nt_mult_samples$gdc_cases.submitter_id[which(
  tcga_nt_mult_samples$gdc_cases.samples.sample_type == "Solid Tissue Normal")]
tcga_matched_normal <- tcga_nt_mult_samples[, which(
  tcga_nt_mult_samples$gdc_cases.submitter_id %in% cases_with_normal)]
dim(tcga_matched_normal) # 58037 1446
info_matched_normal <- colData(tcga_matched_normal)
dim(info_matched_normal) # 1446 864
table(table(tcga_matched_normal$gdc_cases.submitter_id))
# 687 with 2 samples, 3 with 3, 13 with 4, 1 with 5, 1 with 6.
# Total 72 samples.
table(table(tcga_matched_normal$gdc_cases.submitter_id[which(
  tcga_matched_normal$gdc_cases.samples.sample_type == "Primary Tumor")]))
# 687 patients with 1 tumour sample, 3 with 2, 13 with 3, 1 with 4, 1 with 5.
# Total 18 patients with multiple tumour samples, with 54 samples between them.
table(table(tcga_matched_normal$gdc_cases.submitter_id[which(
  tcga_matched_normal$gdc_cases.samples.sample_type == "Solid Tissue Normal")]))
# 705 patients with 1 solid tissue normal, none with >1

# Find and inspect cases with more than one tumour sample
mult_tumour_index <- duplicated(tcga_matched_normal$gdc_cases.submitter_id[
  which(tcga_matched_normal$gdc_cases.samples.sample_type == "Primary Tumor")])
mult_tumour <- unique(tcga_matched_normal$gdc_cases.submitter_id[
  which(tcga_matched_normal$gdc_cases.samples.sample_type == "Primary Tumor")][mult_tumour_index])
length(mult_tumour) # 18

# Previously couldn't find any sample quality-related reasons to choose between multiple samples, 
# but hadn't considered sequencing depth. Will use that now.
grep("read", names(info_matched_normal), value=T)
# read_count_as_reported_by_sra (empty)
# reads_downloaded
# proportion_of_reads_reported_by_sra_downloaded (empty)
# mapped_read_count
# avg_read_length (empty)
grep("count", names(info_matched_normal), value=T)
# Nothing new
grep("cov", names(info_matched_normal), value=T)
# Nothing
grep("depth", names(info_matched_normal), value=T)
# Nothing
tcga_matched_normal_dup <- tcga_matched_normal[
  , which(tcga_matched_normal$gdc_cases.submitter_id %in% mult_tumour)
  ]
table(table(tcga_matched_normal_dup$reads_downloaded))
table(table(tcga_matched_normal_dup$mapped_read_count))
# Both have 72 unique values
cbind(tcga_matched_normal_dup$gdc_cases.submitter_id, 
      tcga_matched_normal_dup$gdc_cases.samples.sample_type, 
      tcga_matched_normal_dup$gdc_cases.samples.is_ffpe, 
      tcga_matched_normal_dup$reads_downloaded / 1e7, 
      tcga_matched_normal_dup$mapped_read_count / 1e7)
# reads_downloaded and mapped_read_count highly correlated. Only two cases with 
# a notable difference. For one there is another tumour sample from the same 
# case with higher values for both. For the other, choosing the sample with the 
# highest mapped_read_count will select a sample with reads_downloaded very 
# close to the highest (which has low mapped_read_count and mismatch between 
# them so not a good choice).
# Expected FFPE samples to have lower counts, but they generally don't. Don't 
# think there's any need to second guess this, so just take the tumour sample 
# with the highest mapped_read_count for each case.
table(table(tcga_matched_normal$mapped_read_count))
# All unique, so can use mapped_read_count as a unique identifier for samples 
# to keep/discard.
discard <- numeric()
for (i in 1:18) {
  temp <- tcga_matched_normal_dup[
    , which(tcga_matched_normal_dup$gdc_cases.submitter_id == mult_tumour[i] & 
              tcga_matched_normal_dup$gdc_cases.samples.sample_type == "Primary Tumor")
    ]
  discard <- c(discard, temp$mapped_read_count[
    which(temp$mapped_read_count != max(temp$mapped_read_count))
    ])
}
tcga_matched_normal <- tcga_matched_normal[, -which(
  tcga_matched_normal$mapped_read_count %in% discard
)]
dim(tcga_matched_normal) # 58037 1410
info_matched_normal <- colData(tcga_matched_normal)
dim(info_matched_normal) # 1410 864

table(table(tcga_matched_normal$gdc_cases.submitter_id))
# 705 with 2 samples
table(table(tcga_matched_normal$gdc_cases.submitter_id[which(
  tcga_matched_normal$gdc_cases.samples.sample_type == "Primary Tumor")]))
# 705 patients with 1 tumour sample
table(table(tcga_matched_normal$gdc_cases.submitter_id[which(
  tcga_matched_normal$gdc_cases.samples.sample_type == "Solid Tissue Normal")]))
# 705 patients with 1 solid tissue normal


# Now have data for 705 matched tumour-normal pairs. Need to find cancers with enough 
# samples for analysis.
sort(table(info_matched_normal$gdc_cases.project.project_id[which(
  info_matched_normal$gdc_cases.samples.sample_type == "Solid Tissue Normal")]), dec=T)
project_ids <- names(sort(table(
  info_matched_normal$gdc_cases.project.project_id[
    which(info_matched_normal$gdc_cases.samples.sample_type == "Solid Tissue Normal")
    ]), dec=T))
info_matched_normal$gdc_cases.project.name[match(project_ids, info_matched_normal$gdc_cases.project.project_id)]
# BRCA 112 pairs breast invasive carcinoma
# KIRC  72 pairs kidney renal clear cell carcinoma
# THCA  59 pairs thyroid carcinoma
# LUAD  58 pairs lung adenocarcinoma
# LIHC  50 pairs liver hepatocellular carcinoma
# LUSC  49 pairs lung squamous cell carcinoma
# PRAD  46 pairs prostate adenocarcinoma
# HNSC  43 pairs head and neck squamous cell carcinoma
# COAD  42 pairs colon adenocarcinoma
# STAD  34 pairs stomach adenocarcinoma
# KIRP  32 pairs kidney renal paillary cell carcinoma
# KICH  25 pairs kidney chromophobe
# UCEC  22 pairs uterine corpus endometrial carcinoma
# BLCA  19 pairs baldder urothelial carcinoma
# ESCA  13 pairs esophogeal carcinoma
# CHOL   9 pairs cholangiocarcinoma
# READ   9 pairs rectum adenocarcinoma
# PAAD   4 pairs pancreatic adenocarcinoma
# CESC   3 pairs cervical squamous cell carcinoma and endocervical adenocarcinoma
# PCPG   3 pairs pheochromocytoma and paraganglioma
# THYM   2 pairs thymona

## Save paired dataset ####
saveRDS(tcga_matched_normal, file=here(
  "recount data/TCGA", "tcga_matched_pairs.rds"
))

