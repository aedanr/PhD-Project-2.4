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
dim(rse_gene)
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
dim(tcga_normal_tumour)
# 58037 10815

# Keep only cases that have multiple samples
mult_samples <- unique(tcga_normal_tumour$gdc_cases.submitter_id[
  duplicated(tcga_normal_tumour$gdc_cases.submitter_id)])
tcga_nt_mult_samples <- tcga_normal_tumour[, which(
  tcga_normal_tumour$gdc_cases.submitter_id %in% mult_samples)]
dim(tcga_nt_mult_samples)
# 58037 1622
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
dim(tcga_matched_normal)
# 58037 1446
info_matched_normal <- colData(tcga_matched_normal)
dim(info_matched_normal)
# 1446 864
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
length(mult_tumour)
# 18

# Look for reasons to choose one tumour sample out of multiple options
# Date matching between tumour and normal?
grep("date", names(info_matched_normal), value=T)
table(table(tcga_matched_normal$gdc_metadata_files.created_datetime.analysis[which(
  tcga_matched_normal$gdc_cases.submitter_id %in% mult_tumour)]))
# gdc_metadata_files.created_datetime.analysis
# gdc_cases.updated_datetime
# gdc_cases.demographic.updated_datetime
# gdc_cases.samples.portions.creation_datetime
cbind(tcga_matched_normal$gdc_cases.submitter_id, 
      tcga_matched_normal$gdc_cases.samples.sample_type, 
      tcga_matched_normal$gdc_metadata_files.created_datetime.analysis)[
        which(tcga_matched_normal$gdc_cases.submitter_id %in% mult_tumour), ]
# Doesn't seem to be any relation betweengdc_metadata_files.created_datetime.analysis and N/T samples
cbind(tcga_matched_normal$gdc_cases.submitter_id, 
      tcga_matched_normal$gdc_cases.samples.sample_type, 
      tcga_matched_normal$gdc_cases.updated_datetime)[
        which(tcga_matched_normal$gdc_cases.submitter_id %in% mult_tumour), ]
# gdc_cases.updated_datetime is the same for all samples from same patient
cbind(tcga_matched_normal$gdc_cases.submitter_id, 
      tcga_matched_normal$gdc_cases.samples.sample_type, 
      tcga_matched_normal$gdc_cases.demographic.updated_datetime)[
        which(tcga_matched_normal$gdc_cases.submitter_id %in% mult_tumour), ]
# gdc_cases.demographic.updated_datetime is the same for all samples from same patient
cbind(tcga_matched_normal$gdc_cases.submitter_id, 
      tcga_matched_normal$gdc_cases.samples.sample_type, 
      tcga_matched_normal$gdc_cases.samples.portions.creation_datetime)[
        which(tcga_matched_normal$gdc_cases.submitter_id %in% mult_tumour), ]
# Doesn't seem to be any relation gdc_cases.samples.portions.creation_datetime and N/T samples.
# Any other fields that could potentially be used to logically choose one tumour sample?
fields <- character(0)
for (i in names(info_matched_normal)) {
  if (class(info_matched_normal[,i]) != "list") {
    if (
      length(table(info_matched_normal[which(
      tcga_matched_normal$gdc_cases.submitter_id %in% mult_tumour),i])) > 1 & 
      length(table(info_matched_normal[which(
        tcga_matched_normal$gdc_cases.submitter_id %in% mult_tumour),i])) < 72
      ) {
      fields <- c(fields, i)
    }
  }
}
fields
fields <- c(
  "gdc_cases.project.primary_site", "gdc_cases.tissue_source_site.bcr_id", 
  "gdc_cases.tissue_source_site.name", "gdc_cases.diagnoses.primary_diagnosis", 
  "gdc_cases.diagnoses.tumor_stage", "gdc_cases.diagnoses.diagnosis_id", 
  "gdc_cases.diagnoses.tissue_or_organ_of_origin", "gdc_cases.samples.intermediate_dimension", 
  "gdc_cases.samples.sample_id", "gdc_cases.samples.is_ffpe", 
  "gdc_cases.samples.pathology_report_uuid", "gdc_cases.samples.days_to_collection", 
  "gdc_cases.samples.initial_weight", "gdc_cases.samples.portions.portion_id", 
  "gdc_cases.samples.portions.analytes.a260_a280_ratio", "cgc_file_disease_type", 
  "cgc_sample_is_ffpe", "cgc_sample_tissue_source_site", "cgc_sample_days_to_collection", 
  "cgc_sample_portion", "cgc_sample_id", "cgc_case_new_tumor_event", "cgc_case_tumor_status", 
  "cgc_case_primary_site", "cgc_case_histological_diagnosis", "cgc_case_pathologic_stage", 
  "cgc_case_batch_number", "cgc_portion_portion_number", "cgc_portion_weight", 
  "cgc_portion_analyte", "cgc_portion_file", "cgc_portion_id", "cgc_slide_percent_tumor_nuclei", 
  "cgc_slide_percent_normal_cells", "cgc_slide_percent_stromal_cells", "cgc_slide_section_location", 
  "cgc_slide_id", "cgc_slide_percent_tumor_cells", "xml_bcr_patient_barcode", 
  "xml_tumor_tissue_site", "xml_residual_tumor", "xml_stage_event_pathologic_stage", 
  "xml_stage_event_gleason_grading", "xml_diagnosis", "xml_diagnosis_subtype", "xml_tumor_levels"
)
cbind(tcga_matched_normal$gdc_cases.submitter_id, 
      tcga_matched_normal$gdc_cases.samples.sample_type, 
      info_matched_normal[, fields[1]])[
        which(tcga_matched_normal$gdc_cases.submitter_id %in% mult_tumour), ]
fields <- fields[c(8,9,10,14,15,17,20,21,28,29,33,34,35,38)]
cbind(tcga_matched_normal$gdc_cases.submitter_id, 
      tcga_matched_normal$gdc_cases.samples.sample_type, 
      info_matched_normal[, fields[1]])[
        which(tcga_matched_normal$gdc_cases.submitter_id %in% mult_tumour), ]; fields[1]
# 3 gdc_cases.samples.is_ffpe
# 5 gdc_cases.samples.portions.analytes.a260_a280_ratio
# 9 cgc_portion_portion_number
# 11 cgc_slide_percent_tumor_nuclei
# 12 cgc_slide_percent_normal_cells
# 13 cgc_slide_percent_stromal_cells
# 14 cgc_slide_percent_tumor_cells
fields <- fields[c(3,5,9,11:14)]
dup.info <- cbind(tcga_matched_normal$gdc_cases.submitter_id, 
                  tcga_matched_normal$gdc_cases.samples.sample_type, 
                  info_matched_normal[, fields])[
                    which(tcga_matched_normal$gdc_cases.submitter_id %in% mult_tumour), ]
names(dup.info) <- c("ID", "Sample type", "FFPE", "A260/280", "Portion #", "% tumour nuclei", 
                     "Slide % normal cells", "Slide % stromal cells", "Slide % tumour cells")
options(showHeadLines=50, showTailLines=50)
dup.info[which(dup.info$ID == mult_tumour[1]), ]
# All have at least one non-FFPE tumour sample, so start by removing FFPE samples
tcga_matched_normal <- tcga_matched_normal[, -which(
  tcga_matched_normal$gdc_cases.submitter_id %in% mult_tumour & 
    tcga_matched_normal$gdc_cases.samples.is_ffpe == TRUE)]
dim(tcga_matched_normal)
# 58037 1426
info_matched_normal <- colData(tcga_matched_normal)
dim(info_matched_normal)
# 1426 864

table(table(tcga_matched_normal$gdc_cases.submitter_id))
# 690 with 2 samples, 14 with 3, 1 with 4.
# Total 46 samples.
table(table(tcga_matched_normal$gdc_cases.submitter_id[which(
  tcga_matched_normal$gdc_cases.samples.sample_type == "Primary Tumor")]))
# 690 patients with 1 tumour sample, 14 with 2, 1 with 3.
# Total 15 patients with multiple tumour samples, with 31 samples between them.
table(table(tcga_matched_normal$gdc_cases.submitter_id[which(
  tcga_matched_normal$gdc_cases.samples.sample_type == "Solid Tissue Normal")]))
# 705 patients with 1 solid tissue normal, none with >1
mult_tumour_index <- duplicated(tcga_matched_normal$gdc_cases.submitter_id[
  which(tcga_matched_normal$gdc_cases.samples.sample_type == "Primary Tumor")])
mult_tumour <- unique(tcga_matched_normal$gdc_cases.submitter_id[
  which(tcga_matched_normal$gdc_cases.samples.sample_type == "Primary Tumor")][mult_tumour_index])
length(mult_tumour)
# 15
fields <- character(0)
for (i in names(info_matched_normal)) {
  if (class(info_matched_normal[,i]) != "list") {
    if (
      length(table(info_matched_normal[which(
        tcga_matched_normal$gdc_cases.submitter_id %in% mult_tumour & 
        tcga_matched_normal$gdc_cases.samples.sample_type == "Primary Tumor"),i])) > 1 & 
      length(table(info_matched_normal[which(
        tcga_matched_normal$gdc_cases.submitter_id %in% mult_tumour & 
        tcga_matched_normal$gdc_cases.samples.sample_type == "Primary Tumor"),i])) < 31
    ) {
      fields <- c(fields, i)
    }
  }
}
length(fields)
fields <- c(
  "gdc_metadata_files.created_datetime.analysis", "gdc_metadata_files.file_size.analysis", 
  "gdc_metadata_files.file_size.experiment", "gdc_metadata_files.file_size.run", 
  "gdc_cases.updated_datetime", "gdc_cases.tissue_source_site.bcr_id", 
  "gdc_cases.samples.sample_id", "gdc_cases.samples.pathology_report_uuid", 
  "gdc_cases.samples.shortest_dimension", "gdc_cases.samples.initial_weight", 
  "gdc_cases.samples.longest_dimension", "gdc_cases.samples.portions.creation_datetime", 
  "gdc_cases.samples.portions.weight", "gdc_cases.samples.portions.analytes.aliquots.aliquot_id", 
  "gdc_cases.samples.portions.analytes.aliquots.concentration", "cgc_sample_shortest_dimension", 
  "cgc_sample_intermediate_dimension", "cgc_sample_initial_weight", 
  "cgc_sample_id", "cgc_sample_longest_dimension", 
  "cgc_case_sample", "cgc_portion_weight", 
  "cgc_portion_slide", "cgc_portion_analyte", 
  "cgc_slide_percent_normal_cells", "cgc_slide_id", 
  "cgc_slide_percent_tumor_cells", "xml_bcr_patient_barcode", 
  "xml_residual_tumor"
)
cbind(tcga_matched_normal$gdc_cases.submitter_id, 
      tcga_matched_normal$gdc_cases.samples.sample_type, 
      info_matched_normal[, fields[1]])[
        which(tcga_matched_normal$gdc_cases.submitter_id %in% mult_tumour), ]
# gdc_cases.samples.sample_id shows that all remaining tumour samples for each 
# patient come from the same sample, as does gdc_cases.samples.pathology_report_uuid 
# for pathology reports. Same for weight, dimensions, cgc_sample_id, cgc_portion_slide, 
# cgc_portion_analyte.
# Don't seem to be any fields that can be used to logically choose which tumour samples 
# to keep, so will just keep first instance for each. Need a unique instance identifier 
# to use to get indices of instances to remove.
table(table(info_matched_normal$bigwig_file))
# 1426 different entries
bigwig_to_delete <- tcga_matched_normal$bigwig_file[
  which(tcga_matched_normal$gdc_cases.samples.sample_type == "Primary Tumor")][mult_tumour_index]
indices_to_delete <- which(tcga_matched_normal$bigwig_file %in% bigwig_to_delete)
identical(unique(tcga_matched_normal$gdc_cases.submitter_id[indices_to_delete]), mult_tumour)
# TRUE
table(tcga_matched_normal$gdc_cases.samples.sample_type[indices_to_delete])
# All "Primary Tumor"
tcga_matched_normal <- tcga_matched_normal[, -indices_to_delete]
dim(tcga_matched_normal)
# 58037 1410
info_matched_normal <- colData(tcga_matched_normal)
dim(info_matched_normal)
# 1410 864

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

## Import data and break down by project and specific diagnoses ####
library(here)
library(recount)
data <- readRDS(here("recount data/TCGA", "tcga_matched_pairs.rds"))

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






