get_info_barcodes <- function (barcodes) {
  IDs <- strsplit(c(barcodes), "-")
  IDs <- plyr::ldply(IDs, rbind)
  colnames(IDs) <- c("project", "tss", "participant", "sample",
                     "portion", "plate", "center")
  cols <- c("project", "tss", "participant")
  IDs$patient <- apply(IDs[, cols], 1, paste, collapse = "-")
  barcode <- barcodes
  IDs <- cbind(IDs, barcode)
  condition <- gsub("11+[[:alpha:]]", "normal", as.character(IDs$sample))
  condition <- gsub("01+[[:alpha:]]", "cancer", condition)
  IDs$condition <- condition
  return(IDs)
}

# Normalized values
prad_norm <- suppressMessages(curatedTCGAData::curatedTCGAData(diseaseCode = "PRAD", assays = "RNASeq2GeneNorm", dry.run = FALSE))

# Differential expression analysis
tprad_norm <- TCGAutils::splitAssays(prad_norm, "01")
tprad_norm_controls <- TCGAutils::splitAssays(prad_norm, "11") # for differential expression analysis

tprad_norm <- Reduce("cbind", lapply(experiments(tprad_norm), function(x) assays(x)[[1]]))
tprad_norm_controls <- Reduce("cbind", lapply(experiments(tprad_norm_controls), function(x) assays(x)[[1]]))
tprad_norm_controls_id <- match(substr(colnames(tprad_norm_controls), 1, 12), substr(colnames(tprad_norm), 1, 12))

tprad_norm_tumors <- tprad_norm[,tprad_norm_controls_id]

tprad_norm_deg <- TCGAbiolinks::TCGAanalyze_DEA(mat1 = tprad_norm_controls, 
                                                mat2 = tprad_norm_tumors, 
                                                Cond1type = "Normal", 
                                                Cond2type = "Tumor",
                                                fdr.cut = 0.01,
                                                logFC.cut = 1,
                                                method = "glmLRT")

# Batch variables
prad_norm_batch <- get_info_barcodes(colnames(tprad_norm))[c("tss", "plate")]
rownames(prad_norm_batch) <- colnames(tprad_norm)

IDs <- get_info_barcodes(colnames(tprad_norm))

# Identify near-zero variance genes
zero_var <- caret::nearZeroVar(t(log2(tprad_norm + 1)))

# Get subtypes
prad_subtype_all <- TCGAbiolinks::TCGAquery_subtype(tumor = "PRAD")
prad_subtype_all$id <- colnames(tprad_norm)[match(prad_subtype_all$sample, substr(colnames(tprad_norm), 1, 15))]

prad_clinical <- prad_subtype_all[c("id", "Race", "Age", "PSA_preop", "Tumor_cellularity_pathology")]

prad_subtype <- prad_norm@colData[c("Subtype")]

# Format Gleason scores into categories
prad_subtype$Gleason_score_primary <- prad_norm@colData$patient.stage_event.gleason_grading.primary_pattern
prad_subtype$Gleason_score_secondary <- prad_norm@colData$patient.stage_event.gleason_grading.secondary_pattern
prad_subtype$Gleason_score <- paste(prad_subtype$Gleason_score_primary, prad_subtype$Gleason_score_secondary, sep = "+")
rownames(prad_subtype) <- colnames(tprad_norm)[match(rownames(prad_subtype), substr(colnames(tprad_norm), 1, 12))]
prad_subtype$Gleason_sum <- prad_subtype$Gleason_score_primary + prad_subtype$Gleason_score_secondary
prad_subtype$Gleason_category <- NA
prad_subtype$Gleason_category[prad_subtype$Gleason_sum < 7] <- "<7"
prad_subtype$Gleason_category[prad_subtype$Gleason_sum > 7] <- ">7"
prad_subtype$Gleason_category[prad_subtype$Gleason_sum == 7 & prad_subtype$Gleason_score_primary == "3"] <- "3+4"
prad_subtype$Gleason_category[prad_subtype$Gleason_sum == 7 & prad_subtype$Gleason_score_primary == "4"] <- "4+3"

# OTP genes
disease_id <- c("EFO_0000673", "EFO_0001663")
assoc_score_fields <- paste(paste("&fields=", c('disease.efo_info.label',
                                                'disease.efo_info.therapeutic_area',
                                                'target.gene_info.symbol',
                                                'association_score.overall',
                                                'disease.id',
                                                'association_score.datatypes'), sep=''), collapse = "")
disease_otp_list <- COPS::retrieveDiseaseGenesOT(c(disease_id), assoc_score_fields)
disease_otp <- plyr::rbind.fill(disease_otp_list)[,-c(10:13)]

otp_gene_filter <- intersect(rownames(tprad_norm), 
                             setdiff(disease_otp$target.gene_info.symbol[disease_otp$association_score.overall > OTP_CUTOFF], 
                                     rownames(tprad_norm)[zero_var]))

combined_gene_filter <- intersect(setdiff(rownames(tprad_norm_deg), rownames(tprad_norm)[zero_var]), otp_gene_filter)

# Survival data
as.data.frame(prad_norm@colData[c("days_to_psa", "patient.days_to_first_biochemical_recurrence", "patient.clinical_cqcf.days_to_preop_psa")])

temp <- prad_norm@colData[c("patientID", "patient.days_to_first_biochemical_recurrence", "days_to_last_followup", 
                            "patient.biochemical_recurrence", "patient.clinical_cqcf.psa_result_preop", "years_to_birth", 
                            "vital_status", "race", "pathology_N_stage", "pathology_T_stage")]
temp$Gleason_category <- prad_subtype$Gleason_category
temp$Subtype <- prad_subtype$Subtype

temp <- temp[!is.na(temp$patient.biochemical_recurrence), ]
# Remove rows with mismatched recurrence information
temp <- temp[!(temp$patient.biochemical_recurrence == "yes" & is.na(temp$patient.days_to_first_biochemical_recurrence)), ]
temp <- temp[!(temp$patient.biochemical_recurrence == "no" & !is.na(temp$patient.days_to_first_biochemical_recurrence)), ]
# Remove rows with deaths
temp <- temp[temp$vital_status != 1, ]
# Remove rows with missing psa, age
temp <- temp[!is.na(temp$patient.clinical_cqcf.psa_result_preop), ]
temp <- temp[!is.na(temp$years_to_birth), ]

prad_survival <- COPS:::survival_preprocess(temp, 
                                            event_time_name = "patient.days_to_first_biochemical_recurrence",
                                            follow_up_time_name = "days_to_last_followup",
                                            event_field_name = "patient.biochemical_recurrence",
                                            event_name = "yes", 
                                            event_time_cutoff = 2000)#Inf)
prad_survival$ID <- colnames(tprad_norm)[match(prad_survival$patientID, substr(colnames(tprad_norm), 1, 12))]


