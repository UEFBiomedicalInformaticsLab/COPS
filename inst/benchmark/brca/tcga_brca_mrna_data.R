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
brca_norm <- suppressMessages(curatedTCGAData::curatedTCGAData(diseaseCode = "BRCA", assays = "RNASeq2GeneNorm", dry.run = FALSE))

# Differential expression analysis
tbrca_norm <- TCGAutils::splitAssays(brca_norm, "01") # tumours
tbrca_norm_controls <- TCGAutils::splitAssays(brca_norm, "11") # healthy controls

tbrca_norm <- Reduce("cbind", lapply(experiments(tbrca_norm), function(x) assays(x)[[1]]))
tbrca_norm_controls <- Reduce("cbind", lapply(experiments(tbrca_norm_controls), function(x) assays(x)[[1]]))
tbrca_norm_controls_id <- match(substr(colnames(tbrca_norm_controls), 1, 12), substr(colnames(tbrca_norm), 1, 12))

tbrca_norm_tumors <- tbrca_norm[,tbrca_norm_controls_id]

tbrca_norm_deg <- TCGAbiolinks::TCGAanalyze_DEA(mat1 = tbrca_norm_controls, 
                                                mat2 = tbrca_norm_tumors, 
                                                Cond1type = "Normal", 
                                                Cond2type = "Tumor",
                                                fdr.cut = 0.01,
                                                logFC.cut = 1,
                                                method = "glmLRT")

# Get subtypes
brca_subtypes <- TCGAbiolinks::TCGAquery_subtype(tumor = "BRCA")

# Filter out Normal-like subtype and samples missing subtype info
subtyped_tumor_ind <- match(brca_subtypes$patient[!brca_subtypes$BRCA_Subtype_PAM50 %in% c("Normal", "NA") & !is.na(brca_subtypes$BRCA_Subtype_PAM50)],
                            substr(colnames(tbrca_norm), 1, 12))
subtyped_tumor_ind <- subtyped_tumor_ind[!is.na(subtyped_tumor_ind)]
tbrca_norm <- tbrca_norm[,subtyped_tumor_ind]

# Batch variables
brca_norm_batch <- get_info_barcodes(colnames(tbrca_norm))[c("tss", "plate")]
rownames(brca_norm_batch) <- colnames(tbrca_norm)

IDs <- get_info_barcodes(colnames(tbrca_norm))

# Identify near-zero variance genes
zero_var <- caret::nearZeroVar(t(log2(tbrca_norm + 1)))

# Create subtype matrix that matches data rows
brca_norm_subtypes_all <- brca_subtypes[match(substr(colnames(tbrca_norm), 1, 12), 
                                              brca_subtypes$patient), 
                                        c("pathologic_stage", "BRCA_Pathology", "BRCA_Subtype_PAM50")]
rownames(brca_norm_subtypes_all) <- colnames(tbrca_norm)

# OTP genes
disease_id <- "EFO_0000305" # BRCA
assoc_score_fields <- paste(paste("&fields=", c('disease.efo_info.label',
                                               'disease.efo_info.therapeutic_area',
                                               'target.gene_info.symbol',
                                               'association_score.overall',
                                               'disease.id',
                                               'association_score.datatypes'), sep=''), collapse = "")
disease_otp <- COPS::retrieveDiseaseGenesOT(c(disease_id), assoc_score_fields)[[1]][,-c(10:13)]
otp_gene_filter <- intersect(rownames(tbrca_norm), 
                             setdiff(disease_otp$target.gene_info.symbol[disease_otp$association_score.datatypes.genetic_association > OTP_CUTOFF], 
                                     rownames(tbrca_norm)[zero_var]))

combined_gene_filter <- intersect(setdiff(rownames(tbrca_norm_deg), rownames(tbrca_norm)[zero_var]), otp_gene_filter)

# Make data frame containing data for survivalanalysis
dat_survival <- COPS:::survival_preprocess(brca_subtypes, 
                                           event_time_name = "days_to_death",
                                           follow_up_time_name = "days_to_last_followup",
                                           event_field_name = "vital_status",
                                           event_name = "Dead", 
                                           event_time_cutoff = 3000)
dat_survival$age <- dat_survival$age_at_initial_pathologic_diagnosis
dat_survival$stage <- dat_survival$pathologic_stage
dat_survival$ID <- colnames(tbrca_norm)[match(dat_survival$patient, substr(colnames(tbrca_norm), 1, 12))]
