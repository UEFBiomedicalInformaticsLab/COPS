library(parallel)

source("load_config.R")
source("brca/brca_default_parameters.R")

## Load data
# Variables:
# tbrca_norm, tbrca_norm_deg, brca_norm_batch, brca_norm_subtypes_all, 
# dat_survival, zero_var, otp_gene_filter
source("brca/tcga_brca_mrna_data.R") 

dat_list <- list(expr = as.data.frame(t(tbrca_norm)))
dat_list$expr$id <- colnames(tbrca_norm)

cv_index <- COPS::cv_fold(dat_list, nfolds = NFOLDS, nruns = NRUNS)

write.csv(cv_index$expr, paste0(path_intermediate_results, "/brca/sc3/cv_index.csv"))

if (file.exists(paste0(path_intermediate_results, "/brca/sc3/clusters.csv.gz"))) {
  #Delete file if it exists
  file.remove(paste0(path_intermediate_results, "/brca/sc3/clusters.csv.gz"))
}
if (file.exists(paste0(path_intermediate_results, "/brca/sc3/metrics.csv.gz"))) {
  #Delete file if it exists
  file.remove(paste0(path_intermediate_results, "/brca/sc3/metrics.csv.gz"))
}
if (file.exists(paste0(path_intermediate_results, "/brca/sc3/chisq_pval.csv.gz"))) {
  #Delete file if it exists
  file.remove(paste0(path_intermediate_results, "/brca/sc3/chisq_pval.csv.gz"))
}
if (file.exists(paste0(path_intermediate_results, "/brca/sc3/bassoc.csv.gz"))) {
  #Delete file if it exists
  file.remove(paste0(path_intermediate_results, "/brca/sc3/bassoc.csv.gz"))
}
if (file.exists(paste0(path_intermediate_results, "/brca/sc3/sassoc.csv.gz"))) {
  #Delete file if it exists
  file.remove(paste0(path_intermediate_results, "/brca/sc3/sassoc.csv.gz"))
}
