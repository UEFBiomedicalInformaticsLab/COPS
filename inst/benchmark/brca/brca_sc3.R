library(parallel)

source("load_config.R")
source("brca/brca_default_parameters.R")

## Load data
# Variables:
# tbrca_norm, tbrca_norm_deg, brca_norm_batch, brca_norm_subtypes_all, 
# dat_survival, zero_var, otp_gene_filter
source("brca/tcga_brca_mrna_data.R") 
MEs <- read.csv(paste0(path_intermediate_results, "/brca/wgcna/module_eigen_genes_deg.csv"), row.names = 1, header = TRUE)

#var_filter <- apply(log2(tbrca_norm + 1), 1, var)

# Run pipeline for all non-zero variance genes
res_sc3 <- COPS::dimred_clusteval_pipeline(list(expr = log2(tbrca_norm[combined_gene_filter,] + 1)), 
                                           batch_label = brca_norm_batch, 
                                           subtype_label = brca_norm_subtypes_all, 
                                           parallel = 1, nruns = 10, 
                                           nfolds = NFOLDS, dimred_methods = "none",
                                           cluster_methods = "SC3",
                                           distance_metric = "euclidean", 
                                           n_clusters = NCLUSTERS,
                                           survival_data = dat_survival,
                                           module_eigs = MEs,
                                           module_cor_threshold = 0.25)

scores_sc3 <- COPS::clusteval_scoring(res_sc3, wsum = (NMI.BRCA_Subtype_PAM50 + 1 - NMI.tss) / 2, summarise = SUMMARISE)
write.csv(scores_sc3$all, paste0(path_intermediate_results, "/brca/sc3/scores.csv"))
write.csv(res_sc3$clusters, gzfile(paste0(path_intermediate_results, "/brca/sc3/clusters.csv.gz")))
