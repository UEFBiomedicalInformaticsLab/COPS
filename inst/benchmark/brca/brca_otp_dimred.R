library(parallel)

source("load_config.R")
source("brca/brca_default_parameters.R")

## Load data
# Variables:
# tbrca_norm, tbrca_norm_deg, brca_norm_batch, brca_norm_subtypes_all, 
# dat_survival, zero_var, otp_gene_filter
source("brca/tcga_brca_mrna_data.R") 
MEs <- read.csv(paste0(path_intermediate_results, "/brca/wgcna/module_eigen_genes_deg.csv"), row.names = 1, header = TRUE)

DIMRED_METHODS <- c("pca", "tsne")

# UMAP 10 neighbours
res_dimred_otp <- COPS::dimred_clusteval_pipeline(list(expr = log2(tbrca_norm[combined_gene_filter,] + 1)), 
                                                  batch_label = brca_norm_batch, 
                                                  subtype_label = brca_norm_subtypes_all, 
                                                  parallel = PARALLEL, nruns = NRUNS, 
                                                  nfolds = NFOLDS, dimred_methods = DIMRED_METHODS,
                                                  pca_dims = PCA_DIMS,
                                                  umap_dims = UMAP_DIMS,
                                                  tsne_perplexities = TSNE_PERPLEXITY,
                                                  cluster_methods = DR_CLUST_METHODS,
                                                  hierarchical_linkage = DR_LINKAGES, 
                                                  gmm_modelNames = GMM_MODEL,
                                                  distance_metric = DR_METRIC, 
                                                  n_clusters = NCLUSTERS,
                                                  survival_data = dat_survival,
                                                  module_eigs = MEs,
                                                  module_cor_threshold = 0.25)

scores_dimred_otp <- COPS::clusteval_scoring(res_dimred_otp, wsum = (NMI.BRCA_Subtype_PAM50 + 1 - NMI.tss) / 2, summarise = SUMMARISE)
write.csv(scores_dimred_otp$all, paste0(path_intermediate_results, "/brca/dimred_otp/scores.csv"))
write.csv(res_dimred_otp$clusters, gzfile(paste0(path_intermediate_results, "/brca/dimred_otp/clusters.csv.gz")))