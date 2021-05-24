library(parallel)

source("load_config.R")
source("brca/brca_default_parameters.R")

## Load data
# Variables:
# tbrca_norm, tbrca_norm_deg, brca_norm_batch, brca_norm_subtypes_all, 
# dat_survival, zero_var, otp_gene_filter
source("brca/tcga_brca_mrna_data.R") 
MEs <- read.csv(paste0(path_intermediate_results, "/brca/wgcna/module_eigen_genes_deg.csv"), row.names = 1, header = TRUE)

brca_degs <- setdiff(rownames(tbrca_norm_deg), rownames(tbrca_norm)[zero_var])

DIMRED_METHODS <- c("umap")

# UMAP 10 neighbours
res_umap10_otp <- COPS::dimred_clusteval_pipeline(list(expr = log2(tbrca_norm[brca_degs,] + 1)), 
                                                  batch_label = brca_norm_batch, 
                                                  subtype_label = brca_norm_subtypes_all, 
                                                  parallel = PARALLEL, nruns = NRUNS, 
                                                  nfolds = NFOLDS, dimred_methods = DIMRED_METHODS,
                                                  pca_dims = PCA_DIMS,
                                                  umap_dims = UMAP_DIMS,
                                                  umap_neighbors = 10, 
                                                  tsne_perplexities = TSNE_PERPLEXITY,
                                                  cluster_methods = DR_CLUST_METHODS,
                                                  hierarchical_linkage = DR_LINKAGES, 
                                                  gmm_modelNames = GMM_MODEL,
                                                  distance_metric = DR_METRIC, 
                                                  n_clusters = NCLUSTERS,
                                                  survival_data = dat_survival,
                                                  module_eigs = MEs,
                                                  module_cor_threshold = 0.25)

scores_umap10_otp <- COPS::clusteval_scoring(res_umap10_otp, wsum = (NMI.BRCA_Subtype_PAM50 + 1 - NMI.tss) / 2, summarise = SUMMARISE)
write.csv(scores_umap10_otp$all, paste0(path_intermediate_results, "/brca/dimred_deg/dimred_umap_10n/scores.csv"))
write.csv(res_umap10_otp$clusters, gzfile(paste0(path_intermediate_results, "/brca/dimred_deg/dimred_umap_10n/clusters.csv.gz")))

# UMAP 20 neighbours
res_umap20_otp <- COPS::dimred_clusteval_pipeline(list(expr = log2(tbrca_norm[brca_degs,] + 1)), 
                                                  batch_label = brca_norm_batch, 
                                                  subtype_label = brca_norm_subtypes_all, 
                                                  parallel = PARALLEL, nruns = NRUNS, 
                                                  nfolds = NFOLDS, dimred_methods = DIMRED_METHODS,
                                                  pca_dims = PCA_DIMS,
                                                  umap_dims = UMAP_DIMS,
                                                  umap_neighbors = 20, 
                                                  tsne_perplexities = TSNE_PERPLEXITY,
                                                  cluster_methods = DR_CLUST_METHODS,
                                                  hierarchical_linkage = DR_LINKAGES, 
                                                  gmm_modelNames = GMM_MODEL,
                                                  distance_metric = DR_METRIC, 
                                                  n_clusters = NCLUSTERS,
                                                  survival_data = dat_survival,
                                                  module_eigs = MEs,
                                                  module_cor_threshold = 0.25)

scores_umap20_otp <- COPS::clusteval_scoring(res_umap20_otp, wsum = (NMI.BRCA_Subtype_PAM50 + 1 - NMI.tss) / 2, summarise = SUMMARISE)
write.csv(scores_umap20_otp$all, paste0(path_intermediate_results, "/brca/dimred_deg/dimred_umap_20n/scores.csv"))
write.csv(res_umap20_otp$clusters, gzfile(paste0(path_intermediate_results, "/brca/dimred_deg/dimred_umap_20n/clusters.csv.gz")))

# UMAP 20 neighbours
res_umap30_otp <- COPS::dimred_clusteval_pipeline(list(expr = log2(tbrca_norm[brca_degs,] + 1)), 
                                                  batch_label = brca_norm_batch, 
                                                  subtype_label = brca_norm_subtypes_all, 
                                                  parallel = PARALLEL, nruns = NRUNS, 
                                                  nfolds = NFOLDS, dimred_methods = DIMRED_METHODS,
                                                  pca_dims = PCA_DIMS,
                                                  umap_dims = UMAP_DIMS,
                                                  umap_neighbors = 10, 
                                                  tsne_perplexities = TSNE_PERPLEXITY,
                                                  cluster_methods = DR_CLUST_METHODS,
                                                  hierarchical_linkage = DR_LINKAGES, 
                                                  gmm_modelNames = GMM_MODEL,
                                                  distance_metric = DR_METRIC, 
                                                  n_clusters = NCLUSTERS,
                                                  survival_data = dat_survival,
                                                  module_eigs = MEs,
                                                  module_cor_threshold = 0.25)

scores_umap30_otp <- COPS::clusteval_scoring(res_umap30_otp, wsum = (NMI.BRCA_Subtype_PAM50 + 1 - NMI.tss) / 2, summarise = SUMMARISE)
write.csv(scores_umap30_otp$all, paste0(path_intermediate_results, "/brca/dimred_deg/dimred_umap_30n/scores.csv"))
write.csv(res_umap30_otp$clusters, gzfile(paste0(path_intermediate_results, "/brca/dimred_deg/dimred_umap_30n/clusters.csv.gz")))


