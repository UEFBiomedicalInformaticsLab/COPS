library(parallel)

source("load_config.R")
source("prad/prad_default_parameters.R")

## Load data
source("prad/tcga_prad_mrna_data.R") 
MEs <- read.csv(paste0(path_intermediate_results, "/prad/wgcna/module_eigen_genes_deg.csv"), row.names = 1, header = TRUE)

prad_degs <- setdiff(rownames(tprad_norm_deg), rownames(tprad_norm_deg)[zero_var])

DIMRED_METHODS <- c("umap")

# UMAP 10 neighbours
res_umap10_otp <- COPS::dimred_clusteval_pipeline(list(expr = log2(tprad_norm[prad_degs,] + 1)), 
                                                  batch_label = prad_norm_batch, 
                                                  subtype_label = prad_subtype, 
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
                                                  survival_data = prad_survival,
                                                  survival_covariate_names = SURVIVAL_COVARIATES,
                                                  module_eigs = MEs,
                                                  module_cor_threshold = 0.25)

scores_umap10_otp <- COPS::clusteval_scoring(res_umap10_otp, wsum = (NMI.Gleason_category + 1 - NMI.tss) / 2, summarise = SUMMARISE)
write.csv(scores_umap10_otp$all, paste0(path_intermediate_results, "/prad/dimred_deg/dimred_umap_10n/scores.csv"))
write.csv(res_umap10_otp$clusters, gzfile(paste0(path_intermediate_results, "/prad/dimred_deg/dimred_umap_10n/clusters.csv.gz")))

# UMAP 20 neighbours
res_umap20_otp <- COPS::dimred_clusteval_pipeline(list(expr = log2(tprad_norm[prad_degs,] + 1)), 
                                                  batch_label = prad_norm_batch, 
                                                  subtype_label = prad_subtype, 
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
                                                  survival_data = prad_survival,
                                                  survival_covariate_names = SURVIVAL_COVARIATES,
                                                  module_eigs = MEs,
                                                  module_cor_threshold = 0.25)

scores_umap20_otp <- COPS::clusteval_scoring(res_umap20_otp, wsum = (NMI.Gleason_category + 1 - NMI.tss) / 2, summarise = SUMMARISE)
write.csv(scores_umap20_otp$all, paste0(path_intermediate_results, "/prad/dimred_deg/dimred_umap_20n/scores.csv"))
write.csv(res_umap20_otp$clusters, gzfile(paste0(path_intermediate_results, "/prad/dimred_deg/dimred_umap_20n/clusters.csv.gz")))

# UMAP 20 neighbours
res_umap30_otp <- COPS::dimred_clusteval_pipeline(list(expr = log2(tprad_norm[prad_degs,] + 1)), 
                                                  batch_label = prad_norm_batch, 
                                                  subtype_label = prad_subtype, 
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
                                                  survival_data = prad_survival,
                                                  survival_covariate_names = SURVIVAL_COVARIATES,
                                                  module_eigs = MEs,
                                                  module_cor_threshold = 0.25)

scores_umap30_otp <- COPS::clusteval_scoring(res_umap30_otp, wsum = (NMI.Gleason_category + 1 - NMI.tss) / 2, summarise = SUMMARISE)
write.csv(scores_umap30_otp$all, paste0(path_intermediate_results, "/prad/dimred_deg/dimred_umap_30n/scores.csv"))
write.csv(res_umap30_otp$clusters, gzfile(paste0(path_intermediate_results, "/prad/dimred_deg/dimred_umap_30n/clusters.csv.gz")))


