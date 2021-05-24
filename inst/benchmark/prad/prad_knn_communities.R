library(parallel)

source("load_config.R")
source("prad/prad_default_parameters.R")

## Load data
source("prad/tcga_prad_mrna_data.R") 
MEs <- read.csv(paste0(path_intermediate_results, "/prad/wgcna/module_eigen_genes_deg.csv"), row.names = 1, header = TRUE)

# Run pipeline for all non-zero variance genes
res_knn_com <- COPS::dimred_clusteval_pipeline(list(expr = log2(tprad_norm[combined_gene_filter,] + 1)), 
                                           batch_label = prad_norm_batch, 
                                           subtype_label = prad_subtype, 
                                           parallel = PARALLEL, nruns = NRUNS, 
                                           nfolds = NFOLDS, dimred_methods = "none",
                                           cluster_methods = "knn_communities",
                                           distance_metric = "euclidean", 
                                           n_clusters = NCLUSTERS,
                                           survival_data = prad_survival,
                                           survival_covariate_names = SURVIVAL_COVARIATES,
                                           module_eigs = MEs,
                                           module_cor_threshold = 0.25)

scores_knn_com <- COPS::clusteval_scoring(res_knn_com, wsum = (NMI.Gleason_category + 1 - NMI.tss) / 2, summarise = SUMMARISE)
write.csv(scores_knn_com$all, paste0(path_intermediate_results, "/prad/knn_communities/scores.csv"))
write.csv(res_knn_com$clusters, gzfile(paste0(path_intermediate_results, "/prad/knn_communities/clusters.csv.gz")))
