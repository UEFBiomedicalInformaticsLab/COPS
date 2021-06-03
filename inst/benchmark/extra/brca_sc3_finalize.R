library(parallel)

source("load_config.R")
source("brca/brca_default_parameters.R")

## Load data
# Variables:
# tbrca_norm, tbrca_norm_deg, brca_norm_batch, brca_norm_subtypes_all, 
# dat_survival, zero_var, otp_gene_filter
source("brca/tcga_brca_mrna_data.R") 
MEs <- read.csv(paste0(path_intermediate_results, "/brca/wgcna/module_eigen_genes_deg.csv"), row.names = 1, header = TRUE)

args <- commandArgs(trailingOnly = TRUE)

path <- args[1]



res_rwr_pw <- COPS::dimred_clusteval_pipeline(dat_list,
                                              batch_label = brca_norm_batch, 
                                              subtype_label = brca_norm_subtypes_all, 
                                              parallel = PARALLEL, nruns = NRUNS, 
                                              nfolds = NFOLDS, dimred_methods = c("none"),
                                              cluster_methods = OG_CLUST_METHODS,
                                              hierarchical_linkage = OG_LINKAGES, 
                                              distance_metric = OG_METRIC, 
                                              correlation_method = COR_METHOD, 
                                              n_clusters = NCLUSTERS,
                                              survival_data = dat_survival,
                                              module_eigs = MEs,
                                              module_cor_threshold = 0.25)

clustering_evaluation <- function(dat,
                                  batch_label_names = NULL,
                                  subtype_label_names = NULL,
                                  n_clusters = 2:5,
                                  cluster_methods = c("hierarchical","diana","kmeans"),
                                  distance_metric = "euclidean",
                                  correlation_method = "spearman",
                                  hierarchical_linkage = "complete",
                                  kmeans_num_init = 100,
                                  kmeans_max_iters = 100,
                                  kmeans_tol = 0.0001,
                                  gmm_modelNames = NULL,  
                                  gmm_shrinkage = 0.01, 
                                  knn_neighbours = 30, 
                                  knn_jaccard = TRUE, 
                                  ...) {}

clusters <- read.csv(paste0(path, "cluster.csv.gz"), row.names = 1, header = TRUE)


# Stability

# Survival

# Modules







scores_rwr_pw <- COPS::clusteval_scoring(res_rwr_pw, wsum = (NMI.BRCA_Subtype_PAM50 + 1 - NMI.tss) / 2, summarise = SUMMARISE)
write.csv(scores_rwr_pw$all, paste0(path, score_file))
write.csv(res_rwr_pw$clusters, gzfile(paste0(path, clust_file)))



