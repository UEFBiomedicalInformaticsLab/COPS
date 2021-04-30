library(parallel)

source("load_config.R")
source("prad/prad_default_parameters.R")

## Load data
source("prad/tcga_prad_mrna_data.R") 
MEs <- read.csv(paste0(path_intermediate_results, "/prad/wgcna/module_eigen_genes_deg.csv"), row.names = 1, header = TRUE)

args <- commandArgs(trailingOnly = TRUE)

path <- args[1]

score_file <- args[2]
clust_file <- args[3]

dat_list <- list()
for (i in 4:length(args)) {
  dat_list[[strsplit(args[[i]], ".", fixed = TRUE)[[1]][1]]] <- read.csv(paste0(path, args[i]), header = TRUE, row.names = 1)
}

names(dat_list) <- paste0(names(dat_list), "_RWR")

for (i in 1:length(dat_list)) {
  colnames(dat_list[[i]]) <- gsub(".", "-", colnames(dat_list[[i]]), fixed = TRUE)
  dat_list[[i]] <- dat_list[[i]][apply(dat_list[[i]], 1, var) > 0,]
}

res_rwr_pw <- COPS::dimred_clusteval_pipeline(dat_list,
                                              batch_label = prad_norm_batch, 
                                              subtype_label = prad_subtype, 
                                              parallel = PARALLEL, nruns = NRUNS, 
                                              nfolds = NFOLDS, dimred_methods = c("none"),
                                              cluster_methods = OG_CLUST_METHODS,
                                              hierarchical_linkage = OG_LINKAGES, 
                                              distance_metric = OG_METRIC, 
                                              correlation_method = COR_METHOD, 
                                              n_clusters = NCLUSTERS,
                                              survival_data = prad_survival,
                                              survival_covariate_names = SURVIVAL_COVARIATES,
                                              module_eigs = MEs)

scores_rwr_pw <- COPS::clusteval_scoring(res_rwr_pw, wsum = (NMI.Gleason_category + 1 - NMI.tss) / 2, summarise = SUMMARISE)
write.csv(scores_rwr_pw$all, paste0(path, score_file))
write.csv(res_rwr_pw$clusters, gzfile(paste0(path, clust_file)))

