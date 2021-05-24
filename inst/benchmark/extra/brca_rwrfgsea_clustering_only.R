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

scores_rwr_pw <- COPS::clusteval_scoring(res_rwr_pw, wsum = (NMI.BRCA_Subtype_PAM50 + 1 - NMI.tss) / 2, summarise = SUMMARISE)
write.csv(scores_rwr_pw$all, paste0(path, score_file))
write.csv(res_rwr_pw$clusters, gzfile(paste0(path, clust_file)))

