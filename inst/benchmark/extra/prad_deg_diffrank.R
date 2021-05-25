library(parallel)

source("load_config.R")
source("prad/prad_default_parameters.R")

## Load data
source("prad/tcga_prad_mrna_data.R") 
MEs <- read.csv(paste0(path_intermediate_results, "/prad/wgcna/module_eigen_genes_deg.csv"), row.names = 1, header = TRUE)

prad_degs <- setdiff(rownames(tprad_norm_deg), rownames(tprad_norm_deg)[zero_var])

## Load annotations
source("prad/prad_default_annotations.R")

res_diffrank_otp <- COPS::dimred_clusteval_pipeline(list(expr = tprad_norm[prad_degs,]),
                                                pathway_enrichment_method = "DiffRank",
                                                min.size = 5, max.size = 200,
                                                gene_set_list = list_db_annots,
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
                                                module_eigs = MEs,
                                                module_cor_threshold = 0.25)

scores_diffrank_otp <- COPS::clusteval_scoring(res_diffrank_otp, wsum = (NMI.Gleason_category + 1 - NMI.tss) / 2, summarise = SUMMARISE)
write.csv(scores_diffrank_otp$all, paste0(path_intermediate_results, "/prad/diffrank_deg/scores.csv"))
write.csv(res_diffrank_otp$clusters, gzfile(paste0(path_intermediate_results, "/prad/diffrank_deg/clusters.csv.gz")))
