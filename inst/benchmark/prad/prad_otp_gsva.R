library(parallel)

source("load_config.R")
source("prad/prad_default_parameters.R")

## Load data
source("prad/tcga_prad_mrna_data.R") 
MEs <- read.csv(paste0(path_intermediate_results, "/prad/wgcna/module_eigen_genes_deg.csv"), row.names = 1, header = TRUE)

## Load annotations
source("prad/prad_default_annotations.R")

res_gsva_otp <- COPS::dimred_clusteval_pipeline(list(expr = log2(tprad_norm[combined_gene_filter,]+1)),
                                                pathway_enrichment_method = "GSVA",
                                                db_annots = db_annots, 
                                                min.size = 5, max.size = 200,
                                                key_name = "SYMBOL",
                                                gs_subcats = GENE_SETS, 
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

scores_gsva_otp <- COPS::clusteval_scoring(res_gsva_otp, wsum = (NMI.Gleason_category + 1 - NMI.tss) / 2, summarise = SUMMARISE)
write.csv(scores_gsva_otp$all, paste0(path_intermediate_results, "/prad/gsva_otp/scores.csv"))
write.csv(res_gsva_otp$clusters, gzfile(paste0(path_intermediate_results, "/prad/gsva_otp/clusters.csv.gz")))
