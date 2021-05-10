library(parallel)

source("load_config.R")
source("brca/brca_default_parameters.R")

## Load data
# Variables:
# tbrca_norm, tbrca_norm_deg, brca_norm_batch, brca_norm_subtypes_all, 
# dat_survival, zero_var, otp_gene_filter
source("brca/tcga_brca_mrna_data.R") 
MEs <- read.csv(paste0(path_intermediate_results, "/brca/wgcna/module_eigen_genes_deg.csv"), row.names = 1, header = TRUE)

## Load annotations
source("brca/brca_default_annotations.R")

res_diffrank_otp <- COPS::dimred_clusteval_pipeline(list(expr = tbrca_norm[combined_gene_filter,]),
                                                pathway_enrichment_method = "DiffRank",
                                                db_annots = db_annots, 
                                                min.size = 5, max.size = 200,
                                                key_name = "SYMBOL",
                                                gs_subcats = GENE_SETS, 
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

scores_diffrank_otp <- COPS::clusteval_scoring(res_diffrank_otp, wsum = (NMI.BRCA_Subtype_PAM50 + 1 - NMI.tss) / 2, summarise = SUMMARISE)
write.csv(scores_diffrank_otp$all, paste0(path_intermediate_results, "/brca/diffrank_otp/scores.csv"))
write.csv(res_diffrank_otp$clusters, gzfile(paste0(path_intermediate_results, "/brca/diffrank_otp/clusters.csv.gz")))
