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

res_gsva_otp <- COPS::dimred_clusteval_pipeline(list(expr = log2(tbrca_norm[combined_gene_filter,]+1)),
                                                pathway_enrichment_method = "GSVA",
                                                min.size = 5, max.size = 200,
                                                gene_set_list = list_db_annots,
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

scores_gsva_otp <- COPS::clusteval_scoring(res_gsva_otp, wsum = (NMI.BRCA_Subtype_PAM50 + 1 - NMI.tss) / 2, summarise = SUMMARISE)
write.csv(scores_gsva_otp$all, paste0(path_intermediate_results, "/brca/gsva_otp/scores.csv"))
write.csv(res_gsva_otp$clusters, gzfile(paste0(path_intermediate_results, "/brca/gsva_otp/clusters.csv.gz")))
