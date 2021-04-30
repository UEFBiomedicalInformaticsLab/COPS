library(parallel)

source("load_config.R")
source("brca/brca_default_parameters.R")

## Load data
# Variables:
# tbrca_norm, tbrca_norm_deg, brca_norm_batch, brca_norm_subtypes_all, 
# dat_survival, zero_var, otp_gene_filter
source("brca/tcga_brca_mrna_data.R") 

## Load annotations
source("brca/brca_default_annotations.R")

# Create a gene co-expression network using WGCNA
WGCNA::enableWGCNAThreads(PARALLEL)

gene_correlation <- cor(t(tbrca_norm[-zero_var,]), method = "spearman")

brca_gcn <- WGCNA::signumAdjacencyFunction(gene_correlation, threshold = RWR_GCN_THRESHOLD)
brca_gcn <- igraph::graph_from_adjacency_matrix(brca_gcn, mode = "undirected", weighted = NULL)

# Separate DEGs
gene_filter_up <- intersect(setdiff(rownames(tbrca_norm_deg)[tbrca_norm_deg$logFC > 0], 
                                    rownames(tbrca_norm)[zero_var]), 
                            otp_gene_filter)
gene_filter_down <- intersect(setdiff(rownames(tbrca_norm_deg)[tbrca_norm_deg$logFC < 0], 
                                      rownames(tbrca_norm)[zero_var]), 
                              otp_gene_filter)

# Run RWR
rwr.up <- COPS::rwr_wrapper(tbrca_norm[gene_filter_up,], 
                            brca_gcn, 
                            rwr_seed_size = RWR_SEEDS, 
                            rwr_restart_probability = RWR_RESTART_PROB, 
                            parallel = PARALLEL)
rwr.down <- COPS::rwr_wrapper(tbrca_norm[gene_filter_down,], 
                              brca_gcn, 
                            rwr_seed_size = RWR_SEEDS, 
                            rwr_restart_probability = RWR_RESTART_PROB,
                            parallel = PARALLEL)

write.csv(as.matrix(rwr.up), gzfile(paste0(path_intermediate_results, "/brca/rwr/gcn_rwr_otp/rwr.up.csv.gz")))
write.csv(as.matrix(rwr.down), gzfile(paste0(path_intermediate_results, "/brca/rwr/gcn_rwr_otp/rwr.down.csv.gz")))

res <- COPS::fgsea_wrapper((rwr.up - rwr.down), list_db_annots, rwr_cutoff = 0, parallel = PARALLEL)

out <- list(KEGG_RWR = t(res[,grep("^KEGG", colnames(res))]), 
            GO_RWR = t(res[,grep("^GO", colnames(res))]), 
            REACTOME_RWR = t(res[,grep("^REACTOME", colnames(res))]))

write.csv(out$KEGG, gzfile(paste0(path_intermediate_results, "/brca/rwr/gcn_rwr_otp/KEGG.csv.gz")))
write.csv(out$GO, gzfile(paste0(path_intermediate_results, "/brca/rwr/gcn_rwr_otp/GO.csv.gz")))
write.csv(out$REACTOME, gzfile(paste0(path_intermediate_results, "/brca/rwr/gcn_rwr_otp/REACTOME.csv.gz")))
