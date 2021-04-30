library(parallel)

source("load_config.R")
source("prad/prad_default_parameters.R")

## Load data
source("prad/tcga_prad_mrna_data.R") 

## Load annotations
source("prad/prad_default_annotations.R")

# Create a gene co-expression network using WGCNA
WGCNA::enableWGCNAThreads(PARALLEL)

gene_correlation <- cor(t(tprad_norm[-zero_var,]), method = "spearman")

prad_gcn <- WGCNA::signumAdjacencyFunction(gene_correlation, threshold = RWR_GCN_THRESHOLD)
prad_gcn <- igraph::graph_from_adjacency_matrix(prad_gcn, mode = "undirected", weighted = NULL)

# Separate DEGs
gene_filter_up <- intersect(setdiff(rownames(tprad_norm_deg)[tprad_norm_deg$logFC > 0], 
                                    rownames(tprad_norm)[zero_var]), 
                            otp_gene_filter)
gene_filter_down <- intersect(setdiff(rownames(tprad_norm_deg)[tprad_norm_deg$logFC < 0], 
                                      rownames(tprad_norm)[zero_var]), 
                              otp_gene_filter)

# Run RWR
rwr.up <- COPS::rwr_wrapper(tprad_norm[gene_filter_up,], 
                            prad_gcn, 
                            rwr_seed_size = RWR_SEEDS, 
                            rwr_restart_probability = RWR_RESTART_PROB, 
                            parallel = PARALLEL)
rwr.down <- COPS::rwr_wrapper(tprad_norm[gene_filter_down,], 
                              prad_gcn, 
                            rwr_seed_size = RWR_SEEDS, 
                            rwr_restart_probability = RWR_RESTART_PROB,
                            parallel = PARALLEL)

write.csv(as.matrix(rwr.up), gzfile(paste0(path_intermediate_results, "/prad/rwr/gcn_rwr_otp/rwr.up.csv.gz")))
write.csv(as.matrix(rwr.down), gzfile(paste0(path_intermediate_results, "/prad/rwr/gcn_rwr_otp/rwr.down.csv.gz")))

res <- COPS::fgsea_wrapper((rwr.up - rwr.down), list_db_annots, rwr_cutoff = 0, parallel = PARALLEL)

out <- list(KEGG_RWR = t(res[,grep("^KEGG", colnames(res))]), 
            GO_RWR = t(res[,grep("^GO", colnames(res))]), 
            REACTOME_RWR = t(res[,grep("^REACTOME", colnames(res))]))

write.csv(out$KEGG, gzfile(paste0(path_intermediate_results, "/prad/rwr/gcn_rwr_otp/KEGG.csv.gz")))
write.csv(out$GO, gzfile(paste0(path_intermediate_results, "/prad/rwr/gcn_rwr_otp/GO.csv.gz")))
write.csv(out$REACTOME, gzfile(paste0(path_intermediate_results, "/prad/rwr/gcn_rwr_otp/REACTOME.csv.gz")))
