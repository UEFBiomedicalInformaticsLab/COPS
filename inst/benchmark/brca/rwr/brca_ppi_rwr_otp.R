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

# Load PPI (https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-019-09186-x/MediaObjects/41467_2019_9186_MOESM4_ESM.xlsx)
PPI <- read.table("~/tcga/PPI_entrez_ensemble.txt", header = TRUE)
PPI <- igraph::graph_from_data_frame(PPI[,1:2], directed = FALSE)

PPI_symbols <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, igraph::V(PPI)$name, "SYMBOL", "ENTREZID")
ambiguous_ind <- PPI_symbols %in% PPI_symbols[duplicated(PPI_symbols)]
PPI_symbols[ambiguous_ind] <- igraph::V(PPI)$name[ambiguous_ind]
igraph::V(PPI)$name <- PPI_symbols

# Separate DEGs
gene_filter_up <- intersect(setdiff(rownames(tbrca_norm_deg)[tbrca_norm_deg$logFC > 0], 
                                    rownames(tbrca_norm)[zero_var]), 
                            otp_gene_filter)
gene_filter_down <- intersect(setdiff(rownames(tbrca_norm_deg)[tbrca_norm_deg$logFC < 0], 
                                      rownames(tbrca_norm)[zero_var]), 
                              otp_gene_filter)

# Run RWR
rwr.up <- COPS::rwr_wrapper(tbrca_norm[gene_filter_up,], 
                            PPI, 
                            rwr_seed_size = RWR_SEEDS, 
                            rwr_restart_probability = RWR_RESTART_PROB, 
                            parallel = PARALLEL)
rwr.down <- COPS::rwr_wrapper(tbrca_norm[gene_filter_down,], 
                              PPI, 
                              rwr_seed_size = RWR_SEEDS, 
                              rwr_restart_probability = RWR_RESTART_PROB,
                              parallel = PARALLEL)

write.csv(as.matrix(rwr.up), gzfile(paste0(path_intermediate_results, "/brca/rwr/ppi_rwr_otp/rwr.up.csv.gz")))
write.csv(as.matrix(rwr.down), gzfile(paste0(path_intermediate_results, "/brca/rwr/ppi_rwr_otp/rwr.down.csv.gz")))

res <- COPS::fgsea_wrapper((rwr.up - rwr.down), list_db_annots, rwr_cutoff = 0, parallel = PARALLEL)

out <- list(KEGG = t(res[,grep("^KEGG", colnames(res))]), 
            GO = t(res[,grep("^GO", colnames(res))]), 
            REACTOME = t(res[,grep("^REACTOME", colnames(res))]),
            HALLMARK = t(res[,grep("^HALLMARK", colnames(res))]))

write.csv(out$KEGG, gzfile(paste0(path_intermediate_results, "/brca/rwr/ppi_rwr_otp/KEGG.csv.gz")))
write.csv(out$GO, gzfile(paste0(path_intermediate_results, "/brca/rwr/ppi_rwr_otp/GO.csv.gz")))
write.csv(out$REACTOME, gzfile(paste0(path_intermediate_results, "/brca/rwr/ppi_rwr_otp/REACTOME.csv.gz")))
write.csv(out$HALLMARK, gzfile(paste0(path_intermediate_results, "/brca/rwr/ppi_rwr_otp/HALLMARK.csv.gz")))