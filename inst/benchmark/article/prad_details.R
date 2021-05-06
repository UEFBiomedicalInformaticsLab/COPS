## PRAD details
source("load_config.R")
plot_scale <- 1.6

library(parallel)
library(ggplot2)

# need to run pareto first
if(!exists("prad_pres")) {
  source("article/plot_data.R")
}

# Load clusterings
source("prad/prad_load_clusterings.R")

# Load other data
source("prad/prad_default_parameters.R")
source("prad/tcga_prad_mrna_data.R")
source("prad/prad_default_annotations.R")

# Select a Pareto optimal result with the highest survival relevance
prad_pres_final_selection <- prad_pres
prad_pres_final_selection$Approach <- prad_pres_final_selection$Approach_unmodified
prad_pres_final_selection <- prad_pres_final_selection[which.max(prad_pres_final_selection$SurvivalPValue_score),]

# This script was made for analysing pathway based clustering results. 
# Details of DR based results could be generated with gene expression. 
# Regardless, the best result in the default configuration should be GSVA. 

# Load pathway based clusters and activity profiles
prad_selected_clusters <- list()
prad_pw_enrichment <- list()
if (any(grepl("^DiffRank", prad_pres_final_selection$Approach))) {
  prad_selected_clusters$diffrank <- plyr::join(prad_pres_final_selection[c("datname", "drname", "k", "m")], 
                                                prad_clust_pw_otp_diffrank[prad_clust_pw_otp_diffrank$run == 1 & 
                                                                             prad_clust_pw_otp_diffrank$fold == 6, ], 
                                                by = c("datname", "drname", "k", "m"), 
                                                type = "inner")
  prad_diffrank <- COPS::DiffRank(tprad_norm[combined_gene_filter,], list_db_annots, parallel = PARALLEL)
  prad_pw_enrichment <- c(prad_pw_enrichment, prad_diffrank)
}
if (any(grepl("^GSVA", prad_pres_final_selection$Approach))) {
  prad_selected_clusters$gsva <- plyr::join(prad_pres_final_selection[c("datname", "drname", "k", "m")], 
                                            prad_clust_pw_otp_gsva[prad_clust_pw_otp_gsva$run == 1 & 
                                                                     prad_clust_pw_otp_gsva$fold == 6, ], 
                                            by = c("datname", "drname", "k", "m"), 
                                            type = "inner")
  prad_gsva <- COPS::genes_to_pathways(log2(tprad_norm[combined_gene_filter,] + 1), "GSVA", list_db_annots, parallel = 1)
  names(prad_gsva) <- gsub("PW", "GSVA", names(prad_gsva))
  prad_pw_enrichment <- c(prad_pw_enrichment, prad_gsva)
}
if (any(grepl("^GCN", prad_pres_final_selection$Approach))) {
  prad_selected_clusters$rwr_gcn <- plyr::join(prad_pres_final_selection[c("datname", "drname", "k", "m")], 
                                               prad_clust_pw_otp_gcn_rwr_fgsea[prad_clust_pw_otp_gcn_rwr_fgsea$run == 1 & 
                                                                               prad_clust_pw_otp_gcn_rwr_fgsea$fold == 6, ], 
                                               by = c("datname", "drname", "k", "m"), 
                                               type = "inner")
  prad_selected_clusters$rwr_gcn$datname <- gsub("RWR", "GCN_RWR_FGSEA", prad_selected_clusters$rwr_gcn$datname)
  prad_gcn_rwr_kegg <- read.csv(paste0(path_intermediate_results,"/prad/rwr/gcn_rwr_otp/KEGG.csv.gz"), header = TRUE, row.names = 1)
  prad_gcn_rwr_go <- read.csv(paste0(path_intermediate_results,"/prad/rwr/gcn_rwr_otp/GO.csv.gz"), header = TRUE, row.names = 1)
  prad_gcn_rwr_reactome <- read.csv(paste0(path_intermediate_results,"/prad/rwr/gcn_rwr_otp/REACTOME.csv.gz"), header = TRUE, row.names = 1)
  prad_gcn <- list(KEGG_GCN_RWR_FGSEA = prad_gcn_rwr_kegg,
                   GO_GCN_RWR_FGSEA = prad_gcn_rwr_go,
                   REACTOME_GCN_RWR_FGSEA = prad_gcn_rwr_reactome)
  
  prad_pw_enrichment <- c(prad_pw_enrichment, prad_gcn)
}
if (any(grepl("^PPI", prad_pres_final_selection$Approach))) {
  prad_selected_clusters$rwr_ppi <- plyr::join(prad_pres_final_selection[c("datname", "drname", "k", "m")], 
                                               prad_clust_pw_otp_ppi_rwr_fgsea[prad_clust_pw_otp_ppi_rwr_fgsea$run == 1 & 
                                                                               prad_clust_pw_otp_ppi_rwr_fgsea$fold == 6, ], 
                                               by = c("datname", "drname", "k", "m"), 
                                               type = "inner")
  prad_selected_clusters$rwr_ppi$datname <- gsub("RWR", "PPI_RWR_FGSEA", prad_selected_clusters$rwr_ppi$datname)
  prad_ppi_rwr_kegg <- read.csv(paste0(path_intermediate_results,"/prad/rwr/ppi_rwr_otp/KEGG.csv.gz"), header = TRUE, row.names = 1)
  prad_ppi_rwr_go <- read.csv(paste0(path_intermediate_results,"/prad/rwr/ppi_rwr_otp/GO.csv.gz"), header = TRUE, row.names = 1)
  prad_ppi_rwr_reactome <- read.csv(paste0(path_intermediate_results,"/prad/rwr/ppi_rwr_otp/REACTOME.csv.gz"), header = TRUE, row.names = 1)
  prad_ppi <- list(KEGG_PPI_RWR_FGSEA = prad_ppi_rwr_kegg,
                   GO_PPI_RWR_FGSEA = prad_ppi_rwr_go,
                   REACTOME_PPI_RWR_FGSEA = prad_ppi_rwr_reactome)
  
  prad_pw_enrichment <- c(prad_pw_enrichment, prad_ppi)
}

prad_selected_clusters_split <- lapply(prad_selected_clusters, function(x) split(x, x[c("datname", "m", "k")], drop = TRUE))
prad_selected_clusters_split <- do.call(c, prad_selected_clusters_split)
prad_selected_clusters_split <- lapply(prad_selected_clusters, rename_methods)
length(prad_selected_clusters_split) # Should be 1

kw_test_p <- list()
anova_p <- list()
correlation_pw <- list()
median_pw <- list()
mean_pw <- list()
for (i in names(prad_selected_clusters_split)) {
  pw_i <- prad_selected_clusters_split[[i]]$datname[1]
  id_i <- prad_selected_clusters_split[[i]]$id
  clust_i <- prad_selected_clusters_split[[i]]$cluster
  pw_dat_i <- data.frame(t(prad_pw_enrichment[[pw_i]]))
  pw_names <- colnames(pw_dat_i)
  pw_dat_i$cluster <- clust_i[match(gsub(".", "-", rownames(pw_dat_i), fixed = TRUE), id_i)]
  
  kw_test_p[[i]] <- c()
  anova_p[[i]] <- c()
  for (j in pw_names) {
    anova_p[[i]][j] <- summary(aov(as.formula(paste(j, "~ cluster")), pw_dat_i))[[1]][["Pr(>F)"]][1]
    kw_test_p[[i]][j] <- kruskal.test(as.formula(paste(j, "~ cluster")), pw_dat_i)$p.value
    clust_cor_j <- c()
    for (k in unique(clust_i)) {
      clust_cor_j[k] <- cor(pw_dat_i$clust == k, pw_dat_i[[j]], method = "spearman")
    }
    correlation_pw[[i]] <- rbind(correlation_pw[[i]], data.frame(t(clust_cor_j)))
    rownames(correlation_pw[[i]])[nrow(correlation_pw[[i]])] <- j
    median_pw[[i]] <- rbind(median_pw[[i]], data.frame(t(tapply(pw_dat_i[[j]], pw_dat_i$cluster, median, na.rm = TRUE))))
    rownames(median_pw[[i]])[nrow(median_pw[[i]])] <- j
    mean_pw[[i]] <- rbind(mean_pw[[i]], data.frame(t(tapply(pw_dat_i[[j]], pw_dat_i$cluster, mean, na.rm = TRUE))))
    rownames(mean_pw[[i]])[nrow(mean_pw[[i]])] <- j
  }
}

kw_test_p_adj <- lapply(kw_test_p, p.adjust, method = "BH")
anova_p_adj <- lapply(anova_p, p.adjust, method = "BH")

# From here on the script only supports one clustering result (which is selected at the start)

# Correlation based plot
pw_order <- hclust(dist(as.matrix(correlation_pw[[1]])), method = "average")$order
pw_correlation <- reshape::melt(as.matrix(correlation_pw[[1]]))
colnames(pw_correlation) <- c("pathway", "cluster", "correlation")
pw_correlation$pathway <- factor(pw_correlation$pathway, levels = rownames(correlation_pw[[1]])[pw_order])

pw_correlation <- pw_correlation[pw_correlation$pathway %in% names(kw_test_p_adj[[1]])[kw_test_p_adj[[1]] < 1e-50], ]

pw_correlation$cluster <- factor(gsub("X", "", pw_correlation$cluster))
pathway_cor_plot <- ggplot(pw_correlation, aes(cluster, pathway, fill = correlation)) + geom_tile() + 
  scale_fill_distiller(palette = "RdBu", limits = c(-1,1)) + 
  theme_bw() + #coord_fixed() + 
  ylab("Most informative pathways")

# Median based plot
pw_order <- hclust(dist(as.matrix(median_pw[[1]])), method = "average")$order
pw_median <- reshape::melt(as.matrix(median_pw[[1]]))
colnames(pw_median) <- c("pathway", "cluster", "median")
pw_median$pathway <- factor(pw_median$pathway, levels = rownames(median_pw[[1]])[pw_order])

pw_median <- pw_median[pw_median$pathway %in% names(kw_test_p_adj[[1]])[kw_test_p_adj[[1]] < 1e-50], ]

pw_median$cluster <- factor(gsub("X", "", pw_median$cluster))
ggplot(pw_median, aes(cluster, pathway, fill = median)) + geom_tile() + 
  colorspace::scale_fill_continuous_divergingx(palette = "RdBu", mid = 0, rev = TRUE) + 
  theme_bw() + #coord_fixed() + 
  ylab("Most informative pathways")

# Mean
pw_order <- hclust(dist(as.matrix(mean_pw[[1]])), method = "average")$order
pw_mean <- reshape::melt(as.matrix(mean_pw[[1]]))
colnames(pw_mean) <- c("pathway", "cluster", "mean")
pw_mean$pathway <- factor(pw_mean$pathway, levels = rownames(mean_pw[[1]])[pw_order])

pw_mean <- pw_mean[pw_mean$pathway %in% names(anova_p_adj[[1]])[anova_p_adj[[1]] < 1e-50], ]

pw_mean$cluster <- factor(gsub("X", "", pw_mean$cluster))
ggplot(pw_mean, aes(cluster, pathway, fill = mean)) + geom_tile() + 
  colorspace::scale_fill_continuous_divergingx(palette = "RdBu", mid = 0, rev = TRUE) + 
  theme_bw() + #coord_fixed() + 
  ylab("Most informative pathways")

### Modules
module_eigen_genes <- read.csv(paste0(path_intermediate_results,"/prad/wgcna/module_eigen_genes_deg.csv"), header = TRUE, row.names = 1)

colnames(module_eigen_genes) <- paste0("module", 1:ncol(module_eigen_genes))

module_data <- t(module_eigen_genes)
module_names <- rownames(module_data)
module_clustering_data_list_selected <- lapply(prad_selected_clusters_split, 
                                               function(x) data.frame(x, t(module_data[,match(x$id, colnames(module_data))])))
module_correlation <- function(dat) {
  # Create matrix of one vs. all correlations for each cluster
  clust_cor <- lapply(module_names, function(x) sapply(unique(dat$cluster), function(y) cor(dat[[x]], dat$cluster == y, method = "spearman")))
  clust_cor_mat <- do.call("rbind", args = clust_cor)
  colnames(clust_cor_mat) <- 1:ncol(clust_cor_mat)
  rownames(clust_cor_mat) <- module_names
  module_order <- stats::hclust(dist(clust_cor_mat), method = "average")$order
  clust_cor_mat <- reshape::melt(clust_cor_mat)
  colnames(clust_cor_mat) <- c("WGCNA_module", "cluster", "cor")
  clust_cor_mat$WGCNA_module <- factor(clust_cor_mat$WGCNA_module, levels = module_names[module_order])
  clust_cor_mat$cluster <- factor(clust_cor_mat$cluster)
  clust_cor_mat$Method <- dat$Method[1]
  return(clust_cor_mat)
}
module_correlations_selected <- lapply(module_clustering_data_list_selected, module_correlation)
module_correlations_selected <- plyr::rbind.fill(module_correlations_selected)

module_correlations_selected$correlation <- module_correlations_selected$cor
module_cor_plot <- ggplot(module_correlations_selected, aes(cluster, WGCNA_module, fill = correlation)) + geom_tile() + 
  scale_fill_distiller(palette = "RdBu", limits = c(-1,1)) + 
  theme_bw() + #coord_fixed() + 
  ylab("WGCNA module / closest pathway")

## Survival
#library(survminer)
prad_survival_plot_data <- prad_survival
prad_survival_plot_data$Cluster <- prad_selected_clusters_split$gsva$cluster[match(prad_survival_plot_data$ID, prad_selected_clusters_split$gsva$id)]
prad_survival_plot_data$Cluster <- factor(prad_survival_plot_data$Cluster)

prad_survival_plot_data[["N_stage"]] <- prad_survival_plot_data$pathology_N_stage
prad_survival_plot_data[["preop_PSA"]] <- prad_survival_plot_data$patient.clinical_cqcf.psa_result_preop

prad_coxph <- survival::coxph(survival::Surv(time, event) ~ N_stage + preop_PSA + Cluster, 
                              data = prad_survival_plot_data)
cox_hazard_ratios <- survminer::ggforest(prad_coxph, data = prad_survival_plot_data, fontsize = 0.8)

survival::survdiff(survival::Surv(time, event) ~ Cluster, data = prad_survival_plot_data)

prad_km <- survminer::surv_fit(survival::Surv(time, event) ~ Cluster, 
                               data = prad_survival_plot_data)
km_plot <- survminer::ggsurvplot(prad_km, data = prad_survival_plot_data, risk.table = TRUE, 
                                 palette = "Dark2", 
                                 #color = "Cluster", 
                                 pval = TRUE, pval.method = TRUE)

km_plot_grob <- ggplotGrob(km_plot$plot + labs(tag = "A"))
km_table_grob <- ggplotGrob(km_plot$table)
km_grob_max_widths <- grid::unit.pmax(km_plot_grob$widths[2:5], km_table_grob$widths[2:5])
km_plot_grob$widths[2:5] <- as.list(km_grob_max_widths)
km_table_grob$widths[2:5] <- as.list(km_grob_max_widths)

## Combined plots
pathway_cor_plot_grob <- ggplotGrob(pathway_cor_plot + labs(tag = "B") + theme(axis.text.y = element_text(size = 8)))
module_cor_plot_grob <- ggplotGrob(module_cor_plot + labs(tag = "D") + 
                                     scale_fill_distiller(palette = "PuOr", limits = c(-1,1)) + 
                                     ylab("WGCNA module"))

if (save_plots_pdf) pdf(paste0(path_plots,"/prad_best_clusters.pdf"), width = plot_scale * 8, height = plot_scale * 8)
if (save_plots_svg) svg(paste0(path_plots,"/prad_best_clusters.svg"), width = plot_scale * 10, height = plot_scale * 10)
gridExtra::grid.arrange(km_plot_grob, 
                        km_table_grob, 
                        cox_hazard_ratios + labs(tag = "C") +  theme(plot.tag = element_text(face = "plain")), 
                        pathway_cor_plot_grob, 
                        module_cor_plot_grob, 
                        layout_matrix = matrix(c(1,2,3, 1,2,3, 4,4,3, 4,4,5), 3, 4), 
                        heights = c(3,1,4), 
                        widths = c(1,1,1,1))
if (save_plots_svg) dev.off()
if (save_plots_pdf) dev.off()

### Pathway genes
# Pathway similarity
pw_filter <-  names(kw_test_p_adj[[1]])[kw_test_p_adj[[1]] < 1e-50]
pw_indicators <- plyr::rbind.fill.matrix(lapply(list_db_annots[pw_filter], function(x) t(sapply(x, function(y) 1))))
pw_indicators[is.na(pw_indicators)] <- 0
pw_indicators <- pw_indicators[, combined_gene_filter[combined_gene_filter %in% colnames(pw_indicators)]]
pw_gene_expression <- tprad_norm[colnames(pw_indicators),]

clusters <- clust_i
names(clusters) <- colnames(tprad_norm)

pw_gene_kw <- apply(pw_gene_expression, 1, function(x) kruskal.test(x ~ clust_i)$p.value)
pw_gene_kw[is.na(pw_gene_kw)] <- 1


if (save_plots_pdf) {
  pheatmap::pheatmap(log2(pw_gene_expression[, order(clusters)] + 1), 
                     annotation_col = data.frame(Cluster = factor(clusters[order(clusters)])), 
                     clustering_method = "average", 
                     clustering_distance_cols = "correlation", 
                     clustering_distance_rows = "correlation", 
                     cluster_cols = FALSE, labels_col = "",
                     main = "Expression of genes in significant pathways for PRAD",
                     filename = paste0(path_plots, "/prad_pw_genes.pdf"), width = plot_scale * 6, height = plot_scale * 6)
} 
if (save_plots_svg) {
  pheatmap::pheatmap(log2(pw_gene_expression[, order(clusters)] + 1), 
                     annotation_col = data.frame(Cluster = factor(clusters[order(clusters)])), 
                     clustering_method = "average", 
                     clustering_distance_cols = "correlation", 
                     clustering_distance_rows = "correlation", 
                     cluster_cols = FALSE, labels_col = "",
                     main = "Expression of genes in significant pathways for PRAD",
                     filename = paste0(path_plots, "/prad_pw_genes.svg"), width = plot_scale * 6, height = plot_scale * 6)
}
if (!save_plots_pdf & !save_plots_svg) {
  pheatmap::pheatmap(log2(pw_gene_expression[, order(clusters)] + 1), 
                     annotation_col = data.frame(Cluster = factor(clusters[order(clusters)])), 
                     clustering_method = "average", 
                     clustering_distance_cols = "correlation", 
                     clustering_distance_rows = "correlation", 
                     cluster_cols = FALSE, labels_col = "",
                     main = "Expression of genes in significant pathways for PRAD")
}
