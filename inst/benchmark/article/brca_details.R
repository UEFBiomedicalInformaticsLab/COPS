## BRCA details
source("load_config.R")
plot_scale <- 1.6

library(parallel)
library(ggplot2)

# need to run pareto first
if(!exists("brca_pres")) {
  source("article/plot_data.R")
}

# Load clustering
source("brca/brca_load_clusterings.R")

# Load other data
source("brca/brca_default_parameters.R")
source("brca/tcga_brca_mrna_data.R")
source("brca/brca_default_annotations.R")

# Select a Pareto optimal result with the highest survival relevance
brca_pres_final_selection <- brca_pres
brca_pres_final_selection$Approach <- brca_pres_final_selection$Approach_unmodified
brca_pres_final_selection <- brca_pres_final_selection[brca_pres_final_selection$Approach == "GCN RWR-FGSEA",]
brca_pres_final_selection <- brca_pres_final_selection[which.max(brca_pres_final_selection$SurvivalPValue_score),]

# This script was made for analysing pathway based clustering results. 
# Details of DR based results could be generated with gene expression. 
# Regardless, the best result in the default configuration should be GCN RWR-FGSEA. 

# Load pathway based clusters and activity profiles
brca_selected_clusters <- list()
brca_pw_enrichment <- list()
if (any(grepl("^DiffRank", brca_pres_final_selection$Approach))) {
  brca_selected_clusters$diffrank <- plyr::join(brca_pres_final_selection[c("datname", "drname", "k", "m")], 
                                                brca_clust_pw_otp_diffrank[brca_clust_pw_otp_diffrank$run == 1 & 
                                                                           brca_clust_pw_otp_diffrank$fold == 6, ], 
                                                by = c("datname", "drname", "k", "m"),
                                                type = "inner")
  brca_diffrank <- COPS::DiffRank(tbrca_norm[combined_gene_filter,], list_db_annots, parallel = PARALLEL)
  brca_pw_enrichment <- c(brca_pw_enrichment, brca_diffrank)
}
if (any(grepl("^GSVA", brca_pres_final_selection$Approach))) {
  brca_selected_clusters$gsva <- plyr::join(brca_pres_final_selection[c("datname", "drname", "k", "m")], 
                                            brca_clust_pw_otp_gsva[brca_clust_pw_otp_gsva$run == 1 & 
                                                                   brca_clust_pw_otp_gsva$fold == 6, ], 
                                            by = c("datname", "drname", "k", "m"),
                                            type = "inner")
  brca_gsva <- COPS::genes_to_pathways(log2(tbrca_norm[combined_gene_filter,] + 1), "GSVA", list_db_annots, parallel = PARALLEL)
  names(brca_gsva) <- gsub("PW", "GSVA", names(brca_gsva))
  brca_pw_enrichment <- c(brca_pw_enrichment, brca_gsva)
}
if (any(grepl("^GCN", brca_pres_final_selection$Approach))) {
  brca_selected_clusters$rwr_gcn <- plyr::join(brca_pres_final_selection[c("datname", "drname", "k", "m")], 
                                               brca_clust_pw_otp_gcn_rwr_fgsea[brca_clust_pw_otp_gcn_rwr_fgsea$run == 1 & 
                                                                               brca_clust_pw_otp_gcn_rwr_fgsea$fold == 6, ], 
                                               by = c("datname", "drname", "k", "m"),
                                               type = "inner")
  brca_selected_clusters$rwr_gcn$datname <- gsub("RWR", "GCN_RWR_FGSEA", brca_selected_clusters$rwr_gcn$datname)
  brca_gcn_rwr_kegg <- read.csv(paste0(path_intermediate_results,"/brca/rwr/gcn_rwr_otp/KEGG.csv.gz"), header = TRUE, row.names = 1)
  brca_gcn_rwr_go <- read.csv(paste0(path_intermediate_results,"/brca/rwr/gcn_rwr_otp/GO.csv.gz"), header = TRUE, row.names = 1)
  brca_gcn_rwr_reactome <- read.csv(paste0(path_intermediate_results,"/brca/rwr/gcn_rwr_otp/REACTOME.csv.gz"), header = TRUE, row.names = 1)
  brca_gcn <- list(KEGG_GCN_RWR_FGSEA = brca_gcn_rwr_kegg,
                   GO_GCN_RWR_FGSEA = brca_gcn_rwr_go,
                   REACTOME_GCN_RWR_FGSEA = brca_gcn_rwr_reactome)
  
  brca_pw_enrichment <- c(brca_pw_enrichment, brca_gcn)
}
if (any(grepl("^PPI", brca_pres_final_selection$Approach))) {
  brca_selected_clusters$rwr_ppi <- plyr::join(brca_pres_final_selection[c("datname", "drname", "k", "m")], 
                                               brca_clust_pw_otp_ppi_rwr_fgsea[brca_clust_pw_otp_ppi_rwr_fgsea$run == 1 & 
                                                                               brca_clust_pw_otp_ppi_rwr_fgsea$fold == 6, ], 
                                               by = c("datname", "drname", "k", "m"),
                                               type = "inner")
  brca_selected_clusters$rwr_ppi$datname <- gsub("RWR", "PPI_RWR_FGSEA", brca_selected_clusters$rwr_ppi$datname)
  brca_ppi_rwr_kegg <- read.csv(paste0(path_intermediate_results,"/brca/rwr/ppi_rwr_otp/KEGG.csv.gz"), header = TRUE, row.names = 1)
  brca_ppi_rwr_go <- read.csv(paste0(path_intermediate_results,"/brca/rwr/ppi_rwr_otp/GO.csv.gz"), header = TRUE, row.names = 1)
  brca_ppi_rwr_reactome <- read.csv(paste0(path_intermediate_results,"/brca/rwr/ppi_rwr_otp/REACTOME.csv.gz"), header = TRUE, row.names = 1)
  brca_ppi <- list(KEGG_PPI_RWR_FGSEA = brca_ppi_rwr_kegg,
                   GO_PPI_RWR_FGSEA = brca_ppi_rwr_go,
                   REACTOME_PPI_RWR_FGSEA = brca_ppi_rwr_reactome)
  
  brca_pw_enrichment <- c(brca_pw_enrichment, brca_ppi)
}

brca_selected_clusters_split <- lapply(brca_selected_clusters, function(x) split(x, x[c("datname", "m", "k")], drop = TRUE))
brca_selected_clusters_split <- do.call(c, brca_selected_clusters_split)
brca_selected_clusters_split <- lapply(brca_selected_clusters, rename_methods)
length(brca_selected_clusters_split) # Should be 1

kw_test_p <- list()
anova_p <- list()
correlation_pw <- list()
median_pw <- list()
mean_pw <- list()
for (i in names(brca_selected_clusters_split)) {
  pw_i <- brca_selected_clusters_split[[i]]$datname[1]
  id_i <- brca_selected_clusters_split[[i]]$id
  clust_i <- brca_selected_clusters_split[[i]]$cluster
  pw_dat_i <- data.frame(t(brca_pw_enrichment[[pw_i]]))
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
    # Scaling specifically for RWR-FGSEA which can have very large differences in scale
    median_pw[[i]] <- rbind(median_pw[[i]], 
                            data.frame(t(tapply(pw_dat_i[[j]], 
                                                pw_dat_i$cluster, 
                                                function(x) median(scale(x, 
                                                                         center = FALSE, 
                                                                         scale = TRUE), 
                                                                   na.rm = TRUE)))))
    rownames(median_pw[[i]])[nrow(median_pw[[i]])] <- j
    mean_pw[[i]] <- rbind(mean_pw[[i]], 
                          data.frame(t(tapply(pw_dat_i[[j]], 
                                              pw_dat_i$cluster, 
                                              mean, 
                                              na.rm = TRUE))))
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
pw_filter <-  names(kw_test_p_adj[[1]])[kw_test_p_adj[[1]] < 1e-50]
pw_order <- hclust(dist(as.matrix(median_pw[[1]][pw_filter, ])), method = "average")$order
pw_median <- reshape::melt(as.matrix(median_pw[[1]][pw_filter,]))
colnames(pw_median) <- c("pathway", "cluster", "median")
pw_median$pathway <- factor(pw_median$pathway, levels = pw_filter[pw_order])

pw_median$cluster <- factor(gsub("X", "", pw_median$cluster))
pathway_median_plot <- ggplot(pw_median, aes(cluster, pathway, fill = median)) + geom_tile() + 
  scale_fill_distiller(palette = "RdBu") + 
  theme_bw() + #coord_fixed() + 
  ylab("Most informative pathways")

# Mean
pw_order <- hclust(dist(as.matrix(mean_pw[[1]])), method = "average")$order
pw_mean <- reshape::melt(as.matrix(mean_pw[[1]]))
colnames(pw_mean) <- c("pathway", "cluster", "mean")
pw_mean$pathway <- factor(pw_mean$pathway, levels = rownames(mean_pw[[1]])[pw_order])

pw_mean <- pw_mean[pw_mean$pathway %in% names(kw_test_p_adj[[1]])[kw_test_p_adj[[1]] < 1e-50], ]

pw_mean$cluster <- factor(gsub("X", "", pw_mean$cluster))
pathway_mean_plot <- ggplot(pw_mean, aes(cluster, pathway, fill = mean)) + geom_tile() + 
  colorspace::scale_fill_continuous_divergingx(palette = "RdBu", mid = 0, rev = TRUE) + 
  theme_bw() + #coord_fixed() + 
  ylab("Most informative pathways")

### Modules
module_eigen_genes <- read.csv(paste0(path_intermediate_results,"/brca/wgcna/module_eigen_genes_deg.csv"), header = TRUE, row.names = 1)

colnames(module_eigen_genes) <- paste0("module", 1:ncol(module_eigen_genes))

module_data <- t(module_eigen_genes)
module_names <- rownames(module_data)
module_clustering_data_list_selected <- lapply(brca_selected_clusters_split, 
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

### Survival
#library(survminer)
brca_survival_plot_data <- dat_survival
brca_survival_plot_data$Cluster <- brca_selected_clusters_split[[1]]$cluster[match(brca_survival_plot_data$ID, 
                                                                                      brca_selected_clusters_split[[1]]$id)]
brca_survival_plot_data$Cluster <- factor(brca_survival_plot_data$Cluster)

brca_survival_plot_data$Stage <- brca_survival_plot_data$pathologic_stage
brca_survival_plot_data$Age <- brca_survival_plot_data$age_at_initial_pathologic_diagnosis

brca_coxph <- survival::coxph(survival::Surv(time, event) ~ Stage + Age + Cluster, 
                              data = brca_survival_plot_data)
cox_hazard_ratios <- survminer::ggforest(brca_coxph, data = brca_survival_plot_data, fontsize = 0.8)

survival::survdiff(survival::Surv(time, event) ~ Cluster, data = brca_survival_plot_data)

brca_km <- survminer::surv_fit(survival::Surv(time, event) ~ Cluster, 
                               data = brca_survival_plot_data)
km_plot <- survminer::ggsurvplot(brca_km, data = brca_survival_plot_data, risk.table = TRUE, 
                                 palette = "Dark2", 
                                 #color = "Cluster", 
                                 pval = TRUE, pval.method = TRUE)

km_plot_grob <- ggplotGrob(km_plot$plot + labs(tag = "A"))
km_table_grob <- ggplotGrob(km_plot$table)
km_grob_max_widths <- grid::unit.pmax(km_plot_grob$widths[2:5], km_table_grob$widths[2:5])
km_plot_grob$widths[2:5] <- as.list(km_grob_max_widths)
km_table_grob$widths[2:5] <- as.list(km_grob_max_widths)

## Combined plots
pathway_cor_plot_grob <- ggplotGrob(pathway_median_plot + labs(tag = "B"))
module_cor_plot_grob <- ggplotGrob(module_cor_plot + labs(tag = "D") + 
                                     scale_fill_distiller(palette = "PuOr", limits = c(-1,1)) + 
                                     ylab("WGCNA module"))

if (save_plots_pdf) pdf(paste0(path_plots,"/brca_best_clusters.pdf"), width = plot_scale * 8, height = plot_scale * 8)
if (save_plots_svg) svg(paste0(path_plots,"/brca_best_clusters.svg"), width = plot_scale * 10, height = plot_scale * 10)
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

# Pathway gene expression
pw_filter <-  names(kw_test_p_adj[[1]])[kw_test_p_adj[[1]] < 1e-50]
pw_indicators <- plyr::rbind.fill.matrix(lapply(list_db_annots[pw_filter], function(x) t(sapply(x, function(y) 1))))
pw_indicators[is.na(pw_indicators)] <- 0
pw_gene_expression <- tbrca_norm[colnames(pw_indicators)[colnames(pw_indicators) %in% rownames(tbrca_norm)],]

clusters <- clust_i
names(clusters) <- colnames(tbrca_norm)
clusters <- data.frame(Cluster = factor(clusters[order(clusters)]))
pw_gene_kw <- apply(pw_gene_expression, 1, function(x) kruskal.test(x ~ clust_i)$p.value)
pw_gene_kw[is.na(pw_gene_kw)] <- 1

clusters$PAM50 <- brca_subtypes$BRCA_Subtype_PAM50[match(substr(rownames(clusters), 1, 12), brca_subtypes$patient)]

if (save_plots_pdf) {
  pheatmap::pheatmap(log2(pw_gene_expression[pw_gene_kw < 1e-50, rownames(clusters)] + 1), 
                     annotation_col = clusters, clustering_method = "average", 
                     clustering_distance_cols = "correlation", clustering_distance_rows = "correlation", cluster_cols = FALSE, 
                     labels_col = "", main = "Most informative pathway gene expression in BRCA tumours",
                     filename = paste0(path_plots, "/brca_pw_genes.pdf"), width = plot_scale * 6, height = plot_scale * 10)
} 
if (save_plots_svg) {
  pheatmap::pheatmap(log2(pw_gene_expression[pw_gene_kw < 1e-50, rownames(clusters)] + 1), 
                     annotation_col = clusters, clustering_method = "average", 
                     clustering_distance_cols = "correlation", clustering_distance_rows = "correlation", cluster_cols = FALSE, 
                     labels_col = "", main = "Most informative pathway gene expression in BRCA tumours",
                     filename = paste0(path_plots, "/brca_pw_genes.svg"), width = plot_scale * 6, height = plot_scale * 10)
}
if (!save_plots_pdf & !save_plots_svg) {
  pheatmap::pheatmap(log2(pw_gene_expression[pw_gene_kw < 1e-50, rownames(clusters)] + 1), 
                     annotation_col = clusters, clustering_method = "average", 
                     clustering_distance_cols = "correlation", clustering_distance_rows = "correlation", cluster_cols = FALSE, 
                     labels_col = "", main = "Most informative pathway gene expression in BRCA tumours")
}

