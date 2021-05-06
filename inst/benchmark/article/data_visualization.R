source("load_config.R")

library(ggplot2)

plot_scale <- 1.6
supplementary_plot_scale <- plot_scale

## Visualization

triple_viz <- function(expr, category, category_label, tsne_perplexity = 45, umap_neighbors = 20) {
  expr_pca <- FactoMineR::PCA(expr, scale.unit = FALSE, ncp = 2, graph = FALSE)
  expr_pca_dat <- as.data.frame(expr_pca$ind$coord)
  expr_pca_dat$category <- category
  eig_percentages <- expr_pca$eig[,"percentage of variance"]
  eig_percentages <- as.character(signif(eig_percentages, 3))
  p1 <- ggplot(expr_pca_dat, aes(Dim.1, Dim.2, color = category)) + geom_point(shape = "+", size = 3) + 
    theme_bw() + scale_color_brewer(palette = "Dark2") + 
    labs(x = paste0("PC1 (", eig_percentages[1], "%)"), y = paste0("PC2 (", eig_percentages[2], "%)"), color = category_label) +
    ggtitle("PCA")
  
  expr_tsne <- Rtsne::Rtsne(expr,
                            dims = 2,
                            perplexity = tsne_perplexity,
                            initial_dims = min(50, dim(expr)[2]),
                            check_duplicates = FALSE,
                            pca = TRUE,
                            partial_pca = TRUE,
                            verbose = FALSE)$Y
  expr_tsne <- as.data.frame(expr_tsne)
  expr_tsne$category <- category
  p2 <- ggplot(expr_tsne, aes(V1, V2, color = category)) + geom_point(shape = "+", size = 3) + 
    theme_bw() + scale_color_brewer(palette = "Dark2") + 
    labs(x = "Z1", y = "Z2", color = category_label) +
    ggtitle("t-SNE")
  
  expr_umap <- uwot::umap(expr, n_neighbors = umap_neighbors, n_components = 2, pca = min(50, dim(expr)[2]), verbose = FALSE, init = "normlaplacian")
  expr_umap <- data.frame(Dim.1 = expr_umap[,1], Dim.2 = expr_umap[,2])
  expr_umap$category <- category
  p3 <- ggplot(expr_umap, aes(Dim.1, Dim.2, color = category)) + geom_point(shape = "+", size = 3) + 
    theme_bw() + scale_color_brewer(palette = "Dark2") + 
    labs(x = "Z1", y = "Z2", color = category_label) + 
    ggtitle("UMAP")
  
  return(list(PCA = p1, tSNE = p2, UMAP = p3))
}

# BRCA
source("brca/brca_default_parameters.R")
source("brca/tcga_brca_mrna_data.R")

brca_gene_filter <- rownames(tbrca_norm)[-zero_var]

set.seed(0)
brca_plots <- triple_viz(t(log2(tbrca_norm[brca_gene_filter,]+1)), brca_norm_subtypes_all$BRCA_Subtype_PAM50, "BRCA Subtype")
brca_legend <- cowplot::get_legend(brca_plots$PCA + theme(legend.position = "bottom") + guides(color = guide_legend(nrow = 1)))

if (save_plots_pdf) pdf(paste0(path_plots, "/brca_visualizations.pdf"), 
                               width = supplementary_plot_scale * 6, 
                               height = supplementary_plot_scale * 2.25)
if (save_plots_svg) svg(paste0(path_plots, "/brca_visualizations.svg"), 
                        width = supplementary_plot_scale * 6, 
                        height = supplementary_plot_scale * 2.25)
gridExtra::grid.arrange(grobs = c(lapply(brca_plots, function(x) x + theme(legend.position = "none")), list(brca_legend)), 
                        nrow = 2, ncol = 3, layout_matrix = matrix(c(1,4,2,4,3,4), 2, 3), widths = c(1,1,1), heights = c(1,0.1))
if (save_plots_svg) dev.off()
if (save_plots_pdf) dev.off()

# PRAD
source("prad/prad_default_parameters.R")
source("prad/tcga_prad_mrna_data.R")

prad_gene_filter <- rownames(tprad_norm)[-zero_var]

set.seed(0)
prad_plots <- triple_viz(t(log2(tprad_norm[prad_gene_filter,]+1)), prad_subtype$Gleason_category, "Gleason score category")
prad_legend <- cowplot::get_legend(prad_plots$PCA + theme(legend.position = "bottom") + guides(color = guide_legend(nrow = 1)))

if (save_plots_pdf) pdf(paste0(path_plots,"/prad_visualizations.pdf"), 
                               width = supplementary_plot_scale * 6, 
                               height = supplementary_plot_scale * 2.25)
if (save_plots_svg) svg(paste0(path_plots, "/prad_visualizations.svg"), 
                        width = supplementary_plot_scale * 6, 
                        height = supplementary_plot_scale * 2.25)
gridExtra::grid.arrange(grobs = c(lapply(prad_plots, function(x) x + theme(legend.position = "none")), list(prad_legend)), 
                        nrow = 2, ncol = 3, layout_matrix = matrix(c(1,4,2,4,3,4), 2, 3), widths = c(1,1,1), heights = c(1,0.1))
if (save_plots_svg) dev.off()
if (save_plots_pdf) dev.off()



