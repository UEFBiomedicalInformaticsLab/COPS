source("load_config.R")

plot_scale <- 1.6

plot_point_size <- 1.5

layout_heights <- c(3.8,4.2,1.3)

source("load_config.R")
source("brca/brca_load_results.R")
brca_scores <- brca_scores[brca_scores$Clustering != "HC_single",]

temp <- brca_scores <- brca_scores[brca_scores$Approach %in% c("DR*", "DR"),]
temp <- temp[temp$k %in% 3:6, ]
# Manual renaming to keep details
temp$drname <- gsub("^pca", "PCA ", temp$drname)
temp$drname <- gsub("^tsne", "t-SNE ", temp$drname)
temp$drname <- gsub("^umap", "UMAP ", temp$drname)
temp$drname <- gsub("$", "d", temp$drname)
temp$drname <- gsub("_", "d ", temp$drname)
temp$drname <- gsub("nd", "n", temp$drname)
temp$drname[grep("t-SNE", temp$drname)] <- gsub("d", "p", temp$drname[grep("t-SNE", temp$drname)])
brca_scores <- temp
brca_scores$Embedding <- brca_scores$drname

source("prad/prad_load_results.R")
prad_scores <- prad_scores[prad_scores$Clustering != "HC_single",]

temp <- prad_scores <- prad_scores[prad_scores$Approach %in% c("DR*", "DR"),]
temp <- temp[temp$k %in% 3:6, ]
# Manual renaming to keep details
temp$drname <- gsub("^pca", "PCA ", temp$drname)
temp$drname <- gsub("^tsne", "t-SNE ", temp$drname)
temp$drname <- gsub("^umap", "UMAP ", temp$drname)
temp$drname <- gsub("$", "d", temp$drname)
temp$drname <- gsub("_", "d ", temp$drname)
temp$drname <- gsub("nd", "n", temp$drname)
temp$drname[grep("t-SNE", temp$drname)] <- gsub("d", "p", temp$drname[grep("t-SNE", temp$drname)])
prad_scores <- temp
prad_scores$Embedding <- prad_scores$drname

library(ggplot2)
library(cowplot)

brca_extra_silhouette_umap <- ggplot(brca_scores[grep("^UMAP", brca_scores$Embedding),], 
                                     aes(x = k, y = Silhouette, group = Method, color = Embedding)) + geom_line() + 
  theme_bw() + scale_colour_brewer(palette = "Set3") + geom_point(size = plot_point_size) + facet_grid(Approach ~ Clustering, scales = "fixed") + 
  theme(axis.text.x = element_blank(), axis.title.x = element_blank()) + 
  theme(legend.position="bottom") + 
  guides(shape = guide_legend(ncol = 3, title.position = "top"), color = guide_legend(ncol = 3, title.position = "top"))
brca_extra_silhouette_umap_legend <-  cowplot::get_legend(brca_extra_silhouette_umap)
brca_extra_silhouette_umap <- brca_extra_silhouette_umap + theme(legend.position = "none")

prad_extra_silhouette_umap <- ggplot(prad_scores[grep("^UMAP", prad_scores$Embedding),], 
                                     aes(x = k, y = Silhouette, group = Method, color = Embedding)) + geom_line() + 
  theme_bw() + scale_colour_brewer(palette = "Set3") + geom_point(size = plot_point_size) + facet_grid(Approach ~ Clustering, scales = "fixed") + 
  theme(legend.position="bottom") + 
  guides(shape = guide_legend(ncol = 3, title.position = "top"), color = guide_legend(ncol = 3, title.position = "top"))
prad_extra_silhouette_umap <- prad_extra_silhouette_umap + theme(legend.position = "none")

brca_extra_silhouette_pca <- ggplot(brca_scores[grep("^PCA", brca_scores$Embedding),], 
                                     aes(x = k, y = Silhouette, group = Method, color = Embedding)) + geom_line() + 
  theme_bw() + scale_colour_brewer(palette = "Dark2") + geom_point(size = plot_point_size) + facet_grid(Approach ~ Clustering, scales = "fixed") + 
  theme(axis.text.x = element_blank(), axis.title.x = element_blank()) + 
  theme(legend.position="bottom") + 
  guides(shape = guide_legend(ncol = 3, title.position = "top"), color = guide_legend(ncol = 3, title.position = "top"))
brca_extra_silhouette_pca_legend <-  cowplot::get_legend(brca_extra_silhouette_pca)
brca_extra_silhouette_pca <- brca_extra_silhouette_pca + theme(legend.position = "none")

prad_extra_silhouette_pca <- ggplot(prad_scores[grep("^PCA", prad_scores$Embedding),], 
                                     aes(x = k, y = Silhouette, group = Method, color = Embedding)) + geom_line() + 
  theme_bw() + scale_colour_brewer(palette = "Dark2") + geom_point(size = plot_point_size) + facet_grid(Approach ~ Clustering, scales = "fixed") + 
  theme(legend.position="bottom") + 
  guides(shape = guide_legend(ncol = 3, title.position = "top"), color = guide_legend(ncol = 3, title.position = "top"))
prad_extra_silhouette_pca <- prad_extra_silhouette_pca + theme(legend.position = "none")


if (save_plots_pdf) pdf(paste0(path_plots, "/silhouette_more_umap.pdf"), width = plot_scale * 4, height = plot_scale * 5)
if (save_plots_svg) svg(paste0(path_plots, "/silhouette_more_umap.svg"), width = plot_scale * 4, height = plot_scale * 5)
gridExtra::grid.arrange(brca_extra_silhouette_umap + ggtitle("BRCA"), prad_extra_silhouette_umap + ggtitle("PRAD"), 
                        brca_extra_silhouette_umap_legend, 
                        nrow = 3, heights = layout_heights)
if (save_plots_svg) dev.off()
if (save_plots_pdf) dev.off()

if (save_plots_pdf) pdf(paste0(path_plots, "/silhouette_more_pca.pdf"), width = plot_scale * 4, height = plot_scale * 5)
if (save_plots_svg) svg(paste0(path_plots, "/silhouette_more_pca.svg"), width = plot_scale * 4, height = plot_scale * 5)
gridExtra::grid.arrange(brca_extra_silhouette_pca + ggtitle("BRCA"), prad_extra_silhouette_pca + ggtitle("PRAD"), 
                        brca_extra_silhouette_pca_legend, 
                        nrow = 3, heights = layout_heights)
if (save_plots_svg) dev.off()
if (save_plots_pdf) dev.off()














