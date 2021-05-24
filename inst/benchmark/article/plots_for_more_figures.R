source("load_config.R")

plot_scale <- 1.6

plot_point_size <- 1.5

layout_heights <- c(3.8,4.2,1.3)

source("brca/brca_load_results.R")
brca_scores <- brca_scores[brca_scores$Clustering != "HC_single",]

temp <- brca_scores <- brca_scores[brca_scores$Approach %in% c("DR*", "DR"),]
#temp <- temp[temp$k %in% 3:6, ]
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
#temp <- temp[temp$k %in% 3:6, ]
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

# UMAP

brca_umap_scores <- brca_scores[grepl("^UMAP", brca_scores$Embedding),]
brca_umap_scores$Dimensions <- sapply(strsplit(brca_umap_scores$Embedding, " "), function(x) x[2])
brca_umap_scores$Neighbours <- sapply(strsplit(brca_umap_scores$Embedding, " "), function(x) x[3])
brca_umap_scores$Method <- paste(brca_umap_scores$Approach, brca_umap_scores$Clustering, sep = "+")
brca_umap_scores$Dimensions <- factor(brca_umap_scores$Dimensions, levels = paste0(2:10,"d"))

brca_umap_scores$Unsupervised_wsum <- brca_umap_scores$Silhouette + 
  brca_umap_scores$ClusteringStabilityJaccard + 
  brca_umap_scores$Module_score + 
  brca_umap_scores$SurvivalPValue_score + 
  1 - brca_umap_scores$NMI.tss

ggplot(brca_umap_scores, aes(x = Dimensions, y = Unsupervised_wsum, group = Method, color = Clustering, shape = Approach)) + geom_line() + 
  theme_bw() + scale_colour_brewer(palette = "Dark2") + geom_point(size = plot_point_size) + facet_grid(Neighbours ~ k, scales = "fixed") + 
  #theme(axis.text.x = element_blank(), axis.title.x = element_blank()) + 
  theme(legend.position="bottom") + 
  guides(shape = guide_legend(ncol = 3, title.position = "top"), color = guide_legend(ncol = 3, title.position = "top"))

brca_umap_scores[order(brca_umap_scores$Unsupervised_wsum, decreasing = TRUE)[1:10],]
brca_umap_overall <- plyr::ddply(brca_umap_scores, c("Dimensions", "Neighbours"), function(x) data.frame(overall_wsum = mean(x$Unsupervised_wsum)))
brca_umap_overall[order(brca_umap_overall$overall_wsum, decreasing = TRUE),]

temp <- brca_umap_scores[order(brca_umap_scores$Unsupervised_wsum, decreasing = TRUE),]
temp[temp$ClusteringStabilityJaccard > 0.8,][1:10,]

prad_umap_scores <- prad_scores[grepl("^UMAP", prad_scores$Embedding),]
prad_umap_scores$Dimensions <- sapply(strsplit(prad_umap_scores$Embedding, " "), function(x) x[2])
prad_umap_scores$Neighbours <- sapply(strsplit(prad_umap_scores$Embedding, " "), function(x) x[3])
prad_umap_scores$Method <- paste(prad_umap_scores$Approach, prad_umap_scores$Clustering, sep = "+")
prad_umap_scores$Dimensions <- factor(prad_umap_scores$Dimensions, levels = paste0(2:10,"d"))

prad_umap_scores$Unsupervised_wsum <- prad_umap_scores$Silhouette + 
  prad_umap_scores$ClusteringStabilityJaccard + 
  prad_umap_scores$Module_score + 
  prad_umap_scores$SurvivalPValue_score + 
  1 - prad_umap_scores$NMI.tss

ggplot(prad_umap_scores, aes(x = Dimensions, y = Unsupervised_wsum, group = Method, color = Clustering, shape = Approach)) + geom_line() + 
  theme_bw() + scale_colour_brewer(palette = "Dark2") + geom_point(size = plot_point_size) + facet_grid(Neighbours ~ k, scales = "fixed") + 
  #theme(axis.text.x = element_blank(), axis.title.x = element_blank()) + 
  theme(legend.position="bottom") + 
  guides(shape = guide_legend(ncol = 3, title.position = "top"), color = guide_legend(ncol = 3, title.position = "top"))

prad_umap_scores[order(prad_umap_scores$Unsupervised_wsum, decreasing = TRUE)[1:10],]
prad_umap_overall <- plyr::ddply(prad_umap_scores, c("Dimensions", "Neighbours"), function(x) data.frame(overall_wsum = mean(x$Unsupervised_wsum)))
prad_umap_overall[order(prad_umap_overall$overall_wsum, decreasing = TRUE),]

# PCA

brca_pca_scores <- brca_scores[grepl("^PCA", brca_scores$Embedding),]
brca_pca_scores$Dimensions <- sapply(strsplit(brca_pca_scores$Embedding, " "), function(x) x[2])
brca_pca_scores$Method <- paste(brca_pca_scores$Approach, brca_pca_scores$Clustering, sep = "+")
brca_pca_scores$Dimensions <- factor(brca_pca_scores$Dimensions, levels = paste0(2:10,"d"))

brca_pca_scores$Unsupervised_wsum <- brca_pca_scores$Silhouette + 
  brca_pca_scores$ClusteringStabilityJaccard + 
  brca_pca_scores$Module_score + 
  brca_pca_scores$SurvivalPValue_score + 
  1 - brca_pca_scores$NMI.tss

ggplot(brca_pca_scores, aes(x = k, y = Unsupervised_wsum, group = Method, color = Clustering, shape = Approach)) + geom_line() + 
  theme_bw() + scale_colour_brewer(palette = "Dark2") + geom_point(size = plot_point_size) + facet_grid(. ~ Dimensions, scales = "fixed") + 
  #theme(axis.text.x = element_blank(), axis.title.x = element_blank()) + 
  theme(legend.position="bottom") + 
  guides(shape = guide_legend(ncol = 3, title.position = "top"), color = guide_legend(ncol = 3, title.position = "top"))

brca_pca_scores[order(brca_pca_scores$Unsupervised_wsum, decreasing = TRUE)[1:10],]
brca_pca_overall <- plyr::ddply(brca_pca_scores, c("Dimensions"), function(x) data.frame(overall_wsum = mean(x$Unsupervised_wsum)))
brca_pca_overall[order(brca_pca_overall$overall_wsum, decreasing = TRUE),]


prad_pca_scores <- prad_scores[grepl("^PCA", prad_scores$Embedding),]
prad_pca_scores$Dimensions <- sapply(strsplit(prad_pca_scores$Embedding, " "), function(x) x[2])
prad_pca_scores$Method <- paste(prad_pca_scores$Approach, prad_pca_scores$Clustering, sep = "+")
prad_pca_scores$Dimensions <- factor(prad_pca_scores$Dimensions, levels = paste0(2:10,"d"))

prad_pca_scores$Unsupervised_wsum <- prad_pca_scores$Silhouette + 
  prad_pca_scores$ClusteringStabilityJaccard + 
  prad_pca_scores$Module_score + 
  prad_pca_scores$SurvivalPValue_score + 
  1 - prad_pca_scores$NMI.tss

ggplot(prad_pca_scores, aes(x = k, y = Unsupervised_wsum, group = Method, color = Clustering, shape = Approach)) + geom_line() + 
  theme_bw() + scale_colour_brewer(palette = "Dark2") + geom_point(size = plot_point_size) + facet_grid(. ~ Dimensions, scales = "fixed") + 
  #theme(axis.text.x = element_blank(), axis.title.x = element_blank()) + 
  theme(legend.position="bottom") + 
  guides(shape = guide_legend(ncol = 3, title.position = "top"), color = guide_legend(ncol = 3, title.position = "top"))

prad_pca_scores[order(prad_pca_scores$Unsupervised_wsum, decreasing = TRUE)[1:10],]
prad_pca_overall <- plyr::ddply(prad_pca_scores, c("Dimensions"), function(x) data.frame(overall_wsum = mean(x$Unsupervised_wsum)))
prad_pca_overall[order(prad_pca_overall$overall_wsum, decreasing = TRUE),]


## Silhouette scores for UMAP and PCA
brca_extra_silhouette_umap <- ggplot(brca_scores[grep("^UMAP.*20n", brca_scores$drname),], 
                                     aes(x = k, y = Silhouette, group = Method, color = Embedding)) + geom_line() + 
  theme_bw() + scale_colour_brewer(palette = "Set3") + geom_point(size = plot_point_size) + facet_grid(Approach ~ Clustering, scales = "fixed") + 
  theme(axis.text.x = element_blank(), axis.title.x = element_blank()) + 
  theme(legend.position="bottom") + 
  guides(shape = guide_legend(ncol = 3, title.position = "top"), color = guide_legend(ncol = 3, title.position = "top"))
brca_extra_silhouette_umap_legend <-  cowplot::get_legend(brca_extra_silhouette_umap)
brca_extra_silhouette_umap <- brca_extra_silhouette_umap + theme(legend.position = "none")

prad_extra_silhouette_umap <- ggplot(prad_scores[grep("^UMAP.*20n", prad_scores$drname),], 
                                     aes(x = k, y = Silhouette, group = Method, color = Embedding)) + geom_line() + 
  theme_bw() + scale_colour_brewer(palette = "Set3") + geom_point(size = plot_point_size) + facet_grid(Approach ~ Clustering, scales = "fixed") + 
  theme(legend.position="bottom") + 
  guides(shape = guide_legend(ncol = 3, title.position = "top"), color = guide_legend(ncol = 3, title.position = "top"))
prad_extra_silhouette_umap <- prad_extra_silhouette_umap + theme(legend.position = "none")

brca_extra_silhouette_pca <- ggplot(brca_scores[grep("^PCA", brca_scores$drname),], 
                                     aes(x = k, y = Silhouette, group = Method, color = Embedding)) + geom_line() + 
  theme_bw() + scale_colour_brewer(palette = "Set3") + geom_point(size = plot_point_size) + facet_grid(Approach ~ Clustering, scales = "fixed") + 
  theme(axis.text.x = element_blank(), axis.title.x = element_blank()) + 
  theme(legend.position="bottom") + 
  guides(shape = guide_legend(ncol = 3, title.position = "top"), color = guide_legend(ncol = 3, title.position = "top"))
brca_extra_silhouette_pca_legend <-  cowplot::get_legend(brca_extra_silhouette_pca)
brca_extra_silhouette_pca <- brca_extra_silhouette_pca + theme(legend.position = "none")

prad_extra_silhouette_pca <- ggplot(prad_scores[grep("^PCA", prad_scores$drname),], 
                                     aes(x = k, y = Silhouette, group = Method, color = Embedding)) + geom_line() + 
  theme_bw() + scale_colour_brewer(palette = "Set3") + geom_point(size = plot_point_size) + facet_grid(Approach ~ Clustering, scales = "fixed") + 
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














