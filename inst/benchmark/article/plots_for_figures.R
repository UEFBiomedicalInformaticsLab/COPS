library(ggplot2)
plot_scale <- 1.6

plot_point_size <- 1.5

layout_heights <- c(3.8,4.2,1.3)

source("load_config.R")
source("article/plot_data.R")

### Breast cancer

brca_p1 <- ggplot(brca_scores, aes(x = k, y = Silhouette, group = Method, color = Embedding, shape = Clustering)) + geom_line() + 
  theme_bw() + scale_colour_brewer(palette = "Dark2") + geom_point(size = plot_point_size) + facet_wrap(Approach ~ ., scales = "fixed", ncol = 3) + 
  theme(axis.text.x = element_blank(), axis.title.x = element_blank()) + ylim(0.05, 0.7) + 
  theme(legend.position="bottom") + 
  guides(shape = guide_legend(ncol = 3, title.position = "top"), color = guide_legend(ncol = 3, title.position = "top"))
brca_p1_legend <-  cowplot::get_legend(brca_p1)
brca_p1 <- brca_p1 + theme(legend.position = "none")

brca_p2 <- ggplot(brca_scores, aes(x = k, y = cNMI, group = Method, color = Embedding, shape = Clustering)) + geom_line() + 
  theme_bw() + scale_colour_brewer(palette = "Dark2") + geom_point(size = plot_point_size) + facet_wrap(Approach ~ ., scales = "fixed", ncol = 3) + 
  theme(axis.text.x = element_blank(), axis.title.x = element_blank()) + 
  theme(legend.position = "none")

brca_p3 <- ggplot(brca_scores, aes(x = k, y = Module_score, group = Method, color = Embedding, shape = Clustering)) + geom_line() + 
  theme_bw() + scale_colour_brewer(palette = "Dark2") + geom_point(size = plot_point_size) + facet_wrap(Approach ~ ., scales = "fixed", ncol = 3) + 
  theme(axis.text.x = element_blank(), axis.title.x = element_blank()) + ylim(0, 1) + 
  theme(legend.position = "none")

brca_p4 <- ggplot(brca_scores, aes(x = k, y = ARI.BRCA_Subtype_PAM50, group = Method, color = Embedding, shape = Clustering)) + geom_line() + 
  theme_bw() + scale_colour_brewer(palette = "Dark2") + geom_point(size = plot_point_size) + facet_wrap(Approach ~ ., scales = "fixed", ncol = 3) + 
  theme(axis.text.x = element_blank(), axis.title.x = element_blank()) + 
  theme(legend.position = "none") + ylab("PAM50-subtype ARI")

brca_p5 <- ggplot(brca_scores, aes(x = k, y = ClusteringStabilityJaccard, group = Method, color = Embedding, shape = Clustering)) + geom_line() + 
  theme_bw() + scale_colour_brewer(palette = "Dark2") + geom_point(size = plot_point_size) + facet_wrap(Approach ~ ., scales = "fixed", ncol = 3) + 
  theme(axis.text.x = element_blank(), axis.title.x = element_blank()) + 
  theme(legend.position = "none") + ylab("Clustering Stability JI")

brca_p6 <- ggplot(brca_scores, aes(x = k, y = ClusteringStabilityARI, group = Method, color = Embedding, shape = Clustering)) + geom_line() + 
  theme_bw() + scale_colour_brewer(palette = "Dark2") + geom_point(size = plot_point_size) + facet_wrap(Approach ~ ., scales = "fixed", ncol = 3) + 
  theme(axis.text.x = element_blank(), axis.title.x = element_blank()) + 
  theme(legend.position = "none")  + ylab("Clustering Stability ARI")

brca_p7 <- ggplot(brca_scores, aes(x = k, y = NMI.BRCA_Subtype_PAM50, group = Method, color = Embedding, shape = Clustering)) + geom_line() + 
  theme_bw() + scale_colour_brewer(palette = "Dark2") + geom_point(size = plot_point_size) + facet_wrap(Approach ~ ., scales = "fixed", ncol = 3) + 
  theme(axis.text.x = element_blank(), axis.title.x = element_blank()) + 
  theme(legend.position = "none")  + ylab("PAM50-subtype NMI")

brca_p8 <- ggplot(brca_scores, aes(x = k, y = NMI.tss, group = Method, color = Embedding, shape = Clustering)) + geom_line() + 
  theme_bw() + scale_colour_brewer(palette = "Dark2") + geom_point(size = plot_point_size) + facet_wrap(Approach ~ ., scales = "fixed", ncol = 3) + 
  theme(axis.text.x = element_blank(), axis.title.x = element_blank()) + 
  theme(legend.position = "none")  + ylab("TSS NMI")

brca_p9 <- ggplot(brca_scores, aes(x = k, y = NMI.plate, group = Method, color = Embedding, shape = Clustering)) + geom_line() + 
  theme_bw() + scale_colour_brewer(palette = "Dark2") + geom_point(size = plot_point_size) + facet_wrap(Approach ~ ., scales = "fixed", ncol = 3) + 
  theme(axis.text.x = element_blank(), axis.title.x = element_blank()) + 
  theme(legend.position = "none")  + ylab("Plate NMI")

brca_bp_quantiles <- plyr::ddply(brca_scores_survival_sampled, c("Approach", "Embedding", "Clustering", "k"), 
                                 function(x) quantile(x$SurvivalPValue, probs = c(0, 0.025, 0.25, 0.5, 0.75, 0.975, 1), na.rm = TRUE))
colnames(brca_bp_quantiles)[5:11] <- c("Q0", "Q0025", "Q025", "Q05", "Q075", "Q0975", "Q1")
brca_bp_quantiles$IQR <- log(brca_bp_quantiles$Q075) - log(brca_bp_quantiles$Q025) 
brca_bp_quantiles$ymax <- apply(cbind(exp(log(brca_bp_quantiles$Q075) + brca_bp_quantiles$IQR * 1.5), brca_bp_quantiles$Q1), 1, min)
brca_bp_quantiles$ymin <- apply(cbind(exp(log(brca_bp_quantiles$Q025) - brca_bp_quantiles$IQR * 1.5), brca_bp_quantiles$Q0), 1, max)

brca_survival_pw <- ggplot(brca_bp_quantiles[brca_bp_quantiles$Embedding %in% c("GO", "KEGG", "REACTOME"),], 
                        aes(x = k, fill = Embedding)) + 
  geom_boxplot(aes(lower = Q025, upper = Q075, middle = Q05, ymin = ymin, ymax = ymax), 
               outlier.shape = NA, stat = "identity", lwd = 0.25) + 
  theme_bw() + scale_fill_brewer(palette = "Dark2") + 
  theme(legend.position = "bottom") + 
  facet_grid(Clustering ~ Approach, scales = "fixed") + 
  scale_y_continuous(trans = scales::trans_new("reverse_log", function(x) -log(x), function(y) exp(-y), breaks = scales::log_breaks()), 
                     limits = c(1.0, 5e-8))
brca_survival_pw_legend <-  cowplot::get_legend(brca_survival_pw)
brca_survival_pw <- brca_survival_pw + theme(legend.position = "none")

brca_survival_dr1 <- ggplot(brca_bp_quantiles[!brca_bp_quantiles$Embedding %in% c("GO", "KEGG", "REACTOME") &
                                                brca_bp_quantiles$Clustering %in% c("diana", "GMM", "kmeans"),], 
                           aes(x = k, fill = Embedding)) + 
  geom_boxplot(aes(lower = Q025, upper = Q075, middle = Q05, ymin = ymin, ymax = ymax), 
               outlier.shape = NA, stat = "identity", lwd = 0.25) + 
  theme_bw() + scale_fill_brewer(palette = "Set1") + 
  theme(legend.position = "bottom") + 
  facet_grid(Clustering ~ Approach, scales = "fixed") + 
  scale_y_continuous(trans = scales::trans_new("reverse_log", function(x) -log(x), function(y) exp(-y), breaks = scales::log_breaks()), 
                     limits = c(1.0, 5e-8))
brca_survival_dr1_legend <-  cowplot::get_legend(brca_survival_dr1)
brca_survival_dr1 <- brca_survival_dr1 + theme(legend.position = "none")

brca_survival_dr2 <- ggplot(brca_bp_quantiles[!brca_bp_quantiles$Embedding %in% c("GO", "KEGG", "REACTOME") &
                                                !brca_bp_quantiles$Clustering %in% c("diana", "GMM", "kmeans"),], 
                           aes(x = k, fill = Embedding)) + 
  geom_boxplot(aes(lower = Q025, upper = Q075, middle = Q05, ymin = ymin, ymax = ymax), 
               outlier.shape = NA, stat = "identity", lwd = 0.25) + 
  theme_bw() + scale_fill_brewer(palette = "Set1") + 
  theme(legend.position = "bottom") + 
  facet_grid(Clustering ~ Approach, scales = "fixed") + 
  scale_y_continuous(trans = scales::trans_new("reverse_log", function(x) -log(x), function(y) exp(-y), breaks = scales::log_breaks()), 
                     limits = c(1.0, 5e-8))
brca_survival_dr2 <- brca_survival_dr2 + theme(legend.position = "none")

### Prostate cancer

prad_p1 <- ggplot(prad_scores, aes(x = k, y = Silhouette, group = Method, color = Embedding, shape = Clustering)) + geom_line() + 
  theme_bw() + scale_colour_brewer(palette = "Dark2") + geom_point(size = plot_point_size) + facet_wrap(Approach ~ ., scales = "fixed", ncol = 3) + 
  ylim(0.05, 0.75) + 
  theme(legend.position="bottom")
prad_p1_legend <-  cowplot::get_legend(prad_p1)
prad_p1 <- prad_p1 + theme(legend.position = "none")

prad_p2 <- ggplot(prad_scores, aes(x = k, y = cNMI, group = Method, color = Embedding, shape = Clustering)) + geom_line() + 
  theme_bw() + scale_colour_brewer(palette = "Dark2") + geom_point(size = plot_point_size) + facet_wrap(Approach ~ ., scales = "fixed", ncol = 3) + 
  theme(legend.position = "none")

prad_p3 <- ggplot(prad_scores, aes(x = k, y = Module_score, group = Method, color = Embedding, shape = Clustering)) + geom_line() + 
  theme_bw() + scale_colour_brewer(palette = "Dark2") + geom_point(size = plot_point_size) + facet_wrap(Approach ~ ., scales = "fixed", ncol = 3) + 
  ylim(0, 1) + 
  theme(legend.position = "none")

prad_p4 <- ggplot(prad_scores, aes(x = k, y = ARI.Gleason_category, group = Method, color = Embedding, shape = Clustering)) + geom_line() + 
  theme_bw() + scale_colour_brewer(palette = "Dark2") + geom_point(size = plot_point_size) + facet_wrap(Approach ~ ., scales = "fixed", ncol = 3) + 
  theme(legend.position = "none") + ylab("Gleason cat. ARI")

prad_p5 <- ggplot(prad_scores, aes(x = k, y = ClusteringStabilityJaccard, group = Method, color = Embedding, shape = Clustering)) + geom_line() + 
  theme_bw() + scale_colour_brewer(palette = "Dark2") + geom_point(size = plot_point_size) + facet_wrap(Approach ~ ., scales = "fixed", ncol = 3) + 
  theme(legend.position = "none") + ylab("Clustering Stability JI")

prad_p6 <- ggplot(prad_scores, aes(x = k, y = ClusteringStabilityARI, group = Method, color = Embedding, shape = Clustering)) + geom_line() + 
  theme_bw() + scale_colour_brewer(palette = "Dark2") + geom_point(size = plot_point_size) + facet_wrap(Approach ~ ., scales = "fixed", ncol = 3) + 
  theme(legend.position = "none") + ylab("Clustering Stability ARI")

prad_p7 <- ggplot(prad_scores, aes(x = k, y = NMI.Gleason_category, group = Method, color = Embedding, shape = Clustering)) + geom_line() + 
  theme_bw() + scale_colour_brewer(palette = "Dark2") + geom_point(size = plot_point_size) + facet_wrap(Approach ~ ., scales = "fixed", ncol = 3) + 
  theme(legend.position = "none") + ylab("Gleason category NMI")

prad_p8 <- ggplot(prad_scores, aes(x = k, y = NMI.tss, group = Method, color = Embedding, shape = Clustering)) + geom_line() + 
  theme_bw() + scale_colour_brewer(palette = "Dark2") + geom_point(size = plot_point_size) + facet_wrap(Approach ~ ., scales = "fixed", ncol = 3) + 
  theme(legend.position = "none") + ylab("TSS NMI")

prad_p9 <- ggplot(prad_scores, aes(x = k, y = NMI.plate, group = Method, color = Embedding, shape = Clustering)) + geom_line() + 
  theme_bw() + scale_colour_brewer(palette = "Dark2") + geom_point(size = plot_point_size) + facet_wrap(Approach ~ ., scales = "fixed", ncol = 3) + 
  theme(legend.position = "none") + ylab("Plate NMI")

prad_bp_quantiles <- plyr::ddply(prad_scores_survival_sampled, c("Approach", "Embedding", "Clustering", "k"), 
                                 function(x) quantile(x$SurvivalPValue, probs = c(0, 0.025, 0.25, 0.5, 0.75, 0.975, 1), na.rm = TRUE))
colnames(prad_bp_quantiles)[5:11] <- c("Q0", "Q0025", "Q025", "Q05", "Q075", "Q0975", "Q1")
prad_bp_quantiles$IQR <- log(prad_bp_quantiles$Q075) - log(prad_bp_quantiles$Q025) 
prad_bp_quantiles$ymax <- apply(cbind(exp(log(prad_bp_quantiles$Q075) + prad_bp_quantiles$IQR * 1.5), prad_bp_quantiles$Q1), 1, min)
prad_bp_quantiles$ymin <- apply(cbind(exp(log(prad_bp_quantiles$Q025) - prad_bp_quantiles$IQR * 1.5), prad_bp_quantiles$Q0), 1, max)


prad_survival_pw <- ggplot(prad_bp_quantiles[prad_bp_quantiles$Embedding %in% c("GO", "KEGG", "REACTOME"),], 
              aes(x = k, fill = Embedding)) + 
  geom_boxplot(aes(lower = Q025, upper = Q075, middle = Q05, ymin = ymin, ymax = ymax), 
               outlier.shape = NA, stat = "identity", lwd = 0.25) + 
  theme_bw() + scale_fill_brewer(palette = "Dark2") + 
  theme(legend.position = "bottom") + 
  facet_grid(Clustering ~ Approach, scales = "fixed") + 
  scale_y_continuous(trans = scales::trans_new("reverse_log", function(x) -log(x), function(y) exp(-y), breaks = scales::log_breaks()), 
                     limits = c(1.0, 5e-6))
prad_survival_pw_legend <-  cowplot::get_legend(prad_survival_pw)
prad_survival_pw <- prad_survival_pw + theme(legend.position = "none")

prad_survival_dr1 <- ggplot(prad_bp_quantiles[!prad_bp_quantiles$Embedding %in% c("GO", "KEGG", "REACTOME") &
                                                prad_bp_quantiles$Clustering %in% c("diana", "GMM", "kmeans"),], 
                            aes(x = k, fill = Embedding)) + 
  geom_boxplot(aes(lower = Q025, upper = Q075, middle = Q05, ymin = ymin, ymax = ymax), 
               outlier.shape = NA, stat = "identity", lwd = 0.25) + 
  theme_bw() + scale_fill_brewer(palette = "Set1") + 
  theme(legend.position = "bottom") + 
  facet_grid(Clustering ~ Approach, scales = "fixed") + 
  scale_y_continuous(trans = scales::trans_new("reverse_log", function(x) -log(x), function(y) exp(-y), breaks = scales::log_breaks()), 
                     limits = c(1.0, 5e-6))
prad_survival_dr1_legend <-  cowplot::get_legend(prad_survival_dr1)
prad_survival_dr1 <- prad_survival_dr1 + theme(legend.position = "none")

prad_survival_dr2 <- ggplot(prad_bp_quantiles[!prad_bp_quantiles$Embedding %in% c("GO", "KEGG", "REACTOME") &
                                                !prad_bp_quantiles$Clustering %in% c("diana", "GMM", "kmeans"),], 
                            aes(x = k, fill = Embedding)) + 
  geom_boxplot(aes(lower = Q025, upper = Q075, middle = Q05, ymin = ymin, ymax = ymax), 
               outlier.shape = NA, stat = "identity", lwd = 0.25) + 
  theme_bw() + scale_fill_brewer(palette = "Set1") + 
  theme(legend.position = "bottom") + 
  facet_grid(Clustering ~ Approach, scales = "fixed") + 
  scale_y_continuous(trans = scales::trans_new("reverse_log", function(x) -log(x), function(y) exp(-y), breaks = scales::log_breaks()), 
                     limits = c(1.0, 5e-6))
prad_survival_dr2 <- prad_survival_dr2 + theme(legend.position = "none")

### Combined plots

if (save_plots_pdf) pdf(paste0(path_plots, "/silhouette.pdf"), width = plot_scale * 4, height = plot_scale * 5)
if (save_plots_svg) svg(paste0(path_plots, "/silhouette.svg"), width = plot_scale * 4, height = plot_scale * 5)
gridExtra::grid.arrange(brca_p1 + ggtitle("BRCA"), prad_p1 + ggtitle("PRAD"), 
                        brca_p1_legend, 
                        nrow = 3, heights = layout_heights)
if (save_plots_svg) dev.off()
if (save_plots_pdf) dev.off()

if (save_plots_pdf) pdf(paste0(path_plots, "/stability_ji.pdf"), width = plot_scale * 4, height = plot_scale * 5)
if (save_plots_svg) svg(paste0(path_plots, "/stability_ji.svg"), width = plot_scale * 4, height = plot_scale * 5)
gridExtra::grid.arrange(brca_p5 + ggtitle("BRCA") + ylab("Stability J"), prad_p5 + ggtitle("PRAD") + ylab("Stability J"), brca_p1_legend, 
                        nrow = 3, heights = layout_heights)
if (save_plots_svg) dev.off()
if (save_plots_pdf) dev.off()

if (save_plots_pdf) pdf(paste0(path_plots, "/cnmi.pdf"), width = plot_scale * 4, height = plot_scale * 5)
if (save_plots_svg) svg(paste0(path_plots, "/cnmi.svg"), width = plot_scale * 4, height = plot_scale * 5)
gridExtra::grid.arrange(brca_p2 + ggtitle("BRCA"), prad_p2 + ggtitle("PRAD"), brca_p1_legend, 
                        nrow = 3, heights = layout_heights)
if (save_plots_svg) dev.off()
if (save_plots_pdf) dev.off()

if (save_plots_pdf) pdf(paste0(path_plots, "/modules.pdf"), width = plot_scale * 4, height = plot_scale * 5)
if (save_plots_svg) svg(paste0(path_plots, "/modules.svg"), width = plot_scale * 4, height = plot_scale * 5)
gridExtra::grid.arrange(brca_p3 + ggtitle("BRCA"), prad_p3 + ggtitle("PRAD"), brca_p1_legend, 
                        nrow = 3, heights = layout_heights)
if (save_plots_svg) dev.off()
if (save_plots_pdf) dev.off()

if (save_plots_pdf) pdf(paste0(path_plots, "/survival_pw.pdf"), width = plot_scale * 4, height = plot_scale * 5)
if (save_plots_svg) svg(paste0(path_plots, "/survival_pw.svg"), width = plot_scale * 4, height = plot_scale * 5)
gridExtra::grid.arrange(brca_survival_pw + ggtitle("BRCA BK-CL") + ylab("Survival p-value"), 
                        prad_survival_pw + ggtitle("PRAD BK-CL") + ylab("Survival p-value"), 
                        brca_survival_pw_legend, nrow = 3, heights = c(4,4,0.5))#, widths = c(4,4.5,3))
if (save_plots_svg) dev.off()
if (save_plots_pdf) dev.off()

if (save_plots_pdf) pdf(paste0(path_plots, "/survival_dr.pdf"), width = plot_scale * 4, height = plot_scale * 5)
if (save_plots_svg) svg(paste0(path_plots, "/survival_dr.svg"), width = plot_scale * 4, height = plot_scale * 5)
gridExtra::grid.arrange(brca_survival_dr2 + ggtitle("BRCA DR HC") + ylab("Survival p-value"), 
                        brca_survival_dr1 + ggtitle("BRCA DR other") + theme(axis.text.y = element_blank()), 
                        prad_survival_dr2 + ggtitle("PRAD DR HC") + ylab("Survival p-value"), 
                        prad_survival_dr1 + ggtitle("PRAD DR other") + theme(axis.text.y = element_blank()), 
                        brca_survival_dr1_legend, nrow = 3, ncol = 2, 
                        layout_matrix = matrix(c(1,3,5,2,4,5), 3, 2), heights = c(4,4,0.5), widths = c(4.15,3.85))
if (save_plots_svg) dev.off()
if (save_plots_pdf) dev.off()

#########################################
# Supplementary
supplementary_plot_scale <- plot_scale

if (save_plots_pdf) pdf(paste0(path_plots, "/subtype_ari.pdf"), width = supplementary_plot_scale * 4, height = supplementary_plot_scale * 5)
if (save_plots_svg) svg(paste0(path_plots, "/subtype_ari.svg"), width = supplementary_plot_scale * 4, height = supplementary_plot_scale * 5)
gridExtra::grid.arrange(brca_p4 + ggtitle("BRCA"), prad_p4 + ggtitle("PRAD"), brca_p1_legend, 
                        nrow = 3, heights = layout_heights)
if (save_plots_svg) dev.off()
if (save_plots_pdf) dev.off()

if (save_plots_pdf) pdf(paste0(path_plots, "/stability_ari.pdf"), width = supplementary_plot_scale * 4, height = supplementary_plot_scale * 5)
if (save_plots_svg) svg(paste0(path_plots, "/stability_ari.svg"), width = supplementary_plot_scale * 4, height = supplementary_plot_scale * 5)
gridExtra::grid.arrange(brca_p6 + ggtitle("BRCA"), prad_p6 + ggtitle("PRAD"), brca_p1_legend, 
                        nrow = 3, heights = layout_heights)
if (save_plots_svg) dev.off()
if (save_plots_pdf) dev.off()


if (save_plots_pdf) pdf(paste0(path_plots, "/subtype_nmi.pdf"), width = supplementary_plot_scale * 4, height = supplementary_plot_scale * 5)
if (save_plots_svg) svg(paste0(path_plots, "/subtype_nmi.svg"), width = supplementary_plot_scale * 4, height = supplementary_plot_scale * 5)
gridExtra::grid.arrange(brca_p7 + ggtitle("BRCA"), prad_p7 + ggtitle("PRAD"), brca_p1_legend, 
                        nrow = 3, heights = layout_heights)
if (save_plots_svg) dev.off()
if (save_plots_pdf) dev.off()


if (save_plots_pdf) pdf(paste0(path_plots, "/batch_nmi.pdf"), width = supplementary_plot_scale * 4, height = supplementary_plot_scale * 5)
if (save_plots_svg) svg(paste0(path_plots, "/batch_nmi.svg"), width = supplementary_plot_scale * 4, height = supplementary_plot_scale * 5)
gridExtra::grid.arrange(brca_p8 + ggtitle("BRCA"), prad_p8 + ggtitle("PRAD"), brca_p1_legend, 
                        nrow = 3, heights = layout_heights)
if (save_plots_svg) dev.off()
if (save_plots_pdf) dev.off()

if (save_plots_pdf) pdf(paste0(path_plots, "/plate_nmi.pdf"), width = supplementary_plot_scale * 4, height = supplementary_plot_scale * 5)
if (save_plots_svg) svg(paste0(path_plots, "/plate_nmi.svg"), width = supplementary_plot_scale * 4, height = supplementary_plot_scale * 5)
gridExtra::grid.arrange(brca_p9 + ggtitle("BRCA"), prad_p9 + ggtitle("PRAD"), brca_p1_legend, 
                        nrow = 3, heights = layout_heights)
if (save_plots_svg) dev.off()
if (save_plots_pdf) dev.off()

### Pareto
# BRCA

#plot_palette <- "Dark2"
plot_palette <- "Set3"

silhouette_limits <- c(0, 0.9)
cnmi_limits <- c(0.475, 0.775)
survival_limits <- c(0, 1)
module_limits <- c(0, 1)
stability_limits <- c(0.7, 1)

plot11 <- ggplot(brca_pres, aes(Silhouette, cNMI)) + 
  geom_point(aes(color = Approach, shape = Clustering, size = k)) + 
  theme_bw() +
  scale_colour_brewer(palette = plot_palette) + theme(legend.position = "none") 
plot21 <- ggplot(brca_pres, aes(Silhouette, SurvivalPValue_score)) + 
  geom_point(aes(color = Approach, shape = Clustering, size = k)) + 
  theme_bw() + 
  scale_colour_brewer(palette = plot_palette) + 
  theme(legend.position = "none") 
plot31 <- ggplot(brca_pres, aes(Silhouette, Module_score)) + 
  geom_point(aes(color = Approach, shape = Clustering, size = k)) + 
  theme_bw() + 
  scale_colour_brewer(palette = plot_palette) + theme(legend.position = "none")
plot41 <- ggplot(brca_pres, aes(Silhouette, ClusteringStabilityJaccard)) + 
  geom_point(aes(color = Approach, shape = Clustering, size = k)) + 
  theme_bw() + 
  scale_colour_brewer(palette = plot_palette) + theme(legend.position = "none")

plot22 <- ggplot(brca_pres, aes(cNMI, SurvivalPValue_score)) + 
  geom_point(aes(color = Approach, shape = Clustering, size = k)) + 
  theme_bw() + 
  scale_colour_brewer(palette = plot_palette) +
  theme(legend.position = "none") 
plot32 <- ggplot(brca_pres, aes(cNMI, Module_score)) + 
  geom_point(aes(color = Approach, shape = Clustering, size = k)) + 
  theme_bw() + 
  scale_colour_brewer(palette = plot_palette) + 
  theme(legend.position = "none") 
plot42 <- ggplot(brca_pres, aes(cNMI, ClusteringStabilityJaccard)) + 
  geom_point(aes(color = Approach, shape = Clustering, size = k)) + 
  theme_bw() + 
  scale_colour_brewer(palette = plot_palette) + 
  theme(legend.position = "none") 

plot33 <- ggplot(brca_pres, aes(SurvivalPValue_score, Module_score)) + 
  geom_point(aes(color = Approach, shape = Clustering, size = k)) + 
  theme_bw() + 
  scale_colour_brewer(palette = plot_palette) + 
  theme(legend.position = "none") 
plot43 <- ggplot(brca_pres, aes(SurvivalPValue_score, ClusteringStabilityJaccard)) + 
  geom_point(aes(color = Approach, shape = Clustering, size = k)) + 
  theme_bw() + 
  scale_colour_brewer(palette = plot_palette) + 
  theme(legend.position = "none") 

plot44 <- ggplot(brca_pres, aes(Module_score, ClusteringStabilityJaccard)) + 
  geom_point(aes(color = Approach, shape = Clustering, size = k)) + 
  theme_bw() + 
  scale_colour_brewer(palette = plot_palette) + 
  theme(legend.position = "none") 

legend_plot <- ggplot(brca_pres, aes(Silhouette, cNMI)) + 
  geom_point(aes(color = Approach, shape = Clustering, size = k)) + 
  scale_colour_brewer(palette = plot_palette) + theme_bw() + 
  theme(legend.box = "horizontal") + 
  guides(shape = guide_legend(ncol = 1, order = 2), 
         color = guide_legend(ncol = 2, order = 1),
         size = guide_legend(ncol = 1, order = 3))
pareto_legend <- cowplot::get_legend(legend_plot)

layout <- rbind(c(1,2,3,4), c(NA,5,6,7), c(NA,NA,8,9), c(11,11,11,10)) 

if (save_plots_pdf) pdf(paste0(path_plots, "/brca_pareto.pdf"), width = plot_scale * 7, height = plot_scale * 6)
if (save_plots_svg) svg(paste0(path_plots, "/brca_pareto.svg"), width = plot_scale * 6, height = plot_scale * 6)
gridExtra::grid.arrange(plot11, plot21, plot31, plot41, plot22, plot32, plot42, plot33, plot43, plot44, pareto_legend, layout_matrix = layout)
if (save_plots_svg) dev.off()
if (save_plots_pdf) dev.off()

# PRAD

#plot_palette <- scale_colour_brewer(palette = "Dark2")
plot_palette <- scale_colour_brewer(palette = "Set3")
#plot_palette <- scale_color_viridis_d()

silhouette_limits <- c(0, 0.9)
cnmi_limits <- c(0.475, 0.775)
survival_limits <- c(0, 1)
module_limits <- c(0, 1)
stability_limits <- c(0.7, 1)

plot11 <- ggplot(prad_pres, aes(Silhouette, cNMI)) + 
  geom_point(aes(color = Approach, shape = Clustering, size = k)) + 
  theme_bw() +
  plot_palette + theme(legend.position = "none") 
plot21 <- ggplot(prad_pres, aes(Silhouette, SurvivalPValue_score)) + 
  geom_point(aes(color = Approach, shape = Clustering, size = k)) + 
  theme_bw() + 
  plot_palette + 
  theme(legend.position = "none") 
plot31 <- ggplot(prad_pres, aes(Silhouette, Module_score)) + 
  geom_point(aes(color = Approach, shape = Clustering, size = k)) + 
  theme_bw() + 
  plot_palette + theme(legend.position = "none")
plot41 <- ggplot(prad_pres, aes(Silhouette, ClusteringStabilityJaccard)) + 
  geom_point(aes(color = Approach, shape = Clustering, size = k)) + 
  theme_bw() + 
  plot_palette + theme(legend.position = "none")

plot22 <- ggplot(prad_pres, aes(cNMI, SurvivalPValue_score)) + 
  geom_point(aes(color = Approach, shape = Clustering, size = k)) + 
  theme_bw() + 
  plot_palette +
  theme(legend.position = "none") 
plot32 <- ggplot(prad_pres, aes(cNMI, Module_score)) + 
  geom_point(aes(color = Approach, shape = Clustering, size = k)) + 
  theme_bw() + 
  plot_palette + 
  theme(legend.position = "none") 
plot42 <- ggplot(prad_pres, aes(cNMI, ClusteringStabilityJaccard)) + 
  geom_point(aes(color = Approach, shape = Clustering, size = k)) + 
  theme_bw() + 
  plot_palette + 
  theme(legend.position = "none") 

plot33 <- ggplot(prad_pres, aes(SurvivalPValue_score, Module_score)) + 
  geom_point(aes(color = Approach, shape = Clustering, size = k)) + 
  theme_bw() + 
  plot_palette + 
  theme(legend.position = "none") 
plot43 <- ggplot(prad_pres, aes(SurvivalPValue_score, ClusteringStabilityJaccard)) + 
  geom_point(aes(color = Approach, shape = Clustering, size = k)) + 
  theme_bw() + 
  plot_palette + 
  theme(legend.position = "none") 

plot44 <- ggplot(prad_pres, aes(Module_score, ClusteringStabilityJaccard)) + 
  geom_point(aes(color = Approach, shape = Clustering, size = k)) + 
  theme_bw() + 
  plot_palette + 
  theme(legend.position = "none") 

legend_plot <- ggplot(prad_pres, aes(Silhouette, cNMI)) + 
  geom_point(aes(color = Approach, shape = Clustering, size = k)) + 
  plot_palette + theme_bw() + 
  theme(legend.box = "horizontal") + 
  guides(shape = guide_legend(ncol = 1, order = 2), 
         color = guide_legend(ncol = 2, order = 1),
         size = guide_legend(ncol = 1, order = 3))
pareto_legend <- cowplot::get_legend(legend_plot)

layout <- rbind(c(1,2,3,4), c(NA,5,6,7), c(NA,NA,8,9), c(11,11,11,10)) 

if (save_plots_pdf) pdf(paste0(path_plots, "/prad_pareto.pdf"), width = plot_scale * 7, height = plot_scale * 6)
if (save_plots_svg) svg(paste0(path_plots, "/prad_pareto.svg"), width = plot_scale * 6, height = plot_scale * 6)
gridExtra::grid.arrange(plot11, plot21, plot31, plot41, plot22, plot32, plot42, plot33, plot43, plot44, pareto_legend, layout_matrix = layout)
if (save_plots_svg) dev.off()
if (save_plots_pdf) dev.off()



