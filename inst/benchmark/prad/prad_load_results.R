#source("load_config.R")

# Read scores from files
prad_res_dimred_pca <- read.csv(paste0(path_intermediate_results, "/prad/dimred/dimred_pca/scores.csv"), header = TRUE, row.names = 1)
prad_res_dimred_pca$Approach <- "DR*"
prad_res_dimred_tsne <- read.csv(paste0(path_intermediate_results, "/prad/dimred/dimred_tsne/scores.csv"), header = TRUE, row.names = 1)
prad_res_dimred_tsne$Approach <- "DR*"
prad_res_dimred_umap_10n <- read.csv(paste0(path_intermediate_results, "/prad/dimred/dimred_umap_10n/scores.csv"), header = TRUE, row.names = 1)
prad_res_dimred_umap_10n$Approach <- "DR*"
prad_res_dimred_umap_10n$drname <- paste0(prad_res_dimred_umap_10n$drname, "_10n")
prad_res_dimred_umap_20n <- read.csv(paste0(path_intermediate_results, "/prad/dimred/dimred_umap_20n/scores.csv"), header = TRUE, row.names = 1)
prad_res_dimred_umap_20n$Approach <- "DR*"
prad_res_dimred_umap_20n$drname <- paste0(prad_res_dimred_umap_20n$drname, "_20n")
prad_res_dimred_umap_30n <- read.csv(paste0(path_intermediate_results, "/prad/dimred/dimred_umap_30n/scores.csv"), header = TRUE, row.names = 1)
prad_res_dimred_umap_30n$Approach <- "DR*"
prad_res_dimred_umap_30n$drname <- paste0(prad_res_dimred_umap_30n$drname, "_30n")

prad_res_dimred_otp <- read.csv(paste0(path_intermediate_results, "/prad/dimred_otp/scores.csv"), header = TRUE, row.names = 1)
prad_res_dimred_otp$Approach <- "DR"
prad_res_dimred_otp_umap10n <- read.csv(paste0(path_intermediate_results, "/prad/dimred_otp/dimred_umap_10n/scores.csv"), header = TRUE, row.names = 1)
prad_res_dimred_otp_umap10n$Approach <- "DR"
prad_res_dimred_otp_umap10n$drname <- paste0(prad_res_dimred_otp_umap10n$drname, "_10n")
prad_res_dimred_otp_umap20n <- read.csv(paste0(path_intermediate_results, "/prad/dimred_otp/dimred_umap_20n/scores.csv"), header = TRUE, row.names = 1)
prad_res_dimred_otp_umap20n$Approach <- "DR"
prad_res_dimred_otp_umap20n$drname <- paste0(prad_res_dimred_otp_umap20n$drname, "_20n")
prad_res_dimred_otp_umap30n <- read.csv(paste0(path_intermediate_results, "/prad/dimred_otp/dimred_umap_30n/scores.csv"), header = TRUE, row.names = 1)
prad_res_dimred_otp_umap30n$Approach <- "DR"
prad_res_dimred_otp_umap30n$drname <- paste0(prad_res_dimred_otp_umap30n$drname, "_30n")

prad_res_pw_otp_diffrank <- read.csv(paste0(path_intermediate_results, "/prad/diffrank_otp/scores.csv"), header = TRUE, row.names = 1)
prad_res_pw_otp_diffrank$Approach <- "DiffRank"
prad_res_pw_otp_gsva <- read.csv(paste0(path_intermediate_results, "/prad/gsva_otp/scores.csv"), header = TRUE, row.names = 1)
prad_res_pw_otp_gsva$Approach <- "GSVA"

prad_res_pw_otp_gcn_rwr_fgsea <- read.csv(paste0(path_intermediate_results, "/prad/rwr/gcn_rwr_otp/rwr.pw.scores.csv"), header = TRUE, row.names = 1)
prad_res_pw_otp_gcn_rwr_fgsea$Approach <- "GCN RWR-FGSEA"
prad_res_pw_otp_ppi_rwr_fgsea <- read.csv(paste0(path_intermediate_results, "/prad/rwr/ppi_rwr_otp/rwr.pw.scores.csv"), header = TRUE, row.names = 1)
prad_res_pw_otp_ppi_rwr_fgsea$Approach <- "PPI RWR-FGSEA"

# VAE results
prad_res_vae <- read.csv("results/prad_vae_scores.csv", header = TRUE, row.names = 1)
prad_res_vae$Approach <- "DR*"
prad_res_vae$drname <- "VAE"
prad_res_vae$run <- prad_res_vae$run + 1
prad_res_vae$fold <- prad_res_vae$fold + 1

prad_scores <- plyr::rbind.fill(prad_res_dimred_pca, 
                                prad_res_dimred_tsne, 
                                prad_res_dimred_umap_10n, 
                                prad_res_dimred_umap_20n, 
                                prad_res_dimred_umap_30n,
                                prad_res_dimred_otp,
                                prad_res_dimred_otp_umap10n,
                                prad_res_dimred_otp_umap20n,
                                prad_res_dimred_otp_umap30n,
                                prad_res_pw_otp_diffrank,
                                prad_res_pw_otp_gsva,
                                prad_res_pw_otp_gcn_rwr_fgsea,
                                prad_res_pw_otp_ppi_rwr_fgsea,
                                prad_res_vae)


source("article/plot_renaming.R")

prad_scores <- rename_methods(prad_scores)
#prad_scores$Approach <- factor(prad_scores$Approach, unique(prad_scores$Approach)[c(1,3,5,2,4,6)])
prad_scores$cNMI <- (prad_scores$NMI.Gleason_category + 1 - prad_scores$NMI.tss) / 2

summary_by <- c("datname", "drname", "k", "m", "Approach", "Method", "Embedding", "Clustering")

# Separate survival data for box-plots
prad_scores_survival_sampled <- prad_scores[c(summary_by, "SurvivalPValue")]

# Prepare averages
summary_function <- function(x) {
  out <- x[1, summary_by]
  for (i in colnames(x)[!colnames(x) %in% summary_by]) {
    out[[i]] <- mean(x[[i]], na.rm = TRUE)
  }
  return(out)
}

prad_scores <- plyr::ddply(prad_scores, summary_by, summary_function)
prad_scores$fold <- NULL
prad_scores$run <- NULL
