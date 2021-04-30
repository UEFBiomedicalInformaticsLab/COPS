source("load_config.R")

# Read scores from files
#brca_res_dimred <- read.csv("~/tcga/brca/intermediary_files/dimred/scores.csv", header = TRUE, row.names = 1)

brca_res_dimred_pca <- read.csv(paste0(path_intermediate_results, "/brca/dimred/dimred_pca/scores.csv"), header = TRUE, row.names = 1)
brca_res_dimred_pca$Approach <- "DR*"
brca_res_dimred_tsne <- read.csv(paste0(path_intermediate_results, "/brca/dimred/dimred_tsne/scores.csv"), header = TRUE, row.names = 1)
brca_res_dimred_tsne$Approach <- "DR*"
brca_res_dimred_umap_10n <- read.csv(paste0(path_intermediate_results, "/brca/dimred/dimred_umap_10n/scores.csv"), header = TRUE, row.names = 1)
brca_res_dimred_umap_10n$Approach <- "DR*"
brca_res_dimred_umap_10n$drname <- paste0(brca_res_dimred_umap_10n$drname, "_10n")
brca_res_dimred_umap_20n <- read.csv(paste0(path_intermediate_results, "/brca/dimred/dimred_umap_20n/scores.csv"), header = TRUE, row.names = 1)
brca_res_dimred_umap_20n$Approach <- "DR*"
brca_res_dimred_umap_20n$drname <- paste0(brca_res_dimred_umap_20n$drname, "_20n")
brca_res_dimred_umap_30n <- read.csv(paste0(path_intermediate_results, "/brca/dimred/dimred_umap_30n/scores.csv"), header = TRUE, row.names = 1)
brca_res_dimred_umap_30n$Approach <- "DR*"
brca_res_dimred_umap_30n$drname <- paste0(brca_res_dimred_umap_30n$drname, "_30n")

brca_res_dimred_otp <- read.csv(paste0(path_intermediate_results, "/brca/dimred_otp/scores.csv"), header = TRUE, row.names = 1)
brca_res_dimred_otp$Approach <- "DR"
brca_res_dimred_otp_umap10n <- read.csv(paste0(path_intermediate_results, "/brca/dimred_otp/dimred_umap_10n/scores.csv"), header = TRUE, row.names = 1)
brca_res_dimred_otp_umap10n$Approach <- "DR"
brca_res_dimred_otp_umap10n$drname <- paste0(brca_res_dimred_otp_umap10n$drname, "_10n")
brca_res_dimred_otp_umap20n <- read.csv(paste0(path_intermediate_results, "/brca/dimred_otp/dimred_umap_20n/scores.csv"), header = TRUE, row.names = 1)
brca_res_dimred_otp_umap20n$Approach <- "DR"
brca_res_dimred_otp_umap20n$drname <- paste0(brca_res_dimred_otp_umap20n$drname, "_20n")
brca_res_dimred_otp_umap30n <- read.csv(paste0(path_intermediate_results, "/brca/dimred_otp/dimred_umap_30n/scores.csv"), header = TRUE, row.names = 1)
brca_res_dimred_otp_umap30n$Approach <- "DR"
brca_res_dimred_otp_umap30n$drname <- paste0(brca_res_dimred_otp_umap30n$drname, "_30n")

brca_res_pw_otp_diffrank <- read.csv(paste0(path_intermediate_results, "/brca/diffrank_otp/scores.csv"), header = TRUE, row.names = 1)
brca_res_pw_otp_diffrank$Approach <- "DiffRank"
brca_res_pw_otp_gsva <- read.csv(paste0(path_intermediate_results, "/brca/gsva_otp/scores.csv"), header = TRUE, row.names = 1)
brca_res_pw_otp_gsva$Approach <- "GSVA"

brca_res_pw_otp_gcn_rwr_fgsea <- read.csv(paste0(path_intermediate_results, "/brca/rwr/gcn_rwr_otp/rwr.pw.scores.csv"), header = TRUE, row.names = 1)
brca_res_pw_otp_gcn_rwr_fgsea$Approach <- "GCN RWR-FGSEA"
brca_res_pw_otp_ppi_rwr_fgsea <- read.csv(paste0(path_intermediate_results, "/brca/rwr/ppi_rwr_otp/rwr.pw.scores.csv"), header = TRUE, row.names = 1)
brca_res_pw_otp_ppi_rwr_fgsea$Approach <- "PPI RWR-FGSEA"

# VAE results
#brca_res_vae <- read.csv("~/tcga/brca/intermediary_files/vae/scores.csv", header = TRUE, row.names = 1)
#brca_res_vae$Approach <- "DR"
#brca_res_vae$drname <- "VAE"
#brca_res_vae$run <- brca_res_vae$run + 1
#brca_res_vae$fold <- brca_res_vae$fold + 1

brca_scores <- plyr::rbind.fill(brca_res_dimred_pca, 
                                brca_res_dimred_tsne, 
                                brca_res_dimred_umap_10n, 
                                brca_res_dimred_umap_20n, 
                                brca_res_dimred_umap_30n,
                                brca_res_dimred_otp,
                                brca_res_dimred_otp_umap10n,
                                brca_res_dimred_otp_umap20n,
                                brca_res_dimred_otp_umap30n,
                                brca_res_pw_otp_diffrank,
                                brca_res_pw_otp_gsva,
                                brca_res_pw_otp_gcn_rwr_fgsea,
                                brca_res_pw_otp_ppi_rwr_fgsea)


source("article/plot_renaming.R")

brca_scores <- rename_methods(brca_scores)
#brca_scores$Approach <- factor(brca_scores$Approach, unique(brca_scores$Approach)[c(1,3,5,2,4,6)])
brca_scores$cNMI <- (brca_scores$NMI.BRCA_Subtype_PAM50 + 1 - brca_scores$NMI.tss) / 2

summary_by <- c("datname", "drname", "k", "m", "Approach", "Method", "Embedding", "Clustering")

# Separate survival data for box-plots
brca_scores_survival_sampled <- brca_scores[c(summary_by, "SurvivalPValue")]

# Prepare averages
summary_function <- function(x) {
  out <- x[1, summary_by]
  for (i in colnames(x)[!colnames(x) %in% summary_by]) {
    out[[i]] <- mean(x[[i]], na.rm = TRUE)
  }
  return(out)
}

brca_scores <- plyr::ddply(brca_scores, summary_by, summary_function)
brca_scores$fold <- NULL
brca_scores$run <- NULL
