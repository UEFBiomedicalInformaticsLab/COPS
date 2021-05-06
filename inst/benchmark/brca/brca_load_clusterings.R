# Read clustering results from intermediary files

#brca_clust_dimred <- read.csv("~/tcga/brca/intermediary_files/dimred/clusters.csv.gz", header = TRUE, row.names = 1)

brca_clust_dimred_pca <- read.csv(paste0(path_intermediate_results, "/brca/dimred/dimred_pca/clusters.csv.gz"), header = TRUE, row.names = 1)
brca_clust_dimred_pca$Approach <- "DR*"
brca_clust_dimred_tsne <- read.csv(paste0(path_intermediate_results, "/brca/dimred/dimred_tsne/clusters.csv.gz"), header = TRUE, row.names = 1)
brca_clust_dimred_tsne$Approach <- "DR*"
brca_clust_dimred_umap_10n <- read.csv(paste0(path_intermediate_results, "/brca/dimred/dimred_umap_10n/clusters.csv.gz"), header = TRUE, row.names = 1)
brca_clust_dimred_umap_10n$Approach <- "DR*"
brca_clust_dimred_umap_10n$drname <- paste0(brca_clust_dimred_umap_10n$drname, "_10n")
brca_clust_dimred_umap_20n <- read.csv(paste0(path_intermediate_results, "/brca/dimred/dimred_umap_20n/clusters.csv.gz"), header = TRUE, row.names = 1)
brca_clust_dimred_umap_20n$Approach <- "DR*"
brca_clust_dimred_umap_20n$drname <- paste0(brca_clust_dimred_umap_20n$drname, "_20n")
brca_clust_dimred_umap_30n <- read.csv(paste0(path_intermediate_results, "/brca/dimred/dimred_umap_30n/clusters.csv.gz"), header = TRUE, row.names = 1)
brca_clust_dimred_umap_30n$Approach <- "DR*"
brca_clust_dimred_umap_30n$drname <- paste0(brca_clust_dimred_umap_30n$drname, "_30n")

brca_clust_dimred_otp <- read.csv(paste0(path_intermediate_results, "/brca/dimred_otp/clusters.csv.gz"), header = TRUE, row.names = 1)
brca_clust_dimred_otp$Approach <- "DR"
brca_clust_dimred_otp_umap10n <- read.csv(paste0(path_intermediate_results, "/brca/dimred_otp/dimred_umap_10n/clusters.csv.gz"), header = TRUE, row.names = 1)
brca_clust_dimred_otp_umap10n$Approach <- "DR"
brca_clust_dimred_otp_umap10n$drname <- paste0(brca_clust_dimred_otp_umap10n$drname, "_10n")
brca_clust_dimred_otp_umap20n <- read.csv(paste0(path_intermediate_results, "/brca/dimred_otp/dimred_umap_20n/clusters.csv.gz"), header = TRUE, row.names = 1)
brca_clust_dimred_otp_umap20n$Approach <- "DR"
brca_clust_dimred_otp_umap20n$drname <- paste0(brca_clust_dimred_otp_umap20n$drname, "_20n")
brca_clust_dimred_otp_umap30n <- read.csv(paste0(path_intermediate_results, "/brca/dimred_otp/dimred_umap_30n/clusters.csv.gz"), header = TRUE, row.names = 1)
brca_clust_dimred_otp_umap30n$Approach <- "DR"
brca_clust_dimred_otp_umap30n$drname <- paste0(brca_clust_dimred_otp_umap30n$drname, "_30n")

brca_clust_pw_otp_diffrank <- read.csv(paste0(path_intermediate_results, "/brca/diffrank_otp/clusters.csv.gz"), header = TRUE, row.names = 1)
brca_clust_pw_otp_diffrank$Approach <- "DiffRank"
brca_clust_pw_otp_gsva <- read.csv(paste0(path_intermediate_results, "/brca/gsva_otp/clusters.csv.gz"), header = TRUE, row.names = 1)
brca_clust_pw_otp_gsva$Approach <- "GSVA"

brca_clust_pw_otp_gcn_rwr_fgsea <- read.csv(paste0(path_intermediate_results, "/brca/rwr/gcn_rwr_otp/rwr.pw.clusters.csv.gz"), header = TRUE, row.names = 1)
brca_clust_pw_otp_gcn_rwr_fgsea$Approach <- "GCN RWR-FGSEA"
brca_clust_pw_otp_ppi_rwr_fgsea <- read.csv(paste0(path_intermediate_results, "/brca/rwr/ppi_rwr_otp/rwr.pw.clusters.csv.gz"), header = TRUE, row.names = 1)
brca_clust_pw_otp_ppi_rwr_fgsea$Approach <- "PPI RWR-FGSEA"

# VAE clusters
#brca_clust_vae <- read.csv("~/tcga/brca/intermediary_files/vae/clusters.csv.gz", header = TRUE, row.names = 1)
#brca_clust_vae$Approach <- "DR"
#brca_clust_vae$drname <- "VAE"
#brca_clust_vae$run <- brca_clust_vae$run + 1
#brca_clust_vae$fold <- brca_clust_vae$fold + 1
