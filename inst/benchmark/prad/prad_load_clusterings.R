# Read clustering results from intermediary files

#prad_clust_dimred <- read.csv("~/tcga/prad/intermediary_files/dimred/clusters.csv.gz", header = TRUE, row.names = 1)

prad_clust_dimred_pca <- read.csv(paste0(path_intermediate_results, "/prad/dimred/dimred_pca/clusters.csv.gz"), header = TRUE, row.names = 1)
prad_clust_dimred_pca$Approach <- "DR*"
prad_clust_dimred_tsne <- read.csv(paste0(path_intermediate_results, "/prad/dimred/dimred_tsne/clusters.csv.gz"), header = TRUE, row.names = 1)
prad_clust_dimred_tsne$Approach <- "DR*"
prad_clust_dimred_umap_10n <- read.csv(paste0(path_intermediate_results, "/prad/dimred/dimred_umap_10n/clusters.csv.gz"), header = TRUE, row.names = 1)
prad_clust_dimred_umap_10n$Approach <- "DR*"
prad_clust_dimred_umap_10n$drname <- paste0(prad_clust_dimred_umap_10n$drname, "_10n")
prad_clust_dimred_umap_20n <- read.csv(paste0(path_intermediate_results, "/prad/dimred/dimred_umap_20n/clusters.csv.gz"), header = TRUE, row.names = 1)
prad_clust_dimred_umap_20n$Approach <- "DR*"
prad_clust_dimred_umap_20n$drname <- paste0(prad_clust_dimred_umap_20n$drname, "_20n")
prad_clust_dimred_umap_30n <- read.csv(paste0(path_intermediate_results, "/prad/dimred/dimred_umap_30n/clusters.csv.gz"), header = TRUE, row.names = 1)
prad_clust_dimred_umap_30n$Approach <- "DR*"
prad_clust_dimred_umap_30n$drname <- paste0(prad_clust_dimred_umap_30n$drname, "_30n")

prad_clust_dimred_otp <- read.csv(paste0(path_intermediate_results, "/prad/dimred_otp/clusters.csv.gz"), header = TRUE, row.names = 1)
prad_clust_dimred_otp$Approach <- "DR"
prad_clust_dimred_otp_umap10n <- read.csv(paste0(path_intermediate_results, "/prad/dimred_otp/dimred_umap_10n/clusters.csv.gz"), header = TRUE, row.names = 1)
prad_clust_dimred_otp_umap10n$Approach <- "DR"
prad_clust_dimred_otp_umap10n$drname <- paste0(prad_clust_dimred_otp_umap10n$drname, "_10n")
prad_clust_dimred_otp_umap20n <- read.csv(paste0(path_intermediate_results, "/prad/dimred_otp/dimred_umap_20n/clusters.csv.gz"), header = TRUE, row.names = 1)
prad_clust_dimred_otp_umap20n$Approach <- "DR"
prad_clust_dimred_otp_umap20n$drname <- paste0(prad_clust_dimred_otp_umap20n$drname, "_20n")
prad_clust_dimred_otp_umap30n <- read.csv(paste0(path_intermediate_results, "/prad/dimred_otp/dimred_umap_30n/clusters.csv.gz"), header = TRUE, row.names = 1)
prad_clust_dimred_otp_umap30n$Approach <- "DR"
prad_clust_dimred_otp_umap30n$drname <- paste0(prad_clust_dimred_otp_umap30n$drname, "_30n")

prad_clust_pw_otp_diffrank <- read.csv(paste0(path_intermediate_results, "/prad/diffrank_otp/clusters.csv.gz"), header = TRUE, row.names = 1)
prad_clust_pw_otp_diffrank$Approach <- "DiffRank"
prad_clust_pw_otp_gsva <- read.csv(paste0(path_intermediate_results, "/prad/gsva_otp/clusters.csv.gz"), header = TRUE, row.names = 1)
prad_clust_pw_otp_gsva$Approach <- "GSVA"

prad_clust_pw_otp_gcn_rwr_fgsea <- read.csv(paste0(path_intermediate_results, "/prad/rwr/gcn_rwr_otp/rwr.pw.clusters.csv.gz"), header = TRUE, row.names = 1)
prad_clust_pw_otp_gcn_rwr_fgsea$Approach <- "GCN RWR-FGSEA"
prad_clust_pw_otp_ppi_rwr_fgsea <- read.csv(paste0(path_intermediate_results, "/prad/rwr/ppi_rwr_otp/rwr.pw.clusters.csv.gz"), header = TRUE, row.names = 1)
prad_clust_pw_otp_ppi_rwr_fgsea$Approach <- "PPI RWR-FGSEA"

# VAE clusters
#prad_clust_vae <- read.csv("~/tcga/prad/intermediary_files/vae/clusters.csv.gz", header = TRUE, row.names = 1)
#prad_clust_vae$Approach <- "DR"
#prad_clust_vae$drname <- "VAE"
#prad_clust_vae$run <- prad_clust_vae$run + 1
#prad_clust_vae$fold <- prad_clust_vae$fold + 1
