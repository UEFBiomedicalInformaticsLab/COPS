#!/bin/bash
# This runs all the R scripts to generate results included in the patient stratification article
# Scripts and results were put into seperate files in order to conserve memory

source config.sh

mkdir -p $PATH_PLOTS

echo TCGA breast cancer
# Create folders for intermediate results
mkdir -p $PATH_INTERMEDIATES
mkdir -p $PATH_INTERMEDIATES/brca/wgcna
mkdir -p $PATH_INTERMEDIATES/brca/dimred/dimred_umap_10n
mkdir -p $PATH_INTERMEDIATES/brca/dimred/dimred_umap_20n
mkdir -p $PATH_INTERMEDIATES/brca/dimred/dimred_umap_30n
mkdir -p $PATH_INTERMEDIATES/brca/dimred/dimred_pca
mkdir -p $PATH_INTERMEDIATES/brca/dimred/dimred_tsne
mkdir -p $PATH_INTERMEDIATES/brca/dimred_otp
mkdir -p $PATH_INTERMEDIATES/brca/dimred_otp/dimred_umap_10n
mkdir -p $PATH_INTERMEDIATES/brca/dimred_otp/dimred_umap_20n
mkdir -p $PATH_INTERMEDIATES/brca/dimred_otp/dimred_umap_30n
mkdir -p $PATH_INTERMEDIATES/brca/diffrank_otp
mkdir -p $PATH_INTERMEDIATES/brca/gsva_otp
mkdir -p $PATH_INTERMEDIATES/brca/rwr/gcn_rwr_otp
mkdir -p $PATH_INTERMEDIATES/brca/rwr/ppi_rwr_otp

echo Module detection
Rscript brca/brca_gene_modules.R --no-restore --no-save

echo Clustering
echo No gene filters DR
echo PCA
Rscript brca/brca_dimred_pca.R --no-restore --no-save
echo t-SNE
Rscript brca/brca_dimred_tsne.R --no-restore --no-save
echo UMAP 10 neighbours
Rscript brca/brca_dimred_umap_10n.R --no-restore --no-save
echo UMAP 20 neighbours
Rscript brca/brca_dimred_umap_20n.R --no-restore --no-save
echo UMAP 30 neighbours
Rscript brca/brca_dimred_umap_30n.R --no-restore --no-save
echo DEG and OTP intersection DR
echo PCA and t-SNE
Rscript brca/brca_otp_dimred.R --no-restore --no-save
echo UMAP
Rscript brca/brca_otp_dimred_umap.R --no-restore --no-save
echo DEG and OTP intersection GSVA
Rscript brca/brca_otp_gsva.R --no-restore --no-save
echo DEG and OTP intersection DiffRank
Rscript brca/brca_otp_diffrank.R --no-restore --no-save
echo DEG and OTP intersection RWRFGSEA on GCN
Rscript brca/rwr/brca_gcn_rwr_otp.R --no-restore --no-save
Rscript --no-restore --no-save brca/rwr/rwr_clustering_only.R \
  $PATH_INTERMEDIATES/brca/rwr/gcn_rwr_otp/ \
  rwr.pw.scores.csv rwr.pw.clusters.csv.gz \
  KEGG.csv.gz GO.csv.gz REACTOME.csv.gz
echo DEG and OTP intersection RWRFGSEA on PPI
Rscript brca/rwr/brca_ppi_rwr_otp.R --no-restore --no-save
Rscript --no-restore --no-save brca/rwr/rwr_clustering_only.R \
  $PATH_INTERMEDIATES/brca/rwr/ppi_rwr_otp/ \
  rwr.pw.scores.csv rwr.pw.clusters.csv.gz \
  KEGG.csv.gz GO.csv.gz REACTOME.csv.gz

echo TCGA prostate cancer
# Create folders for intermediate results
mkdir -p $PATH_INTERMEDIATES
mkdir -p $PATH_INTERMEDIATES/prad/wgcna
mkdir -p $PATH_INTERMEDIATES/prad/dimred/dimred_umap_10n
mkdir -p $PATH_INTERMEDIATES/prad/dimred/dimred_umap_20n
mkdir -p $PATH_INTERMEDIATES/prad/dimred/dimred_umap_30n
mkdir -p $PATH_INTERMEDIATES/prad/dimred/dimred_pca
mkdir -p $PATH_INTERMEDIATES/prad/dimred/dimred_tsne
mkdir -p $PATH_INTERMEDIATES/prad/dimred_otp
mkdir -p $PATH_INTERMEDIATES/prad/dimred_otp/dimred_umap_10n
mkdir -p $PATH_INTERMEDIATES/prad/dimred_otp/dimred_umap_20n
mkdir -p $PATH_INTERMEDIATES/prad/dimred_otp/dimred_umap_30n
mkdir -p $PATH_INTERMEDIATES/prad/diffrank_otp
mkdir -p $PATH_INTERMEDIATES/prad/gsva_otp
mkdir -p $PATH_INTERMEDIATES/prad/rwr/gcn_rwr_otp
mkdir -p $PATH_INTERMEDIATES/prad/rwr/ppi_rwr_otp

echo Module detection
Rscript prad/prad_gene_modules.R --no-restore --no-save

echo Clustering
echo No gene filters DR
echo PCA
Rscript prad/prad_dimred_pca.R --no-restore --no-save
echo t-SNE
Rscript prad/prad_dimred_tsne.R --no-restore --no-save
echo UMAP 10 neighbours
Rscript prad/prad_dimred_umap_10n.R --no-restore --no-save
echo UMAP 20 neighbours
Rscript prad/prad_dimred_umap_20n.R --no-restore --no-save
echo UMAP 30 neighbours
Rscript prad/prad_dimred_umap_30n.R --no-restore --no-save
echo DEG and OTP intersection DR
echo PCA and t-SNE
Rscript prad/prad_otp_dimred.R --no-restore --no-save
echo UMAP
Rscript prad/prad_otp_dimred_umap.R --no-restore --no-save
echo DEG and OTP intersection GSVA
Rscript prad/prad_otp_gsva.R --no-restore --no-save
echo DEG and OTP intersection DiffRank
Rscript prad/prad_otp_diffrank.R --no-restore --no-save
echo DEG and OTP intersection RWRFGSEA on GCN
Rscript prad/rwr/prad_gcn_rwr_otp.R --no-restore --no-save
Rscript --no-restore --no-save prad/rwr/rwr_clustering_only.R \
  $PATH_INTERMEDIATES/prad/rwr/gcn_rwr_otp/ \
  rwr.pw.scores.csv rwr.pw.clusters.csv.gz \
  KEGG.csv.gz GO.csv.gz REACTOME.csv.gz
echo DEG and OTP intersection RWRFGSEA on PPI
Rscript prad/rwr/prad_ppi_rwr_otp.R --no-restore --no-save
Rscript --no-restore --no-save prad/rwr/rwr_clustering_only.R \
  $PATH_INTERMEDIATES/prad/rwr/ppi_rwr_otp/ \
  rwr.pw.scores.csv rwr.pw.clusters.csv.gz \
  KEGG.csv.gz GO.csv.gz REACTOME.csv.gz


echo Plots
Rscript article/plots_for_figures.R --no-restore --no-save
Rscript article/data_visualization.R --no-restore --no-save
Rscript article/plots_for_more_figures.R --no-restore --no-save






