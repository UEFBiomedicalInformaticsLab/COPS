#!/bin/bash
# This runs all the R scripts to generate results included in the patient stratification article
# Scripts and results were put into seperate files in order to conserve memory

source config.sh

mkdir -p $PATH_PLOTS

echo TCGA breast cancer
# Create folders for intermediate results

mkdir -p $PATH_INTERMEDIATES/brca/knn_communities
mkdir -p $PATH_INTERMEDIATES/brca/spectral
mkdir -p $PATH_INTERMEDIATES/brca/sc3

#echo KNN community detection
#Rscript brca/brca_knn_communities.R --no-restore --no-save

#echo KNN Spectral clustering
#Rscript brca/brca_spectral.R --no-restore --no-save

#echo DR k-means consensus \(SC3\)
#Rscript brca/brca_sc3.R --no-restore --no-save

echo TCGA prostate cancer
mkdir -p $PATH_INTERMEDIATES/prad/knn_communities
mkdir -p $PATH_INTERMEDIATES/prad/spectral
mkdir -p $PATH_INTERMEDIATES/prad/sc3
#echo KNN community detection
#Rscript prad/prad_knn_communities.R --no-restore --no-save

echo KNN Spectral clustering
Rscript prad/prad_spectral.R --no-restore --no-save

## DEG gene filter instead of DEG OTP intersection

mkdir -p $PATH_INTERMEDIATES/brca/dimred_deg
mkdir -p $PATH_INTERMEDIATES/brca/dimred_deg/dimred_umap_10n
mkdir -p $PATH_INTERMEDIATES/brca/dimred_deg/dimred_umap_20n
mkdir -p $PATH_INTERMEDIATES/brca/dimred_deg/dimred_umap_30n
mkdir -p $PATH_INTERMEDIATES/brca/diffrank_deg
mkdir -p $PATH_INTERMEDIATES/brca/gsva_deg

echo TCGA breast cancer DEG filter
echo DR
Rscript extra/brca_deg_dimred.R --no-restore --no-save
echo DR UMAP
Rscript extra/brca_deg_dimred_umap.R --no-restore --no-save
echo DiffRank
Rscript extra/brca_deg_diffrank.R --no-restore --no-save
echo GSVA
Rscript extra/brca_deg_gsva.R --no-restore --no-save


mkdir -p $PATH_INTERMEDIATES/prad/dimred_deg
mkdir -p $PATH_INTERMEDIATES/prad/dimred_deg/dimred_umap_10n
mkdir -p $PATH_INTERMEDIATES/prad/dimred_deg/dimred_umap_20n
mkdir -p $PATH_INTERMEDIATES/prad/dimred_deg/dimred_umap_30n
mkdir -p $PATH_INTERMEDIATES/prad/diffrank_deg
mkdir -p $PATH_INTERMEDIATES/prad/gsva_deg

echo TCGA prostate cancer DEG filter
echo DR
Rscript extra/prad_deg_dimred.R --no-restore --no-save
echo DR UMAP
Rscript extra/prad_deg_dimred_umap.R --no-restore --no-save
echo DiffRank
Rscript extra/prad_deg_diffrank.R --no-restore --no-save
echo GSVA
Rscript extra/prad_deg_gsva.R --no-restore --no-save