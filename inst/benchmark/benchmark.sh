# This runs all the R scripts to generate results included in the patient stratification article
# TODO: create all the necessary folders in bash or R 
#       e.g. dir.create(..., recursive = TRUE, showWarnings = FALSE).
# TODO: read folders from configuration

source config.sh

## TCGA breast cancer
# No gene filters DR
Rscript brca/brca_dimred.R --no-restore --no-save
# DEG and OTP intersection DR
Rscript brca/brca_otp_dimred.R --no-restore --no-save
# DEG and OTP intersection GSVA
Rscript brca/brca_otp_gsva.R --no-restore --no-save
# DEG and OTP intersection DiffRank
Rscript brca/brca_otp_diffrank.R --no-restore --no-save
# DEG and OTP intersection RWRFGSEA on GCN
Rscript brca/rwr/brca_rwr_gcn_fixedn.R --no-restore --no-save
Rscript --no-restore --no-save brca/rwr/rwr_clustering_only.R \
  ~/tcga/brca/intermediary_files/rwr/gcn_fixedn/ \
  rwr.pw.scores.csv rwr.pw.clusters.csv.gz \
  KEGG.csv.gz GO.csv.gz REACTOME.csv.gz
# DEG and OTP intersection RWRFGSEA on PPI
Rscript brca/rwr/brca_rwr_ppi_fixedn.R --no-restore --no-save
Rscript --no-restore --no-save brca/rwr/rwr_clustering_only.R \
  ~/tcga/brca/intermediary_files/rwr/ppi_fixedn/ \
  rwr.pw.scores.csv rwr.pw.clusters.csv.gz \
  KEGG.csv.gz GO.csv.gz REACTOME.csv.gz

## PRAD breast cancer
# No gene filters DR
Rscript prad/prad_dimred.R --no-restore --no-save
# DEG and OTP intersection DR
Rscript prad/prad_otp_dimred.R --no-restore --no-save
# DEG and OTP intersection GSVA
Rscript prad/prad_otp_gsva.R --no-restore --no-save
# DEG and OTP intersection DiffRank
Rscript prad/prad_otp_diffrank.R --no-restore --no-save
# DEG and OTP intersection RWRFGSEA on GCN
Rscript prad/rwr/prad_rwr_gcn_fixedn.R --no-restore --no-save
Rscript --no-restore --no-save prad/rwr/rwr_clustering_only.R \
  ~/tcga/prad/intermediary_files/rwr/gcn_fixedn/ \
  rwr.pw.scores.csv rwr.pw.clusters.csv.gz \
  KEGG.csv.gz GO.csv.gz REACTOME.csv.gz
# DEG and OTP intersection RWRFGSEA on PPI
Rscript prad/rwr/prad_rwr_ppi_fixedn.R --no-restore --no-save
Rscript --no-restore --no-save prad/rwr/rwr_clustering_only.R \
  ~/tcga/prad/intermediary_files/rwr/ppi_fixedn/ \
  rwr.pw.scores.csv rwr.pw.clusters.csv.gz \
  KEGG.csv.gz GO.csv.gz REACTOME.csv.gz



















