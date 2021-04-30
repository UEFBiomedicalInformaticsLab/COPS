#!/bin/bash
## General configuration
# Environmental variables used so that they can be read in R

# Path to save intermediate results to
export PATH_INTERMEDIATES="$HOME/tcga_test/intermediary_files"
# Specify location to save clustering results to
export PATH_PLOTS="$HOME/tcga_test/plot"

# Save plots in .pdf files
export SAVE_PDF="TRUE"

# Save plots in .svg files
export SAVE_SVG="FALSE"