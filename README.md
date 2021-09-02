# Benchmark branch of COPS

This branch of the repository includes the R code used in performing a benchmark analysis of dimensionality reduction and pathway enrichment based clustering results in TCGA BRCA and PRAD RNA-Seq data for https://doi.org/10.1093/bib/bbab314. 

This branch also includes the specific version of COPS used in the analysis.

## Details

Code for benchmark is included in the folder [inst/benchmark](inst/benchmark). 
If you install this package it will be copied to your R library into `COPS/benchmark`.

The benchmark is split into multiple scripts which are run by in the main bash script `benchmark.sh`.
Clustering results, metrics and figures will be written on disk in the locations specified in `config.sh`.

## Requirements

The default configuration requires 64 GB RAM (32 GB with a large swap-file should run just fine). 
The exact requirements depend on the number of runs and folds used in resampling (repeated CV) as well as the number of threads which can be configured.

The results in the article were compiled with msigdbr v7.2.1 and curatedTCGAData v1.12.0.

## Installation

```R
install.packages("devtools")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install()
BiocManager::install(c("Rgraphviz", "graph", "supraHex", "STRINGdb", "GSVA", 
                       "fgsea", "biomaRt", "AnnotationDbi", "org.Hs.eg.db",
                       "SC3", "Spectrum"))
devtools::install_github("theislab/kBet")
devtools::install_github("cran/clusteval")
devtools::install_github("vittoriofortino84/COPS@benchmark")

# Downgrade msigdbr
devtools::install_github("igordot/msigdbr@v7.2.1")

# Packages used in the benchmark
install.packages(c("caret"))
BiocManager::install(c("curatedTCGAData", "TCGAbiolinks", "WGCNA", "survminer"))

# Packages used to plot the figures
install.packages(c("ggplot2", "gridExtra", "cowplot", "pheatmap", "rPref"))
```

## Results

[inst/benchmark/results](inst/benchmark/results)
