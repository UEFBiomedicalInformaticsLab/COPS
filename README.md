# COPS - Clustering algorithms for Omics-driven Patient Stratification

COPS is an R-package for clustering patients based on omics data. The COPS package provides a large suite of feature reduction and clustering algorithms as well as a robust evaluation system with multiple metrics. COPS implements clustering stability analysis with cross-validation. 

## Installation

```R
# Tested on Rocky Linux 8.7, R 4.2.1 and Bioconductor 3.15
install.packages("devtools")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install()
bioc_dependencies <- c("ANF", "AnnotationDbi", "biomaRt", "dnet", "fgsea", 
                       "graph", "GSVA", "iClusterPlus", "MOFA2", "org.Hs.eg.db", 
                       "Rgraphviz", "ROntoTools", "STRINGdb", "supraHex", 
                       "ComplexHeatmap")
BiocManager::install(bioc_dependencies)
devtools::install_github("cran/clusteval")
devtools::install_github("cran/clusterCrit")
devtools::install_github("theislab/kBet")
devtools::install_github("UEFBiomedicalInformaticsLab/COPS")
```
## Usage
The single-omic functionalities of COPS are demonstrated in a vignette where single-view clustering is performed on a psoriasis-related RNA-Seq dataset.
* [Psoriasis RNA-Seq vignette](https://htmlpreview.github.io/?https://github.com/UEFBiomedicalInformaticsLab/COPS/blob/master/vignettes/Psoriasis.html)

Multi-omic methods are demonstrated in another vignette based on the TCGA breast-cancer dataset.
* [Breast cancer multi-omic vignette](https://htmlpreview.github.io/?https://github.com/UEFBiomedicalInformaticsLab/COPS/blob/master/vignettes/Breast_Cancer.html)
