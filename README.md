# COPS - Clustering algorithms for Omics-driven Patient Stratification

The COPS package provides a suite of feature reduction and clustering analysis tools for disease subtype discovery from 
omics data. The tools are wrapped in a convenient pipeline. 

## Installation

```R
# Tested on CentOS Linux 7, R 4.2.0 and Bioconductor 3.15
install.packages("devtools")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install()
BiocManager::install(c("ANF", "AnnotationDbi", "biomaRt", "dnet", "fgsea", "graph", "GSVA", "iClusterPlus", "MOFA2", "org.Hs.eg.db", "Rgraphviz", "ROntoTools", "STRINGdb", "supraHex"))
devtools::install_github("cran/clusteval")
devtools::install_github("theislab/kBet")
devtools::install_github("cran/clusteval")
devtools::install_github("vittoriofortino84/COPS")
```

## Currently implemented
### Knowledge driven feature extraction
* GSVA
* DiffRank
* RWR-FGSEA
* Pathway Induced Kernel (PIK)
* Betweenness Weighted Kernel (BWK)
* Pathway-based MultiOmic Graph Kernel (PAMOGK)
### Dimensionality reduction techniques
* PCA
* t-SNE
* UMAP
### Joint dimensionality reduction techniques
* MOFA2
* IntNMF
### Single-omic clustering algorithms
* k-means++
* agglomerative hierarchical clustering
* DIANA
* GMM
* kernelized k-means
### Multi-omics clustering algorithms
* ANF
* iCluster
* IntNMF
* Multiple Kernel K-Means with Matrix induced Regularization (MKKM-MR)
### Performance metrics
* external: survival analysis, cNMI, ARI, NMI
* internal: silhouette, clustering stability, gene module correlation
* batch effect estimation: chi squared test rejection rate, ARI, NMI, kBET, PCA

## Vignette
[link](https://htmlpreview.github.io/?https://github.com/vittoriofortino84/COPS/blob/master/vignettes/Introduction.html)
