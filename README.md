# COPS - Clustering algorithms for Omics-driven Patient Stratification

The COPS package provides a suite of feature reduction and clustering analysis tools for disease subtype discovery from 
omics data. The tools are wrapped in a convenient pipeline. 

## Installation

```R
# Tested on linux 5.4.105-1, R 4.0.4 and Bioconductor 3.12
install.packages("devtools")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install()
BiocManager::install(c("Rgraphviz", "graph", "supraHex", "STRINGdb", "GSVA", 
                       "fgsea", "biomaRt", "AnnotationDbi", "org.Hs.eg.db"))
devtools::install_github("theislab/kBet")
devtools::install_github("vittoriofortino84/COPS")
```

## Currently implemented
### Knowledge driven feature extraction
* GSVA
* DiffRank
* RWR-FGSEA
### Dimensionality reduction techniques
* PCA
* t-SNE
* UMAP
### Single-omic clustering algorithms
* k-means
* hierarchical, diana
* GMM
### Multi-omics clustering algoritmhs
### Performance metrics
* external: survival analysis, cNMI, ARI, NMI
* internal: silhouette, clustering stability, gene module correlation
* batch effect estimation: chi squared test rejection rate, ARI, NMI, kBET, PCA

## Vignette
[link](https://htmlpreview.github.io/?https://github.com/vittoriofortino84/COPS/blob/master/vignettes/Introduction.html)
