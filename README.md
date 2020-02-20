# COPS - Clustering algorithms for Omics-driven Patient Stratification

The COPS package provides a suite of feature reduction and clustering analysis tools for disease subtype discovery from 
omics data. The tools are wrapped in a convenient pipeline. 

## Currently implemented
### Knowledge driven feature extraction
* pathway and network analysis
### Dimensionality reduction techniques
* PCA
* t-SNE
* UMAP
### Single-omic clustering algorithms
* k-means, pam
* hierarchical, diana, agnes 
* sota
### Multi-omics clustering algoritmhs
### Performance metrics
* external: clustering stability
* internal: silhouette, Dunn, connectivity
* batch effect: chi squared test rejection rate, DSC, kBET, PCA

## Vignette
[link](https://htmlpreview.github.io/?https://github.com/vittoriofortino84/COPS/blob/master/vignettes/Introduction.html)
