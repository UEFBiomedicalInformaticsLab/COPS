# COPS - Clustering algorithms for Omics-driven Patient Stratification

COPS is an R-package for clustering patients based on omics data. The COPS package provides a large suite of feature reduction and clustering algorithms as well as a robust evaluation system with multiple metrics. COPS implements clustering stability analysis with cross-validation. 

## Installation

```R
# Tested with R 4.4.1 and Bioconductor 3.20
install.packages("remotes")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install()
bioc_dependencies <- c(
  "ANF", 
  "AnnotationDbi", 
  "biomaRt", 
  "fgsea", 
  "GSVA", 
  "iClusterPlus", 
  "MOFA2", 
  "org.Hs.eg.db", 
  "ROntoTools", 
  "STRINGdb", 
  "ComplexHeatmap")
BiocManager::install(bioc_dependencies)
remotes::install_github("theislab/kBet")
remotes::install_github("trintala/IntNMF")
remotes::install_github("UEFBiomedicalInformaticsLab/COPS")
# Some packages can cause issues when built from source
if (Sys.info()['sysname'] == "Windows") {
  install.packages("CVXR")
} else {
  tryCatch(
    install.packages("CVXR"), 
    error = function(e) warning(
      paste0(
        "Installation of the CVXR-package failed due to: ", 
        e, 
        "\nCVXR is required for multiple kernel learning, ", 
        "but is not necessary for other methods and can be ", 
        "omitted."
      )
    )
  )
}
```

### Optional packages
```R
# For additional internal metrics
remotes::install_github("cran/clusterCrit")
# Additional clustering method: spectral clustering
install.packages("Spectrum")
# Additional clustering method: SC3
BiocManager::install("SC3")
# Optional high-performance library for convex optimization for 
# multiple-kernel learning approaches
install.packages("Rmosek")
```
## Available methods
COPS enables clustering analysis via a pipeline that includes many options for feature extraction (optional) and clustering algorithms. COPS includes several methods for both single- and multi-omic feature extraction and clustering. Furthermore, the feature extraction methods can be divided into data- and biology-driven, where the latter aim to integrate prior knowledge about biological processes. 

### Biology-driven feature extraction methods

The extraction methods include algorithms that transform molecular data (e.g., gene-expression) to pathway activation scores for each sample and pathway, as well as kernels that represent the molecular data as a kernel matrix computed via kernel functions that take the interactions within a pathway and their structure into account. 

| Method    | Integrated knowledge  | Output | Description  | R-package | Reference   |
|-----------|:----------------------|:-------|:-------------|:----------|:------------|
| GSVA      | Feature-sets corresponding to pathways | Pathway activity matrix | Gene Set Variation Analysis summarises gene-expression at the level of pathways by using non-parametric gene statistics. | GSVA | https://doi.org/10.1186/1471-2105-14-7    |
| DiffRank  | Feature-sets corresponding to pathways | Pathway activity matrix | A ranking based approach to score pathway activity in individual samples.    | COPS | https://doi.org/10.1186/s12920-018-0449-4 |
| RWR-FGSEA | One feature-network, e.g., gene-gene-network, and pathway-feature sets | Pathway activity matrix | Combines Random Walk with Restart and Fast Gene Set Enrichment Analysis to first propagate dysregulated gene information through a network to generate gene-statistics which used to identify sample specific pathway enrichment.   | COPS | https://doi.org/10.1093/bib/bbab314       |
| BWK       | Pathway feature-network | Pathway kernel | Betweenness Weighted Kernel uses the linear kernel, but weighs features based on betweenness centrality within a pathway network. | COPS |  |
| PAMOGK    | Pathway feature-network | Pathway kernel | PAthway (Multi-Omic) Kernel uses RWR to smooth features and is defined as a sum of inner products between features along all shortest paths in a network. | COPS | https://doi.org/10.1093/bioinformatics/btaa655 |
| PIK       | Pathway feature-network | Pathway kernel | Pathway Induced Kernel uses the normalized Laplacian matrices of a pathway network to define the kernel. | COPS | https://doi.org/10.1038/s41540-019-0086-3 |

#### Pathway network generation methods

The pathway networks used for kernel-based biology-driven can be generated with different methods which can be combined with any of the pathway integrating kernels. COPS includes functions to:
* Use the KEGG pathway networks included in the ROntoTools-package. 
* Use a given gene-gene network (e.g., PPI) and pathway gene-sets (e.g., gene ontology) to induce pathway networks which are sub-networks of the full network. 

### Data-driven feature extraction methods

Data-driven feature extraction in COPS is done either via dimensionality reduction or kernels. For either type of method there are a couple of options. 

| Method    | Description  | R-package | Reference    |
|-----------|:-------------|:----------|:-------------|
| PCA       | Principal Componenent Analysis in this package refers to the use of Singular Value Decomposition for dimensionality reduction. The data is projected onto eigenvectors that correspond to the orthogonal directions maximizing the variance of the projection, yielding principal components that can be used as a low-rank approximation of the original data. | FactoMineR |  |
| t-SNE     | T-distributed Stochastic Neighborhood Embedding has been used as a dimensionality reduction tool for clustering, however it has been shown to have poor stability and is therefore not recommended. T-SNE produces non-linear embeddings. | Rtsne | http://jmlr.org/papers/v9/vandermaaten08a.html |
| UMAP      | Uniform Manifold Aproximation and Projection is a faster and more stable alternative to t-SNE. Similarly to t-SNE it produces non-linear embeddings that can have better clustering performance than linear embeddings such as principal components. | uwot | https://doi.org/10.48550/arXiv.1802.03426 |
| Linear kernel | The linear kernel is defined as the standard inner product between two vectors: ```sum(x*y)```. | COPS |  |
| Gaussian kernel (a.k.a RBF)| The Gaussian kernel is defined as ```exp(-gamma*sum((x-y)^2))```, where gamma is a parameter defining the smoothness of the kernel, with lower values increasing smoothness. | COPS |  |
| Jaccard or Tanimoto kernel | The Jaccard and Tanimoto kernels can be used to define a kernel based on set similarity of binary feature vectors, in which case they are identical. | COPS |  |

### Single-omic clustering algorithms

Single-omic clustering algorithms refer to standard clustering algorithms which assume that the data is represented by a single matrix. 

| Method     | Description  | R-package | Reference    | 
|------------|:-------------|:----------|:-------------| 
| k-means++  | K-means with efficient random seeding. | ClusterR | http://ilpubs.stanford.edu:8090/778/ |
| agglomerative hierarchical clustering | Agglomerative hierarchical clustering works bottom up, i.e., it initializes every sample as it's own node and combines them hierarchically.  | cluster |  |
| DIANA      | DIvisive hierarchical clustering ANAlysis works top to bottom, i.e., it initializes as a single node containing all samples and divides them hierarchically. | cluster |  |
| GMM        | Gaussian Mixture Models model the data distribution as a mixture of multi-variate Gaussian distributions. After fitting the model, clusters are assigned based on maximum likelihood. | mclust | https://doi.org/10.32614/RJ-2016-021 |
| kNNG Louvain | The application of the Louvain community detection algorithm on k-Nearest Neighbor Graphs scaled by shared neighbors is used by Phenograph which was implemented in Python. We have implemented the clustering functionality in R. | COPS | https://doi.org/10.1016/j.cell.2015.05.047 |
| Spectral | Spectral clustering based on eigendecomposition of a graph or kernel matrix. | Spectrum | https://doi.org/10.1093/bioinformatics/btz704 |
| Kernel k-means | Kernel k-means, implemented with the regular k-means alternating update algorithm and random initializations similar to k-means++ or with the relaxed spectral approximation optimum discretized and finalized with alternating udpates. | COPS |  |

### Multi-omic algorithms

Multi-omic algorithms aim to integrate data consisting of different types of omics data in a balanced and holistic way. Many methods combine joint dimensionality reduction and clustering algorithms into a single approach like iCluster and IntNMF while MOFA is a joint dimensionality reduction method that can be combined with any of the single-view clustering methods in COPS. ANF, MKKM-MR and ECMC work on similarity networks or kernels and the latter can be extracted using the knowledge-integrating-kernel-based feature extraction methods implemented COPS. 

| Method         | Description                               | R-package    | Reference  |
|----------------|-------------------------------------------|--------------|------------|
| iClusterPlus   | A Bayesian latent-variable model based clustering approach that generates cluster assignment based on joint inference across data types. Uses a modified Monte Carlo Newton-Raphson algorithm for optimization. | iClusterPlus | https://doi.org/10.1073/pnas.1208949110  |
| iClusterBayes  | A Bayesian latent-variable model based clustering approach that generates cluster assignment based on joint inference across data types. Uses Bayesian inference and the Metropolis-Hastings algorithm to sample from the latent variable posterior distribution. | iClusterPlus | https://doi.org/10.1093/biostatistics/kxx017 |
| IntNMF | An integrative approach for disease subtype classification based on non-negative matrix factorization. | IntNMF | https://doi.org/10.1371/journal.pone.0176278 |
| MOFA2          | Reconstructs a low-dimensional representation of the data using computationally efficient variational inference. Used for joint dimensionality reduction after which COPS applies single-omic clustering algorithms to the shared representation. | MOFA2 | https://doi.org/10.1186/s13059-020-02015-1 |
| ANF            | Based on  patient affinity networks that are calculated for each omic data type and fused into one network that is used for spectral clustering.  | ANF | https://doi.org/10.1016/j.ymeth.2018.05.020 |
| MKKM-MR        | Multiple Kernel K-Means with Matrix-induced Regularization calculates optimal weights for summing kernels such that redundancy is lowered while diversity of the selected kernels is increased.  | COPS | https://doi.org/10.1609/aaai.v30i1.10249 |
| ECMC           | Multi-view clustering with enhanced consensus is a multiple kernel method that aims to decompose the kernels corresponding to each view into a consensus and a disagreement kernel. The consensus kernel alignment between views is optimised and the combined kernel can be clustered with kernelised k-means. | COPS | https://doi.org/10.1186/s12920-017-0306-x |

## Usage
The single-omic functionalities of COPS are demonstrated in a vignette where single-omic clustering is performed on a psoriasis-related RNA-Seq dataset.
* [Psoriasis RNA-Seq vignette](https://htmlpreview.github.io/?https://github.com/UEFBiomedicalInformaticsLab/COPS/blob/master/vignettes/Psoriasis.html)

Multi-omic methods are demonstrated in another vignette based on the TCGA breast-cancer dataset.
* [Breast cancer multi-omic vignette](https://htmlpreview.github.io/?https://github.com/UEFBiomedicalInformaticsLab/COPS/blob/master/vignettes/Breast_Cancer.html)
