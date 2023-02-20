#' Breast cancer multi-omic data
#' 
#' A list of matrics corresponding to RNA-Seq TPM, average gene promoter 
#' methylation M-values and copy number variations from GISTIC 2 for genes 
#' involved in selected KEGG pathways. From the TCGA breast cancer dataset 
#' (https://doi.org/10.1038/nature11412).
#' 
#' @docType data
#' 
#' @usage data(brca_multi_omics)
#' 
#' @format An object of class \code{list}
#' 
#' @keywords datasets
"brca_multi_omics"
#' Pathway networks associated with cancer hallmarks
#' 
#' Associations from Zhang et al. 2020 (https://doi.org/10.3389/fgene.2020.00029). 
#' Networks from the KEGG pathway database (https://doi.org/10.1093/nar/gkac963).
#' 
#' @docType data
#' 
#' @usage data(pathway_networks)
#' 
#' @format An object of class \code{list} containing \code{igraph} objects.
#' 
#' @keywords datasets
"pathway_networks"