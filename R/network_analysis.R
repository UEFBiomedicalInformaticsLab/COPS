#' Topology overlap measure
#' 
#' Computes the unsigned topology overlap measure used in WGCNA. 
#' 
#' @param A Adjacency \code{matrix}.
#' 
#' @export
unsigned_TOM_matrix <- function(A) {
  diag(A) <- 0
  L <- Matrix::crossprod(abs(A))
  k <- Matrix::rowSums(A)
  K <- matrix(k)[,rep_len(1, length(k))] 
  K <- pmin(K, t(K))
  W <- (L + A) / (K + 1 - A)
  diag(W) <- 1
  return(W)
}

#' Gene module eigen genes
#' 
#' Computes module eigen genes by applying PCA to gene modules
#' 
#' @param x Gene expression data with samples on columns.
#' @param modules Gene clustering result.
#' @param npc Number of PCs to extract per module.
#' 
#' @export
gene_module_eigengenes <- function(x, modules, npc = 1) {
  split_data <- split(x, f = modules)
  f <- function(y) {
    FactoMineR::PCA(
      t(y),
      scale.unit = FALSE,
      ncp = npc,
      graph = FALSE)$ind$coord[,1:npc]
  }
  eigen_genes <- lapply(split_data, f)
  out <- Reduce(cbind, eigen_genes)
  colnames(out) <- names(split_data)
  return(out)
}