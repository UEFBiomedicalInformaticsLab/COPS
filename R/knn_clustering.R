#' Generate KNN graph from
#'
#' @param dat data matrix with samples on columns
#' @param k number of neighbours
#' @param jaccard_kernel if TRUE computes the Jaccard similarity of nearest 
#'   neighbors for each sample-pair which results in a weighted graph where the 
#'   weight describes the fraction of mutual nearest neighbors
#'
#' @return \code{igraph} object
#' @export
knn_g <- function(
    dat, 
    k = 30, 
    jaccard_kernel = TRUE
) {
  # Distance metric, for now only Euclidean (Minkowski with p=2)
  dat_dist <- as.matrix(dist(
    t(dat), 
    method = "euclidean", 
    diag = TRUE, 
    upper = TRUE))
  # Rank neighbors by distance
  nn_matrix <- apply(dat_dist, 1, function(x) order(order(x)))
  # Don't count edges to self
  diag(nn_matrix) <- Inf
  # Indicator of k neighbors for each node (columns)
  nn_matrix <- nn_matrix < (k + 2)
  
  if (jaccard_kernel) {
    # Compute Jaccard coefficient matrix between columns
    a_mat <- COPS::jaccard_matrix(nn_matrix)
    # Return as a weighted graph
    weighted = TRUE
  } else {
    # Return unweighted graph if no kernel was applied
    a_mat <- nn_matrix
    weighted = FALSE
  }
  return(igraph::graph_from_adjacency_matrix(
    a_mat, 
    "undirected", 
    weighted = weighted, 
    diag = FALSE))
}

#' @describeIn knn_g Louvain community detection on KNNG
#'
#' @param ... passed to \code{\link{knn_g}}
#'
#' @return cluster assignments
#' @export
knn_communities <- function(...) {
  # Obtain KNN graph
  g <- knn_g(dat, k, jaccard_kernel)
  # Apply Louvain
  g_coms <- igraph::cluster_louvain(g)
  
  out <- g_coms$membership
  names(out) <- g_coms$names
  return(out)
}
