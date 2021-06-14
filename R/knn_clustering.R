#' Title
#'
#' @param dat 
#' @param k 
#' @param jaccard_kernel 
#'
#' @return
#' @export
knn_g <- function(dat, k = 30, jaccard_kernel = TRUE) {
  # Distance metric, for now only Euclidean (Minkowski with p=2)
  dat_dist <- as.matrix(dist(t(dat)), method = "euclidean", diag = TRUE, upper = TRUE)
  # Rank neighbors by distance
  nn_matrix <- apply(dat_dist, 1, function(x) order(order(x)))
  # Don't count edges to self
  diag(nn_matrix) <- Inf
  # Indicator of k neighbors for each node (columns)
  nn_matrix <- nn_matrix < (k + 2)
  
  if (jaccard_kernel) {
    # Compute Jaccard coefficient matrix between columns
    j_matrix <- COPS::jaccard_matrix(nn_matrix)
    # Return as a weighted graph
    return(igraph::graph_from_adjacency_matrix(j_matrix, "undirected", weighted = TRUE, diag = FALSE))
  } else {
    # Return unweighted graph if no kernel was applied
    return(igraph::graph_from_adjacency_matrix(nn_matrix, "undirected", weighted = FALSE, diag = FALSE))
  }
}

#' Title
#'
#' @param dat 
#' @param k 
#' @param jaccard_kernel 
#'
#' @return
#' @export
knn_communities <- function(dat, k = 30, jaccard_kernel = TRUE) {
  # Obtain KNN graph
  g <- knn_g(dat, k, jaccard_kernel)
  # Apply Louvain
  g_coms <- igraph::cluster_louvain(g)
  
  out <- g_coms$membership
  names(out) <- g_coms$names
  return(out)
}
