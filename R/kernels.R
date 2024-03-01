#' Parallel computation of node betweenness
#'
#' @param networks \code{list} of \code{igraph} objects
#' @param parallel number of threads
#'
#' @return \code{list} of betweenness centrality vectors
#' @export
#'
#' @importFrom igraph betweenness
node_betweenness_parallel <- function(networks, parallel = 1) {
  parallel_clust <- setup_parallelization(parallel)
  b_list <- tryCatch(
    foreach(
      i = 1:length(networks), 
      .combine = c, 
      .inorder = FALSE
      ) %dopar% {
        res <- igraph::betweenness(networks[[i]])
        out <- list()
        out[[names(networks)[i]]] <- res
        out
      }, 
    finally = close_parallel_cluster(parallel_clust)
  )
  return(b_list)
}

#' Weighted linear kernel
#'
#' @param x feature matrix
#' @param weights feature weights
#'
#' @return kernel \code{matrix}
#' @export
weighted_linear_kernel <- function(x, weights) {
  weights <- weights[names(weights) %in% rownames(x) & !is.na(names(weights))]
  if (length(weights) > 0) {
    x_weighted <- as.matrix(weights)[,rep(1, ncol(x))] * x[names(weights),]
    out <- t(x_weighted) %*% (x_weighted)
  } else {
    out <- NULL
  }
  return(out)
}

#' Pathway induced kernel (Manica et al. 2019)
#'
#' @param x gene feature matrix
#' @param L pathway graph Laplacian matrix
#' @param rwr_smoothing apply feature smoothing
#' @param rwr_restart_prob restart probability for smoothing
#' @param ... ignored
#'
#' @return kernel \code{matrix}
#' @export
PIK <- function(
    x, 
    L, 
    rwr_smoothing = FALSE, 
    rwr_restart_prob = 0.75, 
    ...
) {
  common_genes <- colnames(L)[colnames(L) %in% colnames(x)]
  if (length(common_genes) > 0) {
    y <- matrix(
      0, 
      nrow = nrow(x), 
      ncol = ncol(L), 
      dimnames = list(
        id = rownames(x), 
        gene = colnames(L)
      )
    )
    y[,common_genes] <- x[,common_genes]
    if (rwr_smoothing) {
      y <- binary_node_attribute_smoothing_from_adjacency(
        y, 
        L, 
        rwr_restart_prob = rwr_restart_prob
      )
    }
    out <- (y) %*% L %*% t(y)
  } else {
    out <- NULL
  }
  return(out)
}

#' @describeIn PIK Use KEGG pathways to form pathway induced kernels
#'
#' @param x gene feature matrix
#' @param gene_key column in \code{\link[org.Hs.eg.db]{org.Hs.eg.db}} that KEGG IDs should be translated to
#' @param ... passed on to \code{\link{PIK}}
#'
#' @return \code{list} of KEGG-based pathway-kernel matrices
#' @export
#'
#' @importFrom ROntoTools keggPathwayGraphs setEdgeWeights keggPathwayNames
#' @importFrom igraph graph_from_graphnel as.undirected laplacian_matrix
PIK_KEGG <- function(x, gene_key = "SYMBOL", ...) {
  kegg_pw_net <- KEGG_networks()
  
  if (gene_key == "SYMBOL") {
    for (i in 1:length(kegg_pw_net)) {
      kegg_symbols <-  suppressMessages(
        AnnotationDbi::mapIds(
          org.Hs.eg.db::org.Hs.eg.db, 
          gsub("hsa:", "", names(igraph::V(kegg_pw_net[[i]]))), 
          gene_key, 
          "ENTREZID"
        )
      )
      kegg_pw_net[[i]] <- igraph::set.vertex.attribute(
        kegg_pw_net[[i]] , 
        "name", 
        value = kegg_symbols)
    }
  } else {
    stop(paste("gene_key \"", gene_key, "not supported."))
  }
  
  return(PIK(x, kegg_pw_net, ...))
}

#' Download KEGG pathway graphs
#'
#' @return list of undirected igraph objects
#' @export
KEGG_networks <- function() {
  kpg <- ROntoTools::keggPathwayGraphs("hsa", verbose = FALSE)
  kpg <- ROntoTools::setEdgeWeights(
    kpg, 
    edgeTypeAttr = "subtype", 
    edgeWeightByType = list(
      activation = 1, 
      inhibition = 1, 
      expression = 1, 
      repression = 1
    ), 
    defaultWeight = 1)
  names(kpg) <- ROntoTools::keggPathwayNames("hsa")[names(kpg)]
  
  kegg_pw_net <- lapply(kpg, igraph::graph_from_graphnel)
  kegg_pw_net <- lapply(kegg_pw_net, igraph::as.undirected, mode = "collapse")
  return(kegg_pw_net)
}

#' @describeIn PIK Pathway induced kernel from networks
#'
#' @param x gene feature matrix
#' @param networks list of igraph objects corresponding to pathway graphs
#' @param normalized_laplacian normalize Laplacian for calculations
#' @param parallel number of parallel threads
#' @param ... passed on to \code{\link{PIK}}
#'
#' @return list of kernel matrices
#' @export
#' 
#' @importFrom igraph V get.vertex.attribute delete_vertices laplacian_matrix
PIK_from_networks <- function(x, networks, normalized_laplacian = TRUE, parallel = 1, ...) {
  x_genes <- colnames(x)
  parallel_clust <- setup_parallelization(parallel)
  piks <- tryCatch(
    foreach(
      i = 1:length(networks), 
      .combine = c, 
      .inorder = FALSE,
      .packages = c("igraph")
      ) %dopar% {
        v_names <- igraph::get.vertex.attribute(networks[[i]], "name")
        missing_genes <- v_names[!v_names %in% x_genes]
        missing_genes <- missing_genes[!is.na(missing_genes)]
        missing_nodes <- igraph::V(networks[[i]])[missing_genes]
        net <- igraph::delete_vertices(networks[[i]], missing_nodes)
        L <- igraph::laplacian_matrix(net, normalized = normalized_laplacian)
        L@x[is.nan(L@x)] <- 0
        if (sum(abs(L@x)) > 0) {
         out <- PIK(x, L, ...)
        } else {
         out <- NULL
        }
        out <- list(out)
        names(out) <- names(networks)[i]
        out
      }, 
    finally = close_parallel_cluster(parallel_clust)
  )
  piks <- piks[!sapply(piks, is.null)]
  return(piks)
}

binary_node_attribute_smoothing_from_adjacency <- function(
    x, 
    A, 
    rwr_restart_prob = 0.75
) {
  g <- igraph::graph_from_adjacency_matrix(
    as.matrix(A), 
    mode = "undirected", 
    weighted = TRUE)
  y <- dnet::dRWR(
    g, 
    normalise = "none", 
    setSeeds = t(x), 
    restart = rwr_restart_prob,
    parallel = FALSE)
  return(t(y))
}

#' Extract subnetworks based on pathway genesets
#'
#' @param gene_network \code{igraph} object
#' @param gene_sets  \code{list} of gene sets
#'
#' @return list of kernel matrices
#' 
#' @examples
#' \dontrun{
#' pw_db <- msigdbr::msigdbr(species = "Homo sapiens")
#' pw_db <- dplyr::filter(pw_db, grepl("CP:KEGG", gs_subcat))
#' pw_list <- lapply(split(pw_db, pw_db$gs_name), function(x) x$ensembl_gene)
#' pw_list <- pw_list[which(sapply(pw_list, length) <= 200 & sapply(pw_list, length) >= 5)]
#' ppi_net <- COPS::getHumanPPIfromSTRINGdb(gene_id_mart_column = "ensembl_gene_id")
#' ppi_pws <- COPS::pathway_gene_subnetworks(ppi_net, pw_list)
#' # Check the connected components to see if subnets are properly connected:
#' lapply(ppi_pws, function(x) igraph::components(x)$csize) # Not all are
#' }
#' 
#' @export
pathway_gene_subnetworks <- function(gene_network, gene_sets) {
  out <- list()
  for (pw_i in names(gene_sets)) {
    sub_net_v <- gene_sets[[pw_i]]
    sub_net_v <- sub_net_v[sub_net_v %in% igraph::V(gene_network)$name]
    out[[pw_i]] <- igraph::induced_subgraph(gene_network, sub_net_v, impl = "create_from_scratch")
  }
  return(out)
}

#' @describeIn PIK Extract subnetworks based on pathway genesets and compute PIKs
#'
#' @param x gene feature matrix
#' @param gene_network \code{igraph} object
#' @param gene_sets  \code{list} of gene sets
#' @param ... passed on to \code{\link{PIK}}
#'
#' @return list of kernel matrices
#' @export
PIK_GNGS <- function(x, gene_network, gene_sets, ...) {
  L_pw <- list()
  pw_nets <- pathway_gene_subnetworks(gene_network, gene_sets)
  for (pw_i in names(pw_nets)) {
    L_pw[[pw_i]] <- igraph::laplacian_matrix(pw_nets[[i]], normalized = TRUE)
  }
  return(PIK(x, L_pw, ...))
}

#' Multiple kernel k-means with matrix induced regularization (Liu et al. 2016)
#'
#' @param K_list \code{list} of kernel matrices
#' @param k number of clusters
#' @param lambda \code{numeric} regularization parameter
#' @param tolerance \code{numeric} stopping criterion value
#' @param parallel number of parallel threads used by the SDP solver
#' @param use_mosek if \code{TRUE} the optimization will be run with \code{Rmosek} instead of \code{CVXR},
#'
#' @return a kernel \code{matrix}
#' @export
mkkm_mr <- function(
    K_list, 
    k, 
    lambda, 
    tolerance = 1e-6, 
    parallel = 0, 
    use_mosek = FALSE
) {
  M <- matrix(NA, length(K_list), length(K_list))
  for (i in 1:length(K_list)) {
    for (j in i:length(K_list)) {
      M[i,j] <- sum(K_list[[i]] * K_list[[j]])
      M[j,i] <- M[i,j]
    }
  }
  mu <- rep(1, length(K_list)) / length(K_list)
  objective <- NA
  objective_t <- Inf # at least two iterations (initial soution + 1 iteration)
  first <- TRUE
  while(first | log(objective) - log(objective_t) > tolerance) {
    first <- FALSE
    K <- K_list[[1]] * mu[1]**2
    for (i in 2:length(K_list)) {
      K <- K + K_list[[i]] * mu[i]**2
    }
    H <- mkkm_mr_h_opt(K, k)
    mu <- mkkm_mr_mu_opt(K_list, H, M, lambda, parallel, use_mosek)
    objective <- objective_t
    objective_t <- t(mu) %*% M %*% mu
  }
  K <- K_list[[1]] * mu[1]**2
  for (i in 2:length(K_list)) {
    K <- K + K_list[[i]] * mu[i]**2
  }
  H <- mkkm_mr_h_opt(K, k)
  return(list(K = K, H = H, mu = mu, objective = objective))
}

# min <K, (I - HHT)>
# s.t. HTH = I
# k first eigenvectors is optimal
mkkm_mr_h_opt <- function(
    K, 
    k
) {
  eig_vecs <- eigen(K, symmetric = TRUE)$vectors[,1:k]
  return(eig_vecs)
}

# min muT/2 (2 * Z + lambda * M) mu 
# s.t. sum(mu) = 1, mu_i >= 0
# Z = diag(<K_i, I - HHT>_i)
mkkm_mr_mu_opt <- function(
    K_list, 
    H, 
    M, 
    lambda, 
    parallel = 0, 
    use_mosek = FALSE
) {
  n <- nrow(H)
  HHT <- diag(rep(1, n)) - H %*% t(H)
  Z <- c()
  for (i in 1:length(K_list)) {
    Z[i] <- sum(as.vector(K_list[[i]]) * as.vector(HHT))
  }
  Z <- diag(Z)
  
  if (use_mosek) {
    return(mkkm_mr_mu_opt_mosek(Z, M, lambda, parallel = parallel))
  } else {
    return(mkkm_mr_mu_opt_cvxr(Z, M, lambda, parallel = parallel))
  }
}

mkkm_mr_mu_opt_mosek <- function(
    Z, 
    M, 
    lambda, 
    parallel = 0
) {
  if (!requireNamespace("Rmosek", quietly = TRUE)) {
    stop("Trying to run MKKM-MR with MOSEK, but Rmosek has not been installed.")
  }
  m <- ncol(Z)
  prob <- list(sense = "min")
  prob$iparam <- list(NUM_THREADS = parallel)
  prob$A <- Matrix::Matrix(rep(1, m), nrow = 1, sparse = TRUE)
  prob$bc <- rbind(blc = 1, buc = 1)
  prob$c <- rep(0, m)
  prob$bx <- rbind(blx=rep(0, m),
                   bux=rep(Inf, m))
  
  ni <- (m*m+m)/2 # sparse length
  
  # column index
  l <- rep(0, ni)
  l[cumsum(m:2)+1] <- 1
  l <- cumsum(l) + 1
  # row index
  k <- Reduce("c", lapply(1:m, function(x) x:m))
  
  prob$qobj$i <- k
  prob$qobj$j <- l
  prob$qobj$v <- (2*Z + lambda * M)[cbind(k,l)]
  
  res <- Rmosek::mosek(prob)
  return(res$sol$itr$xx)
}

mkkm_mr_mu_opt_cvxr <- function(
    Z,
    M, 
    lambda, 
    parallel = 0
) {
  if (!requireNamespace("CVXR", quietly = TRUE)) {
    stop("Trying to run MKKM-MR with CVXR, but CVXR has not been installed.")
  }
  parallel <- 0 # parallel solver is not implemented
  m <- ncol(Z)
  mu <- CVXR::Variable(m)
  #objective <- CVXR::Minimize(0.5 * t(mu) %*% (2 * Z + lambda * M) %*% mu)
  objective <- CVXR::Minimize(CVXR::quad_form(mu, 0.5 * (2 * Z + lambda * M)))
  constraints <- list(CVXR::sum_entries(mu) == 1, mu >= 0)
  problem <- CVXR::Problem(objective, constraints)
  solution <- CVXR::solve(problem, parallel = parallel, solver = "SCS")
  mu_sol <- solution$getValue(mu)
  return(mu_sol)
}

mkkm_mr_objective <- function(
    K_list, 
    mu
) {
  obj <- 0
  for (i in 1:(length(K_list)-1)) {
    for (j in (i+1):length(K_list)) {
      obj <- obj + sum(K_list[[i]] * K_list[[j]])
    }
  }
  return(obj)
}

#' Global kernel k-means (Tzortzis & Likas 2008)
#'
#' @param K kernel \code{matrix}
#' @param n_clusters number of clusters
#'
#' @return \code{list} of cluster assignments and k-means objective
#' @export
global_kernel_kmeans <- function(
    K, 
    n_clusters
) {
  k_clusters <- rep(1, nrow(K))
  for (k in 2:n_clusters) {
    res <- list()
    for (i in 1:nrow(K)) {
      k_clusters_init <- k_clusters
      k_clusters_init[i] <- k
      res[[i]] <- kernel_kmeans_algorithm(
        K, 
        k, 
        k_clusters_init, 
        maxiter = Inf)
    }
    best_i <- which.min(sapply(res, function(x) x$E))
    k_clusters <- res[[best_i]]$clusters
  }
  return(res[[best_i]])
}

#' Kernel k-means 
#' 
#' K-means++ with set number of iterations. 
#'
#' @param K kernel \code{matrix}
#' @param n_k number of clusters
#' @param n_initializations number of initializations
#' @param maxiter maximum number of k-means iterations
#' @param parallel number of parallel threads for running different initializations
#' @param ... ignored
#'
#' @return \code{list} of cluster assignments and k-means objective
#' @export
kernel_kmeans <- function(
    K, 
    n_k, 
    n_initializations = 100, 
    maxiter = 1e2, 
    parallel = 1, 
    ...
) {
  out <- list(clusters = NA, E = Inf)
  if(n_initializations > ncol(K)) {
    w1 <- "K-means++ is deterministic for a given initialization."
    w2 <- "Specified initializations exceed the number of data points."
    w3 <- "Setting it to the number of data points ..."
    warning(paste(w1, w2, w3))
    n_initializations <- ncol(K)
  }
  lower_error <- function(x, y) {
    if(x$E <= y$E) {
      return(x)
    } else {
      return(y)
    }
  }
  random_seeds <- sample(1:ncol(K), n_initializations, replace = FALSE)
  parallel_clust <- setup_parallelization(parallel)
  out <- tryCatch(
    foreach(
      ri = random_seeds, 
      .combine = lower_error, 
      .export = c(), 
      .packages = c(), 
      .inorder = FALSE
    ) %dopar% {
      out_i <- kernel_kmeanspp(K = K, n_k = n_k, seed = ri, maxiter = maxiter)
      out_i
    }, 
    finally = close_parallel_cluster(parallel_clust)
  )
  return(out)
}

kernel_kmeanspp <- function(
    K, 
    n_k, 
    seed, 
    maxiter = 1e2
) {
  centroids <- seed
  for (i in 2:n_k) {
    min_similarity <- apply(K[, centroids, drop = FALSE], 1, min)
    min_similarity[centroids] <- Inf
    ci <- which.min(min_similarity)
    centroids <- c(centroids, ci)
  }
  init <- apply(K[,centroids, drop = FALSE], 1, which.max)
  return(kernel_kmeans_algorithm(K = K, n_k = n_k, init, maxiter = maxiter))
}

kernel_kmeans_algorithm <- function(
    K, 
    n_k, 
    init, 
    maxiter = 1e2
) {
  K_clusters <- init
  diff <- TRUE
  it <- 0
  K_dist <- matrix(NA, nrow(K), n_k)
  while (diff) {
    if (it > maxiter) {
      warning("Kernel k-means did not converge.")
      break()
    }
    for (i in 1:nrow(K)) {
      for (j in 1:n_k) {
        j_ind <- K_clusters == j
        K_dist[i,j] <- K[i,i] - 
          2 * mean(K[i, j_ind]) + 
          mean(K[j_ind, j_ind])
      }
    }
    K_clusters_t <- apply(K_dist, 1, which.min)
    diff <- any(K_clusters != K_clusters_t)
    K_clusters <- K_clusters_t
    it <- it + 1
  }
  out <- list(clusters = K_clusters)
  out$E <- 0
  for (j in 1:n_k) {
    j_ind <- K_clusters == j
    out$E <- out$E + sum(K_dist[j_ind, j])
  }
  return(out)
}

kernel_kmeans_spectral_approximation <- function(
    eigen_vectors, 
    k
) {
    init_qr <- qr(t(eigen_vectors[,1:k]), LAPACK = TRUE)
    r_11 <- init_qr$qr[,1:k]
    r_11[lower.tri(r_11)] <- 0
    r_12 <- init_qr$qr[,(k+1):nrow(eigen_vectors)]
    r_11_inv <- solve(r_11)
    r_hat <- r_11_inv %*% cbind(r_11, r_12)
    init <- rep(NA, nrow(eigen_vectors))
    init[init_qr$pivot] <- apply(abs(r_hat), 2, which.max)
    return(init)
}

# Utilities
#' Kernel centering
#'
#' @param x kernel matrix
#'
#' @return centered kernel matrix
#' @export
center_kernel <- function(x) {
  H <- diag(rep(1, times = ncol(x))) - 1 / ncol(x)
  return(H %*% x %*% H)
}

#' Kernel normalization to unit norm
#'
#' @param x  kernel matrix
#'
#' @return normalized kernel matrix
#' @export
normalize_kernel <- function(x) {
  D <- diag(1/sqrt(diag(x)))
  return(D %*% x %*% D)
}

scale_kernel_norm <- function(x) {
  return(x / Matrix::norm(x, "F"))
}

