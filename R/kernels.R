#' Parallel computation of node betweenness
#'
#' @param networks \code{list} of \code{igraph} objects
#' @param parallel number of threads
#' @param pathway_node_betweenness_endpoints If \code{TRUE}, includes path endpoints 
#'   in shortest path ensuring that all non-isolated nodes have betweenness > 0.
#' @param pathway_first_shortest_path If \code{TRUE}, uses only the first found 
#'   path in each breadth first search between node pairs of each network. 
#'
#' @return \code{list} of betweenness centrality vectors
#' @export
#'
#' @importFrom igraph betweenness
node_betweenness_parallel <- function(
    networks, 
    parallel = 1, 
    pathway_node_betweenness_endpoints = TRUE, 
    pathway_first_shortest_path = FALSE
) {
  if (pathway_first_shortest_path) {
    betweenness_fun <- node_betweenness_first_path_only
  } else if (pathway_node_betweenness_endpoints) {
    betweenness_fun <- node_betweenness_with_endpoint
  } else {
    betweenness_fun <- igraph::betweenness
  }
  
  parallel_clust <- setup_parallelization(parallel)
  b_list <- tryCatch(
    foreach(
      i = 1:length(networks), 
      .combine = c, 
      .inorder = FALSE
      ) %dopar% {
        res <- betweenness_fun(networks[[i]])
        out <- list()
        out[[names(networks)[i]]] <- res
        out
      }, 
    finally = close_parallel_cluster(parallel_clust)
  )
  return(b_list)
}

node_betweenness_first_path_only <- function(G) {
  n <- igraph::vcount(G)
  w <- rep(0, n)
  #names(w) <- igraph::as_ids(igraph::V(G))
  for (i in 1:(n-1)) {
    sp <- igraph::shortest_paths(G, i, igraph::V(G)[(i+1):n])
    if (length(sp$vpath) > 0) {
      #ends <- as.character(sapply(sp$vpath, function(x) rev(x)[1]))
      #end_counts <- table(ends)
      plengths <- sapply(sp$vpath, length)
      #ec_reps <- rep(end_counts[ends], plengths)
      sp_unlist <- unlist(lapply(sp$vpath, as.character))
      #sp_sum <- tapply(1 / ec_reps, sp_unlist, sum)
      sp_sum <- table(sp_unlist)
      #if (any(is.na(sp_sum))) stop("grr")
      w[as.integer(names(sp_sum))] <- w[as.integer(names(sp_sum))] + sp_sum
    }
  }
  names(w) <- igraph::get.vertex.attribute(G, "name")
  return(w)
}

node_betweenness_with_endpoint3 <- function(G) {
  n <- igraph::vcount(G)
  w <- rep(0, n)
  for (i in 1:(n-1)) {
    sp <- igraph::all_shortest_paths(G, i, igraph::V(G)[(i+1):n])
    if (length(sp$res) > 0) {
      ends <- as.character(sapply(sp$res, function(x) rev(x)[1]))
      end_counts <- table(ends)
      for (j in 1:length(ends)) {
        w[sp[["res"]][[j]]] <- w[sp[["res"]][[j]]] + 1 / end_counts[ends[j]]
      }
    }
  }
  names(w) <- igraph::get.vertex.attribute(G, "name")
  return(w)
}

node_betweenness_with_endpoint2 <- function(G) {
  n <- igraph::vcount(G)
  w <- rep(0, n)
  #names(w) <- igraph::as_ids(igraph::V(G))
  for (i in 1:(n-1)) {
    sp <- igraph::all_shortest_paths(G, i, igraph::V(G)[(i+1):n])
    if (length(sp$res) > 0) {
      ends <- as.character(sapply(sp$res, function(x) rev(x)[1]))
      end_counts <- table(ends)
      plengths <- sapply(sp$res, length)
      ec_reps <- rep(end_counts[ends], plengths)
      sp_unlist <- unlist(lapply(sp$res, as.character))
      sp_sum <- tapply(1 / ec_reps, sp_unlist, sum)
      if (any(is.na(sp_sum))) stop("grr")
      w[as.integer(names(sp_sum))] <- w[as.integer(names(sp_sum))] + sp_sum
    }
  }
  names(w) <- igraph::get.vertex.attribute(G, "name")
  return(w)
}

node_betweenness_with_endpoint <- function(G) {
  # Get betweenness without endpoints
  w <- igraph::betweenness(G)
  # Endpoints can be added by considering reachability
  # Get reachability via distances
  d <- igraph::distances(G)
  # Don't count self
  diag(d) <- Inf
  r <- apply(!is.infinite(d), 1, sum)
  return(w + r)
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
#' @param parallel number of parallel threads used by the quadratic solver
#' @param use_mosek If \code{TRUE}, the optimization will be run with 
#'   \code{Rmosek} instead of \code{CVXR}. 
#' @param mkkm_mr_maxiter maximum number of iterations, usually the algorithm 
#'   converges soon after 2 iterations
#' @param no_stop If \code{TRUE}, always runs \code{mkkm_mr_maxiter} iterations
#'
#' @return a kernel \code{matrix}
#' @export
mkkm_mr <- function(
    K_list, 
    k, 
    lambda, 
    tolerance = 1e-6, 
    parallel = 0, 
    use_mosek = FALSE, 
    mosek_verbosity = 0L, 
    mkkm_mr_maxiter = 10, 
    no_stop = FALSE
) {
  M <- matrix(NA, length(K_list), length(K_list))
  for (i in 1:length(K_list)) {
    for (j in i:length(K_list)) {
      M[i,j] <- sum(K_list[[i]] * K_list[[j]])
      M[j,i] <- M[i,j]
    }
  }
  mu <- rep(1, length(K_list)) / length(K_list)
  # at least two iterations (initial solution + 1 iteration)
  objective <- c(Inf)
  for(it in 1:mkkm_mr_maxiter) {
    K <- K_list[[1]] * mu[1]**2
    for (i in 2:length(K_list)) {
      K <- K + K_list[[i]] * mu[i]**2
    }
    H <- mkkm_mr_h_opt(K, k)
    K_target <- diag(nrow(H)) - H %*% t(H)
    objective[it+1] <- sum(K*K_target) + lambda / 2 * t(mu) %*% M %*% mu
    
    stop_con <- (objective[it] - objective[it+1]) / objective[it+1] < tolerance
    if (stop_con & !no_stop) break
    
    mu <- mkkm_mr_mu_opt(
      K_list = K_list, 
      H = H, 
      M = M, 
      lambda = lambda, 
      parallel, use_mosek, 
      mosek_verbosity = mosek_verbosity
    )
  }
  K <- K_list[[1]] * mu[1]**2
  for (i in 2:length(K_list)) {
    K <- K + K_list[[i]] * mu[i]**2
  }
  H <- mkkm_mr_h_opt(K, k)
  return(list(K = K, H = H, mu = mu, objective = objective[-1]))
}

# min <K, (I - HHT)>
# s.t. HTH = I
# k first eigenvectors is optimal
mkkm_mr_h_opt <- function(
    K, 
    k
) {
  eig_dec <- eigen(K, symmetric = TRUE)
  eig_vecs <- eig_dec$vectors[,1:k]
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
    use_mosek = FALSE, 
    mosek_verbosity = 0L
) {
  n <- nrow(H)
  HHT <- diag(rep(1, n)) - H %*% t(H)
  Z <- c()
  for (i in 1:length(K_list)) {
    Z[i] <- sum(as.vector(K_list[[i]]) * as.vector(HHT))
  }
  Z <- Matrix::Diagonal(x = Z)
  
  if (use_mosek) {
    return(
      mkkm_mr_mu_opt_mosek(
        Z = Z, 
        M = M, 
        lambda = lambda, 
        parallel = parallel, 
        mosek_verbosity = mosek_verbosity
      )
    )
  } else {
    return(
      mkkm_mr_mu_opt_cvxr(
        Z = Z, 
        M = M, 
        lambda = lambda, 
        parallel = parallel
      )
    )
  }
}

mkkm_mr_mu_opt_mosek <- function(
    Z, 
    M, 
    lambda, 
    parallel = 0, 
    mosek_verbosity = 0L
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
  
  Q <- (2*Z + lambda * M)
  ind <- lower.tri(Q, diag = TRUE)
  indw <- which(ind, arr.ind = TRUE)
  indv <- as.vector(ind)
  prob$qobj$i <- indw[,1]
  prob$qobj$j <- indw[,2]
  prob$qobj$v <- Q@x[indv]
  
  res <- Rmosek::mosek(prob, opts = list(verbose = mosek_verbosity))
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
        max_iter = Inf)
    }
    best_i <- which.min(sapply(res, function(x) x$E))
    k_clusters <- res[[best_i]]$clusters
  }
  return(res[[best_i]])
}

#' Kernel k-means 
#' 
#' Kernel k-means with different algorithm options. Spectral relaxation uses 
#' standard randomly initialized k-means on the eigen vectors of the kernel matrix 
#' while the QR decomposition of the eigen vectors yields a single solution 
#' directly. The last option is to use the kernel matrix to optimize average 
#' distances without utilizing the spectral relaxation. 
#' 
#' The \code{tol} parameter is only used by the spectral relaxation algorithm 
#' which makes use of \code{\link[ClusterR]{KMeans_rcpp}}. Other iterative 
#' algorithms are considered converged only if cluster assignments do not change. 
#'
#' @param K kernel \code{matrix}
#' @param n_k number of clusters
#' @param algorithm one of "spectral", "spectral_qr", or "kernelized"
#' @param spectral_qr_refine refine QR result with kernelized k-means
#' @param kernel_eigen_vectors eigenvectors of the kernel matrix can be pre-computed
#' @param max_iter maximum number of iterations
#' @param num_init number of kmeans++ initializations for kernelized k-means and spectral clustering
#' @param tol delta error convergence threshold for spectral clustering
#' @param parallel number of threads for \code{\link{kernelized_kmeans}}
#' @param ... ignored
#'
#' @return \code{list} of cluster assignments and k-means objective
#' @export
kernel_kmeans <- function(
    K, 
    n_k, 
    algorithm = "spectral_qr", 
    spectral_qr_refine = TRUE, 
    kernel_eigen_vectors = NULL, 
    max_iter = 100,
    num_init = 100, 
    tol = 1e-8, 
    parallel = 1, 
    ...
) {
  if (algorithm == "spectral_qr") {
    if (is.null(kernel_eigen_vectors)) {
      kernel_eigen_vectors <- eigen(K, symmetric = TRUE)$vectors
    }
    res <- kernel_kmeans_spectral_approximation(
      kernel_eigen_vectors[,1:n_k], 
      k = n_k
    )
    if (spectral_qr_refine) {
      res <- kernel_kmeans_algorithm(
        K = K, 
        n_k = n_k, 
        init = res, 
        max_iter = max_iter
      )
    }
    return(res)
  } 
  if (algorithm == "spectral") {
    if (is.null(kernel_eigen_vectors)) {
      kernel_eigen_vectors <- eigen(K, symmetric = TRUE)$vectors
    }
    res <- ClusterR::KMeans_rcpp(
      data = kernel_eigen_vectors[,1:n_k], 
      clusters = n_k, 
      num_init = num_init, 
      max_iters = max_iter,
      tol = tol)
    return(res)
  } 
  if (algorithm == "kernelized") {
    res <- kernelized_kmeans(
      K = K, 
      n_k = n_k, 
      num_init = num_init, 
      max_iter = max_iter, 
      parallel = parallel
    )
    return(res)
  }
  stop(paste("Undefined kernel k-mean algorithm:", algorithm))
}



#' Kernelized k-means 
#' 
#' K-means++ without spectral relaxation. 
#'
#' @param K kernel \code{matrix}
#' @param n_k number of clusters
#' @param num_init number of initializations
#' @param max_iter maximum number of k-means iterations
#' @param parallel number of parallel threads for running different initializations
#' @param ... ignored
#'
#' @return \code{list} of cluster assignments and k-means objective
#' @export
kernelized_kmeans <- function(
    K, 
    n_k, 
    num_init = 100, 
    max_iter = 1e2, 
    parallel = 1, 
    ...
) {
  out <- list(clusters = NA, E = Inf)
  if(num_init > ncol(K)) {
    w1 <- "K-means++ is deterministic for a given initial point."
    w2 <- "Specified initializations exceed the number of data points."
    w3 <- "Setting it to the number of data points ..."
    warning(paste(w1, w2, w3))
    num_init <- ncol(K)
  }
  lower_error <- function(x, y) {
    if(x$E <= y$E) {
      return(x)
    } else {
      return(y)
    }
  }
  random_seeds <- sample(1:ncol(K), num_init, replace = FALSE)
  parallel_clust <- setup_parallelization(parallel)
  out <- tryCatch(
    foreach(
      ri = random_seeds, 
      .combine = lower_error, 
      .export = c(), 
      .packages = c(), 
      .inorder = FALSE
    ) %dopar% {
      kernel_kmeanspp(
        K = K, 
        n_k = n_k, 
        seed = ri, 
        max_iter = max_iter
      )
    }, 
    finally = close_parallel_cluster(parallel_clust)
  )
  return(out)
}

kernel_kmeanspp <- function(
    K, 
    n_k, 
    seed, 
    max_iter = 1e2
) {
  centroids <- seed
  for (i in 2:n_k) {
    min_similarity <- apply(K[, centroids, drop = FALSE], 1, min)
    min_similarity[centroids] <- Inf
    ci <- which.min(min_similarity)
    centroids <- c(centroids, ci)
  }
  init <- apply(K[,centroids, drop = FALSE], 1, which.max)
  return(
    kernel_kmeans_algorithm(
      K = K, 
      n_k = n_k, 
      init, 
      max_iter = max_iter
    )
  )
}

kernel_kmeans_algorithm <- function(
    K, 
    n_k, 
    init, 
    max_iter = 1e2, 
    ...
) {
  K_clusters <- init
  diff <- TRUE
  it <- 0
  K_dist <- matrix(NA, nrow(K), n_k)
  while (diff) {
    if (it > max_iter) {
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

#' Multi-omic kernels
#'
#' Process a set of omics observations into kernels
#'
#' @param dat_list List of input \code{data.frame}s for input.
#' @param data_is_kernels If \code{TRUE}, input data is assumed to be kernel matrices. 
#'   Otherwise kernels are computed based on input data and the 
#'   \code{kernels} parameter. 
#' @param kernels Character vector of kernel names to use for different views. 
#'   See details. 
#' @param kernels_center Logical vector specifying which kernels should be 
#'   centered. Repeated for each view if length 1.
#' @param kernels_normalize Logical vector specifying which kernels should be 
#'   normalized Repeated for each view if length 1.
#' @param kernels_scale_norm Logical vector specifying which kernels should be 
#'   scaled to unit F-norm. Repeated for each view if length 1.
#' @param kernel_gammas Numeric vector specifying gamma for the gaussian kernel. 
#' @param pathway_networks List of \code{igraph} objects containing pathway 
#'   networks. Required for pathway kernels. 
#' @param pathway_node_betweenness_endpoints see \code{\link{node_betweenness_parallel}}
#' @param pathway_first_shortest_path see \code{\link{node_betweenness_parallel}}
#' @param kernel_rwr_restart Restart probability for RWR, applies to both RWR-BWK and PAMOGK. 
#' @param kernel_rwr_seeds Seed selection strategy for RWR, one of: 
#'   "discrete", "continuous", or "threshold". Applies to both RWR-BWK and PAMOGK. 
#'   See details below. 
#' @param kernel_rwr_seed_under_threshold z-score threshold for under-expressed, applies to both RWR-BWK and PAMOGK. 
#' @param kernel_rwr_seed_over_threshold z-score threshold for over-expressed, applies to both RWR-BWK and PAMOGK. 
#' @param kernel_rwr_verbose See \code{\link[dnet]{dRWR}}, applies to both RWR-BWK and PAMOGK.
#' @param gene_id_list If data has been pre-processed by the \code{COPS} pipeline, 
#'   the genes of each omic need to be provided as a list. 
#' @param zero_var_removal If set, removes all zero variance features from each omic.
#' @param mvc_threads Number of threads to use for supported operations. 
#' @param preprocess_data If \code{TRUE}, applies \code{\link{data_preprocess}}.
#' @param pathway_rwr_parallelization parallelizes pathway network RWR (experimental)
#' @param ... Extra arguments are ignored. 
#'
#' @return \code{list} of kernels
#' 
#' Supported kernels: 
#' \itemize{
#'   \item "linear" - Linear kernel based on standard dot product. 
#'   \item "gaussian", "rbf" - Gaussian kernel, a.k.a. radial basis function.
#'   \item "jaccard" - Kernel based on Jaccard index. Used for binary features. 
#'   \item "tanimoto" - For now, this is identical to "jaccard".
#'   \item "BWK" - Betweenness Weighted Kernel. Uses pathway networks to compute 
#'     betweenness centrality which is used to weight features in a linear 
#'     pathway kernels. 
#'   \item "RWR-BWK" - BWK with RWR and z-score based seeding similar to PAMOGK.
#'   \item "PAMOGK" - PAthway Multi-Omics Graph Kernel (Tepeli et al. 2021). 
#'     Uses z-scores, RWR and shortest paths in pathway networks to create 
#'     pathway kernels. 
#'   \item "PIK" - Pathway Induced Kernel (Manica et al. 2019). Uses pathway 
#'     network adjacency matrices (specifically normalized Laplacians) to define 
#'     pathway kernels. 
#' }
#' Please note that for pathway kernels, the input data must always be mapped to 
#' genes and that the names must match with the gene names in the pathways. 
#' The default set of pathways is KEGG molecular pathways with gene symbols. 
#' 
#' PAMOGK RWR seed weight options:
#' \itemize{
#'   \item "discrete" - 1 if |z| > t, 0 otherwise. 
#'   \item "continuous" - z
#'   \item "threshold" - z if |z| > t, 0 otherwise
#' }
#' Regardless of the option, the seeds are divided into two sets based on the 
#' sign of the z-score. Each set has a separate smoothing step and the end 
#' result is two different kernels per pathway per omic. This impacts the 
#' RWR label smoothing by changing the initial distribution. 
#' 
#' @export
get_multi_omic_kernels <- function(
    dat_list, 
    data_is_kernels = FALSE, 
    kernels = rep_len("linear", length(dat_list)),
    kernels_center = TRUE,
    kernels_normalize = TRUE,
    kernels_scale_norm = FALSE,
    kernel_gammas = rep_len(0.5, length(dat_list)),
    pathway_networks = NULL,
    pathway_node_betweenness_endpoints = TRUE, 
    pathway_first_shortest_path = FALSE, 
    kernel_rwr_restart = 0.7,
    kernel_rwr_seeds = "discrete", 
    kernel_rwr_seed_under_threshold = qnorm(0.025), 
    kernel_rwr_seed_over_threshold = qnorm(0.975), 
    kernel_rwr_dnet = TRUE, 
    kernel_rwr_verbose = FALSE, 
    gene_id_list = NULL, 
    zero_var_removal = TRUE, 
    mvc_threads = 1, 
    preprocess_data = TRUE, 
    pathway_rwr_parallelization = FALSE, 
    ...
) {
  if (preprocess_data) {
    dat_processed <- data_preprocess(dat_list)
    dat_list <- dat_processed[["dat_list"]]
    gene_id_list <- dat_processed[["gene_id_list"]]
    
    for (j in 1:length(dat_list)) {
      sel <- grep("^dim[0-9]+$", colnames(dat_list[[j]]))
      
      if (data_is_kernels & length(sel) > nrow(dat_list[[j]])) {
        stop("Input kernels are not square!")
      }
      dat_list[[j]] <- as.matrix(as.data.frame(dat_list[[j]])[,sel])
    }
  }
  if (zero_var_removal & !data_is_kernels) {
    dat_list <- lapply(dat_list, function(x) x[,apply(x, 2, var) > 0])
  }
  # Centering, normalization and scaling within subsets
  if (length(kernels_center) != length(dat_list)) {
    kernels_center <- rep_len(kernels_center, length(dat_list))
  }
  if (length(kernels_normalize) != length(dat_list)) {
    kernels_normalize <- rep_len(kernels_normalize, length(dat_list))
  }
  if (length(kernels_scale_norm) != length(dat_list)) {
    kernels_scale_norm <- rep_len(kernels_scale_norm, length(dat_list))
  }
  if (data_is_kernels) {
    multi_omic_kernels <- dat_list
    multi_omic_kernels[kernels_center] <- lapply(
      multi_omic_kernels[kernels_center],
      center_kernel)
    multi_omic_kernels[kernels_normalize] <- lapply(
      multi_omic_kernels[kernels_normalize],
      normalize_kernel)
    multi_omic_kernels[kernels_scale_norm] <- lapply(
      multi_omic_kernels[kernels_scale_norm],
      scale_kernel_norm)
    return(multi_omic_kernels)
  }
  # Pathway-based kernels need pathway networks
  if (any(kernels %in% c("PIK", "BWK", "PAMOGK")) &
      is.null(pathway_networks)) {
    w1 <- "No pathway networks specified for pathway kernel."
    w2 <- "Defaulting to KEGG pathways mapped to gene symbols."
    warning(paste(w1, w2))
    kegg_nets <- KEGG_networks()
    for (net_i in 1:length(kegg_nets)) {
      kegg_entrez <- gsub("hsa:", "", names(igraph::V(kegg_nets[[net_i]])))
      kegg_symbol <- AnnotationDbi::mapIds(
        org.Hs.eg.db::org.Hs.eg.db, 
        kegg_entrez, 
        "SYMBOL", 
        "ENTREZID")
      kegg_nets[[net_i]] <- igraph::set.vertex.attribute(
        kegg_nets[[net_i]], 
        "name", 
        value = kegg_symbol)
    }
    pathway_networks <- kegg_nets
  }
  if (any(kernels %in% c("BWK", "PAMOGK"))) {
    nw_weights <- node_betweenness_parallel(
      pathway_networks, 
      mvc_threads, 
      pathway_node_betweenness_endpoints = pathway_node_betweenness_endpoints, 
      pathway_first_shortest_path = pathway_first_shortest_path)
    nw_weights <- lapply(nw_weights, sqrt)
  }
  # Construct kernels
  multi_omic_kernels <- list()
  for (i in 1:length(dat_list)) {
    if (kernels[i] == "linear") {
      temp <- dat_list[[i]] %*% t(dat_list[[i]])
      if (kernels_center[i]) temp <- center_kernel(temp)
      if (kernels_normalize[i]) temp <- normalize_kernel(temp)
      if (kernels_scale_norm[i]) temp <- scale_kernel_norm(temp)
      temp <- list(temp)
      names(temp) <- names(dat_list)[i]
      multi_omic_kernels <- c(multi_omic_kernels, temp)
    } else if (kernels[i] %in% c("gaussian", "rbf")) {
      temp <- exp(- kernel_gammas[i] * as.matrix(dist(dat_list[[i]]))**2)
      temp <- list(temp)
      names(temp) <- names(dat_list)[i]
      multi_omic_kernels <- c(multi_omic_kernels, temp)
    } else if (kernels[i] %in% c("jaccard", "tanimoto")) {
      temp <- jaccard_matrix(t(dat_list[[i]]))
      temp[is.nan(temp)] <- 0
      diag(temp) <- 1
      if (kernels_center[i]) temp <- center_kernel(temp)
      if (kernels_normalize[i]) temp <- normalize_kernel(temp)
      if (kernels_scale_norm[i]) temp <- scale_kernel_norm(temp)
      temp <- list(temp)
      names(temp) <- names(dat_list)[i]
      multi_omic_kernels <- c(multi_omic_kernels, temp)
    } else if (kernels[i] %in% c("PIK", "BWK", "PAMOGK")) {
      gene_col_ind <- as.integer(gsub("^dim", "", colnames(dat_list[[i]])))
      temp <- dat_list[[i]]
      colnames(temp) <- gene_id_list[[i]][gene_col_ind]
      if (kernels[i] == "PIK") {
        temp <- scale(temp, scale = TRUE) # z-scores
        temp <- PIK_from_networks(temp, pathway_networks, parallel = mvc_threads)
        names(temp) <- paste0(names(dat_list)[i], "_", names(temp))
        temp <- lapply(temp, as.matrix)
        if (kernels_normalize[i]) {
          temp <- lapply(temp, normalize_kernel)
          temp <- lapply(temp, function(x) {x[is.na(x)]  <- 0;return(x)})
        }
        if (kernels_scale_norm[i]) temp <- lapply(temp, scale_kernel_norm)
        multi_omic_kernels <- c(multi_omic_kernels, temp)
      } else if (kernels[i] == "BWK") {
        temp <- t(temp)
        temp <- lapply(nw_weights, function(w) weighted_linear_kernel(temp, w))
        names(temp) <- paste0(names(dat_list)[i], "_", names(temp))
        temp <- temp[!sapply(temp, is.null)]
        temp <- lapply(temp, as.matrix)
        temp <- temp[which(sapply(temp, function(x) var(as.vector(x))) > 0)]
        if (kernels_center[i]) temp <- lapply(temp, center_kernel)
        if (kernels_normalize[i]) {
          temp <- lapply(temp, normalize_kernel)
          temp <- lapply(temp, function(x) {x[is.na(x)]  <- 0;return(x)})
        }
        if (kernels_scale_norm[i]) temp <- lapply(temp, scale_kernel_norm)
        multi_omic_kernels <- c(multi_omic_kernels, temp)
      } else if (kernels[i] %in% c("RWR-BWK", "PAMOGK")) {
        temp <- scale(temp, scale = TRUE) # z-scores
        if (mvc_threads > 1) {
          rwr_threads <- mvc_threads
        } else {
          rwr_threads <- NULL
        }
        
        if (kernel_rwr_seeds == "discrete") {
          seed_up <- t(temp) > kernel_rwr_seed_over_threshold
          seed_dn <- t(temp) < kernel_rwr_seed_under_threshold
        } else if (kernel_rwr_seeds == "continuous") {
          seed_up <- t(temp)
          seed_dn <- -t(temp)
          seed_up[seed_up < 0] <- 0
          seed_dn[seed_dn < 0] <- 0
        } else if (kernel_rwr_seeds == "threshold") {
          seed_up <- seed_dn <- t(temp)
          seed_up[seed_up < kernel_rwr_seed_over_threshold] <- 0
          seed_dn[seed_dn > kernel_rwr_seed_under_threshold] <- 0
          seed_dn <- -seed_dn
        } else {
          stop(paste0(
            "kernel_rwr_seeds option '", 
            kernel_rwr_seeds, 
            "' not recognized. "
          ))
        }
        up_gene_ind <- seed_up > 0
        dn_gene_ind <- seed_dn > 0
        
        if (kernel_rwr_dnet) {
          if (!requireNamespace("dnet", quietly = TRUE)) {
            stop("Please install the dnet-package or set the 'kernel_rwr_dnet' option to FALSE.")
          }
          if (kernel_rwr_verbose) {
            rwr_message_wrapper <- function(x) return(x)
          } else {
            rwr_message_wrapper <- function(x) return(suppressMessages(suppressWarnings(x)))
          }
          if (pathway_rwr_parallelization) {
            parallel_clust <- setup_parallelization(mvc_threads)
            multi_omic_kernels_j <- tryCatch(
              foreach(
                jnet = pathway_networks, 
                jnw = nw_weights, 
                .combine = c
              ) %dopar% {
                kernels_j <- list()
                pw_gene_ind <- rownames(seed_up) %in% names(igraph::V(jnet))
                if (!any(pw_gene_ind)) next
                any_up_gene <- apply(up_gene_ind[pw_gene_ind, , drop = FALSE], 2, any)
                # Skip pathways where only one sample has seeds
                if (sum(any_up_gene)>1) {
                  k_up <- rwr_message_wrapper(dnet::dRWR(
                    jnet, 
                    normalise = "laplacian",
                    setSeeds = seed_up,
                    restart = kernel_rwr_restart,
                    normalise.affinity.matrix = "none",
                    #parallel = mvc_threads > 1,
                    #multicores = rwr_threads, 
                    verbose = kernel_rwr_verbose
                  ))
                  rownames(k_up) <- names(igraph::V(jnet))
                  colnames(k_up) <- rownames(temp)
                  k_up <- weighted_linear_kernel(as.matrix(k_up), jnw)
                  k_up <- process_kernel(
                    K = k_up, 
                    center = kernels_center[i], 
                    normalize = kernels_normalize[i], 
                    scale = kernels_scale_norm[i]
                  )
                  kernels_j <- c(kernels_j, list(k_up))
                }
                # Same check for down genes
                any_dn_gene <- apply(dn_gene_ind[pw_gene_ind, , drop = FALSE], 2, any)
                if (sum(any_dn_gene)>1) {
                  k_dn <- rwr_message_wrapper(dnet::dRWR(
                    jnet, 
                    normalise = "laplacian",
                    setSeeds = seed_dn,
                    restart = kernel_rwr_restart,
                    normalise.affinity.matrix = "none",
                    #parallel = mvc_threads > 1,
                    #multicores = rwr_threads, 
                    verbose = kernel_rwr_verbose
                  ))
                  rownames(k_dn) <- names(igraph::V(jnet))
                  colnames(k_dn) <- rownames(temp)
                  k_dn <- weighted_linear_kernel(as.matrix(k_dn), jnw)
                  k_dn <- process_kernel(
                    K = k_dn, 
                    center = kernels_center[i], 
                    normalize = kernels_normalize[i], 
                    scale = kernels_scale_norm[i]
                  )
                  kernels_j <- c(kernels_j, list(k_dn))
                }
                kernels_j
              }, 
              finally = close_parallel_cluster(parallel_clust)
            )
            multi_omic_kernels <- c(multi_omic_kernels, multi_omic_kernels_j)
          } else {
            for (j in 1:length(pathway_networks)) {
              jnet <- pathway_networks[[j]]
              jnw <- nw_weights[[j]]
              pw_gene_ind <- rownames(seed_up) %in% names(igraph::V(jnet))
              if (!any(pw_gene_ind)) next
              any_up_gene <- apply(up_gene_ind[pw_gene_ind, , drop = FALSE], 2, any)
              # Skip pathways where only one sample has seeds
              if (sum(any_up_gene)>1) {
                
                k_up <- rwr_message_wrapper(dnet::dRWR(
                  jnet, 
                  normalise = "laplacian",
                  setSeeds = seed_up,
                  restart = kernel_rwr_restart,
                  normalise.affinity.matrix = "none",
                  parallel = mvc_threads > 1,
                  multicores = rwr_threads, 
                  verbose = kernel_rwr_verbose
                ))
                rownames(k_up) <- names(igraph::V(jnet))
                colnames(k_up) <- rownames(temp)
                k_up <- weighted_linear_kernel(as.matrix(k_up), jnw)
                k_up <- process_kernel(
                  K = k_up, 
                  center = kernels_center[i], 
                  normalize = kernels_normalize[i], 
                  scale = kernels_scale_norm[i]
                )
                multi_omic_kernels <- c(multi_omic_kernels, list(k_up))
              }
              # Same check for down genes
              any_dn_gene <- apply(dn_gene_ind[pw_gene_ind, , drop = FALSE], 2, any)
              if (sum(any_dn_gene)>1) {
                k_dn <- rwr_message_wrapper(dnet::dRWR(
                  jnet, 
                  normalise = "laplacian",
                  setSeeds = seed_dn,
                  restart = kernel_rwr_restart,
                  normalise.affinity.matrix = "none",
                  parallel = mvc_threads > 1,
                  multicores = rwr_threads, 
                  verbose = kernel_rwr_verbose
                ))
                rownames(k_dn) <- names(igraph::V(jnet))
                colnames(k_dn) <- rownames(temp)
                k_dn <- weighted_linear_kernel(as.matrix(k_dn), jnw)
                k_dn <- process_kernel(
                  K = k_dn, 
                  center = kernels_center[i], 
                  normalize = kernels_normalize[i], 
                  scale = kernels_scale_norm[i]
                )
                multi_omic_kernels <- c(multi_omic_kernels, list(k_dn))
              }
            }
          }
        } else {
          for(j in 1:length(pathway_networks)) {
            jnet <- pathway_networks[[j]]
            jnw <- nw_weights[[j]]
            kernels_j <- list()
            pw_gene_ind <- rownames(seed_up) %in% names(igraph::V(jnet))
            if (!any(pw_gene_ind)) next
            any_up_gene <- apply(up_gene_ind[pw_gene_ind, , drop = FALSE], 2, any)
            # Skip pathways where only one sample has seeds
            if (sum(any_up_gene)>1) {
              k_up <- seeded_rwr(
                X = seed_up, 
                G = jnet, 
                p = kernel_rwr_restart, 
                graph_normalization = "laplacian", 
                affinity_normalization = TRUE
              )
              rownames(k_up) <- names(igraph::V(jnet))
              colnames(k_up) <- rownames(temp)
              k_up <- weighted_linear_kernel(as.matrix(k_up), jnw)
              k_up <- process_kernel(
                K = k_up, 
                center = kernels_center[i], 
                normalize = kernels_normalize[i], 
                scale = kernels_scale_norm[i]
              )
              multi_omic_kernels <- c(multi_omic_kernels, list(k_up))
            }
            # Same check for down genes
            any_dn_gene <- apply(dn_gene_ind[pw_gene_ind, , drop = FALSE], 2, any)
            if (sum(any_dn_gene)>1) {
              k_dn <- seeded_rwr(
                X = seed_dn, 
                G = jnet, 
                p = kernel_rwr_restart, 
                graph_normalization = "laplacian", 
                affinity_normalization = TRUE
              )
              rownames(k_dn) <- names(igraph::V(jnet))
              colnames(k_dn) <- rownames(temp)
              k_dn <- weighted_linear_kernel(as.matrix(k_dn), jnw)
              k_dn <- process_kernel(
                K = k_dn, 
                center = kernels_center[i], 
                normalize = kernels_normalize[i], 
                scale = kernels_scale_norm[i]
              )
              multi_omic_kernels <- c(multi_omic_kernels, list(k_dn))
            }
          }
        }
      } 
    } else {
      stop(paste0("Kernel \"", kernels[i], "\" is not supported."))
    }
  }
  if (!check_dims(multi_omic_kernels)) stop("Kernel dimensions mis-match.")
  return(multi_omic_kernels)
}

check_dims <- function(x_list) {
  if (length(x_list) <= 1) return(TRUE)
  return(all(sapply(x_list[-1], function(x) identical(dim(x), dim(x_list[[1]])))))
}

process_kernel <- function(
    K, 
    center = TRUE, 
    normalize = TRUE, 
    scale = FALSE
) {
  if (!is.null(K)) {
    if (var(as.vector(K)) > 0) {
      if (center) K <- center_kernel(K)
      if (normalize) {
        K <- normalize_kernel(K)
        K[is.na(K)] <- 0
      }
      if (scale) K <- scale_kernel_norm(K)
      return(K)
    }
  }
  return(NULL)
}

#' Kernel clustering utilities
#' 
#' Centering should be applied before normalization. 
#'
#' @param x kernel matrix
#'
#' @return centered kernel matrix
#' @export
center_kernel <- function(x) {
  H <- diag(rep(1, times = ncol(x))) - 1 / ncol(x)
  return(H %*% x %*% H)
}

#' @describeIn center_kernel Normalize kernel space to unit norm
#'
#' @param x kernel matrix
#'
#' @return normalized kernel matrix
#' @export
normalize_kernel <- function(x) {
  D <- diag(1/sqrt(diag(x)))
  return(D %*% x %*% D)
}

#' @describeIn center_kernel Scale kernel to unit Frobenius norm
#'
#' @param x kernel matrix
#'
#' @return normalized kernel matrix
#' @export
scale_kernel_norm <- function(x) {
  return(x / Matrix::norm(x, "F"))
}
