if (FALSE) {
  ###################################
  # Pathway induced kernel
  
  # Load PPI
  PPI <- read.table("~/tcga/PPI_entrez_ensemble.txt", header = TRUE)
  PPI <- igraph::graph_from_data_frame(PPI[,1:2], directed = FALSE)
  
  PPI_symbols <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, igraph::V(PPI)$name, "SYMBOL", "ENTREZID")
  ambiguous_ind <- PPI_symbols %in% PPI_symbols[duplicated(PPI_symbols)]
  PPI_symbols[ambiguous_ind] <- igraph::V(PPI)$name[ambiguous_ind]
  igraph::V(PPI)$name <- PPI_symbols
  
  db_annots <- msigdbr::msigdbr(species = "Homo sapiens", "C2", "CP:KEGG")
  #db_annots <- dplyr::filter(db_annots, grepl(paste(GENE_SETS, collapse = "|"), gs_subcat))
  list_db_annots <- lapply(split(db_annots, db_annots$gs_name), function(x) unique(x$gene_symbol))
  list_db_annots <- list_db_annots[which(sapply(list_db_annots, length) <= 200 & sapply(list_db_annots, length) >= 5)]
  
  # KEGG pathways (including networks)
  #UPDATE_CACHE = TRUE # Set to true to load updated pathway definitions.
  #print("Loading pathway DB")
  kpg <- ROntoTools::keggPathwayGraphs("hsa", verbose = TRUE)
  #print("Setting edge weights")
  kpg <- ROntoTools::setEdgeWeights(kpg, 
                                    edgeTypeAttr = "subtype", 
                                    edgeWeightByType = list(activation = 1, inhibition = 1, expression = 1, repression = 1), 
                                    defaultWeight = 0)
  #print("Setting pathway names.")
  names(kpg) <- ROntoTools::keggPathwayNames("hsa")[names(kpg)]
  
  kegg_pw_net <- lapply(kpg, igraph::graph_from_graphnel)
  kegg_pw_net <- lapply(kegg_pw_net, igraph::as.undirected, mode = "collapse")
  
  kegg_pw_net_L <- lapply(kegg_pw_net, igraph::laplacian_matrix, normalized = TRUE)
  kegg_pw_net_L <- lapply(kegg_pw_net_L, function(x) {x@x[is.nan(x@x)] <- 0;return(x)})
  #kegg_pw_net_L_eig <- lapply(kegg_pw_net_L, function(x) eigen(as.matrix(x))$values) # all positive or very close to 0 (numeric imprecision)
  
  #kegg_subnets_entrez <- lapply(kegg_pw_net, function(x) suppressMessages(AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, gsub("hsa:", "", igraph::V(x)$name), "SYMBOL", "ENTREZID")))
  # Translate rows to symbols
  for (i in 1:length(kegg_pw_net_L)) {
    kegg_entrez_i <-  suppressMessages(AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, gsub("hsa:", "", rownames(kegg_pw_net_L[[i]])), "SYMBOL", "ENTREZID"))
    rownames(kegg_pw_net_L[[i]]) <- colnames(kegg_pw_net_L[[i]]) <-  kegg_entrez_i
  }
  
  PIK <- function(x, gene_subnet_l) {
    common_genes <- rownames(gene_subnet_l)[rownames(gene_subnet_l) %in% colnames(x)]
    out <- (x[, common_genes, drop = FALSE]) %*% gene_subnet_l[common_genes, common_genes] %*% t(x[, common_genes, drop = FALSE])
    return(out)
  }
  
  prad_mrna_scaled <- scale(t(log2(prad_tumour_mrna+1)))
  
  prad_pik <- lapply(kegg_pw_net_L, function(x) PIK(prad_mrna_scaled, x))
  
  prad_pik <- prad_pik[which(sapply(prad_pik, function(x) var(as.vector(x))) > 0)]
  
  pik_eig <- lapply(prad_pik, function(x) eigen(as.matrix(x))$values)
  all(sapply(prad_pik, function(x) all(diag(as.matrix(x)) > 0)))
}


#' KEGG pathway induced kernels
#'
#' @param x gene expression matrix
#' @param normalized_laplacian 
#' @param gene_key
#'
#' @return
#' @export
#'
#' @importFrom ROntoTools keggPathwayGraphs setEdgeWeights keggPathwayNames
#' @importFrom igraph graph_from_graphnel as.undirected laplacian_matrix
PIK_KEGG <- function(x, normalized_laplacian = TRUE, gene_key = "SYMBOL") {
  kpg <- ROntoTools::keggPathwayGraphs("hsa", verbose = TRUE)
  kpg <- ROntoTools::setEdgeWeights(kpg, 
                                    edgeTypeAttr = "subtype", 
                                    edgeWeightByType = list(activation = 1, inhibition = 1, expression = 1, repression = 1), 
                                    defaultWeight = 0)
  names(kpg) <- ROntoTools::keggPathwayNames("hsa")[names(kpg)]
  
  kegg_pw_net <- lapply(kpg, igraph::graph_from_graphnel)
  kegg_pw_net <- lapply(kegg_pw_net, igraph::as.undirected, mode = "collapse")
  
  kegg_pw_net_L <- lapply(kegg_pw_net, igraph::laplacian_matrix, normalized = normalized_laplacian)
  kegg_pw_net_L <- lapply(kegg_pw_net_L, function(x) {x@x[is.nan(x@x)] <- 0;return(x)})
  if (gene_key != "ENTREZID") {
    # Translate rows names
    for (i in 1:length(kegg_pw_net_L)) {
      kegg_entrez_i <-  suppressMessages(AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, gsub("hsa:", "", rownames(kegg_pw_net_L[[i]])), gene_key, "ENTREZID"))
      rownames(kegg_pw_net_L[[i]]) <- colnames(kegg_pw_net_L[[i]]) <-  kegg_entrez_i
    }
  }
  common_genes <- rownames(gene_subnet_l)[rownames(gene_subnet_l) %in% colnames(x)]
  out <- (x[, common_genes, drop = FALSE]) %*% gene_subnet_l[common_genes, common_genes] %*% t(x[, common_genes, drop = FALSE])
  return(out)
}

PIK_GNGS <- function(x, gene_network, gene_sets, normalize = FALSE) {
  L_pw <- list()
  V_pw <- list()
  
  #igraph::laplacian_matrix(gene_network)
  
  for (pw_i in names(list_db_annots)[grep("^KEGG_", names(list_db_annots))]) {
    sub_net_v <- list_db_annots[[pw_i]]
    sub_net_v <- sub_net_v[sub_net_v %in% igraph::V(gene_network)$name]
    sub_net_v <- sub_net_v[sub_net_v %in% colnames(x)]
    sub_net <- igraph::induced_subgraph(gene_network, sub_net_v, impl = "create_from_scratch")
    L_pw[[pw_i]] <- igraph::laplacian_matrix(sub_net, normalized = TRUE)
    V_pw[[pw_i]] <- sub_net_v # need this for later
  }
  k_pw <- list()
  for (pw_i in names(L_pw)) {
    k_pw[[pw_i]] <- (x[, V_pw[[pw_i]], drop = FALSE]) %*% L_pw[[pw_i]] %*% t(x[, V_pw[[pw_i]], drop = FALSE])
    # Normalize
    if (normalize) {
      d_temp <- c()
      for (i in 1:ncol(k_pw[[pw_i]])) {
        d_temp[i] <- k_pw[[pw_i]][i]
      }
      d_temp <- diag(1/sqrt(d_temp))
      k_pw[[pw_i]] <- d_temp %*% k_pw[[pw_i]] %*% d_temp
    }
  }
}

# Will only work for affinity matrix K (nonnegative) 
# NOT FINISHED
mkkm_mr_postprocessing <- function(K, H) {
  D <- diag(apply(K - diag(diag(K)), 1, sum))
  Z <- 1/sqrt(D) %*% H
}

#' Multiple kernel k-means with matrix induced regularization (Liu et al. 2016)
#'
#' @param K_list 
#' @param k 
#' @param lambda 
#' @param tolerance 
#'
#' @return
#' @export
mkkm_mr <- function(K_list, k, lambda, tolerance = 1e-6, parallel = 0) {
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
    mu <- mkkm_mr_mu_opt(K_list, H, M, lambda, parallel)
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

# min <K (I - HHT)>
# s.t.HTH = I
# k first eigenvectors is optimal
mkkm_mr_h_opt <- function(K, k) {
  #pca <- FactoMineR::PCA(K, scale.unit = FALSE, ncp = k, graph = FALSE)
  #eig_vecs <- sweep(pca$va$coord, 2, sqrt(pca$eig[1:k, 1]), FUN = "/")
  eig_vecs <- eigen(K, symmetric = TRUE)$vectors[,1:k]
  #t(eig_vecs) %*% (eig_vecs) # ~ Identity
  return(eig_vecs)
}

# min muT/2 (2 * Z + lambda * M) mu 
# s.t. sum(mu) = 1
# Z = diag(<K_i, I - HHT>_i)
mkkm_mr_mu_opt <- function(K_list, H, M, lambda, parallel = 0) {
  n <- nrow(H)
  HHT <- diag(rep(1, n)) - H %*% t(H)
  Z <- c()
  for (i in 1:length(K_list)) {
    Z[i] <- sum(as.vector(K_list[[i]]) * as.vector(HHT))
  }
  Z <- diag(Z)
  
  prob <- list(sense = "min")
  prob$iparam <- list(NUM_THREADS = parallel)
  prob$A <- Matrix::Matrix(rep(1, length(K_list)), nrow = 1, sparse = TRUE)
  prob$bc <- rbind(blc = 1, buc = 1)
  prob$c <- rep(0, length(K_list))
  prob$bx <- rbind(blx=rep(0, length(K_list)),
                   bux=rep(Inf, length(K_list)))
  
  m <- ncol(H)
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

mkkm_mr_objective <- function(K_list, mu) {
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
#' @param K 
#' @param n_clusters 
#'
#' @return
#' @export
global_kernel_kmeans <- function(K, n_clusters) {
  k_clusters <- rep(1, nrow(K))
  for (k in 2:n_clusters) {
    res <- list()
    for (i in 1:nrow(K)) {
      k_clusters_init <- k_clusters
      k_clusters_init[i] <- k
      res[[i]] <- kernel_kmeans_algorithm(K, k, k_clusters_init)
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
#' @param K 
#' @param n_k 
#' @param n_initializations 
#' @param ... 
#'
#' @return
#' @export
kernel_kmeans <- function(K, n_k, n_initializations = 100, ...) {
  out <- list(clusters = NA, E = Inf)
  if(n_initializations > ncol(K)) {
    w1 <- "K-means++ is deterministic for a given initialization."
    w2 <- "Specified initializations exceed the number of data points."
    w3 <- "Setting it to the number of data points ..."
    warning(paste(w1, w2, w3))
    n_initializations <- ncol(K)
  }
  ri <- sample(1:ncol(K), n_initializations, replace = FALSE)
  for (i in 1:n_initializations) {
    out_i <- kernel_kmeanspp(K = K, n_k = n_k, seed = ri[i], ...)
    if(out_i$E < out$E) out <- out_i
  }
  return(out)
}

kernel_kmeanspp <- function(K, n_k, seed, ...) {
  centroids <- seed
  for (i in 2:n_k) {
    min_similarity <- apply(K[, centroids, drop = FALSE], 1, min)
    min_similarity[centroids] <- Inf
    ci <- which.min(min_similarity)
    centroids <- c(centroids, ci)
  }
  init <- apply(K[,centroids, drop = FALSE], 1, which.max)
  return(kernel_kmeans_algorithm(K = K, n_k = n_k, init, ...))
}

kernel_kmeans_algorithm <- function(K, n_k, init, maxiter = 1e3) {
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