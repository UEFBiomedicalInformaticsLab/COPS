#' Random walk with restart and FGSEA
#'
#' Runs random walk in a given network starting from seed genes selected from most over/under expressed in the data intersected with 
#' known disease genes. The gene affinities are then processed with FGSEA to yield pathway features. 
#'
#' @param expr a gene expression matrix with samples on columns
#' @param gene_network an \code{igraph.object} with nodes matching to \code{expr} rows
#' @param gene_set_list list of character vectors containing gene sets for FGSEA
#' @param disease_genes a character vector containing gene IDs associated with the target disease
#' @param rwr_gene_lists list containing exactly two character vectors corresponding to gene ids. Used to separate genes for double RWR 
#'   (e.g. up and down regulated seeds) RWR affinities of second seed list will get subtracted from the RWR affinities of the first seed list. 
#' @param rwr_seed_size integer, controls the number of gene seed candidates to intersect with the disease genes
#' @param min_size integer, minimum size of gene sets
#' @param max_size integer, maximum size of gene sets
#' @param parallel integer, number of threads
#' @param verbose TRUE/FALSE
#' @param rwr_restart_probability restart probability used for RWR. See \code{\link[dnet]{dRWR}} for more details, defaults to one sixth of rows.
#' @param rwr_adjacency_normalization method used to normalise the adjacency matrix of the input graph. See \code{\link[dnet]{dRWR}} for more details.
#' @param rwr_affinity_normalization method used to normalise the rwr affinity matrix. See \code{\link[dnet]{dRWR}} for more details.
#' @param fgsea_input_cutoff a cutoff value used to select most visited genes for FGSEA. 
#' @param rwr_ecdf if TRUE, uses \code{\link[COPS]{ecdf_transform}} to transform expression values into empirical cumulative probability density.
#' @param second_seed_list_reverse_order if TRUE, seeds for second RWR are selected from based on lowest expression or ecdf.
#' @param rwr_return_seeds if TRUE, returns seeds in output list (not implemented)
#' @param fgsea_nperm a numeric value determining the number of permutations used in \code{\link[fgsea]{fgseaSimple}}
#' @param logp_scaling if \code{TRUE}, multiplies FGSEA normalized enrichment scores (NESs) with the log p.value to give more weight to the most significantly enriched pathways.
#' @param ... extra arguments are ignored
#'
#' @return matrix of enrichment scores
#' @export
#' 
#' @examples library(parallel)
#' library(COPS)
#' 
#' ad_data <- ad_ge_micro_zscore
#' ad_wgcna_net <- coexpression_network_unweighted(ad_data)
#' kegg_annotations <- msigdbr::msigdbr(
#'     species = "Homo sapiens", 
#'     category = "C2", 
#'     subcategory = "CP:KEGG")
#' list_kegg_annotations <- lapply(
#'     split(kegg_annotations, kegg_annotations$gs_name), 
#'     function(x) x$ensembl_gene) 
#' 
#' ad_rwrfgsea_res <- RWRFGSEA(
#'     ad_data, 
#'     ad_wgcna_net, 
#'     list_kegg_annotations[1:3], 
#'     parallel = 2)
#' 
#' # Separate genes (e.g. based on DEA)
#' set.seed(0)
#' # toy example random fold-change
#' ad_fc <- sample(c(1,-1), nrow(ad_data), replace = TRUE)
#' 
#' ad_gene_lists <- list(
#'     rownames(ad_data)[ad_fc > 0], 
#'     rownames(ad_data)[ad_fc < 0])
#' ad_rwrfgsea_res <- RWRFGSEA(
#'     ad_data, 
#'     ad_wgcna_net, 
#'     list_kegg_annotations[1:3], 
#'     rwr_genelists = ad_gene_lists, 
#'     rwr_ecdf = TRUE, 
#'     second_seed_list_reverse_order = TRUE, 
#'     parallel = 2)
#' 
#' \dontrun{
#' # OTP genes for dermatitis, moderate download size
#' assoc_score_fields <- paste(
#'     paste0(
#'         "&fields=", 
#'         c("target.gene_info.symbol", "association_score.datatypes")
#'     ), 
#'     collapse = ""
#' )
#' disease_otp <- COPS::retrieveDiseaseGenesOT(
#'     "MONDO_0002406", 
#'     assoc_score_fields)[[1]]
#' otp_genes <- disease_otp$target.gene_info.symbol[
#'     disease_otp$association_score.datatypes.genetic_association > 0]
#' 
#' ad_rwrfgsea_res <- RWRFGSEA(
#'     ad_data, 
#'     ad_wgcna_net, 
#'     list_kegg_annotations[1:3], 
#'     disease_genes = otp_genes, 
#'     rwr_genelists = ad_gene_lists, 
#'     rwr_ecdf = TRUE, 
#'     second_seed_list_reverse_order = TRUE, 
#'     parallel = 2)
#' }
#' 
#' @importFrom foreach foreach %dopar%
#' @importFrom plyr aaply
#' @importFrom igraph V
#' @importFrom dnet dRWR
#' @importFrom fgsea fgsea
RWRFGSEA <- function(
    expr, 
    gene_network,
    gene_set_list, 
    disease_genes = NULL, 
    rwr_gene_lists = NULL, 
    rwr_seed_size = NULL, 
    min_size = 5, 
    max_size = 200, 
    parallel = 1,
    verbose = FALSE,
    rwr_restart_probability = 0.75,
    rwr_adjacency_normalization = "laplacian",
    rwr_affinity_normalization = "none",
    fgsea_input_cutoff = 0,
    rwr_ecdf = FALSE,
    second_seed_list_reverse_order = FALSE,
    rwr_return_seeds = FALSE,
    fgsea_nperm = 10000,
    logp_scaling = TRUE,
    ...
) {
  gene_set_list <- gene_set_list[
    sapply(
      gene_set_list, 
      function(x) length(x) >= min_size & length(x) <= max_size)]
  
  if (!is.null(rwr_gene_lists)) {
    if (length(rwr_gene_lists) == 2) {
      expr_list <- list(
        expr[rwr_gene_lists[[1]], ], 
        expr[rwr_gene_lists[[2]], ])
    } else {
      stop("If set, rwr_gene_lists must have length two.")
    }
  } else {
    expr_list <- list(expr)
  }
  
  rwr_results <- list()
  for (i in 1:length(expr_list)) {
    sample_names <- colnames(expr_list[[i]])
    if (!is.null(disease_genes)) {
      expr_list[[i]] <- expr_list[[i]][
        intersect(
          rownames(expr_list[[i]]), 
          disease_genes
        ),
      ]
    }
    
    # Rank disease-genes within each sample
    if (rwr_ecdf) {
      expr_list[[i]] <- ecdf_transform(expr_list[[i]], parallel)
    }
    
    if (second_seed_list_reverse_order & i == 2) {
      expr_list[[i]] <- -expr_list[[i]]
    }
    
    if (is.null(rwr_seed_size)) {
      rwr_seed_size <- nrow(expr_list[[i]]) %/% 6
    }
    
    rwr_results[[i]] <- rwr_wrapper(
      expr_list[[i]], 
      gene_network, 
      rwr_seed_size = rwr_seed_size, 
      parallel = 1, 
      rwr_restart_probability = rwr_restart_probability, 
      rwr_adjacency_normalization = rwr_adjacency_normalization,
      rwr_affinity_normalization = rwr_affinity_normalization)
  }
  
  if (length(rwr_results) == 2) {
    rwr_results <- rwr_results[[1]] - rwr_results[[2]]
  } else {
    rwr_results <- rwr_results[[1]]
  }
  
  rwrfgsea_results <- t(COPS::fgsea_wrapper(
    rwr_results, 
    gene_set_list, 
    parallel = parallel, 
    fgsea_input_cutoff = fgsea_input_cutoff,
    fgsea_nperm = fgsea_nperm, 
    logp_scaling = logp_scaling))
  
  out <- rwrfgsea_results
  
  if (rwr_return_seeds) {
    attributes(out)$rwr_seeds <- NULL # TODO: implement
  }
  
  return(out)
}

#' @describeIn RWRFGSEA Transforms gene-expression values into gene-network node-affinity values
#'
#' @export
rwr_wrapper <- function(
    expr, 
    gene_network, 
    rwr_seed_size = nrow(expr) %/% 6, 
    parallel = 1, 
    rwr_restart_probability = 0.75, 
    rwr_adjacency_normalization = "laplacian",
    rwr_affinity_normalization = "none",
    ...
) {
  gene_ranking <- apply(expr, 2, rank)
  rownames(gene_ranking) <- rownames(expr)
  
  parallel_clust <- setup_parallelization(parallel)
  
  rwr <- tryCatch(
    dnet::dRWR(
      gene_network, 
      setSeeds = gene_ranking > nrow(gene_ranking) - rwr_seed_size, 
      normalise = rwr_adjacency_normalization,
      restart = rwr_restart_probability, 
      normalise.affinity.matrix = rwr_affinity_normalization, 
      parallel = parallel > 1, 
      multicores = parallel, 
      verbose = FALSE
    ), 
    finally = close_parallel_cluster(parallel_clust)
  )
  rownames(rwr) <- names(igraph::V(gene_network))
  colnames(rwr) <- colnames(expr)
  
  out <- rwr
  return(out)
}

#' @describeIn RWRFGSEA Wrapper for fast pre-ranked gene set enrichment analysis (\code{\link[fgsea]{fgseaSimple}})
#'
#' @param data_matrix numeric values used for ranking, assumes samples on columns
#'
#' @export
fgsea_wrapper <- function(
    data_matrix, 
    gene_set_list, 
    fgsea_input_cutoff = 0, 
    parallel = 1, 
    fgsea_nperm = 10000, 
    logp_scaling = TRUE, 
    ...
) {
  parallel_clust <- setup_parallelization(parallel)
  # Use foreach to speed up per sample FGSEA
  res <- tryCatch(
    foreach(
      genes_i = lapply(1:ncol(data_matrix), function(i) data_matrix[,i]), 
      .combine = rbind, #list, #rbind,
      .export = c(),
      .multicombine = TRUE,
      .maxcombine = max(ncol(data_matrix), 2)
    ) %dopar% {
      # Sum up duplicated gene id:s
      genes_i <- tapply(genes_i, names(genes_i), sum)
      genes_i <- genes_i[abs(genes_i) > fgsea_input_cutoff]
      res_i <- fgsea::fgseaSimple(
        gene_set_list, 
        genes_i, 
        nperm = fgsea_nperm, 
        maxSize = min(
          max(sapply(gene_set_list, length)), 
          length(genes_i) - 1), 
        nproc = 1)
      res_out <- rep(NA, length(gene_set_list))
      if (logp_scaling) {
       res_ii <- res_i$NES * (-log10(res_i$pval)) 
      } else {
       res_ii <- res_i$NES
      }
      res_out[match(res_i$pathway, names(gene_set_list))] <- res_ii
      
      res_out
    }, 
    finally = close_parallel_cluster(parallel_clust)
  )
  rownames(res) <- colnames(data_matrix)
  colnames(res) <- names(gene_set_list)
  # Remove pathways with only NA
  res <- res[, apply(res, 2, function(x) !all(is.na(x)))]
  #for (i in 1:nrow(res)) {
  #  res[i, is.na(res[i,])] <- 0
  #}
  res[is.na(res)] <- 0
  return(res)
}

seeded_rwr <- function(
    X, 
    G, 
    p = 0.7, 
    graph_normalization = c("laplacian", "transition", "none"), 
    cholesky = !igraph::is_directed(G), 
    affinity_normalization = TRUE
) {
  graph_normalization <- match.arg(
    graph_normalization, 
    choices = c("laplacian", "transition", "none")
  )
  v_names <- igraph::vertex_attr(G, "name")
  x_names <- rownames(X)
  x_ind <- match(x_names, v_names)
  x_ind_nna <- !is.na(x_ind)
  X_matched <- Matrix::sparseMatrix(
    i = c(), 
    j = c(), 
    x = numeric(), 
    dims = c(length(v_names), ncol(X)), 
    dimnames = list(v_names, colnames(X))
  )
  X_matched[x_ind[x_ind_nna],] <- X[x_ind_nna,]
  
  X_scaled <- Matrix::colScale(X_matched, 1 / Matrix::colSums(X_matched))
  
  A <- igraph::as_adjacency_matrix(G, type = "both", sparse = TRUE)
  if (graph_normalization == "transition") {
    M <- Matrix::colScale(A, 1 / Matrix::colSums(A))
  } else if (graph_normalization == "laplacian") {
    d <- igraph::degree(G)
    M <- Matrix::colScale(A, 1 / sqrt(d))
    M <- Matrix::rowScale(M, 1 / sqrt(d))
  }
  M <- Matrix::.sparseDiagonal(nrow(M), 1) - (1-p) * M
  if (cholesky & graph_normalization != "transition") {
    M <- Matrix::Cholesky(M)
  }
  Y <- Matrix::solve(M, p * X_scaled)
  if (affinity_normalization) {
    Y <- Matrix::colScale(Y, 1 / Matrix::colSums(Y))
  }
  return(Y)
}
