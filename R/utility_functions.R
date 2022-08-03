#' Cross-validation fold permutation
#' 
#' Creates cross-validation folds of data for downstream analysis. 
#'
#' @param dat_list list of data matrices with samples on columns
#' @param nfolds number of cross-validation folds
#' @param nruns number of cross-validation replicates
#' @param stratified_cv if \code{TRUE}, try to maximize separation of batch labels within folds
#' @param mixed_cv if \code{TRUE}, try to minimize separation of batch labels within folds
#' @param ... extra arguments are ignored
#'
#' @return list of data.frames with added columns "fold", "run" and "cv_index" as well as 
#'         duplicated rows of the original data corresponding to different folds.
#' @export
#'
#' @importFrom plyr join
cv_fold <- function(dat_list, 
                    nfolds = 5, 
                    nruns = 2, 
                    stratified_cv = FALSE, 
                    mixed_cv = FALSE,
                    cv_stratification_var = NULL,
                    ...) {
  out <- list()
  for (i in 1:length(dat_list)) {
    folded <- list()
    for (j in 1:nruns) {
      if (!is.null(cv_stratification_var) & (stratified_cv | mixed_cv)) {
        a_ind <- lapply(table(cv_stratification_var), function(x) sample(1:x, x))
        b_ind <- sample(1:length(unique(cv_stratification_var)), length(unique(cv_stratification_var)))
        c_ind <- cumsum(table(cv_stratification_var)[unique(cv_stratification_var)[b_ind]])
        cv_index <- c()
        for (u in 1:length(b_ind)) {
          un <- unique(cv_stratification_var)[b_ind[u]]
          cv_index[cv_stratification_var == un] <- a_ind[[un]] + ifelse(u > 1, c_ind[u-1], 0)
        }
        if (stratified_cv) {
          # Stratified cv folds such that holdout set labels mostly do not match to rest of data
          cv_index <- cv_index %/% -(length(cv_stratification_var) %/% -nfolds) + 1
        } else {
          # Mixed cv folds such that labels are evenly distributed within folds
          cv_index <- cv_index %% nfolds + 1
        }
      } else {
        # Completely random folds
        cv_index <- sample(1:nrow(dat_list[[i]])) %% nfolds + 1
      }
      # Got index, create folds +1 extra "fold" with whole data
      # TODO: The reference fold is the same accross all runs, maybe include it only once? 
      #       Possible incompatibility with external methods. 
      #       Also note that some methods are stochastic. 
      folded[[j]] <- list()
      for (f in 1:(nfolds+1)) {
        # TODO: fix downstream support so that test set can be included too
        #tempfold <- dat_list[[i]][cv_index != f, ]
        #tempfold$fold <- f
        #tempfold$run <- j
        #tempfold$cv_index <- cv_index[cv_index != f]
        folded[[j]][[f]] <- data.table(fold = f, run = j, 
                                       cv_index = cv_index[cv_index != f], 
                                       id = dat_list[[i]]$id[cv_index != f])
      }
      folded[[j]] <- data.table::rbindlist(folded[[j]])
    }
    out[[i]] <- data.table::rbindlist(folded)
  }
  names(out) <- names(dat_list)
  return(out)
}

#' Empirical cumulative density function transformation
#'
#' Estimate each feature (row) distribution using Gaussian kernels with sigma corresponding to sd. 
#' 
#'
#' @param x numerical matrix
#' @param parallel number of threads
#'
#' @return matrix of row-wise ecdf values with dimensions matching input
#' @export
ecdf_transform <- function(x, 
                           parallel = 1) {
  x_sd <- apply(x, 1, sd)
  parallel_clust <- setup_parallelization(parallel)
  score <- tryCatch(foreach(i = 1:ncol(x), 
                   .combine = cbind,
                   .export = c(),
                   .multicombine = TRUE,
                   .maxcombine = max(ncol(x), 2)) %dopar% {
                     score_i <- x * 0
                     for (j in (1:ncol(x))[-i]) {
                       # z-score with respect to each kernel
                       score_i[,j] <- pnorm((x[,j] - x[,i]) / x_sd)
                     }
                     # sum over kernels
                     apply(score_i, 1, sum) / (ncol(x) - 1) 
                   }, finally = close_parallel_cluster(parallel_clust))
  out <- score - 0.5
  colnames(out) <- colnames(x)
  return(out)
}

#' Jaccard index between indicator matrix columns
#'
#' @param x indicator or binary feature matrix
#'
#' @return
#' @export
jaccard_matrix <- function(x) {
  A <- t(x) %*% x
  B <- t(x) %*% (1-x)
  out <- A / (A + B + t(B))
  if (is.null(colnames(A)) | any(duplicated(colnames(A)))) {
    colnames(out) <- rownames(out) <- 1:ncol(out)
  }
  return(out)
}


#' Unweighted gene co-expression network constructor
#' 
#' Generates a gene co-expression network by thresholding gene expression correlations with \code{\link[WGCNA]{signumAdjacencyFunction}}.
#'
#' @param dat gene expression data, samples on columns
#' @param correlation_method correlation method
#' @param cor_threshold numeric threshold used to define edges/links (see \code{\link[WGCNA]{pickHardThreshold}})
#'
#' @return
#' @export
coexpression_network_unweighted <- function(dat, 
                                            correlation_method = "spearman", 
                                            cor_threshold = 0.5) {
  cor_mat <- cor(t(dat), method = correlation_method)
  coexpr_net <- WGCNA::signumAdjacencyFunction(cor_mat, threshold = cor_threshold)
  coexpr_net <- igraph::graph_from_adjacency_matrix(coexpr_net, mode = "undirected", weighted = NULL)
  return(coexpr_net)
}

#' Plot similarity matrix as a heatmap
#' 
#' Order using hierarchical clustering
#'
#' @param sim_mat similarity matrix
#' @param method hclust method
#' @param palette color distiller palette
#' @param limits bounds of similarity
#' @param palette_direction set color direction
#' @param title plot title
#'
#' @return
#' @export
#' @importFrom ggplot2 ggplot aes geom_tile theme_bw coord_fixed ggtitle scale_fill_distiller theme element_blank
plot_similarity_matrix <- function(sim_mat, 
                                   method = "average", 
                                   palette = "RdBu", 
                                   #pos_color = "#6E0700", 
                                   #neg_color = "#05006E", 
                                   #midpoint = 0, 
                                   limits = c(-1,1), 
                                   palette_direction = 1, 
                                   title = NULL) {
  # Remove NAs before hclust
  sim_mat[is.na(sim_mat)] <- 0
  hc <- hclust(as.dist(-limits[2] - sim_mat), method = method)
  dat <- reshape2::melt(sim_mat[hc$order, rev(hc$order)])
  ggplot(dat, aes(Var1, Var2, fill = value)) + geom_tile() + theme_bw() + coord_fixed() + ggtitle(title) + 
    scale_fill_distiller(palette = palette, limits = limits, direction = palette_direction) + 
    #scale_fill_gradientn(colors = c(pos_color, "#FFFFFF", neg_color), 
    #                     values = c(0, (limits[1] + midpoint) / limits[1] + limits[2], 1), 
    #                     limits = limits) + 
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
          axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
}

#' Retrieve disease-gene associations from the Open Targets platforms. 
#' 
#' Utility function to extract relevant disease-gene associations
#' 
#' @param diseases a character vector indicating the disease names. 
#' @param fields a character vector indicating the data types used for the infrence of disease-gene associations.
#' Check the Open Target platforms for more details.
#' 
#' @return a data frame including disease-gene association found for each specified data type.
#' @export
#' @importFrom httr content GET
#' @importFrom jsonlite fromJSON
#' @importFrom AnnotationDbi mapIds
retrieveDiseaseGenesOT <- function(diseases, fields) {
  server <- 'https://platform-api.opentargets.io/v3/platform'
  endpoint_prmtrs <- '/public/association/filter'
  OT_list <- list()
  for(d in diseases) {
    optional_prmtrs <- paste('?size=10000&disease=', d, fields, "&order=association_score.overall", sep="")
    uri <- paste(server, endpoint_prmtrs, optional_prmtrs,sep='')
    get_association_json <- httr::content(httr::GET(uri),'text')
    get_association_usable <- jsonlite::fromJSON(get_association_json, flatten = TRUE)
    OT_score <- get_association_usable$data
    print(dim(OT_score))
    if(length(which(duplicated(OT_score$target.gene_info.symbol) == TRUE)) > 0)
      OT_score = OT_score[-which(duplicated(OT_score$target.gene_info.symbol)),]
    print(dim(OT_score))
    rownames(OT_score) <- OT_score$target.gene_info.symbol
    annot <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, OT_score$target.gene_info.symbol,'ENTREZID','SYMBOL')
    OT_score$entrezId < rep(NA,nrow(OT_score))
    OT_score[names(annot),'entrezId'] <- annot
    OT_list[[length(OT_list)+1]] <- OT_score
  }
  names(OT_list) <- diseases
  return(OT_list)
}

#' Retrieve human protein-protein interaction network from STRINGdb.
#' 
#' Utility function to extract a gene subnetwork from STRINGdb including only the seed genes and their interactions. 
#' 
#' @param gene.diseases a character vector indicating the gene seeds. 
#' @param cutoff a numeric value indicating the cutoff for the edge scores.  
#' @param directed a boolean value indicating the type of grpah to be generated. 
#' 
#' @return an igraph object.
#' @export
#' @importFrom STRINGdb STRINGdb
#' @importFrom igraph graph_from_edgelist graph_from_adjacency_matrix as_adjacency_matrix
#' @importFrom biomaRt useEnsembl getBM
#' @importFrom Matrix rowSums
getHumanPPIfromSTRINGdb <- function(gene.diseases, cutoff = 700, directed = FALSE) {
  string_db <- STRINGdb::STRINGdb$new(version="10", species=9606, score_threshold = cutoff)
  if(directed) {
    gene.stringdb <- string_db$map(my_data_frame = gene.diseases, my_data_frame_id_col_names='target.gene_info.symbol', 
                                   removeUnmappedRows=TRUE)
    string_inter <- string_db$get_interactions(gene.stringdb$STRING_id)
    idx_from <- match(x = string_inter$from, table = gene.stringdb$STRING_id)
    idx_to <- match(x = string_inter$to, table = gene.stringdb$STRING_id)
    print(idx_from)
    print(idx_to)
    print(gene.stringdb$target.gene_info.symbol[idx_to])
    print(gene.stringdb$target.gene_info.symbol[idx_from])
    ppi_network <- data.frame(node1=gene.stringdb$target.gene_info.symbol[idx_from], 
                              node2=gene.stringdb$target.gene_info.symbol[idx_to], 
                              weight=string_inter$combined_score/1000,stringsAsFactors = FALSE)
    g.ppi_network <- igraph::graph_from_edgelist(as.matrix(ppi_network[,1:2]), directed=TRUE)
  } 
  else {
    hs_all_proteins = string_db$get_proteins(); rownames(hs_all_proteins) = hs_all_proteins[,1]
    hs_ex_proteins = hs_all_proteins[which((hs_all_proteins$preferred_name %in% 
                                              rownames(gene.diseases)) == TRUE),]; 
    g = string_db$get_subnetwork(hs_ex_proteins$protein_external_id)
    # create adjacency matrix
    adj_matrix <- igraph::as_adjacency_matrix(g)
    # map gene ids to protein ids
    # get gene/protein ids via Biomart
    #mart <- biomaRt::useMart(host = 'grch37.ensembl.org', 
    #                         biomart='ENSEMBL_MART_ENSEMBL', 
    #                         dataset='hsapiens_gene_ensembl')
    mart <- biomaRt::useEnsembl("ensembl", "hsapiens_gene_ensembl")
    # extract protein ids from the human network
    protein_ids <- sapply(strsplit(rownames(adj_matrix), '\\.'), function(x) x[2])
    # get protein to gene id mappings
    mart_results <- biomaRt::getBM(attributes = c("hgnc_symbol","ensembl_peptide_id"),
                                   filters = "ensembl_peptide_id", values = protein_ids,
                                   mart = mart)
    ### replace protein ids with gene ids
    ix <- match(protein_ids, mart_results$ensembl_peptide_id)
    ix <- ix[!is.na(ix)]
    newnames <- protein_ids
    newnames[match(mart_results[ix,'ensembl_peptide_id'], newnames)] <- mart_results[ix, 'hgnc_symbol']
    rownames(adj_matrix) <- newnames
    colnames(adj_matrix) <- newnames
    ppi <- adj_matrix[!duplicated(newnames), !duplicated(newnames)]
    nullrows <- Matrix::rowSums(ppi)==0
    ppi <- ppi[!nullrows,!nullrows] 
    g.ppi_network <- igraph::graph_from_adjacency_matrix(ppi)
  }
  return(g.ppi_network)
}

#' Gene expression to gene seeds to RWR on PPI into FGSEA scores
#'
#' Conveniently wraps \code{\link{RWRFGSEA}} using 
#' \code{\link{retrieveDiseaseGenesOT}} for disease associated genes, 
#' \code{STRINGdb} human PPI as a network, \code{msigdbr} for gene set annotations and 
#' \code{biomaRt} for matching gene annotations. 
#'
#' @param dat A gene expression matrix with samples on columns.
#' @param disease_id Integer ID in Open Targets platform.
#' @param otp_score Name of association score column in Open Targets.
#' @param otp_cutoff Numeric association score cutoff for Open Targets platform genes.
#' @param ppi_cutoff Numeric PPI link score cutoff.
#' @param pw_min.size Minimum gene set size to use.
#' @param pw_max.size Maximum gene set size to use. 
#' @param dat_gene_key Data gene annotation type.
#' @param gs_subcats Gene set subcategories in \code{\link[msigdbr]{msigdbr}} to retrieve for enrichment analysis. 
#' @param directed_ppi Whether to generate a directed network.
#' @param ... extra arguments are passed on to \code{\link{RWRFGSEA}}
#'
#' @return list of enrichment scores
#' @export
#'
#' @importFrom STRINGdb STRINGdb
#' @importFrom igraph V
#' @importFrom msigdbr msigdbr
#' @importFrom dplyr filter
#' @importFrom biomaRt useEnsembl getBM
expressionToRWFeatures <- function(dat, 
                                   disease_id, 
                                   otp_score = "association_score.datatypes.rna_expression", 
                                   otp_cutoff = 0.0, 
                                   ppi_cutoff = 700, 
                                   pw_min.size = 5,
                                   pw_max.size = 200, 
                                   dat_gene_key = "SYMBOL", 
                                   gs_subcats = c("BP", "MF", "CP:KEGG", "CP:REACTOME"), 
                                   directed_ppi = TRUE,
                                   ...) {
  assoc_score_fields <- paste(paste("&fields=", c('disease.efo_info.label',
                                                  'disease.efo_info.therapeutic_area',
                                                  'target.gene_info.symbol',
                                                  'association_score.overall',
                                                  'disease.id',
                                                  'association_score.datatypes'), sep=''), collapse = "")
  disease_otp <- retrieveDiseaseGenesOT(c(disease_id), assoc_score_fields)[[1]][,-c(10:13)]
  gene.diseases <- disease_otp[which(disease_otp[[otp_score]] > otp_cutoff),]
  
  sdb_res <- STRINGdb::STRINGdb$new(version="10", species=9606, score_threshold = ppi_cutoff)
  gene.network <- sdb_res$get_graph()
  igraph::V(gene.network)$name <- gsub("^9606\\.", "", igraph::V(gene.network)$name)
  
  if (!directed_ppi) {
    gene.network <- igraph::as.undirected(gene.network, mode = "collapse")
  }
  
  # Get pathways information from msigdb (https://www.gsea-msigdb.org/)
  db_annots <- msigdbr::msigdbr(species = "Homo sapiens")
  db_annots <- dplyr::filter(db_annots, grepl(paste(gs_subcats, collapse = "|"), gs_subcat))
  
  # Harmonize gene labels
  mart_attributes <- c("ensembl_gene_id", "ensembl_peptide_id", "entrezgene_id")
  if (dat_gene_key == "SYMBOL") {
    mart_attributes <- c(mart_attributes, "hgnc_symbol")
  }
  mart <- biomaRt::useEnsembl("ensembl", "hsapiens_gene_ensembl")
  mart_results <- biomaRt::getBM(attributes = mart_attributes, mart = mart)
  
  ppi_temp <- mart_results$ensembl_gene_id[match(igraph::V(gene.network)$name, 
                                                 mart_results$ensembl_peptide_id)]
  igraph::V(gene.network)$name[!is.na(ppi_temp)] <- ppi_temp[!is.na(ppi_temp)]
  
  disease.genes <- mart_results$ensembl_gene_id[match(gene.diseases$entrezId, 
                                                      mart_results$entrezgene_id)]
  disease.genes <- disease.genes[!is.na(disease.genes)]
  
  db_annots$ensembl_gene_id <- mart_results$ensembl_gene_id[match(db_annots$entrez_gene, 
                                                                  mart_results$entrezgene_id)]
  db_annots <- db_annots[!is.na(db_annots$ensembl_gene_id),]
  
  if (dat_gene_key == "SYMBOL") {
    rownames_ensembl <- mart_results$ensembl_gene_id[match(rownames(dat), mart_results$hgnc_symbol)]
    dat <- dat[!is.na(rownames_ensembl),]
    rownames(dat) <- rownames_ensembl[!is.na(rownames_ensembl)]
  } else 
    if (dat_gene_key == "ENTREZ") {
      rownames_ensembl <-  mart_results$ensembl_gene_id[match(rownames(dat), mart_results$entrez_gene)]
      dat <- dat[!is.na(rownames_ensembl),]
      rownames(dat) <- rownames_ensembl[!is.na(rownames_ensembl)]
    }
  
  dat <- dat[!duplicated(rownames(dat)),]
  
  list_db_annots <- lapply(split(db_annots, db_annots$gs_name), function(x) x$ensembl_gene_id)
  list_db_annots <- list_db_annots[which(sapply(list_db_annots, length) <= pw_max.size & 
                                           sapply(list_db_annots, length) >= pw_min.size)]
  
  out <- RWRFGSEA(dat, gene.network, list_db_annots, disease.genes, ...)
  return(out)
}

#' Data visualization using PCA, t-SNE and UMAP
#'
#' @param data 
#' @param category 
#' @param category_label 
#' @param tsne_perplexity 
#' @param umap_neighbors 
#'
#' @return
#' @export
#' 
#' @importFrom FactoMineR PCA
#' @importFrom Rtsne Rtsne
#' @importFrom uwot umap
#' @importFrom ggplot2 ggplot aes geom_point scale_color_brewer theme_bw labs ggtitle
triple_viz <- function(data, category, category_label, tsne_perplexity = 45, umap_neighbors = 20, tsne = TRUE) {
  res_pca <- FactoMineR::PCA(data, scale.unit = FALSE, ncp = 2, graph = FALSE)
  res_pca_dat <- as.data.frame(res_pca$ind$coord)
  res_pca_dat <- cbind(res_pca_dat, category)
  colnames(res_pca_dat)[3] <- "category"
  eig_percentages <- res_pca$eig[,"percentage of variance"]
  eig_percentages <- as.character(signif(eig_percentages, 3))
  p1 <- ggplot(res_pca_dat, aes(Dim.1, Dim.2, color = category)) + geom_point(shape = "+", size = 3) + 
    theme_bw() + scale_color_brewer(palette = "Dark2") + 
    labs(x = paste0("PC1 (", eig_percentages[1], "%)"), y = paste0("PC2 (", eig_percentages[2], "%)"), color = category_label) +
    ggtitle("PCA")
  
  if (tsne) {
    res_tsne <- Rtsne::Rtsne(data,
                             dims = 2,
                             perplexity = tsne_perplexity,
                             initial_dims = min(50, dim(data)),
                             check_duplicates = FALSE,
                             pca = TRUE,
                             partial_pca = TRUE,
                             verbose = FALSE)$Y
    res_tsne <- as.data.frame(res_tsne)
    res_tsne <- cbind(res_tsne, category)
    colnames(res_tsne)[3] <- "category"
    p2 <- ggplot(res_tsne, aes(V1, V2, color = category)) + geom_point(shape = "+", size = 3) + 
      theme_bw() + scale_color_brewer(palette = "Dark2") + 
      labs(x = "Z1", y = "Z2", color = category_label) +
      ggtitle("t-SNE")
  } else {
    p2 <- NULL
  }
  
  res_umap <- uwot::umap(data, n_neighbors = umap_neighbors, n_components = 2, pca = min(50, dim(data)), verbose = FALSE, init = "normlaplacian")
  res_umap <- data.frame(Dim.1 = res_umap[,1], Dim.2 = res_umap[,2])
  res_umap <- cbind(res_umap, category)
  colnames(res_umap)[3] <- "category"
  p3 <- ggplot(res_umap, aes(Dim.1, Dim.2, color = category)) + geom_point(shape = "+", size = 3) + 
    theme_bw() + scale_color_brewer(palette = "Dark2") + 
    labs(x = "Z1", y = "Z2", color = category_label) + 
    ggtitle("UMAP")
  
  return(list(PCA = p1, tSNE = p2, UMAP = p3))
}

#' Rbind modification which fills missing columns with NA using base R functions
#'
#' @param a 
#' @param b 
#'
#' @return
#' @export
rbind_fill <- function(a,b) {
  all_cols <- union(colnames(a), colnames(b))
  a_fill <- all_cols[!(all_cols %in% colnames(a))]
  a_fill_mat <- matrix(NA, nrow = nrow(a), ncol = length(a_fill))
  colnames(a_fill_mat) <- a_fill
  a <- cbind(a, a_fill_mat)
  
  b_fill <- all_cols[!(all_cols %in% colnames(b))]
  b_fill_mat <- matrix(NA, nrow = nrow(b), ncol = length(b_fill))
  colnames(b_fill_mat) <- b_fill
  b <- cbind(b, b_fill_mat)
  
  return(rbind(a, b))
}

#' Cbind modification which fills missing rows with NA using base R functions
#'
#' @param a 
#' @param b 
#'
#' @return
#' @export
cbind_fill <- function(a,b) {
  all_rows <- union(rownames(a), rownames(b))
  a_fill <- all_rows[!(all_rows %in% rownames(a))]
  a_fill_mat <- matrix(NA, nrow = length(a_fill), ncol = ncol(a))
  rownames(a_fill_mat) <- a_fill
  colnames(a_fill_mat) <- colnames(a)
  a <- rbind(a, a_fill_mat)
  
  b_fill <- all_rows[!(all_rows %in% rownames(b))]
  b_fill_mat <- matrix(NA, nrow = length(b_fill), ncol = ncol(b))
  rownames(b_fill_mat) <- b_fill
  colnames(b_fill_mat) <- colnames(b)
  b <- rbind(b, b_fill_mat)
  
  return(cbind(a, b))
}

#' Plot p-values in -log10 scale with original labels
#'
#' @param x 
#' @param target 
#' @param x_axis_var 
#' @param color_var 
#' @param palette 
#' @param by 
#' @param facetx 
#' @param facety 
#' @param limits 
#'
#' @return
#' @export
#'
#' @importFrom plyr ddply
#' @importFrom ggplot2 ggplot theme_bw scale_fill_brewer theme scale_y_continuous facet_grid
#' @importFrom scales trans_new log_breaks
plot_pvalues <- function(x, 
                         target, 
                         x_axis_var = NULL, 
                         color_var = NULL, 
                         group_var = NULL,
                         palette = "Dark2",
                         by = c("Approach", "Embedding", "Clustering", "k"), 
                         facetx = NULL, 
                         facety = NULL, 
                         limits = NULL) {
  bp_quantiles <- plyr::ddply(x, by, function(a) quantile(a[[target]], probs = c(0, 0.025, 0.25, 0.5, 0.75, 0.975, 1), na.rm = TRUE))
  colnames(bp_quantiles)[length(by) + 1:7] <- c("Q0", "Q0025", "Q025", "Q05", "Q075", "Q0975", "Q1")
  bp_quantiles$IQR <- log(bp_quantiles$Q075) - log(bp_quantiles$Q025) 
  bp_quantiles$ymax <- apply(cbind(exp(log(bp_quantiles$Q075) + bp_quantiles$IQR * 1.5), bp_quantiles$Q1), 1, min)
  bp_quantiles$ymin <- apply(cbind(exp(log(bp_quantiles$Q025) - bp_quantiles$IQR * 1.5), bp_quantiles$Q0), 1, max)
  if (is.null(limits)) limits <- c(1, 10^floor(min(log10(bp_quantiles$ymin))))
  
  if(!is.null(facetx)) {
    # TODO: check if this works with all configurations
    if (!is.null(facety)) {
      temp_facets <- facet_grid(bp_quantiles[[facetx]] ~ bp_quantiles[[facety]], scales = "fixed")
    } else {
      temp_facets <- facet_grid( ~ bp_quantiles[[facetx]], scales = "fixed")
    }
  } else {
    temp_facets <- NULL
  }
  
  # Deal with string aes in boxplot
  con1 <- !is.null(x_axis_var)
  con2 <- !is.null(color_var)
  con3 <- !is.null(group_var)
  
  if (con1 & con2 & con3) temp_aes <- aes_string(x = x_axis_var, fill = color_var, group = group_var)
  if (con1 & con2 & !con3) temp_aes <- aes_string(x = x_axis_var, fill = color_var)
  if (con1 & !con2 & con3) temp_aes <- aes_string(x = x_axis_var, group = group_var)
  if (con1 & !con2 & !con3) temp_aes <- aes_string(x = x_axis_var)
  if (!con1 & !con2 & con3) temp_aes <- aes_string(group = group_var)
  if (!con1 & con2 & con3) temp_aes <- aes_string(group = group_var, fill = color_var)
  if (!con1 & con2 & !con3) temp_aes <- aes_string(fill = color_var)
  if (!con1 & !con2 & !con3) temp_aes <- aes_string()
  
  temp <- ggplot(bp_quantiles, temp_aes) + 
    geom_boxplot(aes(lower = Q025, upper = Q075, middle = Q05, ymin = ymin, ymax = ymax), 
                 outlier.shape = NA, stat = "identity", lwd = 0.25) + 
    theme_bw() + scale_fill_brewer(palette = "Dark2") + 
    theme(legend.position = "bottom") + 
    scale_y_continuous(trans = scales::trans_new("reverse_log", function(x) -log(x), 
                                                 function(y) exp(-y), breaks = scales::log_breaks()), 
                       limits = limits)
  if (!is.null(temp_facets)) temp <- temp + temp_facets
  return(temp)
}

setup_parallelization <- function(parallel) {
  if (is.null(parallel)) return(NULL)
  if (parallel > 1) {
    parallel_clust <- parallel::makeCluster(parallel)
    doParallel::registerDoParallel(parallel_clust)
    return(parallel_clust)
  }
  foreach::registerDoSEQ()
  return(NULL)
}

close_parallel_cluster <- function(cluster) {
  if(!is.null(cluster)) parallel::stopCluster(cluster)
}

split_by_safe <- function(x, by) {
  if (!is.null(x) & nrow(x) > 0) {
    by <- by[by %in% colnames(x)]
    if (length(by) > 0) {
      if (data.table::is.data.table(x)) {
        x_list <- split(x, by = by)
      } else {
        # probably data.frame
        x_list <- split(x, x[, by, drop = FALSE])
      }
      x_list <- x_list[sapply(x_list, nrow) > 0]
    } else {
      # Nothing to split by
      x_list <- list(x)
    }
  } else {
    x_list <- NULL
  }
  return(x_list)
}


