#' Cross-validation fold permutation
#' 
#' Creates cross-validation folds of data for downstream analysis. 
#'
#' @param dat_list list of datasets, each either data.table or data.frame (samples x features) with an "id" column or expression matrix (genes x samples) with named columns
#' @param nfolds number of cross-validation folds
#' @param nruns number of cross-validation replicates
#' @param stratified_cv if \code{TRUE}, perform stratified sampling for folds
#' @param anti_stratified if \code{TRUE}, maximize separation of batch labels within folds, opposite of stratified sampling
#' @param cv_stratification_var labels used for stratification
#' @param ... extra arguments are ignored
#'
#' @return list of data.frames with added columns "fold", "run" and "cv_index" as well as 
#'         duplicated rows of the original data corresponding to different folds.
#' @export
#'
#' @importFrom plyr join
#' @importFrom data.table data.table rbindlist
cv_fold <- function(dat_list, 
                    nfolds = 5, 
                    nruns = 2, 
                    stratified_cv = FALSE, 
                    anti_stratified = FALSE,
                    cv_stratification_var = NULL,
                    extra_fold = TRUE, 
                    ...) {
  out <- list()
  for (i in 1:length(dat_list)) {
    if (any(c("data.table", "data.frame") %in% class(dat_list[[i]]))) {
      id <- dat_list[[i]]$id
    } else if ("matrix" %in% class(dat_list[[i]])) {
      id <- colnames(dat_list[[i]])
    } else {
      stop("Unrecognized input class.")
    }
    if (is.null(id)) stop("No identifier for samples found.")
    folded <- list()
    for (j in 1:nruns) {
      if (!is.null(cv_stratification_var) & (stratified_cv)) {
        a_ind <- lapply(table(cv_stratification_var), function(x) sample(1:x, x))
        b_ind <- sample(1:length(unique(cv_stratification_var)), length(unique(cv_stratification_var)))
        c_ind <- cumsum(table(cv_stratification_var)[unique(cv_stratification_var)[b_ind]])
        cv_index <- c()
        for (u in 1:length(b_ind)) {
          un <- unique(cv_stratification_var)[b_ind[u]]
          cv_index[cv_stratification_var == un] <- a_ind[[un]] + ifelse(u > 1, c_ind[u-1], 0)
        }
        if (anti_stratified) {
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
        folded[[j]][[f]] <- data.table::data.table(fold = f, run = j, 
                                                   cv_index = cv_index[cv_index != f], 
                                                   id = id[cv_index != f])
      }
      folded[[j]] <- data.table::rbindlist(folded[[j]])
    }
    out[[i]] <- data.table::rbindlist(folded)
  }
  names(out) <- names(dat_list)
  return(out)
}

data_preprocess <- function(dat, verbose = FALSE) {
  if ("list" %in% class(dat)) {
    dat_list <- dat
  } else {
    dat_list <- list(dat)
  }
  if(verbose) print("Processing data sets ..."); flush.console()
  
  # Collect gene names for later. 
  # They are required for pathway-based approaches
  gene_id_list <- lapply(dat_list, rownames)
  
  # Convert data to data.table to optimize memory usage
  for (i in 1:length(dat_list)) {
    id <- colnames(dat_list[[i]])
    dat_list[[i]] <- data.table::as.data.table(t(dat_list[[i]]))
    #data.table::setDT(dat_list[[i]])
    colnames(dat_list[[i]]) <- paste0("dim", 1:ncol(dat_list[[i]]))
    dat_list[[i]]$id <- id
    data.table::setkey(dat_list[[i]], id)
  }
  return(list(dat_list = dat_list, gene_id_list = gene_id_list))
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
#' @param gene_ids a character vector indicating target genes to include. 
#' @param cutoff a numeric value indicating the cutoff for the edge scores.  
#' @param directed a boolean value indicating the type of grpah to be generated. 
#' @param version a character value specifying STRINGdb version to query from 
#' 
#' @return an igraph object.
#' @export
#' @importFrom STRINGdb STRINGdb
#' @importFrom igraph graph_from_edgelist graph_from_adjacency_matrix as_adjacency_matrix
#' @importFrom biomaRt useEnsembl getBM
#' @importFrom Matrix rowSums
getHumanPPIfromSTRINGdb <- function(gene_ids = NULL, cutoff = 700, directed = FALSE, 
                                    version = "11", gene_id_mart_column = "hgnc_symbol") {
  string_db <- STRINGdb::STRINGdb$new(version = version, species=9606, score_threshold = cutoff)
  if(directed) {
    gene.stringdb <- string_db$map(my_data_frame = gene_ids, 
                                   my_data_frame_id_col_names='target.gene_info.symbol', 
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
    g = string_db$get_subnetwork(hs_all_proteins$protein_external_id)
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
    mart_results <- biomaRt::getBM(attributes = c(gene_id_mart_column,"ensembl_peptide_id"),
                                   #filters = "ensembl_peptide_id", values = protein_ids,
                                   mart = mart)
    ### replace protein ids with gene ids
    ix <- match(protein_ids, mart_results$ensembl_peptide_id)
    ix <- ix[!is.na(ix)]
    newnames <- protein_ids
    newnames[match(mart_results[ix,'ensembl_peptide_id'], newnames)] <- mart_results[ix, gene_id_mart_column]
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
#' \code{biomaRt} for matching gene IDs. 
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
                                   directed_ppi = FALSE,
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
#' @importFrom ggplot2 ggtitle
triple_viz <- function(data, category, category_label, tsne_perplexity = 45, umap_neighbors = 20, tsne = FALSE) {
  p1 <- pca_viz(data, category, category_label) + ggtitle("PCA")
  if (tsne) {
    p2 <- tsne_viz(data, category, category_label, tsne_perplexity = tsne_perplexity) + ggtitle("t-SNE")
  } else {
    p2 <- NULL
  }
  p3 <- umap_viz(data, category, category_label, umap_neighbors = umap_neighbors) + ggtitle("UMAP")
  
  return(list(PCA = p1, tSNE = p2, UMAP = p3))
}

#' Data visualization using PCA
#'
#' @param data 
#' @param category 
#' @param category_label 
#'
#' @return
#' @export
#' 
#' @importFrom FactoMineR PCA
#' @importFrom ggplot2 ggplot aes geom_point scale_color_brewer theme_bw labs
pca_viz <- function(data, category, category_label) {
  res_pca <- FactoMineR::PCA(data, scale.unit = FALSE, ncp = 2, graph = FALSE)
  res_pca_dat <- as.data.frame(res_pca$ind$coord)
  res_pca_dat <- cbind(res_pca_dat, category)
  colnames(res_pca_dat)[3] <- "category"
  eig_percentages <- res_pca$eig[,"percentage of variance"]
  eig_percentages <- as.character(signif(eig_percentages, 3))
  p1 <- ggplot(res_pca_dat, aes(Dim.1, Dim.2, color = category)) + geom_point(shape = "+", size = 3) + 
    theme_bw() + scale_color_brewer(palette = "Dark2") + 
    labs(x = paste0("PC1 (", eig_percentages[1], "%)"), y = paste0("PC2 (", eig_percentages[2], "%)"), color = category_label)
  
  return(p1)
}

#' Data visualization using UMAP
#'
#' @param data 
#' @param category 
#' @param category_label 
#' @param umap_neighbors 
#'
#' @return
#' @export
#' 
#' @importFrom uwot umap
#' @importFrom ggplot2 ggplot aes geom_point scale_color_brewer theme_bw labs
umap_viz <- function(data, category, category_label, umap_neighbors = 20) {
  res_umap <- uwot::umap(data, n_neighbors = umap_neighbors, n_components = 2, pca = min(50, dim(data)), verbose = FALSE, init = "normlaplacian")
  res_umap <- data.frame(Dim.1 = res_umap[,1], Dim.2 = res_umap[,2])
  res_umap <- cbind(res_umap, category)
  colnames(res_umap)[3] <- "category"
  p1 <- ggplot(res_umap, aes(Dim.1, Dim.2, color = category)) + geom_point(shape = "+", size = 3) + 
    theme_bw() + scale_color_brewer(palette = "Dark2") + 
    labs(x = "Z1", y = "Z2", color = category_label)
  
  return(p1)
}

#' Data visualization using t-SNE
#'
#' @param data 
#' @param category 
#' @param category_label 
#' @param tsne_perplexity 
#'
#' @return
#' @export
#' 
#' @importFrom Rtsne Rtsne
#' @importFrom ggplot2 ggplot aes geom_point scale_color_brewer theme_bw labs
tsne_viz <- function(data, category, category_label, tsne_perplexity = 45) {
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
  p1 <- ggplot(res_tsne, aes(V1, V2, color = category)) + geom_point(shape = "+", size = 3) + 
    theme_bw() + scale_color_brewer(palette = "Dark2") + 
    labs(x = "Z1", y = "Z2", color = category_label)
  
  return(p1)
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
#' To help manage the scale differences, outliers (log p < Q1 - 1.5*IQR) are omitted. 
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

plot_pairwise_metrics <- function(results, metrics) {
  for (i in 1:(length(metrics)-1)) {
    for (j in (i+1):length(metrics)) {
      
    }
  }
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

nlog10_trans <- scales::trans_new("reverse_log", function(x) -log(x), 
                                  function(y) exp(-y), breaks = scales::log_breaks())

#' Pairwise plots for visual Pareto based multi-objective optimization
#'
#' @param scores 
#' @param plot_palette color codes used for coloring based on \code{color_var}
#' @param metrics 
#' @param color_var 
#' @param shape_var 
#' @param size_var 
#'
#' @return
#' @export
#'
#' @importFrom pals watlington
#' @importFrom ggplot2 ggplot ensym scale_color_manual geom_point theme_bw 
#' @importFrom cowplot get_legend
#' @importFrom gridExtra grid.arrange
#' @importFrom scales trans_new log_breaks
pareto_plot <- function(scores, plot_palette = pals::watlington(16),
                        metrics = c("TrainStabilityJaccard", 
                                    "Silhouette", 
                                    "Smallest_cluster_size"), 
                        metrics_scale = rep("identity", length(metrics)), 
                        color_var = "Transform",
                        shape_var = "Clustering",
                        size_var = "k", 
                        size_range = c(2,6),
                        color_scale = ggplot2::scale_color_manual(values = plot_palette)) {
  if (!"data.frame" %in% class(scores)) {
    if (is.null(scores$all)) {
      stop("Please provide a data.frame of scores")
    } else {
      # Assume that we are dealing with COPS::clusteval_scoring output
      scores <- scores$all
    }
  }
  
  # Convert size_var to numeric to avoid warnings
  size_var_conv <- as.numeric(as.character(scores[[size_var]]))
  if (all(!is.na(size_var_conv))) scores[[size_var]] <- size_var_conv
  
  plot_list <- list()
  # Create list of pairwise metric scatter plots
  for (i in 1:(length(metrics)-1)) {
    for (j in (i+1):length(metrics)) {
      i_name <- metrics[i]
      j_name <- metrics[j]
      i_scale <- switch(metrics_scale[i], identity = "identity", nlog10 = nlog10_trans)
      j_scale <- switch(metrics_scale[j], identity = "identity", nlog10 = nlog10_trans)
      plot_ij <- ggplot(scores, aes(x = !!ggplot2::ensym(i_name), 
                                    y = !!ggplot2::ensym(j_name), 
                                    color = !!ggplot2::ensym(color_var), 
                                    shape = !!ggplot2::ensym(shape_var), 
                                    size = !!ggplot2::ensym(size_var))) + 
        geom_point() + theme_bw() + scale_color_manual(values = plot_palette) + 
        theme(legend.position = "none") + scale_x_continuous(trans = i_scale) + 
        scale_y_continuous(trans = j_scale, position = "right") +
        scale_size(range = size_range)
      plot_list <- c(plot_list, list(plot_ij))
    }
  }
  # Create plot for legend extraction
  legend_plot <- ggplot(scores, aes(x = 1, y = 1, 
                                    color = !!ggplot2::ensym(color_var),
                                    shape = !!ggplot2::ensym(shape_var), 
                                    size = !!ggplot2::ensym(size_var))) + 
    geom_point() + theme_bw() + scale_color_manual(values = plot_palette) +
    theme(legend.box = "horizontal") + scale_size(range = size_range) + 
    guides(shape = guide_legend(ncol = 1, order = 2), 
           color = guide_legend(ncol = 1, order = 1),
           size = guide_legend(ncol = 1, order = 3))
  pareto_legend <- cowplot::get_legend(legend_plot)
  
  # Define layout for comparing metrics
  layout <- matrix(NA, length(metrics) - 1, length(metrics) - 1)
  layout[lower.tri(layout, diag = TRUE)] <- 1:length(plot_list)
  layout <- layout[,(length(metrics)-1):1]
  n <- length(metrics)-1
  m <- floor(n/2)
  layout[1:m,1:m] <- length(plot_list) + 1
  
  plot_list <- c(plot_list, list(pareto_legend))
  
  plot_out <- gridExtra::grid.arrange(grobs = plot_list, layout_matrix = layout)
  
  return(plot_out)
}


#' Reorder factors in scores to organize plots
#'
#' @param x scores from COPS pipeline
#'
#' @return
#' @export
reorder_method_factors <- function(x) {
  gsva_ind <- grepl("GSVA", x$Approach)
  diff_ind <- grepl("DiffRank", x$Approach)
  rwrfgsea_ind <- grepl("RWR-FGSEA", x$Approach)
  other_ind <- !(gsva_ind|diff_ind|rwrfgsea_ind)
  pw_settings <- unique(x$Embedding[!other_ind])
  x$Transform <- factor(x$Transform, c(unique(x$Transform[other_ind]), 
                                     unique(x$Transform[gsva_ind]),
                                     unique(x$Transform[diff_ind]),
                                     unique(x$Transform[rwrfgsea_ind])))
  x$Embedding <- factor(x$Embedding, c(unique(x$Embedding[other_ind]), 
                                       unique(x$Embedding[!other_ind])))
  x$Approach <- factor(x$Approach, c(unique(x$Approach[other_ind]), 
                                     unique(x$Approach[gsva_ind]),
                                     unique(x$Approach[diff_ind]),
                                     unique(x$Approach[rwrfgsea_ind])))
  return(x)
}


#' Renames internal method and score names to something more understandable. 
#'
#' @param x scores from COPS pipeline
#'
#' @return
#' @export
format_scores <- function(x) {
  if (class(x) == "list" & "all" %in% names(x)) {
    out <- list()
    out$all <- format_scores(x$all)
    out$best <- format_scores(x$best)
    return(out)
  }
  # Factor for grouping observations in plots
  x$Method <- paste0(x$datname, "+", x$drname, "+", x$m)
  
  # Approach, either DR or specific PW enrichment method name
  pathway_approaches <- grepl("_RWRFGSEA$|_GSVA$|_DiffRank$", x$datname)
  x$Approach <- NA
  x$Approach[!pathway_approaches] <- "DR"
  x$Approach[grepl("_RWRFGSEA$", x$datname)] <- "RWR-FGSEA"
  x$Approach[grepl("_GSVA$", x$datname)] <- "GSVA"
  x$Approach[grepl("_DiffRank$", x$datname)] <- "DiffRank"
  
  # Embedding, describes the features which are used for clustering
  x$Embedding <- NA
  x$Embedding[pathway_approaches] <- x$datname[pathway_approaches]
  x$Embedding[!pathway_approaches] <- ifelse(!is.na(as.numeric(as.character(x$datname[!pathway_approaches]))), 
                                             "", x$datname[!pathway_approaches])
  x$Embedding[!pathway_approaches] <- paste0(x$drname, "+", x$Embedding[!pathway_approaches])
  # remove "original" tag which is just used to indicate a skipped DR step
  x$Embedding <- gsub("\\+original$", "", x$Embedding)
  x$Embedding <- gsub("\\+$", "", x$Embedding)
  # remove redundant pathway method tags (included in Transform)
  x$Embedding <- gsub("_RWRFGSEA|_GSVA|_DiffRank", "", x$Embedding)
  # format methods and dimension numbers
  x$Embedding <- paste0(gsub("^pca", "PCA, ", x$Embedding), ifelse(grepl("^pca", x$Embedding), "d", ""))
  x$Embedding <- paste0(gsub("^tsne", "t-SNE, ", x$Embedding), ifelse(grepl("^tsne", x$Embedding), "d", ""))
  x$Embedding <- paste0(gsub("^umap", "UMAP, ", x$Embedding), ifelse(grepl("^umap", x$Embedding), "d", ""))
  
  # Transform, same as Embedding except that pathway gene sets are appended with 
  # enrichment method name (used for Pareto plots)
  x$Transform <- NA
  x$Transform[!pathway_approaches] <- x$Embedding
  x$Transform[pathway_approaches] <- gsub("_", " ", x$datname)
  
  # Clustering method
  x$Clustering <- x$m
  x$Clustering <- gsub("model", "GMM", x$Clustering)
  x$Clustering <- gsub("kmeans", "k-means", x$Clustering)
  x$Clustering <- gsub("hierarchical", "HC", x$Clustering)
  x$Clustering <- gsub("_average$", " (average)", x$Clustering)
  x$Clustering <- gsub("_ward$", " (Ward)", x$Clustering)
  x$Clustering <- gsub("_complete$", " (complete)", x$Clustering)
  x$Clustering <- gsub("^diana$", "DIANA", x$Clustering)
  
  # Survival
  colnames(x)[colnames(x) == "cluster_significance"] <- "SurvivalPValue"
  
  # Stability
  colnames(x) <- gsub("^TrainStability", "ClusteringStability", colnames(x))
  colnames(x) <- gsub("^TestStability", "ProjectionClusteringStability", colnames(x))
  
  # Other
  x$k <- factor(x$k)
  
  return(x)
}

multi_view_cv_fold <- function(dat_list, nfolds = 5, nruns = 2, ...) {
  if (!all(Reduce("&", lapply(dat_list[-1], function(x) x$id == dat_list[[1]]$id)))) {
    warning("Colnames in all views do not match.")
    stop("Cross-validation for missing sample views not implemented.")
  }
  return(cv_fold(dat_list = dat_list[1], nfolds = nfolds, nruns = nruns, ...))
}

subset_cv_data <- function(dat_list, cv_index, data_is_kernels = FALSE) {
  dat_i <- list()
  non_data_cols <- list()
  if(data_is_kernels) {
    for (j in 1:length(dat_list)) {
      if (sum(grepl("^dim[0-9]+$", colnames(dat_list[[j]]))) > nrow(dat_list[[j]])) {
        stop("Input kernels are not square!")
      }
      ij_ind <- match(cv_index$id, dat_list[[j]]$id)
      dat_i[[j]] <- as.matrix(as.data.frame(dat_list[[j]])[ij_ind, paste0("dim", ij_ind)])
      
      temp <- merge(cv_index, dat_list[[j]], by = "id")
      sel <- grep("^dim[0-9]+$", colnames(temp))
      if ("data.table" %in% class(temp)) {
        non_data_cols[[j]] <- temp[,-..sel]
      } else {
        non_data_cols[[j]] <- temp[,-sel]
      }
    }
  } else {
    for (j in 1:length(dat_list)) {
      dat_i[[j]] <- merge(cv_index, dat_list[[j]], by = "id")
      sel <- grep("^dim[0-9]+$", colnames(dat_i[[j]]))
      if ("data.table" %in% class(dat_i[[j]])) {
        non_data_cols[[j]] <- dat_i[[j]][,-..sel]
        dat_i[[j]] <- as.matrix(dat_i[[j]][,..sel])
      } else {
        non_data_cols[[j]] <- dat_i[[j]][,-sel]
        dat_i[[j]] <- as.matrix(dat_i[[j]][,sel])
      }
    }
  }
  names(dat_i) <- names(dat_list)
  names(non_data_cols) <- names(dat_list)
  return(list(dat_i = dat_i, non_data_cols = non_data_cols))
}

