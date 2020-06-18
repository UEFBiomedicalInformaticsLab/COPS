#' Transform a gene-level data matrix into path-level information
#' 
#' Utility function to extract pathway-based features from gene expression data using GSVA.
#' 
#' @param dat a numeric matrix representing gene expression profiles. Genes on the rows and samples on the columns.
#' @param study_batch factor indicating the individual studies. 
#' A NULL value implicates that the input data matrix represents the original dataset (without removing the batches).
#' @param min.size a numeric value indicating the minimum size of gene sets included
#' @param max.size a numeric value indicating the maximum size of gene sets included
#' @param parallel a numeric value indicating the number of processors to use when doing the calculations in parallel.
#' @param verbose controls verbosity
#' @param kcdf distribution name for \code{\link[GSVA]{gsva}} empirical distribution kernel
#' @param gs_subcats a character vector indicating msigdbr gene set subcategory names to include in the analysis
#' 
#' @return a list of three data frames:\cr
#'        - \strong{KEGG_PW}: a character variable containing entrez gene ids;\cr
#'        - \strong{GO_PW}: a character variable containing gene symbols;\cr
#'        - \strong{REACTOME_PW}: a numeric variable containing gene-disease scores.
#' @export
#' @importFrom AnnotationDbi mapIds
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @importFrom dplyr filter
#' @importFrom msigdbr msigdbr
#' @importFrom GSVA gsva
fromGeneToPathwayFeatures <- function(dat, study_batch = NULL, 
                                      min.size = 5, max.size = 200, 
                                      parallel = 4,
                                      verbose = FALSE,
                                      key_name = "ENSEMBL",
                                      kcdf = "Gaussian",
                                      gs_subcats = c("BP", "MF", "CP:KEGG", "CP:REACTOME")
                                      #rnaseq = FALSE # not implemented (bioconductor GSVA behind GitHub)
                                      ) {
  # extract pathways information from msigdb (https://www.gsea-msigdb.org/)
  db_annots = msigdbr::msigdbr(species = "Homo sapiens")
  db_annots <- dplyr::filter(db_annots, grepl(paste(gs_subcats, collapse = "|"), gs_subcat))
  list_db_annots <- lapply(split(db_annots, db_annots$gs_name), function(x) x$gene_symbol)
  list_db_annots <- list_db_annots[which(sapply(list_db_annots, length) < max.size)]
  
  ke_pathways <- NULL
  go_pathways <- NULL
  re_pathways <- NULL
  if(!is.null(study_batch)) {
    if(verbose) print("The dataset in input corresponds to the original dataset")
    list_dat <- lapply(unique(study_batch), function(x) dat[,which(x == as.character(study_batch)), drop = FALSE])
    ke_pathways <- Reduce(plyr::rbind.fill.matrix, lapply(list_dat, function(s) {
      if (key_name != "SYMBOL") rownames(s) <- suppressMessages(as.character(AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, rownames(s), "SYMBOL", key_name)))
      t(suppressWarnings(GSVA::gsva(s, list_db_annots[grep("KEGG", names(list_db_annots))], mx.diff=TRUE, 
                 verbose=FALSE, parallel.sz=parallel, min.sz=min.size, max.sz=max.size, kcdf = kcdf)))#, rnaseq = rnaseq))) # later version for rnaseq?
    }))
    go_pathways <- Reduce(plyr::rbind.fill.matrix, lapply(list_dat, function(s) {
      if (key_name != "SYMBOL") rownames(s) <- suppressMessages(as.character(AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, rownames(s), "SYMBOL", key_name)))
      t(suppressWarnings(GSVA::gsva(s, list_db_annots[grep("GO_", names(list_db_annots))], mx.diff=TRUE, 
                 verbose=FALSE, parallel.sz=parallel, min.sz=min.size, max.sz=max.size, kcdf = kcdf)))#, rnaseq = rnaseq))) # later version for rnaseq?
    }))
    re_pathways <- Reduce(plyr::rbind.fill.matrix, lapply(list_dat, function(s) {
      if (key_name != "SYMBOL") rownames(s) <- suppressMessages(as.character(AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, rownames(s), "SYMBOL", key_name)))
      t(suppressWarnings(GSVA::gsva(s, list_db_annots[grep("REACTOME", names(list_db_annots))], mx.diff=TRUE, 
                 verbose=FALSE, parallel.sz=parallel, min.sz=min.size, max.sz=max.size, kcdf = kcdf)))#, rnaseq = rnaseq))) # later version for rnaseq?
    }))
    ke_pathways <- t(ke_pathways)
    go_pathways <- t(go_pathways)
    re_pathways <- t(re_pathways)
    ke_pathways[is.na(ke_pathways)] <- 0
    go_pathways[is.na(go_pathways)] <- 0
    re_pathways[is.na(re_pathways)] <- 0
  }
  else {
    if (key_name != "SYMBOL") rownames(dat) <- suppressMessages(as.character(AnnotationDbi::mapIds(org.Hs.eg.db, rownames(dat), "SYMBOL", key_name)))
    ke_pathways <- suppressWarnings(GSVA::gsva(dat, list_db_annots[grep("KEGG", names(list_db_annots))], mx.diff=TRUE, 
                        verbose=FALSE, parallel.sz=parallel, max.sz = max.size, kcdf = kcdf))#, rnaseq = rnaseq)) # later version for rnaseq?
    go_pathways <- suppressWarnings(GSVA::gsva(dat, list_db_annots[grep("GO", names(list_db_annots))], mx.diff=TRUE, 
                        verbose=FALSE, parallel.sz=parallel, max.sz = max.size, kcdf =  kcdf))#, rnaseq = rnaseq)) # later version for rnaseq?
    re_pathways <- suppressWarnings(GSVA::gsva(dat, list_db_annots[grep("REACTOME", names(list_db_annots))], mx.diff=TRUE, 
                        verbose=FALSE, parallel.sz=parallel, max.sz = max.size, kcdf = kcdf))#, rnaseq = rnaseq)) # later version for rnaseq?
  }
  colnames(ke_pathways) <- colnames(dat)
  colnames(go_pathways) <- colnames(dat)
  colnames(re_pathways) <- colnames(dat)
  
  return(list(KEGG_PW = ke_pathways,
              GO_PW = go_pathways,
              REACTOME_PW = re_pathways))
}

#' Single sample gene set enrichment based on expression rank analysis
#' 
#' DiffRank by Wang et al. BMC Medical Genomics 2019
#' 
#' @param expr gene expression
#' @param list_db_annots list of gene sets
#' 
#' @return 
#' @export
DiffRank <- function(expr, list_db_annots, parallel = 1) {
  if (parallel > 1) {
    parallel_clust <- parallel::makeCluster(parallel)
    doParallel::registerDoParallel(parallel_clust)
  } else {
    foreach::registerDoSEQ()
  }
  ranks <- apply(expr, 2, function(x) order(order(x)))
  out <- matrix(NA, length(list_db_annots), ncol(expr))
  rownames(out) <- names(list_db_annots)
  colnames(out) <- colnames(expr)
  out <- foreach(i = list_db_annots, 
                .combine = rbind,
                .export = c(),
                .multicombine = TRUE,
                .maxcombine = length(list_db_annots)) %dopar% {
    ind <- which(rownames(expr) %in% i)
    apply(ranks[ind,,drop=FALSE] - nrow(ranks)/2, 2, mean) - apply(ranks[-ind,,drop=FALSE] - nrow(ranks)/2, 2, mean)
  }
  if (parallel > 1) parallel::stopCluster(parallel_clust)
  rownames(out) <- names(list_db_annots)
  return(list(KEGG_DiffRank = out[grep("^KEGG", rownames(out)),], 
              GO_DiffRank = out[grep("^GO", rownames(out)),], 
              REACTOME_DiffRank = out[grep("^REACTOME", rownames(out)),]))
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
#' Conveniently wraps \code{\link{fromGeneToNetworksToPathwayFeatures}} using 
#' \code{\link{retrieveDiseaseGenesOT}} for disease associated genes, 
#' \code{STRINGdb} human PPI as a network, \code{msigdbr} for gene set annotations and 
#' \code{biomaRt} for harmonization.
#'
#' @param dat a gene expression matrix with samples on columns
#' @param disease_id integer ID in Open Targets platform
#' @param otp_cutoff numeric association score cutoff for Open Targets platform
#' @param ppi_cutoff numeric PPI link score cutoff
#' @param ... extra arguments are passed on to \code{\link{fromGeneToNetworksToPathwayFeatures}}
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
  
  out <- fromGeneToNetworksToPathwayFeatures(dat, disease.genes, gene.network, list_db_annots, ...)
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
ecdf_transform <- function(x, parallel = 1) {
  if (parallel > 1) {
    parallel_clust <- parallel::makeCluster(parallel)
    doParallel::registerDoParallel(parallel_clust)
  } else {
    # Avoid warnings related to %dopar%
    foreach::registerDoSEQ()
  }
  x_sd <- apply(x, 1, sd)
  score <- foreach(i = 1:ncol(x), 
                   .combine = cbind,
                   .export = c(),
                   .multicombine = TRUE,
                   .maxcombine = ncol(x)) %dopar% {
                     score_i <- x * 0
                     for (j in (1:ncol(x))[-i]) {
                       # z-score with respect to each kernel
                       score_i[,j] <- pnorm((x[,j] - x[,i]) / x_sd)
                     }
                     # sum over kernels
                     apply(score_i, 1, sum) / (ncol(x) - 1) 
                   }
  out <- score - 0.5
  colnames(out) <- colnames(x)
  if (parallel > 1) parallel::stopCluster(parallel_clust)
  return(out)
}

#' Jaccard index between indicator matrix columns
#'
#' @param x 
#'
#' @return
#' @export
jaccard_matrix <- function(x) {
  A <- t(x) %*% x
  B <- t(x) %*% (1-x)
  return(A / (A + B + t(B)))
}


#' Random walk with restart and FGSEA worker
#'
#' Runs random walk in a given network starting from seed genes selected from most over/under expressed in the data intersected with 
#' known disease genes. The gene affinities are then processed with FGSEA to yield pathway features. 
#'
#' @param dat a gene expression matrix with samples on columns
#' @param disease.genes a character vector containing Entrez IDs of genes associated with the target disease
#' @param gene.network an \code{igraph.object} with nodes matching to \code{dat} rows
#' @param list_db_annots list of character vectors containing gene sets for FGSEA
#' @param top.ranked.genes integer, controls the number of gene seed candidates to intersect with the disease genes
#' @param min.size integer, minimum size of gene sets
#' @param max.size integer, maximum size of gene sets
#' @param parallel integer, number of threads
#' @param verbose boolean, verbosity of \code{dnet::dRWR}
#' @param rwr_restart the restart probability used for RWR. See \code{dnet::dRWR} for more details.
#' @param rwr_norm the way to normalise the adjacency matrix of the input graph. See \code{dnet::dRWR} for more details.
#' @param rwr_cutoff the cuoff value to select the most visited genes.  
#' @param ... extra arguments are ignored
#'
#' @return list of enrichment scores
#' @export
#' 
#' @importFrom foreach foreach %dopar%
#' @importFrom plyr aaply
#' @importFrom igraph V
#' @importFrom dnet dRWR
#' @importFrom fgsea fgsea
fromGeneToNetworksToPathwayFeatures <- function(dat, 
                                                disease.genes, # a character vector
                                                gene.network,  # igraph.object
                                                list_db_annots, # list of pathway annotated gene sets
                                                top.ranked.genes = 100, 
                                                parallel = 1,
                                                verbose = FALSE,
                                                rwr_restart = 0.33,
                                                rwr_norm = "none",
                                                rwr_affinity_norm = "quantile",
                                                rwr_cutoff = 0,
                                                kernelized_sorting = FALSE,
                                                reg_up_down = FALSE,
                                                return_seeds = FALSE,
                                                ...) {
  sample_names <- colnames(dat)
  dat <- dat[intersect(rownames(dat), disease.genes),]
  
  # Rank disease-genes within each sample
  if (kernelized_sorting) {
    expr <- ecdf_transform(dat, parallel)
  } else {
    expr <- dat
  }
  
  if (reg_up_down) {
    gene.seeds.up <- apply(expr, 2, function(s) { #, dg, ntop, gn) {
      sorted.genes <- names(s)[order(s, decreasing = TRUE)]
      out <- intersect(sorted.genes, disease.genes)[1:top.ranked.genes]
      # Remove NA when intersection is smaller than number of genes selected
      out <- out[!is.na(out)]
      out <- as.numeric(igraph::V(gene.network)$name %in% out)
      out
    })#, dg = disease.genes, ntop = top.ranked.genes, gn = gene.network)
    rownames(gene.seeds.up) <- igraph::V(gene.network)$name
    
    # Apply random walk (dnet) to extend the set of genes
    rwr.top.genes.up <- dnet::dRWR(gene.network, setSeeds = gene.seeds.up, normalise = rwr_norm,
                                restart = rwr_restart, normalise.affinity.matrix = rwr_affinity_norm, 
                                parallel = parallel > 1, multicores = parallel, verbose = verbose)
    rownames(rwr.top.genes.up) <- igraph::V(gene.network)$name
    
    gene.seeds.down <- apply(expr, 2, function(s) { #, dg, ntop, gn) {
      sorted.genes <- names(s)[order(s, decreasing = FALSE)]
      out <- intersect(sorted.genes, disease.genes)[1:top.ranked.genes]
      # Remove NA when intersection is smaller than number of genes selected
      out <- out[!is.na(out)]
      out <- as.numeric(igraph::V(gene.network)$name %in% out)
      out
    })#, dg = disease.genes, ntop = top.ranked.genes, gn = gene.network)
    rownames(gene.seeds.down) <- igraph::V(gene.network)$name
    
    # Apply random walk (dnet) to extend the set of genes
    rwr.top.genes.down <- dnet::dRWR(gene.network, setSeeds = gene.seeds.down, normalise = rwr_norm,
                                restart = rwr_restart, normalise.affinity.matrix = rwr_affinity_norm, 
                                parallel = parallel > 1, multicores = parallel, verbose = verbose)
    rownames(rwr.top.genes.down) <- igraph::V(gene.network)$name
    
    # Combine
    rwr.top.genes <- rwr.top.genes.up - rwr.top.genes.down
    gene.seeds <- rbind(gene.seeds.up, gene.seeds.down)
  } else {
    gene.seeds <- apply(expr, 2, function(s) { #, dg, ntop, gn) {
      sorted.genes <- names(s)[order(abs(s), decreasing = TRUE)]
      out <- intersect(sorted.genes, disease.genes)[1:top.ranked.genes]
      # Remove NA when intersection is smaller than number of genes selected
      out <- out[!is.na(out)]
      out <- as.numeric(igraph::V(gene.network)$name %in% out)
      out
    })#, dg = disease.genes, ntop = top.ranked.genes, gn = gene.network)
    rownames(gene.seeds) <- igraph::V(gene.network)$name
    
    # Apply random walk (dnet) to extend the set of genes
    rwr.top.genes <- dnet::dRWR(gene.network, setSeeds = gene.seeds, normalise = rwr_norm,
                                restart = rwr_restart, normalise.affinity.matrix = rwr_affinity_norm, 
                                parallel = parallel > 1, multicores = parallel, verbose = verbose)
    rownames(rwr.top.genes) = igraph::V(gene.network)$name
  }
  
  if (parallel > 1) {
    parallel_clust <- parallel::makeCluster(parallel)
    doParallel::registerDoParallel(parallel_clust)
  } else {
    # Avoid warnings related to %dopar%
    foreach::registerDoSEQ()
  }
  
  # Apply FGSEA
  #res <- array(NA, c(ncol(dat), length(list_db_annots)), 
  #             list(id = colnames(dat), pathway = names(list_db_annots)))
  #for (i in 1:ncol(rwr.top.genes)) {
  # Use foreach to speed up per sample FGSEA
  res <- foreach(genes_i = lapply(1:ncol(dat), function(x) rwr.top.genes[,x]), 
                             .combine = rbind,
                             .export = c(),
                             .multicombine = TRUE,
                             .maxcombine = ncol(dat)) %dopar% {
    #genes_i <- rwr.top.genes[,i]
    # Sum up duplicated gene id:s
    genes_i <- tapply(genes_i, names(genes_i), sum)
    genes_i <- genes_i[abs(genes_i) > rwr_cutoff]
    res_i <- fgsea::fgsea(list_db_annots, 
                          genes_i, 
                          nperm=10000, 
                          maxSize=min(max(sapply(list_db_annots, length)), length(genes_i) - 1), 
                          nproc = 1)
    #res[i,match(res_i$pathway, names(list_db_annots))] <- res_i$NES * (-log10(res_i$padj)) 
    res_out <- rep(NA, length(list_db_annots))
    res_out[match(res_i$pathway, names(list_db_annots))] <- res_i$NES * (-log10(res_i$padj)) 
    res_out
  }
  rownames(res) <- sample_names
  colnames(res) <- names(list_db_annots)
  # Remove pathways with only NA
  res <- res[, apply(res, 2, function(x) !all(is.na(x)))]
  res[is.na(res)] <- 0
  
  if (parallel > 1) parallel::stopCluster(parallel_clust)
  
  genes_out <- as.matrix(rwr.top.genes)
  colnames(genes_out) <- sample_names
  
  out <- list(GENE_RWR = genes_out,
              KEGG_RWR = t(res[,grep("^KEGG", colnames(res))]), 
              GO_RWR = t(res[,grep("^GO", colnames(res))]), 
              REACTOME_RWR = t(res[,grep("^REACTOME", colnames(res))]))
  
  if (return_seeds) {
    out$GENE_SEEDS <- gene.seeds
  }
  
  return(out)
}

#' FGSEA wrapper
#'
#' @param ranking_matrix 
#' @param list_db_annots 
#' @param rwr_cutoff 
#' @param parallel 
#'
#' @return
#' @export
fgsea_wrapper <- function(ranking_matrix, list_db_annots, rwr_cutoff = 0, parallel = 1) {
  if (parallel > 1) {
    parallel_clust <- parallel::makeCluster(parallel)
    doParallel::registerDoParallel(parallel_clust)
  } else {
    # Avoid warnings related to %dopar%
    foreach::registerDoSEQ()
  }
  # Use foreach to speed up per sample FGSEA
  res <- foreach(genes_i = lapply(1:ncol(ranking_matrix), function(i) ranking_matrix[,i]), 
                 .combine = rbind,
                 .export = c(),
                 .multicombine = TRUE,
                 .maxcombine = ncol(ranking_matrix)) %dopar% {
                   # Sum up duplicated gene id:s
                   genes_i <- tapply(genes_i, names(genes_i), sum)
                   genes_i <- genes_i[abs(genes_i) > rwr_cutoff]
                   res_i <- fgsea::fgsea(list_db_annots, 
                                         genes_i, 
                                         nperm=10000, 
                                         maxSize=min(max(sapply(list_db_annots, length)), length(genes_i) - 1), 
                                         nproc = 1)
                   #res[i,match(res_i$pathway, names(list_db_annots))] <- res_i$NES * (-log10(res_i$padj)) 
                   res_out <- rep(NA, length(list_db_annots))
                   res_out[match(res_i$pathway, names(list_db_annots))] <- res_i$NES * (-log10(res_i$padj)) 
                   res_out
                 }
  rownames(res) <- colnames(ranking_matrix)
  colnames(res) <- names(list_db_annots)
  # Remove pathways with only NA
  res <- res[, apply(res, 2, function(x) !all(is.na(x)))]
  res[is.na(res)] <- 0
  
  if (parallel > 1) parallel::stopCluster(parallel_clust)
  return(res)
}



#' Visualize gene seed similarity for given settings
#'
#' @param dat 
#' @param disease_targets 
#' @param n_top_expr 
#' @param n_top_target 
#' @param mart_results 
#' @param gene.network 
#' @param plot_title 
#' @param otp_score 
#' @param kernelized_sorting 
#' @param reg_up_down 
#' @param parallel 
#'
#' @return
#' @export
seed_similarity <- function(dat, disease_targets, n_top_expr, n_top_target, 
                            mart_results, gene.network, plot_title,
                            otp_score = "association_score.datatypes.rna_expression",
                            kernelized_sorting = FALSE, 
                            reg_up_down = FALSE,
                            parallel = 1) {
  gene.diseases <- disease_targets[order(disease_targets[[otp_score]], 
                                         decreasing = TRUE),]
  disease.genes <- mart_results$ensembl_gene_id[match(gene.diseases$entrezId, 
                                                      mart_results$entrezgene_id)]
  disease.genes <- disease.genes[!is.na(disease.genes)][1:min(n_top_target, sum(disease_targets[[otp_score]] > 0))]
  
  dat <- dat[intersect(rownames(dat), disease.genes),]
  
  if (parallel > 1) {
    parallel_clust <- parallel::makeCluster(parallel)
    doParallel::registerDoParallel(parallel_clust)
  } else {
    foreach::registerDoSEQ()
  }
  # Rank disease-genes within each sample
  if (kernelized_sorting) {
    gene_sd <- apply(dat, 1, sd)
    gene_expr_score <- foreach(i = 1:ncol(dat), 
                               .combine = cbind,
                               .export = c(),
                               .multicombine = TRUE,
                               .maxcombine = ncol(dat)) %dopar% {
                                 gene_expr_score_i <- dat * 0
                                 for (j in (1:ncol(dat))[-i]) {
                                   # z-score with respect to each kernel
                                   gene_expr_score_i[,j] <- pnorm((dat[,j] - dat[,i]) / gene_sd)
                                 }
                                 # sum over kernels
                                 apply(gene_expr_score_i, 1, sum) / (ncol(dat) - 1) 
                               }
    expr <- gene_expr_score - 0.5
  } else {
    expr <- dat
  }
  
  if (reg_up_down) {
    gene.seeds.up <- apply(expr, 2, function(s) {
      sorted.genes <- names(s)[order(s, decreasing = TRUE)]
      out <- intersect(sorted.genes, disease.genes)[1:n_top_expr]
      # Remove NA when intersection is smaller than number of genes selected
      out <- out[!is.na(out)]
      out <- as.numeric(igraph::V(gene.network)$name %in% out)
      out
    })
    
    gene.seeds.down <- apply(expr, 2, function(s) {
      sorted.genes <- names(s)[order(s, decreasing = FALSE)]
      out <- intersect(sorted.genes, disease.genes)[1:n_top_expr]
      # Remove NA when intersection is smaller than number of genes selected
      out <- out[!is.na(out)]
      out <- as.numeric(igraph::V(gene.network)$name %in% out)
      out
    })
    
    # Combine
    gene.seeds <- rbind(gene.seeds.up, gene.seeds.down)
  } else {
    gene.seeds <- apply(expr, 2, function(s) {
      sorted.genes <- names(s)[order(abs(s), decreasing = TRUE)]
      out <- intersect(sorted.genes, disease.genes)[1:n_top_expr]
      # Remove NA when intersection is smaller than number of genes selected
      out <- out[!is.na(out)]
      out <- as.numeric(igraph::V(gene.network)$name %in% out)
      out
    })
  }
  
  # Fast matrix based Jaccard calculation
  A <- t(gene.seeds) %*% gene.seeds
  B <- t(gene.seeds) %*% (1-gene.seeds)
  seeds_jaccard <- A / (A + B + t(B))
  
  if (parallel > 1) parallel::stopCluster(parallel_clust)
  
  # Set samples with 0 gene seeds to zero
  seeds_jaccard[is.na(seeds_jaccard)] <- 0
  gplots::heatmap.2(seeds_jaccard, symm = TRUE, dist = as.dist, density.info = "none", 
                    trace = "none", main = plot_title, breaks = seq(0,1,0.01), revC = TRUE)
}
