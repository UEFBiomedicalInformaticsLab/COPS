#' Transform a gene-level data matrix into path-level information
#' 
#' This is a utility function which wraps several methods to extract pathway-based features from gene expression data
#' 
#' @param expr a numeric matrix representing gene expression profiles. Genes on the rows and samples on the columns.
#' @param enrichment_method options: "DiffRank", "GSVA", "RWRFGSEA"
#' @param gene_set_list list of gene sets (character vectors) corresponding to pathway annotations or gene sets, 
#'   if \code{NULL} gene sets are retrieved with \code{\link[msigdbr]{msigdbr}}.
#' @param batch_label_pw batch labels for batch-wise enrichment, normal enrichment performed if NULL
#' @param min_size a numeric value indicating the minimum size of gene sets included
#' @param max_size a numeric value indicating the maximum size of gene sets included
#' @param parallel a numeric value indicating the number of processors to use when doing the calculations in parallel.
#' @param verbose controls verbosity
#' @param gene_key_expr if \code{is.null(gene_set_list)} and \code{expr} rownames are not gene symbols, this specifies the column name in 
#'   \code{\link[org.Hs.eg.db]{org.Hs.eg.db}} to translate pathway gene symbols to. The default value results in gene symbol based gene 
#'   sets when \code{is.null(gene_set_list)}. 
#' @param gsva_kcdf distribution name for \code{\link[GSVA]{gsva}} empirical probability distribution kernel.
#' @param gs_subcats if \code{is.null(gene_set_list)}, this character vector indicates \code{\link[msigdbr]{msigdbr}}) gene set subcategory 
#'   names that are included in the analysis.
#' @param ... extra arguments are passed to \code{\link[COPS]{RWRFGSEA}}
#' 
#' @return a list of data.frames corresponding to the transformed features based on the selected gene sets (only supports GO, KEGG and REACTOME at the moment)
#' @export
#' 
#' @examples library(parallel)
#' library(COPS)
#' 
#' ## GSVA example
#' ad_gsva <- genes_to_pathways(as.matrix(ad_ge_micro_zscore), "GSVA", 
#'                              parallel = 2, gene_key_expr = "ENSEMBL", 
#'                              gs_subcats = "CP:KEGG")
#' # batch-wise
#' ad_gsva <- genes_to_pathways(as.matrix(ad_ge_micro_zscore), "GSVA", 
#'                              batch_label_pw = ad_studies, parallel = 2, 
#'                              gene_key_expr = "ENSEMBL", gs_subcats = "CP:KEGG")
#' 
#' ## DiffRank example with batch-wise enrichment
#' ad_diffrank <- genes_to_pathways(ad_ge_micro_zscore, "DiffRank", 
#'                                  parallel = 2, gene_key_expr = "ENSEMBL", 
#'                                  gs_subcats = "CP:KEGG")
#' # batch-wise
#' ad_diffrank <- genes_to_pathways(ad_ge_micro_zscore, "DiffRank", 
#'                                  batch_label_pw = ad_studies, 
#'                                  parallel = 2, gene_key_expr = "ENSEMBL", 
#'                                  gs_subcats = "CP:KEGG")
#' 
#' ## RWRFGSEA example
#' ad_data <- ad_ge_micro_zscore
#' rownames(ad_data) <- as.character(AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, 
#'                                                         rownames(ad_data), "SYMBOL", "ENSEMBL"))
#' ad_wgcna_net <- coexpression_network_unweighted(ad_data)
#' kegg_annotations <- msigdbr::msigdbr(species = "Homo sapiens", 
#'                                      category = "C2", 
#'                                      subcategory = "CP:KEGG")
#' list_kegg_annotations <- lapply(split(kegg_annotations, kegg_annotations$gs_name), 
#'                                 function(x) x$human_gene_symbol) 
#' 
#' ad_rwrfgsea <- genes_to_pathways(ad_data, "RWRFGSEA", gene_set_list = list_kegg_annotations, 
#'                                  gene_network = ad_wgcna_net, parallel = 2)
#' 
#' @importFrom AnnotationDbi mapIds
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @importFrom dplyr filter
#' @importFrom msigdbr msigdbr
#' @importFrom GSVA gsva
genes_to_pathways <- function(expr, 
                              enrichment_method = "GSVA",
                              gene_set_list = NULL,
                              batch_label_pw = NULL, 
                              min_size = 5, 
                              max_size = 200, 
                              parallel = 1,
                              verbose = FALSE,
                              gene_key_expr = "SYMBOL",
                              gs_subcats = c("BP", "MF", "CP:KEGG", "CP:REACTOME"),
                              gsva_kcdf = "Gaussian",
                              ...
) {
  #expr <- t(expr)
  if (is.null(gene_set_list)) {
    # extract pathways information from msigdb (https://www.gsea-msigdb.org/)
    db_annots <- msigdbr::msigdbr(species = "Homo sapiens")
    db_annots <- dplyr::filter(db_annots, grepl(paste(gs_subcats, collapse = "|"), gs_subcat))
    
    if (gene_key_expr != "SYMBOL") {
      db_annots$gene_id <- suppressMessages(as.character(AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, 
                                                                               db_annots$human_gene_symbol, 
                                                                               column = gene_key_expr, 
                                                                               keytype = "SYMBOL")))
    } else {
      db_annots$gene_id <- db_annots$human_gene_symbol
    }
    
    gene_set_list <- lapply(split(db_annots, db_annots$gs_name), function(x) x$gene_id)
  }
  
  gene_set_list <- gene_set_list[which(sapply(gene_set_list, length) <= max_size & sapply(gene_set_list, length) >= min_size)]
  
  if(!is.null(batch_label_pw)) {
    # Batch-wise enrichment analysis
    if (any(table(batch_label_pw) < 10)) {
      warning("Batch-wise pathway enrichment results may be unreliable due to low number of samples in some batches.")
    }
    if (any(table(batch_label_pw) < 2)) {
      stop("Batch-wise pathway enrichment failed, all batches must have at least 2 samples.")
    }
    batch_dat_list <- lapply(unique(batch_label_pw), function(x) expr[, which(x == as.character(batch_label_pw)), drop = FALSE])
    
    # Apply the specified enrichment methods
    if (enrichment_method == "GSVA") {
      enriched_dat <- lapply(batch_dat_list, GSVA::gsva, gset.idx.list = gene_set_list, mx.diff = TRUE, 
                             verbose = verbose, parallel.sz = parallel, kcdf = gsva_kcdf,
                             min.sz = min_size, max.sz = max_size)
    } else if (enrichment_method == "DiffRank") {
      enriched_dat <- lapply(batch_dat_list, DiffRank, gene_set_list = gene_set_list, parallel = parallel)
    } else if (enrichment_method == "RWRFGSEA") {
      enriched_dat <- lapply(batch_dat_list, RWRFGSEA, gene_set_list = gene_set_list, parallel = parallel, verbose = verbose, ...)
    } else {
      stop(paste("Unsupported pathway enrichment method:", enrichment_method))
    }
    enriched_dat <- t(Reduce(plyr::rbind.fill.matrix, lapply(enriched_dat, t)))
  } else {
    # Regular enrichment analysis
    if (enrichment_method == "GSVA") {
      enriched_dat <- suppressWarnings(GSVA::gsva(expr, gset.idx.list = gene_set_list, mx.diff = TRUE, 
                                                  verbose = verbose, parallel.sz = parallel, 
                                                  kcdf = gsva_kcdf,
                                                  min.sz = min_size, max.sz = max_size))#, rnaseq = rnaseq)) # later version for rnaseq?
    } else if (enrichment_method == "DiffRank") {
      enriched_dat <- COPS::DiffRank(expr, gene_set_list, parallel)
    } else if (enrichment_method == "RWRFGSEA") {
      enriched_dat <- RWRFGSEA(expr, gene_set_list = gene_set_list, parallel = parallel, verbose = verbose, ...)
    } else {
      stop(paste("Unsupported pathway enrichment method:", enrichment_method))
    }
  }
  
  # Format output
  if (enrichment_method ==  "GSVA") {
    out <- list()
    out$KEGG_GSVA <- enriched_dat[grep("^KEGG_", rownames(enriched_dat)),, drop = FALSE]
    out$GO_GSVA <- enriched_dat[grep("^GO_|^GOBP_", rownames(enriched_dat)),, drop = FALSE]
    out$REACTOME_GSVA <- enriched_dat[grep("^REACTOME_", rownames(enriched_dat)),, drop = FALSE]
  } else if (enrichment_method ==  "DiffRank") {
    out <- list()
    out$KEGG_DiffRank <- enriched_dat[grep("^KEGG", rownames(enriched_dat)),, drop = FALSE]
    out$GO_DiffRank <- enriched_dat[grep("^GO_|^GOBP_", rownames(enriched_dat)),, drop = FALSE]
    out$REACTOME_DiffRank <- enriched_dat[grep("^REACTOME", rownames(enriched_dat)),, drop = FALSE]
  } else if (enrichment_method ==  "RWRFGSEA") {
    out <- list()
    out$KEGG_RWRFGSEA <- enriched_dat[grep("^KEGG_", rownames(enriched_dat)),, drop = FALSE]
    out$GO_RWRFGSEA <- enriched_dat[grep("^GO_|^GOBP_", rownames(enriched_dat)),, drop = FALSE]
    out$REACTOME_RWRFGSEA <- enriched_dat[grep("^REACTOME_", rownames(enriched_dat)),, drop = FALSE]
  } else {
    # This is never run
    out <- enriched_dat
  }
  
  out <- out[sapply(out, nrow) > 0]
  
  return(out)
}

#' Cross-validation wrapped pathway enrichment
#'
#' This function wraps the pathway enrichment based data transformations for use in a 'cross-validated' clustering pipeline. 
#' Assumes a parallel backend has been registered for \code{\link[foreach]{foreach}} with \code{\link[doParallel]{registerDoParallel}}. 
#' Setting \code{parallel} does nothing. 
#'
#' @param dat_list list of data sets corresponding gene expression matrices where the gene names have been replaced by dimension IDs (dim1, dim2, ...).
#' @param cv_index list of data.frames corresponding to cross-validation fold indicators as produced by \code{\link[COPS]{cv_fold}}
#' @param gene_id_list list of gene name vectors of the corresponding columns in dat_list
#' @param parallel ignored and set to 1 for spawned subprocesses
#' @param ... arguments passed on to \code{\link[COPS]{genes_to_pathways}}
#'
#' @return
#' @export
cv_pathway_enrichment <- function(dat_list, cv_index, gene_id_list, parallel = 1, ...) {
  temp_list <- list()
  for (i in 1:length(cv_index)) {
    temp <- cv_index[[i]]
    datname <- names(cv_index)[i]
    if (is.null(datname)) datname <- i
    temp$datname <- datname
    temp <- split(temp, by = c("run", "fold"))
    temp <- lapply(temp, function(x) as.data.frame(merge(dat_list[[datname]], x, by = "id")))
    temp <- lapply(temp, function(x) list(expr = x, gene_ids = gene_id_list[[i]]))
    temp_list <- c(temp_list, temp)
  }
  
  out <- foreach(i = temp_list, 
                 .combine = c,
                 .export = c("genes_to_pathways"), #"dat_list"),
                 .packages = c("GSVA", "fgsea", "dnet", "msigdbr", "AnnotationDbi", "org.Hs.eg.db")) %dopar% {
                   sel <- grep("^dim[0-9]+$", colnames(i$expr))
                   temp <- i$expr[, sel]
                   colnames(temp) <- i$gene_ids
                   # Parallel = 1 for subprocesses
                   pw_temp <- genes_to_pathways(t(temp), parallel = 1, ...)
                   if(length(pw_temp) == 0) stop("Pathway enrichment failed, empty result.")
                   pw_temp <- lapply(pw_temp, function(x) {
                     temp2 <- as.data.frame(t(x))
                     if (ncol(temp2)>0) {
                       colnames(temp2) <- paste0("dim", 1:ncol(temp2)) # replace pw names with dimX for compatibility
                       return(cbind(i$expr[,-sel], temp2))
                     } else {
                       return(data.frame())
                     }
                   })
                   # It is possible for the enrichment to return 0 enriched features.
                   # Here we remove all cv-folds where this happens, which is
                   # unsafe if reference fold is removed but other folds remain. 
                   # TODO: investigate this issue further
                   pw_temp <- pw_temp[sapply(pw_temp, ncol) > 0] 
                   for (j in 1:length(pw_temp)) {
                     pw_temp[[j]]$datname <- names(pw_temp)[j]
                   }
                   pw_temp
                 }
  return(out)
}

#' @describeIn genes_to_pathways DiffRank by Wang et al. BMC Medical Genomics 2019
#' 
#' @param expr gene expression matrix, samples on columns and genes on rows
#' @param gene_set_list list of gene sets with gene names that correspond to rows in \strong{expr}
#' @param parallel a numeric value indicating the number of processors to use when doing the calculations in parallel.
#' 
#' @return 
#' @export
DiffRank <- function(expr, gene_set_list, parallel = 1) {
  if (parallel > 1) {
    parallel_clust <- parallel::makeCluster(parallel)
    doParallel::registerDoParallel(parallel_clust)
  } else {
    foreach::registerDoSEQ()
  }
  ranks <- apply(expr, 2, function(x) order(order(x)))
  out <- matrix(NA, length(gene_set_list), ncol(expr))
  rownames(out) <- names(gene_set_list)
  colnames(out) <- colnames(expr)
  out <- foreach(i = gene_set_list, 
                .combine = rbind,
                .export = c(),
                .multicombine = TRUE,
                .maxcombine = length(gene_set_list)) %dopar% {
    ind <- which(rownames(expr) %in% i)
    # Compare mean ranks of pw genes vs non-pw genes
    apply(ranks[ind,,drop=FALSE] - nrow(ranks)/2, 2, mean) - apply(ranks[-ind,,drop=FALSE] - nrow(ranks)/2, 2, mean)
  }
  if (parallel > 1) parallel::stopCluster(parallel_clust)
  rownames(out) <- names(gene_set_list)
  out <- out[!apply(out, 1, function(x) all(is.na(x))),]
  
  return(out)
  #return(list(KEGG_DiffRank = out[grep("^KEGG", rownames(out)),], 
  #            GO_DiffRank = out[grep("^GO_|^GOBP_", rownames(out)),], 
  #            REACTOME_DiffRank = out[grep("^REACTOME", rownames(out)),]))
}

