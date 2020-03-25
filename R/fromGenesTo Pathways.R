###########################
####   COPS functions   ###
###########################

### TESTING FUNCTIONS
fgtpw_test <- function() {
  provola = getHumanPPIfromSTRINGdb(gene.diseases, directed = FALSE) # TODO: fix error on directed = FALSE
  
  pso_comb_centered <- Reduce("rbind", lapply(split(as.data.frame(t(pso_comb_unscaled)), batch), 
                                              function(x) scale(x, scale = FALSE)))
  
  datRW <- COPS:::expressionToRWFeatures(pso_comb, "EFO_0000676", otp_cutoff = 0.6) # Psoriasis
  #datRW <- COPS:::expressionToRWFeatures(tbrca, "EFO_0000305") # Breast cancer
  datRW[is.na(datRW)] <- 0 
  # Split by pw annotation source (KEGG, GO, REACTOME)
  datRW_list <- split(as.data.frame(t(datRW)), sapply(strsplit(colnames(datRW), "_"), function(x) x[1]))
  names(datRW_list) <- paste0(names(datRW_list), "_RWR")
  
  datGSVA <- COPS:::fromGeneToPathwayFeatures(pso_comb_unscaled, batch, parallel = 10)
  
  res <- COPS::dimred_clusteval_pipeline(c(list(zscore = pso_comb_scaled), datRW_list, datGSVA), 
                                         batch_label = batch, 
                                         include_original = FALSE, 
                                         parallel = 10, nruns = 10, 
                                         nfolds = 5, dimred_methods = c("pca", "umap"),
                                         cluster_methods = c("hierarchical", "kmeans"),
                                         metric = "euclidean", 
                                         n_clusters = 2:6)
  scores <- COPS::clusteval_scoring(res, w_connectivity = 0, w_dunn = 0, w_silhouette = 0)
  best <- COPS::dimred_cluster(c(list(zscore = pso_comb_scaled), datRW_list, datGSVA), scores$best)
  GGally::ggpairs(as.data.frame(best$embedding), aes(color = factor(batch), alpha = 0.4)) + ggtitle("Batches in best embeddings")
  GGally::ggpairs(as.data.frame(COPS::dimred_cluster(c(list(zscore = pso_comb_scaled), datRW_list, datGSVA), 
                                                     scores$all[5,])$embedding), aes(color = factor(batch), alpha = 0.4))
  
  # Heatmaps of gene seeds
  assoc_score_fields = paste(paste("&fields=", c('disease.efo_info.label',
                                                 'disease.efo_info.therapeutic_area',
                                                 'target.gene_info.symbol',
                                                 'association_score.overall',
                                                 'disease.id',
                                                 'association_score.datatypes'), sep=''), collapse = "")
  disease_otp = retrieveDiseaseGenesOT(c(disease_id), assoc_score_fields)[[1]][,-c(10:13)]
  
  
  mart <- biomaRt::useEnsembl("ensembl", "hsapiens_gene_ensembl", "useast.ensembl.org")
  dis_mart_results <- biomaRt::getBM(attributes = c("ensembl_gene_id", "entrezgene_id"),
                                     filters = "entrezgene_id", 
                                     values = disease_otp$entrezId,
                                     mart = mart)
  otp_cutoff <- 0.6
  gene.diseases = disease_otp[which(disease_otp$association_score.overall > otp_cutoff),]
  disease.genes <- mart_results$ensembl_gene_id[match(gene.diseases$entrezId, 
                                                      mart_results$entrezgene_id)]
  disease.genes <- disease.genes[!is.na(disease.genes)]
  sum(rownames(dat) %in% disease.genes)
  mean(disease.genes %in% rownames(dat))
  
  top.ranked.genes <- nrow(dat) %/% 20
  gene.seeds <- plyr::aaply(dat, 2, function(s) {
    sorted.genes <- names(s)[order(abs(s), decreasing = TRUE)[1:top.ranked.genes]]
    out <- intersect(sorted.genes, disease.genes)
    out <- as.numeric(igraph::V(gene.network)$name %in% out)
    out
  })
  
  
  f_jaccard <- function(x, y) {sum(x & y) / sum(x | y)}
  seeds_jaccard <- array(NA, dim = c(nrow(gene.seeds), nrow(gene.seeds)))
  diag(seeds_jaccard) <- 1
  for (i in 1:(nrow(gene.seeds)-1)) {
    for (j in (i+1):nrow(gene.seeds)) {
      seeds_jaccard[i,j] <- f_jaccard(gene.seeds[i,], gene.seeds[j,])
      seeds_jaccard[j,i] <- seeds_jaccard[i,j]
    }
  }
  # Set samples with 0 gene seeds to zero
  seeds_jaccard[is.na(seeds_jaccard)] <- 0
  heatmap.2(seeds_jaccard, symm = TRUE, dist = as.dist, density.info = "none", trace = "none")
}


#'Transform a gene-level data matrix into path-level information
#'
#'Utility function to extract pathway-based features from gene expression data.
#'
#'@param dat a numeric matrix representing gene expression profiles. Genes on the rows and samples on the columns.
#'@param study_batch factor indicating the individual studies. 
#'A NULL value implicates that the input data matrix represents the original dataset (without removing the batches).
#'@param min.size a numeric value indicating the minimum size of the resulting gene sets.
#'@param max.size a numeric value indicating the maximum size of the resulting gene sets.
#'@param parallel a numeric value indicating the umber of processors to use when doing the calculations in parallel.
#'@param verbose controls verbosity
#'@param kcdf distribution name for \code{\link[GSVA]{gsva}} empirical distribution kernel
#'
#'@return a list of three data frames:\cr
#'        - \strong{KEGG_PW}: a character variable containing entrez gene ids;\cr
#'        - \strong{GO_PW}: a character variable containing gene symbols;\cr
#'        - \strong{REACTOME_PW}: a numeric variable containing gene-disease scores.
#'@export
#'@importFrom AnnotationDbi mapIds
#'@importFrom dplyr filter '%>%'
#'@importFrom msigdbr msigdbr
#'@importFrom GSVA gsva
fromGeneToPathwayFeatures <- function(dat, study_batch = NULL, 
                                      min.size = 5, max.size = 200, 
                                      parallel = 4,
                                      verbose = FALSE,
                                      kcdf = "Gaussian") {
  # extract pathways information from msigdb (https://www.gsea-msigdb.org/)
  db_annots = msigdbr::msigdbr(species = "Homo sapiens") %>% 
                             dplyr::filter(gs_subcat == "BP" | gs_subcat == "MF" | 
                                           gs_subcat == "CP:KEGG" | gs_subcat == "CP:REACTOME")
  list_db_annots <- lapply(split(db_annots, db_annots$gs_name), function(x) x$gene_symbol)
  list_db_annots <- list_db_annots[which(sapply(list_db_annots, length) < max.size)]
  
  ke_pathways <- NULL
  go_pathways <- NULL
  re_pathways <- NULL
  if(!is.null(study_batch)) {
    if(verbose) print("The dataset in input corresponds to the original dataset")
    list_dat <- lapply(unique(study_batch), function(x) dat[,which(x == as.character(study_batch))])
    ke_pathways <- Reduce(plyr::rbind.fill.matrix, lapply(list_dat, function(s) {
      #rownames(s) <- suppressMessages(as.character(AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, rownames(s), "SYMBOL", "ENSEMBL")))
      t(suppressWarnings(GSVA::gsva(s, list_db_annots[grep("KEGG", names(list_db_annots))], mx.diff=TRUE, 
                 verbose=FALSE, parallel.sz=parallel, min.sz=min.size, max.sz=max.size, kcdf = kcdf)))
    }))
    go_pathways <- Reduce(plyr::rbind.fill.matrix, lapply(list_dat, function(s) {
      #rownames(s) <- suppressMessages(as.character(AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, rownames(s), "SYMBOL", "ENSEMBL")))
      t(suppressWarnings(GSVA::gsva(s, list_db_annots[grep("GO_", names(list_db_annots))], mx.diff=TRUE, 
                 verbose=FALSE, parallel.sz=parallel, min.sz=min.size, max.sz=max.size, kcdf = kcdf)))
    }))
    re_pathways <- Reduce(plyr::rbind.fill.matrix, lapply(list_dat, function(s) {
      #rownames(s) <- suppressMessages(as.character(AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, rownames(s), "SYMBOL", "ENSEMBL")))
      t(suppressWarnings(GSVA::gsva(s, list_db_annots[grep("REACTOME", names(list_db_annots))], mx.diff=TRUE, 
                 verbose=FALSE, parallel.sz=parallel, min.sz=min.size, max.sz=max.size, kcdf = kcdf)))
    }))
    ke_pathways <- t(ke_pathways)
    go_pathways <- t(go_pathways)
    re_pathways <- t(re_pathways)
    ke_pathways[is.na(ke_pathways)] <- 0
    go_pathways[is.na(go_pathways)] <- 0
    re_pathways[is.na(re_pathways)] <- 0
  }
  else {
    rownames(dat) <- suppressMessages(as.character(AnnotationDbi::mapIds(org.Hs.eg.db, rownames(dat), "SYMBOL", "ENSEMBL")))
    ke_pathways <- suppressWarnings(gsva(dat, list_ad_m_df[grep("KEGG", names(list_ad_m_df))], mx.diff=TRUE, 
                        verbose=FALSE, parallel.sz=parallel, max.sz = max.size, kcdf = kcdf))
    go_pathways <- suppressWarnings(gsva(dat, list_ad_m_df[grep("GO", names(list_ad_m_df))], mx.diff=TRUE, 
                        verbose=FALSE, parallel.sz=parallel, max.sz = max.size, kcdf =  kcdf))
    re_pathways <- suppressWarnings(gsva(dat, list_ad_m_df[grep("REACTOME", names(list_ad_m_df))], mx.diff=TRUE, 
                        verbose=FALSE, parallel.sz=parallel, max.sz = max.size, kcdf = kcdf))
  }
  colnames(ke_pathways) <- colnames(dat)
  colnames(go_pathways) <- colnames(dat)
  colnames(re_pathways) <- colnames(dat)
  
  return(list(KEGG_PW = ke_pathways,
              GO_PW = go_pathways,
              REACTOME_PW = re_pathways))
}


#'Retrieve disease-gene associations from the Open Targets platforms. 
#'
#'Utility function to extract relevant disease-gene associations
#'
#'@param diseases a character vector indicating the disease names. 
#'@param fields a character vector indicating the data types used for the infrence of disease-gene associations.
#'Check the Open Target platforms for more details.
#' 
#'@return a data frame including disease-gene association found for each specified data type.
#'@export
#'@importFrom httr content
#'@importFrom jsonlite fromJSON
#'@importFrom AnnotationDbi mapIds
retrieveDiseaseGenesOT <- function(diseases, fields) {
  server <- 'https://platform-api.opentargets.io/v3/platform'
  endpoint_prmtrs <- '/public/association/filter'
  OT_list <- list()
  for(d in diseases) {
    optional_prmtrs <- paste('?size=10000&disease=', d, fields, sep="")
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


#'Retrieve a protein-protein interaction network from STRINGdb.
#'
#'Utility function to extract a gene subnetwork from STRINGdb including only the seed genes and their interactions. 
#'@param gene.diseases a character vector indicating the gene seeds. 
#'@param cutoff a numeric value indicating the cutoff for the edge scores.  
#'@param directed a boolean value indicating the type of grpah to be generated. 
#' 
#'@return an igraph object.
#'@export
#'@import STRINGdb
getHumanPPIfromSTRINGdb <- function(gene.diseases, cutoff = 700, directed = FALSE) {
  string_db <- STRINGdb::STRINGdb$new(version="10", species=9606, score_threshold = cutoff)
  if(directed == TRUE) {
    gene.stringdb <- string_db$map(my_data_frame = gene.diseases, my_data_frame_id_col_names='target.gene_info.symbol', 
                                     removeUnmappedRows=T)
    string_inter <- string_db$get_interactions(gene.stringdb$STRING_id)
    idx_from <- match(x = string_inter$from, table = gene.stringdb$STRING_id)
    idx_to <- match(x = string_inter$to, table = gene.stringdb$STRING_id)
    print(idx_from)
    print(idx_to)
    print(gene.stringdb$target.gene_info.symbol[idx_to])
    print(gene.stringdb$target.gene_info.symbol[idx_from])
    ppi_network <- data.frame(node1=gene.stringdb$target.gene_info.symbol[idx_from], 
                              node2=gene.stringdb$target.gene_info.symbol[idx_to], 
                              weight=string_inter$combined_score/1000,stringsAsFactors = F)
    g.ppi_network <- igraph::graph_from_edgelist(as.matrix(ppi_network[,1:2]), directed=T)
  } 
  else {
    hs_all_proteins = string_db_hs$get_proteins(); rownames(hs_all_proteins) = hs_all_proteins[,1]
    hs_ex_proteins = hs_all_proteins[which((hs_all_proteins$preferred_name %in% 
                                              rownames(gene.diseases)) == TRUE),]; 
    g = string_db_hs$get_subnetwork(hs_ex_proteins$protein_external_id)
    # create adjacency matrix
    adj_matrix <- as_adjacency_matrix(g)
    # map gene ids to protein ids
    # get gene/protein ids via Biomart
    mart <- biomaRt::useMart(host = 'grch37.ensembl.org', 
                             biomart='ENSEMBL_MART_ENSEMBL', 
                             dataset='hsapiens_gene_ensembl')
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

expressionToRWFeatures <- function(dat, 
                                   disease_id, 
                                   otp_cutoff = 0.8, 
                                   ppi_cutoff = 700,
                                   pw_max.size = 200,
                                   ...) {
  assoc_score_fields = paste(paste("&fields=", c('disease.efo_info.label',
                                                 'disease.efo_info.therapeutic_area',
                                                 'target.gene_info.symbol',
                                                 'association_score.overall',
                                                 'disease.id',
                                                 'association_score.datatypes'), sep=''), collapse = "")
  disease_otp = retrieveDiseaseGenesOT(c(disease_id), assoc_score_fields)[[1]][,-c(10:13)]
  gene.diseases = disease_otp[which(disease_otp$association_score.overall > otp_cutoff),]
  
  sdb_res <- STRINGdb::STRINGdb$new(version="10", species=9606, score_threshold = ppi_cutoff)
  gene.network <- sdb_res$get_graph()
  igraph::V(gene.network)$name <- gsub("^9606\\.", "", igraph::V(gene.network)$name)
  
  # Get pathways information from msigdb (https://www.gsea-msigdb.org/)
  db_annots = msigdbr::msigdbr(species = "Homo sapiens") %>% 
    dplyr::filter(gs_subcat == "BP" | gs_subcat == "MF" | 
                    gs_subcat == "CP:KEGG" | gs_subcat == "CP:REACTOME")
  
  
  # Harmonize gene labels
  #mart <- biomaRt::useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  mart <- biomaRt::useEnsembl("ensembl", "hsapiens_gene_ensembl", "useast.ensembl.org")
  # It is actually faster to load whole table rather than post a large filter
  mart_results <- biomaRt::getBM(attributes = c("ensembl_gene_id", "ensembl_peptide_id", "entrezgene_id"), 
                                 mart = mart)
  #ppi_mart_results <- biomaRt::getBM(attributes = c("ensembl_gene_id", "ensembl_peptide_id"), #, "hgnc_symbol", "entrezgene_id"
                                 #filters = "ensembl_peptide_id", 
                                 #values = igraph::V(gene.network)$name,
                                 #mart = mart)
  #dis_pw_mart_results <- biomaRt::getBM(attributes = c("ensembl_gene_id", "entrezgene_id"),
                                     #filters = "entrezgene_id", 
                                     #values = unique(c(gene.diseases$entrezId, db_annots$entrez_gene)),
                                     #mart = mart)
  
  ppi_temp <- mart_results$ensembl_gene_id[match(igraph::V(gene.network)$name, 
                                                 mart_results$ensembl_peptide_id)]
  igraph::V(gene.network)$name[!is.na(ppi_temp)] <- ppi_temp[!is.na(ppi_temp)]
  
  disease.genes <- mart_results$ensembl_gene_id[match(gene.diseases$entrezId, 
                                                      mart_results$entrezgene_id)]
  disease.genes <- disease.genes[!is.na(disease.genes)]
  
  db_annots$ensembl_gene_id <- mart_results$ensembl_gene_id[match(db_annots$entrez_gene, 
                                                                  mart_results$entrezgene_id)]
  db_annots <- db_annots[!is.na(db_annots$ensembl_gene_id),]
  
  if (FALSE) { # If using hgnc symbols
    igraph::V(gene.network)$name <- mart_results$hgnc_symbol[match(igraph::V(gene.network)$name, 
                                                             mart_results$ensembl_peptide_id)]
    rownames(dat) <- mart_results$hgnc_symbol[match(rownames(dat), mart_results$ensembl_gene_id)]
  }
  if (FALSE) { # If using ensembl peptide ids
    rownames(dat) <- mart_results$ensembl_peptide_id[match(rownames(dat), mart_results$ensembl_gene_id)]
    disease.genes <- mart_results$ensembl_peptide_id[match(rownames(gene.diseases), mart_results$hgnc_symbol)]
    disease.genes <- disease.genes[!is.na(disease.genes) & !grepl("^\\s*$", disease.genes)]
  }
  
  #dat <- dat[!is.na(rownames(dat)) & !grepl("^\\s*$", rownames(dat)),]
  
  list_db_annots <- lapply(split(db_annots, db_annots$gs_name), function(x) x$ensembl_gene_id)
  list_db_annots <- list_db_annots[which(sapply(list_db_annots, length) < pw_max.size)]
  
  out <- fromGeneToNetworksToPathwayFeatures(dat, disease.genes, gene.network, list_db_annots, ...)
  return(out)
}

#'@param rwr_restart the restart probability used for RWR. See \code{dnet::dRWR} for more details.
#'@param rwr_norm the way to normalise the adjacency matrix of the input graph. See \code{dnet::dRWR} for more details.
#'@param rwr_cutoff the cuoff value to select the most visited genes.  (TO BE COMPLETED)
fromGeneToNetworksToPathwayFeatures <- function(dat, 
                                                disease.genes, # a character vector
                                                gene.network,  # igraph.object
                                                list_db_annots, # list of pathway annotated gene sets
                                                #study_batch = NULL, 
                                                top.ranked.genes = 100,
                                                min.size = 5, 
                                                max.size = 200, 
                                                parallel = 4,
                                                verbose = FALSE,
                                                rwr_restart = 0.33,
                                                rwr_norm = "quantile",
                                                rwr_cutoff = 0,
                                                ...) {
  # rank disease-genes within each sample
  gene.seeds <- plyr::aaply(dat, 2, function(s) {
                      sorted.genes <- names(s)[order(abs(s), decreasing = TRUE)[1:top.ranked.genes]]
                      out <- intersect(sorted.genes, disease.genes)
                      out <- as.numeric(igraph::V(gene.network)$name %in% out)
                      out
  })
  colnames(gene.seeds) <- igraph::V(gene.network)$name
  #gene.seeds <- unique(do.call(c, gene.seeds))
  #gene.seeds <- as.numeric(igraph::V(gene.network)$name %in% gene.seeds)
  #names(gene.seeds) <- igraph::V(gene.network)$name
  # compile the subnetwork from the PPI
  
  # apply random walk (dnet) to extend the set of genes
  rwr.top.genes <- dnet::dRWR(gene.network, setSeeds = t(gene.seeds), normalise = "none",
                              restart = rwr_restart, normalise.affinity.matrix = rwr_norm, 
                              parallel = parallel > 1, multicores = parallel, verbose = FALSE)
  rownames(rwr.top.genes) = igraph::V(gene.network)$name
  
  # apply fgsea
  res <- array(NA, c(nrow(gene.seeds), length(list_db_annots)), 
               list(id = rownames(gene.seeds), pathway = names(list_db_annots)))
  for (i in 1:ncol(rwr.top.genes)) {
    genes_i <- rwr.top.genes[,i]
    # Sum up duplicated gene id:s (alternative splicings)
    genes_i <- tapply(genes_i, names(genes_i), sum)
    genes_i <- genes_i[genes_i > rwr_cutoff]
    res_i <- fgsea::fgsea(list_db_annots, genes_i, nperm=10000, minSize=1, maxSize=200)
    res[i,match(res_i$pathway, names(list_db_annots))] <- res_i$NES * (-log10(res_i$padj)) 
  }
  
  return(res)
}

# RWR-FGSEA pipeline_DM_TG for drug matrix and tg-gates networks  (TO BE COMPLETED)
fromRWRtoFGSEA <- function() {
  #random walk with restart
  probm<-dRWR(gr, setSeeds=mat, normalise.affinity.matrix='quantile',parallel=T,verbose = F)
  probm<-as.matrix(probm)
  colnames(probm)<-names(d_sign)                        #Columns are chemicals
  rownames(probm)<-V(gr)$name
  for (i in 1:dim(probm)[2]){
    cutof_q<-sort(probm[,i],decreasing = T)[ng]
    probm[probm[,i]<cutof_q,i]<-0
  }
  #Fgsea
  q<-apply(probm,2,function(y)sum(y!=0)) 
  probm<-probm[,q>=ng]                    
  probm<-as.matrix(probm)
  colnames(probm)<-names(q)[which(q>=ng)]
  newscore<-ES<-NES<-pval<-matrix(NA,nrow=length(ptw),ncol=ncol(probm),
                                  dimnames=list(names(ptw),colnames(probm))) 
  for(i in 1:ncol(probm)){
    #for each chemical and its sets of genes we do this procedure
    res<-fgsea(pathways=ptw,stats=probm[probm[,i]!=0,i],nperm=10000,minSize=1,maxSize=200,BPPARAM=p)
    #the pathways more than 500 genes are not considered and those with at least one gene will be considered
    ES[res$pathway,i]<-res$ES
    NES[res$pathway,i]<-res$NES
    pval[res$pathway,i]<-res$pval
    newscore[res$pathway,i]<-(res$NES)*(-log10(res$pval))
  }
  return(list(ES=ES,NES=NES,pval=pval,newscore=newscore))
}# end function
