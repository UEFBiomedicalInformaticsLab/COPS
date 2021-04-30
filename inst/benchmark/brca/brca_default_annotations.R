# Load the functional annotations used by GSE scripts
db_annots <- msigdbr::msigdbr(species = "Homo sapiens")
db_annots <- dplyr::filter(db_annots, grepl(paste(GENE_SETS, collapse = "|"), gs_subcat))
list_db_annots <- lapply(split(db_annots, db_annots$gs_name), function(x) x$gene_symbol)
list_db_annots <- list_db_annots[which(sapply(list_db_annots, length) <= 200 & sapply(list_db_annots, length) >= 5)]