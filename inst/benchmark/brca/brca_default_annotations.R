# Load the functional annotations used by GSE scripts
db_annots <- msigdbr::msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG")
db_annots <- rbind(db_annots, msigdbr::msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME"))
db_annots <- rbind(db_annots, msigdbr::msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP"))
db_annots <- rbind(db_annots, msigdbr::msigdbr(species = "Homo sapiens", category = "H"))

list_db_annots <- lapply(split(db_annots, db_annots$gs_name), function(x) x$gene_symbol)
list_db_annots <- list_db_annots[which(sapply(list_db_annots, length) <= 200 & sapply(list_db_annots, length) >= 5)]