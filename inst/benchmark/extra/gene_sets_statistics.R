# Load the functional annotations used by GSE scripts
db_annots <- msigdbr::msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG")
db_annots <- rbind(db_annots, msigdbr::msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME"))
db_annots <- rbind(db_annots, msigdbr::msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP"))
db_annots <- rbind(db_annots, msigdbr::msigdbr(species = "Homo sapiens", category = "H"))

list_db_annots <- lapply(split(db_annots, db_annots$gs_name), function(x) x$gene_symbol)


gs_stats <- data.frame(source = sapply(strsplit(names(list_db_annots), split = "_"), function(x) x[[1]]), 
                       N = sapply(list_db_annots, length))

ggplot(gs_stats, aes(source, N)) + geom_boxplot()

mean(gs_stats$N[gs_stats$source == "HALLMARK"])

filtered_gs_stats <- gs_stats[gs_stats$N <= 200 & gs_stats$N >= 5, ]

mean(gs_stats$N[gs_stats$source == "KEGG"])
mean(gs_stats$N[gs_stats$source == "GO"])
mean(gs_stats$N[gs_stats$source == "REACTOME"])