library(Seurat)
allcells_css = readRDS(file = "data/processed/allcells_css_reannot.RDS")

cell_type_mk = list()
allcells_css = SetIdent(allcells_css, value = "major_ct")
na_filt = !is.na(allcells_css@meta.data$major_ct)
cell_type_mk[["major_all"]] = FindAllMarkers(allcells_css[,na_filt], 
                                         assay = "SCT", test.use = "wilcox",
                                         pseudocount.use = 0.1, logfc.threshold = 0,
                                         max.cells.per.ident = 10000,
                                         min.cells.feature = 10, verbose = T)

cell_type_mk[["major_healthy"]] = FindAllMarkers(allcells_css[,na_filt & allcells_css@meta.data$Condition=="healthy"], 
                                             assay = "SCT", test.use = "wilcox",
                                             pseudocount.use = 0.1, logfc.threshold = 0,
                                             return.thresh = 1,
                                             min.cells.feature = 3, verbose = T)

saveRDS(cell_type_mk, file = "results/integr_allcells/cell_type_mk_major.RDS")