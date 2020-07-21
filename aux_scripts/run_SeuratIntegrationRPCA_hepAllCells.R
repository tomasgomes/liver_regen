# RUN SEURAT INTEGRATION HEPATOCYTES OF ALL DATASETS WITH rPCA
# RUNTIME: FAILS WITHOUT AN ERROR MESSAGE



# load libraries
library(Seurat)

# define working directory
setwd("..") # project head directory

# Load datasets
load(file = "data/processed/hep_donor_list.RData")

# find features variable across all individual datasets
integr_feat = SelectIntegrationFeatures(hep_donor_list, nfeatures = 5000)

# calculate possible missing Pearson residuals for SCTransform
hep_donor_list = PrepSCTIntegration(hep_donor_list, anchor.features = integr_feat, verbose = T)

hep_donor_list <- lapply(X = hep_donor_list, FUN = function(x) {
  x <- ScaleData(x, features = integr_feat, verbose = FALSE)
  x <- RunPCA(x, features = integr_feat, verbose = FALSE)
})

# finding the anchors for integration
all_cell_anchors = FindIntegrationAnchors(hep_donor_list, normalization.method = "SCT", 
                                          anchor.features = integr_feat, reduction = "rpca",
                                          assay = rep("SCT", length(hep_donor_list)), 
                                          dims = 1:50, verbose = T)

# actual data integration
allgenes = rownames(hep_donor_list[[1]]@assays$RNA@counts)
all_cell_integr = IntegrateData(all_cell_anchors, normalization.method = "SCT", dims = 1:20,
                                verbose = T, features.to.integrate = allgenes)

# run PCA and UMAP to see how it looks
all_cell_integr = RunPCA(all_cell_integr, verbose = F)
all_cell_integr = RunUMAP(all_cell_integr, dims = 1:30)

# save output
saveRDS(all_cell_integr, file = "data/processed/scripts_out/SeuratIntegrRPCA_hepAllCells.RDS")