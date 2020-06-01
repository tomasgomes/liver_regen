# RUN SEURAT INTEGRATION ON ALL DATASETS WITH rPCA
# RUNTIME: ~4h



# load libraries
library(Seurat)

# define working directory
setwd("..") # project head directory

# Load datasets
load(file = "data/processed/healthy_srat_indivclust.RData")
load(file = "data/processed/cond_srat_indivclust.RData")

# Merge datasets
all_cell_srat = c(healthy_srat, cond_srat)

# find features variable across all individual datasets
integr_feat = SelectIntegrationFeatures(all_cell_srat, nfeatures = 5000)

# calculate possible missing Pearson residuals for SCTransform
all_cell_srat = PrepSCTIntegration(all_cell_srat, anchor.features = integr_feat, verbose = T)

all_cell_srat <- lapply(X = all_cell_srat, FUN = function(x) {
  x <- ScaleData(x, features = integr_feat, verbose = FALSE)
  x <- RunPCA(x, features = integr_feat, verbose = FALSE)
})

# finding the anchors for integration - in this case 30 dims should be a good default
all_cell_anchors = FindIntegrationAnchors(all_cell_srat, normalization.method = "SCT", 
                                          anchor.features = integr_feat, reduction = "rpca",
                                          assay = rep("SCT", length(all_cell_srat)), 
                                          verbose = T)

# actual data integration
all_cell_integr = IntegrateData(all_cell_anchors, normalization.method = "SCT", verbose = T)

# run PCA and UMAP to see how it looks
all_cell_integr = RunPCA(all_cell_integr, verbose = F)
all_cell_integr = RunUMAP(all_cell_integr, dims = 1:30)

# save output
saveRDS(all_cell_integr, file = "data/processed/scripts_out/SeuratIntegrRPCA_allcells.RDS")
