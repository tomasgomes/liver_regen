# RUN SEURAT INTEGRATION ON ALL DATASETS WITH rPCA
# RUNTIME: ~45min



# load libraries
library(Seurat)

# define working directory
setwd("..") # project head directory

# Load datasets
load(file = "data/processed/healthy_srat_indivclust.RData")

# find features variable across all individual datasets
integr_feat = SelectIntegrationFeatures(healthy_srat, nfeatures = 5000)

# calculate possible missing Pearson residuals for SCTransform
healthy_srat = PrepSCTIntegration(healthy_srat, anchor.features = integr_feat, verbose = T)

healthy_srat <- lapply(X = healthy_srat, FUN = function(x) {
  x <- ScaleData(x, features = integr_feat, verbose = FALSE)
  x <- RunPCA(x, features = integr_feat, verbose = FALSE)
})


# finding the anchors for integration
all_cell_anchors = FindIntegrationAnchors(healthy_srat, normalization.method = "SCT", 
                                          anchor.features = integr_feat, reduction = "rpca",
                                          assay = rep("SCT", length(healthy_srat)), 
                                          dims = 1:50, verbose = T)

# actual data integration
allgenes = rownames(healthy_srat$sc_H1.EKS.Exp1@assays$RNA@counts)
all_cell_integr = IntegrateData(all_cell_anchors, normalization.method = "SCT", dims = 1:50,
                                verbose = T, features.to.integrate = allgenes)

# run PCA and UMAP to see how it looks
all_cell_integr = RunPCA(all_cell_integr, verbose = F)
all_cell_integr = RunUMAP(all_cell_integr, dims = 1:30)

# save output
saveRDS(all_cell_anchors, file = "data/processed/scripts_out/SeuratIntegrRPCA_allhealthy_anchors.RDS")
saveRDS(all_cell_integr, file = "data/processed/scripts_out/SeuratIntegrRPCA_allhealthy.RDS")
