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
int_healthy = SelectIntegrationFeatures(healthy_srat, nfeatures = 5000)
int_cond = SelectIntegrationFeatures(cond_srat, nfeatures = 5000)
int_all = SelectIntegrationFeatures(all_cell_srat, nfeatures = 5000)
## include also top markers from identified healthy clusters
mk_healthy = read.csv("results/integr_healthy/clusters_res0.7_markers.csv", 
                      header = T, row.names = 1, stringsAsFactors = F)
mk_healthy = mk_healthy[mk_healthy$p_val_adj<=0.05 & mk_healthy$avg_logFC>0.2,]
mk_healthy = mk_healthy[order(mk_healthy$avg_logFC, decreasing = T),]
topmk = tapply(mk_healthy$gene, mk_healthy$cluster, function(x) x[1:100])
topmk = Reduce(union, topmk)
## integration features have healthy markers and var genes from healthy, cond and all
integr_feat = unique(c(topmk, int_healthy, int_cond, int_all))
integr_feat = integr_feat[!is.na(integr_feat)]
print(length(integr_feat))

# using ALL POSSIBLE genes - this is because many important genes are not shared between datasets
all_scale = lapply(all_cell_srat, function(x) rownames(x@assays$SCT@scale.data))
integr_feat = Reduce(intersect, all_scale)
print(length(integr_feat))

# calculate possible missing Pearson residuals for SCTransform
all_cell_srat = PrepSCTIntegration(all_cell_srat, anchor.features = integr_feat, verbose = T)

all_cell_srat <- lapply(X = all_cell_srat, FUN = function(x) {
  x <- ScaleData(x, features = integr_feat, verbose = FALSE)
  x <- RunPCA(x, features = integr_feat, verbose = FALSE)
})

# finding the anchors for integration
all_cell_anchors = FindIntegrationAnchors(all_cell_srat, normalization.method = "SCT", 
                                          anchor.features = integr_feat, reduction = "rpca",
                                          assay = rep("SCT", length(all_cell_srat)), 
                                          max.features = 500, scale = F, 
                                          dims = 1:50, verbose = T)

# actual data integration
allgenes = rownames(healthy_srat$sc_H1.EKS.Exp1@assays$RNA@counts)
all_cell_integr = IntegrateData(all_cell_anchors, normalization.method = "SCT", dims = 1:50,
                                verbose = T, features.to.integrate = allgenes)

# run PCA and UMAP to see how it looks
all_cell_integr = RunPCA(all_cell_integr, verbose = F)
all_cell_integr = RunUMAP(all_cell_integr, dims = 1:30)

# save output
#saveRDS(all_cell_anchors, file = "data/processed/scripts_out/SeuratIntegrRPCA_allcells_anchors.RDS")
saveRDS(all_cell_integr, file = "data/processed/scripts_out/SeuratIntegrRPCA_allcells_v2.RDS")