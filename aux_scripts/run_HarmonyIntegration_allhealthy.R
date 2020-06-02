# RUN SEURAT INTEGRATION ON HEALTHY DATASETS
# RUNTIME: ~30min



# load libraries
library(Seurat)
library(harmony)

# define working directory
setwd("..") # project head directory

# Loading/creating the merged dataset
merged_file = "data/processed/healthy_merged_srat.RDS"
if(!file.exists(merged_file)){
  # Load datasets
  load(file = "data/processed/healthy_srat_indivclust.RData")

  # Merge datasets (here by actually creating a single Seurat object)
  ## this is a long step
  all_healthy_srat = merge(healthy_srat[[1]], healthy_srat[2:length(healthy_srat)], 
                           add.cell.ids = names(healthy_srat))
  
  saveRDS(all_healthy_srat, file = merged_file)
} else{
  # Load merged datasets
  all_healthy_srat = readRDS(merged_file)
}

# Re-run the Seurat normalisation for all cells (I don't think this is required)
#all_cell_srat = SCTransform(all_cell_srat, do.correct.umi = T, verbose = F, seed.use = 1,
#                            vars.to.regress = "nCount_RNA", variable.features.rv.th = 1,
#                            return.only.var.genes = F, variable.features.n = NULL)


# Run PCA on merged object
all_healthy_srat = RunPCA(all_healthy_srat, verbose = FALSE, assay = "SCT",
                          features = rownames(all_healthy_srat@assays$SCT@scale.data))

# Run Harmony
all_healthy_srat = RunHarmony(all_healthy_srat, "orig.ident", tau = 30, 
                              plot_convergence = F, assay.use = "SCT")

# Run UMAP on Harmony dimentions
all_healthy_srat = RunUMAP(all_healthy_srat, reduction = "harmony", dims = 1:10)

# save output
saveRDS(all_healthy_srat, file = "data/processed/scripts_out/HarmIntegr_allhealthy.RDS")


