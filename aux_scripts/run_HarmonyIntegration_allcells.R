# RUN SEURAT INTEGRATION ON ALL DATASETS
# RUNTIME: ~2h



# load libraries
library(Seurat)
library(harmony)

# define working directory
setwd("..") # project head directory

# Loading/creating the merged dataset
merged_file = "data/processed/healthy_cond_merged_srat.RDS"
if(!file.exists(merged_file)){
  # Load datasets
  load(file = "data/processed/healthy_srat_indivclust.RData")
  load(file = "data/processed/cond_srat_indivclust.RData")
  
  # Merge datasets (here by actually creating a single Seurat object)
  ## this is a long step
  all_cell_l = c(healthy_srat, cond_srat)
  all_cell_srat = merge(all_cell_l[[1]], all_cell_l[[2]])
  for(n in names(all_cell_l)){
    all_cell_srat = merge(all_cell_srat, all_cell_l[[n]])
  }
  
  saveRDS(all_cell_srat, file = merged_file)
} else{
  # Load merged datasets
  all_cell_srat = readRDS(merged_file)
}

# Re-run the Seurat normalisation for all cells (I don't think this is required)
#all_cell_srat = SCTransform(all_cell_srat, do.correct.umi = T, verbose = F, seed.use = 1,
#                            vars.to.regress = "nCount_RNA", variable.features.rv.th = 1,
#                            return.only.var.genes = F, variable.features.n = NULL)


# Run PCA on merged object
all_cell_srat = RunPCA(all_cell_srat, features = rownames(all_cell_srat@assays$SCT@scale.data), 
                       verbose = FALSE, assay = "SCT")

# Run Harmony
all_cell_srat = RunHarmony(all_cell_srat, "orig.ident", tau = 30, plot_convergence = F, assay.use = "SCT")

# Run UMAP on Harmony dimentions
all_cell_srat = RunUMAP(all_cell_srat, reduction = "harmony", dims = 1:30)

# save output
saveRDS(all_cell_srat, file = "data/processed/scripts_out/HarmIntegr_allcells.RDS")


