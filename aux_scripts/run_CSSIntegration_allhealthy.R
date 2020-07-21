# RUN CSS INTEGRATION ON ALL DATASETS
# RUNTIME: ~2h



# load libraries
library(Seurat)
library(simspec)

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

# scale data with these genes
## using all genes to avoid bias with variable genes across samples
integr_feat = rownames(all_healthy_srat@assays$SCT@scale.data)
VariableFeatures(all_healthy_srat) = integr_feat
all_healthy_srat = ScaleData(all_healthy_srat, features = integr_feat, use.umi = T,
                             do.scale = F, verbose = FALSE)
all_healthy_srat = RunPCA(all_healthy_srat, verbose = FALSE, assay = "SCT", npcs = 50,
                          features = integr_feat)

# Run CSS
all_healthy_srat = cluster_sim_spectrum(all_healthy_srat, label_tag="unique_name",
                                        redo_pca = T, dims_use = 1:30,
                                        var_genes = integr_feat)

# Run UMAP on CSS dimentions
all_healthy_srat = RunUMAP(all_healthy_srat, reduction="css", 
                           dims = 1:ncol(Embeddings(all_healthy_srat,"css")), 
                           reduction.name="umap_css", reduction.key="UMACSS_")

# save output
saveRDS(all_healthy_srat, file = "data/processed/scripts_out/CSSIntegr_allhealthy.RDS")