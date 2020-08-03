# RUN CSS INTEGRATION ON ALL DATASETS
# RUNTIME: ~2h



# load libraries
library(Seurat)
library(simspec)

# define working directory
setwd("..") # project head directory

# Load datasets
load(file = "data/processed/healthy_srat_indivclust.RData")
load(file = "data/processed/cond_srat_indivclust.RData")

# Loading/creating the merged dataset
merged_file = "data/processed/healthy_cond_merged_srat.RDS"
if(!file.exists(merged_file)){
  # Merge datasets (here by actually creating a single Seurat object)
  ## this is a long step
  all_cell_l = c(healthy_srat, cond_srat)
  
  all_cell_srat = merge(all_cell_l[[1]], all_cell_l[2:length(all_cell_l)], 
                        add.cell.ids = names(all_cell_l))
  
  saveRDS(all_cell_srat, file = merged_file)
} else{
  # Load merged datasets
  all_cell_srat = readRDS(merged_file)
}

# find features variable across all individual datasets
int_healthy = SelectIntegrationFeatures(healthy_srat, nfeatures = 3000)
int_cond = SelectIntegrationFeatures(cond_srat, nfeatures = 3000)
## include also top markers from identified healthy clusters
#mk_healthy = read.csv("results/integr_healthy/clusters_res0.7_markers.csv", 
#                      header = T, row.names = 1, stringsAsFactors = F)
mk_healthy = read.csv("results/integr_healthy/css_allPops_markers.csv", 
                      header = T, stringsAsFactors = F)
mk_healthy = mk_healthy[mk_healthy$p_val_adj<=0.05 & mk_healthy$avg_logFC>0.2,]
mk_healthy = mk_healthy[order(mk_healthy$avg_logFC, decreasing = T),]
topmk = tapply(mk_healthy$gene, mk_healthy$cluster, function(x) x[1:100])
topmk = Reduce(union, topmk)
## integration features have healthy markers and var genes from healthy, cond and all
integr_feat = unique(c(topmk, int_healthy, int_cond))
integr_feat = integr_feat[!is.na(integr_feat)]
integr_feat = integr_feat[integr_feat %in% rownames(all_cell_srat@assays$SCT@data)]
print(length(integr_feat))

# scale data with these genes
VariableFeatures(all_cell_srat) = integr_feat
all_cell_srat = ScaleData(all_cell_srat, features = integr_feat, 
                          #use.umi = T, vars.to.regress = c("nCount_SCT", "Donor"),
                          do.scale = F, verbose = FALSE)
all_cell_srat = RunPCA(all_cell_srat, verbose = FALSE, assay = "SCT", npcs = 60,
                       features = rownames(all_cell_srat@assays$SCT@scale.data))

# Run CSS
all_cell_srat = cluster_sim_spectrum(all_cell_srat, label_tag="unique_name",
                                     redo_pca = T, dims_use = 1:50, 
                                     var_genes = integr_feat)

# Run UMAP on CSS dimentions
all_cell_srat = RunUMAP(all_cell_srat, reduction="css", dims = 1:ncol(Embeddings(all_cell_srat,"css")), 
                        reduction.name="umap_css", reduction.key="UMAPCSS_")

# save output
#saveRDS(all_cell_srat, file = "data/processed/scripts_out/CSSIntegr_allcells_corrCount_pc50.RDS")
saveRDS(all_cell_srat, file = "data/processed/scripts_out/CSSIntegr_allcells_pc50_cssmk.RDS")



