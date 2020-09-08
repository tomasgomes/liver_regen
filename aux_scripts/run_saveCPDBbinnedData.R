# SAVE DATA FOR CELLPHONEDB
## duration: 6hrs


# setup
setwd("../")
library(Seurat)
library(data.table)



# general vars
nbins = 20
cpdb_path = "results/cell_comm/CellPhoneDB_runs"
mincells = 10



# load data
allcells_css = readRDS(file = "data/processed/allcells_css.RDS")
end_cells = readRDS(file = "results/zonation_cond/end_cells_zonation.RDS")
hep_cells = readRDS(file = "results/zonation_cond/hep_cells_zonation.RDS")



# define pseudotime bins 
## hepatocytes
df_h = data.table("zonation_pt" = hep_cells$healthy@meta.data$zonation_pt)
df_h$zonation_int = cut(df_h$zonation_pt, nbins)

df_e = data.table("zonation_pt" = hep_cells$embolised@meta.data$zonation_pt)
df_e = df_h[df_e, on="zonation_pt", roll=Inf, rollends = T]

df_r = data.table("zonation_pt" = hep_cells$regenerating@meta.data$zonation_pt)
df_r = df_h[df_r, on="zonation_pt", roll=Inf, rollends = T]

hep_int_l = list("healthy" = df_h, "embolised" = df_e, "regenerating" = df_r)

## endothelial cells
df_h = data.table("zonation_pt" = end_cells$healthy@meta.data$zonation_pt)
df_h$zonation_int = cut(df_h$zonation_pt, nbins)

df_e = data.table("zonation_pt" = end_cells$embolised@meta.data$zonation_pt)
df_e = df_h[df_e, on="zonation_pt", roll=Inf, rollends = T]

df_r = data.table("zonation_pt" = end_cells$regenerating@meta.data$zonation_pt)
df_r = df_h[df_r, on="zonation_pt", roll=Inf, rollends = T]

end_int_l = list("healthy" = df_h, "embolised" = df_e, "regenerating" = df_r)



# save data
cond_cells = list()
for(cond in unique(allcells_css@meta.data$Condition)){
  print(cond)
  cond_cells[[cond]] = allcells_css[,allcells_css@meta.data$Condition==cond]
  
  # remove doublet cells
  cond_cells[[cond]] = cond_cells[[cond]][,cond_cells[[cond]]@meta.data$allcells_major!="Doublets"]
  
  # remove contaminating hepatocytes
  cond_cells[[cond]] = cond_cells[[cond]][,colnames(cond_cells[[cond]]) %in% colnames(hep_cells[[cond]]) | cond_cells[[cond]]@meta.data$allcells_major!="Hepatocytes"]
  
  # remove contaminating endothelial cells
  cond_cells[[cond]] = cond_cells[[cond]][,colnames(cond_cells[[cond]]) %in% colnames(end_cells[[cond]]) | !startsWith(cond_cells[[cond]]@meta.data$allcells_simp,"LSEC")]
  
  # add zonation labels
  meta = data.frame(row.names = rownames(cond_cells[[cond]]@meta.data),
                    "Cell" = rownames(cond_cells[[cond]]@meta.data), 
                    "cell_type" = as.character(cond_cells[[cond]]@meta.data$allcells_simp),
                    stringsAsFactors = F)
  
  hep_zones = hep_int_l[[cond]]$zonation_int
  bnames = paste0("Hepatocytes_bin", 1:nbins)
  names(bnames) = levels(hep_zones)
  levels(hep_zones) = bnames
  hep_zones = as.character(hep_zones)
  names(hep_zones) = colnames(hep_cells[[cond]])
  meta[names(hep_zones),"cell_type"] = hep_zones
  
  end_zones = end_int_l[[cond]]$zonation_int
  bnames = paste0("LSEC_bin", 1:nbins)
  names(bnames) = levels(end_zones)
  levels(end_zones) = bnames
  end_zones = as.character(end_zones)
  names(end_zones) = colnames(end_cells[[cond]])
  meta[names(end_zones),"cell_type"] = end_zones
  
  # filter meta
  tabct = table(meta$cell_type)
  remct = names(tabct[tabct<mincells])
  meta = meta[!(meta$cell_type %in% remct),]
  
  write.table(meta, file = paste0(cpdb_path, "/", cond, "_meta_names_bins.txt"), 
              sep = "\t", col.names = T, row.names = F, quote = F)
  
  # FILTER, normalise and save
  cond_cells[[cond]] = cond_cells[[cond]][,rownames(meta)]
  cond_cells[[cond]] = suppressWarnings(SCTransform(cond_cells[[cond]], do.correct.umi = T, 
                                                    verbose = F, 
                                                    vars.to.regress=c("unique_name",
                                                                      "nCount_RNA"),
                                                    variable.features.rv.th = 1, seed.use = 1,
                                                    return.only.var.genes = F, 
                                                    variable.features.n = NULL))
  
  dat = cbind(rownames(cond_cells[[cond]]@assays$SCT@data),
              Matrix::as.matrix(cond_cells[[cond]]@assays$SCT@data))
  colnames(dat)[1] = "Gene"
  write.table(dat, file = paste0(cpdb_path, "/", cond, "_exp_norm_bins.txt"), 
              sep = "\t", col.names = T, row.names = F, quote = F)
}


