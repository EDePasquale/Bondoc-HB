# Load libraries
library(Seurat)
library(plyr)
library(dplyr)

# Read in background object
setwd("/data/GI-Informatics/DePasquale/Projects/Bondoc_HB/Resolution_0.5background_no17/Harmony_Background/")
M_back <- readRDS("Seurat_Log_res0.5_backsub.rds")

# Add cluster names
cluster_names<-c("Hepatocyte 1", # 0
                 "Hepatocyte 2", # 1
                 "Hepatocyte 3", # 2
                 "Hepatocyte 4", # 3
                 "Kupffer cells", # 4
                 "Hepatocyte 5", # 5
                 "LSEC 1", # 6
                 "Hepatocyte 6", # 7
                 "Hepatocyte 7", # 8
                 "LSEC 2", # 9
                 "Cycling", # 10
                 "Cholangiocytes", # 11
                 "Stellate cells", # 12
                 "Hepatocyte 8", #13
                 "Immune") # 14
M_back@meta.data[["cluster_names"]] <- mapvalues(M_back@meta.data[["subcluster_imm"]], from=c(0:14), to=cluster_names)
Idents(M_back) = M_back$cluster_names

pdf(file = "TSNE_Names.pdf", width = 10, height = 9)
par(mar=c(2, 2, 2, 2))
  DimPlot(M_back, reduction="tsne", group.by = "cluster_names", label=T, raster=F, repel=T)
dev.off()

pdf(file = "UMAP_Names.pdf", width = 10, height = 9)
par(mar=c(2, 2, 2, 2))
  DimPlot(M_back, reduction="umap", group.by = "cluster_names", label=T, raster=F, repel=T)
dev.off()


# Add cluster names redu
cluster_names<-c("Hepatocyte", # 0
                 "Hepatocyte", # 1
                 "Hepatocyte", # 2
                 "Hepatocyte", # 3
                 "Kupffer cells", # 4
                 "Hepatocyte", # 5
                 "LSEC", # 6
                 "Hepatocyte", # 7
                 "Hepatocyte", # 8
                 "LSEC", # 9
                 "Cycling", # 10
                 "Cholangiocytes", # 11
                 "Stellate cells", # 12
                 "Hepatocyte", #13
                 "Immune") # 14
M_back@meta.data[["cluster_names_redu"]] <- mapvalues(M_back@meta.data[["subcluster_imm"]], from=c(0:14), to=cluster_names)
Idents(M_back) = M_back$cluster_names_redu

pdf(file = "TSNE_Names_Redu.pdf", width = 10, height = 9)
par(mar=c(2, 2, 2, 2))
  DimPlot(M_back, reduction="tsne", group.by = "cluster_names_redu", label=T, raster=F, repel=T)
dev.off()

pdf(file = "UMAP_Names_Redu.pdf", width = 10, height = 9)
par(mar=c(2, 2, 2, 2))
  DimPlot(M_back, reduction="umap", group.by = "cluster_names_redu", label=T, raster=F, repel=T)
dev.off()

saveRDS(M_back, file = "Seurat_Log_res0.5_backsub_names.rds")


##############################################
#                                            #
# Where do these cells lie on the full UMAP? #
#                                            #
##############################################

# Read in full Seurat object (no Lung background)
setwd("/data/GI-Informatics/DePasquale/Projects/Bondoc_HB/Resolution_0.5noLung")
M_full <- readRDS("Seurat_Log_res0.5_meta.rds")

# Transfer all meta data with left_join
M_full_meta=cbind(cell=row.names(M_full@meta.data), M_full@meta.data)
M_back_meta=cbind(cell=row.names(M_back@meta.data), M_back@meta.data)
#M_full_meta[,1]=gsub("_[^_]+$", "", M_full_meta[,1])
#M_back_meta[,1]=gsub("_[^_]+$", "", M_back_meta[,1])
new_meta=left_join(M_full_meta, M_back_meta, by=c("cell"="cell"))

# remove duplicate name mismatches because of trimming
#new_meta=new_meta[which(new_meta$orig.ident.x==new_meta$orig.ident.y),]

print(colnames(M_full@meta.data))
print(colnames(new_meta))

# Add new meta data to new Seurat object
M_full@meta.data<-cbind(M_full@meta.data, new_meta[,22:23])

# Save object
saveRDS(M_full, file = "Seurat_Liver_25_backsub_backGroups.rds")

# Plot Seurat TSNE
pdf(file = "TSNE_Liver_25_Seurat_backsub_backGroups.pdf", width = 10, height = 8)
par(mar=c(2, 2, 2, 2))
  DimPlot(M_full, reduction="tsne", group.by = "cluster_names", label=T, raster=F, shuffle=T)
dev.off()

# Plot Seurat UMAP
pdf(file = "UMAP_Liver_25_Seurat_backsub_backGroups.pdf", width = 10, height = 8)
par(mar=c(2, 2, 2, 2))
  DimPlot(M_full, reduction="umap", group.by = "cluster_names", label=T, raster=F, shuffle=T)
dev.off()


