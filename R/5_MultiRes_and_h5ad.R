##############################
#                            #
# Prepare h5ad Seurat object #
#       Bondoc - HB          #
#      June 15, 2023         #
#                            #
##############################

# Load libraries
library(Seurat)
#BiocManager::install("zellkonverter")
library(zellkonverter) # package for easy Seurat to anndata conversion

# Pull in object with background groups added
setwd("/Volumes/GI-Informatics/DePasquale/Projects/Bondoc_HB/Resolution_0.5noLung/")
M<-readRDS("Seurat_Liver_25_backsub_backGroups.rds")

# Cluster at various resolutions and save as metadata
DefaultAssay(M)<-"integrated"
for(i in 1:10){
  M <- FindClusters(M, resolution = i/10)
  M@meta.data[[paste0("Resolution_", i/10)]]<-M@active.ident
  assign(paste0("p", i), DimPlot(M, group.by = paste0("Resolution_", i/10 ), raster=F))
}

# Remove unnecessary and redundant metadata slots
M@meta.data<-M@meta.data[,c(1,2,3,6,9,10,11,12,14,17,19,21,23,24,26,28,30,32,33)]

# Make UMAP plot with all resolutions
pdf(file = "All_Resolutions_UMAP.pdf", width = 48, height = 32)
par(mar=c(4, 4, 4, 4))
  p1+p2+p3+p4+p5+p6+p7+p8+p9+p10
dev.off()

# Save new Seurat object
saveRDS(M, file = "Seurat_Liver_25_backsub_backGroups_multiRes.rds")

# Convert to h5ad and save
M_SCE=as.SingleCellExperiment(M, assay ="RNA")
writeH5AD(M_SCE, file='Bondoc_HB_MultiRes_RNA.h5ad')
