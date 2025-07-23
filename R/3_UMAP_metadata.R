##############################
#                            #
#        Color UMAP by       #
#        new meta data       #
#         Bondoc — HB        #
#       Erica DePasquale     #
#         26 Apr 2023        #
#                            #
##############################


# Load libraries
library(Seurat)
library(plyr)
library(gdata)

# Read in Seurat object
#setwd("/Users/dep1kb/Documents/Dated_Work_Folders/2023_04_26/")
#setwd("/Volumes/GI-Informatics/DePasquale/Projects/Bondoc_HB/Resolution_0.5background/")
setwd("/Volumes/GI-Informatics/DePasquale/Projects/Bondoc_HB/Resolution_0.5noLung/")
M=readRDS("Seurat_Log_res0.5.rds")
#M=readRDS("Seurat_Log_res0.5_meta.rds")

# Read in metadata file
metadata=read.table("/Volumes/GI-Informatics/DePasquale/Projects/Bondoc_HB/Resolution_0.5noLung/sample_metadata.txt", sep="\t", header=T)
metadata=metadata[metadata$Sample %in% unique(M@meta.data[["orig.ident"]]),]

# Add metadata to the Seurat object
M@meta.data[["SampleType"]] <- mapvalues(M@meta.data[["orig.ident"]], from=metadata$Sample, to=metadata$SampleType)
M@meta.data[["Gender"]] <- mapvalues(M@meta.data[["orig.ident"]], from=metadata$Sample, to=metadata$Gender)
M@meta.data[["Tissue"]] <- mapvalues(M@meta.data[["orig.ident"]], from=metadata$Sample, to=metadata$Tissue)
M@meta.data[["Diagnosis"]] <- mapvalues(M@meta.data[["orig.ident"]], from=metadata$Sample, to=metadata$Diagnosis)


# Plot UMAPs
pdf(file = "UMAP_metadata_Sample.pdf", width = 10, height = 8)
par(mar=c(2, 2, 2, 2))
print(DimPlot(object = M, reduction = "umap", pt.size=0.1, group.by="orig.ident", raster=F, shuffle=T))
dev.off()

pdf(file = "UMAP_metadata_SampleType.pdf", width = 9, height = 8)
par(mar=c(2, 2, 2, 2))
print(DimPlot(object = M, reduction = "umap", pt.size=0.1, group.by="SampleType", raster=F, shuffle=T))
dev.off()

pdf(file = "UMAP_metadata_Gender.pdf", width = 9, height = 8)
par(mar=c(2, 2, 2, 2))
print(DimPlot(object = M, reduction = "umap", pt.size=0.1, group.by="Gender", raster=F, shuffle=T))
dev.off()

pdf(file = "UMAP_metadata_Tissue.pdf", width = 9, height = 8)
par(mar=c(2, 2, 2, 2))
print(DimPlot(object = M, reduction = "umap", pt.size=0.1, group.by="Tissue", raster=F, shuffle=T))
dev.off()

pdf(file = "UMAP_metadata_Diagnosis.pdf", width = 9, height = 8)
par(mar=c(2, 2, 2, 2))
print(DimPlot(object = M, reduction = "umap", pt.size=0.1, group.by="Diagnosis", raster=F, shuffle=T))
dev.off()

pdf(file = "UMAP_metadata_Tissue_Split.pdf", width = 25, height = 8)
par(mar=c(2, 2, 2, 2))
  DimPlot(object = M, reduction = "umap", pt.size=0.1, raster=F, shuffle=T, split.by="Tissue")
dev.off()

pdf(file = "UMAP_metadata_SampleType_Split.pdf", width = 20, height = 8)
par(mar=c(2, 2, 2, 2))
  DimPlot(object = M, reduction = "umap", pt.size=0.1, raster=F, shuffle=T, split.by="SampleType")
dev.off()


# Plot TSNEs
pdf(file = "TSNE_metadata_Sample.pdf", width = 10, height = 8)
par(mar=c(2, 2, 2, 2))
print(DimPlot(object = M, reduction = "tsne", pt.size=0.1, group.by="orig.ident", raster=F, shuffle=T))
dev.off()

pdf(file = "TSNE_metadata_SampleType.pdf", width = 9, height = 8)
par(mar=c(2, 2, 2, 2))
print(DimPlot(object = M, reduction = "tsne", pt.size=0.1, group.by="SampleType", raster=F, shuffle=T))
dev.off()

pdf(file = "TSNE_metadata_Gender.pdf", width = 9, height = 8)
par(mar=c(2, 2, 2, 2))
print(DimPlot(object = M, reduction = "tsne", pt.size=0.1, group.by="Gender", raster=F, shuffle=T))
dev.off()

pdf(file = "TSNE_metadata_Tissue.pdf", width = 9, height = 8)
par(mar=c(2, 2, 2, 2))
print(DimPlot(object = M, reduction = "tsne", pt.size=0.1, group.by="Tissue", raster=F, shuffle=T))
dev.off()

pdf(file = "TSNE_metadata_Diagnosis.pdf", width = 9, height = 8)
par(mar=c(2, 2, 2, 2))
print(DimPlot(object = M, reduction = "tsne", pt.size=0.1, group.by="Diagnosis", raster=F, shuffle=T))
dev.off()


saveRDS(M, file = "Seurat_Log_res0.5_meta.rds")
