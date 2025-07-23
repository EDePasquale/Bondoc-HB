#################################
#                               #
#        Sub-clustering         #
#          Bondoc - HB          #
#       Erica DePasquale        #
#          18 May 2023          #
#                               #
#################################

# Load libraries
library(Seurat)
library(gdata)
library(plyr)


############################################
#                 <><><><>                 #
# Load in all background for subclustering #
#        the immune/cycling cluster        #
#                 <><><><>                 #
############################################

setwd("/Volumes/GI-Informatics/DePasquale/Projects/Bondoc_HB/Resolution_0.5background_no17/")
M=readRDS("Seurat_Log_res0.5.rds")

# Subset to only hep
M_imm=subset(x = M, idents = c(10)) 

# Subcluster (or only rerun DR)
DefaultAssay(M_imm)<-"integrated"
#M_imm <- FindVariableFeatures(M_imm, selection.method = "vst", nfeatures = 2000)
#M_imm <- ScaleData(object = M_imm)
M_imm <- RunPCA(M_imm, npcs = 20, verbose = FALSE)
M_imm <- RunUMAP(M_imm, reduction = "pca", dims = 1:20)
M_imm <- RunTSNE(M_imm, reduction = "pca", dims = 1:20)
M_imm <- FindNeighbors(M_imm, reduction = "pca", dims = 1:20)
M_imm <- FindClusters(M_imm, resolution = 0.05)
# M_imm <- FindMultiModalNeighbors(M_imm, reduction.list = list("pca", "lsi"), dims.list = list(1:30, 2:50)) # playing with https://satijalab.org/seurat/articles/weighted_nearest_neighbor_analysis.html#wnn-analysis-of-10x-multiome-rna-atac-1
# M_imm <- RunUMAP(M_imm, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
# M_imm <- FindClusters(M_imm, graph.name = "wsnn", algorithm = 3, resolution=0.5)

dir.create("Harmony_Background")
setwd("Harmony_Background")

# Visualization
p1 <- DimPlot(M_imm, reduction = "umap", group.by = "orig.ident", raster=FALSE)
p2 <- DimPlot(M_imm, reduction = "umap", label = TRUE, repel = TRUE, raster=FALSE)
p3 <- DimPlot(M_imm, reduction = "tsne", group.by = "orig.ident", raster=FALSE)
p4 <- DimPlot(M_imm, reduction = "tsne", label = TRUE, repel = TRUE, raster=FALSE)

pdf(file = "UMAP_sample_and_cluster.pdf", width = 18, height = 8.5)
par(mar=c(4, 4, 4, 4))
p1 + p2
dev.off()

pdf(file = "UMAP_split_by_sample.pdf", width = 80, height = 8.5)
par(mar=c(4, 4, 4, 4))
DimPlot(M_imm, reduction = "umap", split.by = "orig.ident")
dev.off()

pdf(file = "TSNE_sample_and_cluster.pdf", width = 18, height = 8.5)
par(mar=c(4, 4, 4, 4))
p3 + p4
dev.off()

pdf(file = "TSNE_split_by_sample.pdf", width = 80, height = 8.5)
par(mar=c(4, 4, 4, 4))
DimPlot(M_imm, reduction = "tsne", split.by = "orig.ident")
dev.off()

#############

# Plot Seurat UMAP
pdf(file = "UMAP_Liver_25_Seurat_Groups.pdf", width = 7, height = 6)
par(mar=c(2, 2, 2, 2))
DimPlot(object = M_imm, reduction = "umap",label= TRUE, raster=FALSE)
dev.off()

# Plot Seurat TSNE
pdf(file = "TSNE_Liver_25_Seurat_Groups.pdf", width = 7, height = 6)
par(mar=c(2, 2, 2, 2))
DimPlot(object = M_imm, reduction = "tsne",label= TRUE, raster=FALSE)
dev.off()

#saveRDS(M_imm, file = "Seurat_Log_res0.5_backsub.rds")


#Add new split clusters to full Seurat object
temp=cbind(as.numeric(as.character(M@active.ident)),as.numeric(as.character(M@active.ident)))
row.names(temp)=M@assays[["RNA"]]@data@Dimnames[[2]]
new_clusters=mapvalues(as.numeric(as.character(M_imm@active.ident)), from=c(0,1), to=c(10,14))
temp[names(M_imm@active.ident),2]<-new_clusters
M@meta.data[["subcluster_imm"]] <- temp[,2]

# Plot New Seurat UMAP
pdf(file = "UMAP_Liver_25_Seurat_Groups_New_Clusters.pdf", width = 7, height = 6)
par(mar=c(2, 2, 2, 2))
  DimPlot(M, reduction="umap", group.by = "subcluster_imm", label=T, raster=FALSE)
dev.off()

# Plot New Seurat TSNE
pdf(file = "TSNE_Liver_25_Seurat_Groups_New_Clusters.pdf", width = 7, height = 6)
par(mar=c(2, 2, 2, 2))
  DimPlot(M, reduction="tsne", group.by = "subcluster_imm", label=T, raster=FALSE)
dev.off()

##############

# Save new object
saveRDS(M_imm, file = "Seurat_Log_res0.5_backsub_immONLY.rds")
saveRDS(M, file = "Seurat_Log_res0.5_backsub.rds")

# setwd("/Volumes/GI-Informatics/DePasquale/Projects/Huppert_CZI/WNNA_attempt1/Subcluster_9/Subcluster_hep_withCluster_Normal_0.1")
# M_imm=readRDS("Seurat_CZI_GEX_SCTv2_res0.25_ATAC_subcluster_hep_ONLY.rds")
# M=readRDS("Seurat_CZI_GEX_SCTv2_res0.25_ATAC_subcluster_hep.rds")

