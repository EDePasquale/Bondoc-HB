#################
#               #
#  Single Cell  #
#    Figure 2   #
#   8 Nov 2024  #
#               #
#################

# Load libraries
library(Seurat)
library(plyr)
library(ggplot2)
library(scales)
library(tidyr)
library(dplyr)
library(ggpubr)

# Load in object
M2 <- readRDS("/data/GI-Informatics/DePasquale/Projects/Bondoc_HB/Resolution_0.5noLung/Seurat_Liver_25_backsub_backGroups_multiRes_BTonly_subcluster.rds")

setwd("/data/GI-Informatics/DePasquale/Projects/Bondoc_HB/Resolution_0.5noLung/Prelim_Figs/")


###################
#                 #
# All Hepatocytes #
#                 #
###################

# Subset to Hepatocyte and plot with higher resolution
Idents(M2)<-"subcluster"
M_Hep<-subset(M2, idents=c("Cycling Hepatocyte", "Hepatocyte"))
DimPlot(M_Hep, raster=F, label=T, repel=T)

DefaultAssay(M_Hep)<-"integrated"
M_Hep <- RunPCA(M_Hep, npcs = 30, verbose = FALSE)
M_Hep <- RunUMAP(M_Hep, reduction = "pca", dims = 1:30)
M_Hep <- FindNeighbors(M_Hep, reduction = "pca", dims = 1:30)
M_Hep <- FindClusters(M_Hep, resolution = 0.5)
DimPlot(M_Hep, raster=F, label=T, repel=T)

DefaultAssay(M_Hep)<-"RNA"
# FeaturePlot(M_Hep, features="CYP2E1", order=T, min.cutoff = 0, raster=F) #
# FeaturePlot(M_Hep, features="GLUL", order=T, min.cutoff = 0, raster=F)
# FeaturePlot(M_Hep, features="GLS2", order=T, min.cutoff = 0, raster=F) #
# FeaturePlot(M_Hep, features="ARG1", order=T, min.cutoff = 0, raster=F)
FeaturePlot(M_Hep, features=c("HAL", "SDS", "ALDOB"), order=T, min.cutoff = 0, raster=F)
FeaturePlot(M_Hep, features=c("CYP1A2", "CYP3A4", "ADH1B"), order=T, min.cutoff = 0, raster=F)
# Not helpful with background/tumor differences


############
Idents(M_Hep)<-"SampleType"
M_Hep_B<-subset(M_Hep, idents="Background")
Idents(M_Hep_B)<-"seurat_clusters"

DefaultAssay(M_Hep_B)<-"integrated"
M_Hep_B <- FindVariableFeatures(M_Hep_B)
M_Hep_B <- RunPCA(M_Hep_B, npcs = 30, verbose = FALSE)
M_Hep_B <- RunUMAP(M_Hep_B, reduction = "pca", dims = 1:30)
DimPlot(M_Hep_B, raster=F, label=T, repel=T)

DefaultAssay(M_Hep_B)<-"RNA"
FeaturePlot(M_Hep_B, features=c("GLS2", "HAL", "SDS", "ALDOB"), order=T, min.cutoff = 0, raster=F)
FeaturePlot(M_Hep_B, features=c("GLUL", "CYP1A2", "CYP3A4", "ADH1B"), order=T, min.cutoff = 0, raster=F)
DimPlot(M_Hep_B, group.by = "orig.ident")
FeaturePlot(M_Hep_B, features=c("nCount_RNA", "nFeature_RNA"), order=T, min.cutoff = 0, raster=F) # no clear effects from this
FeaturePlot(M_Hep_B, features=c("percent.mt"), order=T, min.cutoff = 0, raster=F) # big differences here
FeaturePlot(M_Hep_B, features=c("ALB"), order=T, min.cutoff = 0, raster=F) # not much difference here

FeaturePlot(M_Hep_B, reduction="pca", features=c("percent.mt"), order=T, min.cutoff = 0, raster=F, dims=c(1,3)) # 1 associated with mito
FeaturePlot(M_Hep_B, reduction="pca", features=c("percent.mt"), order=T, min.cutoff = 0, raster=F, dims=c(2,3)) # 2 also associated with mito

M_Hep_B <- RunUMAP(M_Hep_B, reduction = "pca", dims = 3:30) # run without mito PCs
DimPlot(M_Hep_B, raster=F, label=T, repel=T)
FeaturePlot(M_Hep_B, features=c("GLS2", "HAL", "SDS", "ALDOB"), order=T, min.cutoff = 0, raster=F)
FeaturePlot(M_Hep_B, features=c("GLUL", "CYP1A2", "CYP3A4", "ADH1B"), order=T, min.cutoff = 0, raster=F)
FeaturePlot(M_Hep_B, features=c("GLS2", "CYP1A2"), order=T, min.cutoff = 0, raster=F) # these contradict. likely something other than zonation is controlling the shape of the UMAP still

####
# Remove trouble samples

Idents(M_Hep_B)<-"orig.ident"
A=unique(M_Hep_B@active.ident)
M_Hep_B_new=subset(M_Hep_B, ident=setdiff(A, c("HB17B", "HB31B")))

DefaultAssay(M_Hep_B_new)<-"integrated"
M_Hep_B_new <- FindVariableFeatures(M_Hep_B_new)
M_Hep_B_new <- RunPCA(M_Hep_B_new, npcs = 30, verbose = FALSE)
M_Hep_B_new <- RunUMAP(M_Hep_B_new, reduction = "pca", dims = 1:30)
DimPlot(M_Hep_B_new, raster=F, label=T, repel=T)

DefaultAssay(M_Hep_B_new)<-"RNA"
FeaturePlot(M_Hep_B_new, features=c("GLS2", "HAL", "SDS", "ALDOB"), order=T, min.cutoff = 0, raster=F)
FeaturePlot(M_Hep_B_new, features=c("GLUL", "CYP1A2", "CYP3A4", "ADH1B"), order=T, min.cutoff = 0, raster=F)
DimPlot(M_Hep_B_new, group.by = "orig.ident")
FeaturePlot(M_Hep_B_new, features=c("nCount_RNA", "nFeature_RNA"), order=T, min.cutoff = 0, raster=F) # no clear effects from this
FeaturePlot(M_Hep_B_new, features=c("percent.mt"), order=T, min.cutoff = 0, raster=F) # big differences here
FeaturePlot(M_Hep_B_new, features=c("ALB"), order=T, min.cutoff = 0, raster=F) # not much difference here

FeaturePlot(M_Hep_B_new, reduction="pca", features=c("percent.mt"), order=T, min.cutoff = 0, raster=F, dims=c(1,3)) # 1 associated with mito
FeaturePlot(M_Hep_B_new, reduction="pca", features=c("percent.mt"), order=T, min.cutoff = 0, raster=F, dims=c(2,3)) # 2 also associated with mito

M_Hep_B_new <- RunUMAP(M_Hep_B_new, reduction = "pca", dims = 3:30) # run without mito PCs
DimPlot(M_Hep_B_new, raster=F, label=T, repel=T)
FeaturePlot(M_Hep_B_new, features=c("GLS2", "HAL", "SDS", "ALDOB"), order=T, min.cutoff = 0, raster=F)
FeaturePlot(M_Hep_B_new, features=c("GLUL", "CYP1A2", "CYP3A4", "ADH1B"), order=T, min.cutoff = 0, raster=F)
FeaturePlot(M_Hep_B_new, features=c("GLS2", "CYP1A2"), order=T, min.cutoff = 0, raster=F) # these contradict. likely something other than zonation is controlling the shape of the UMAP still



#########
Idents(M_Hep)<-"SampleType"
M_Hep_T<-subset(M_Hep, idents="Tumor")
Idents(M_Hep_T)<-"seurat_clusters"

DefaultAssay(M_Hep_T)<-"integrated"
M_Hep_T <- FindVariableFeatures(M_Hep_T)
M_Hep_T <- RunPCA(M_Hep_T, npcs = 30, verbose = FALSE)
M_Hep_T <- RunUMAP(M_Hep_T, reduction = "pca", dims = 1:30)
DimPlot(M_Hep_T, raster=F, label=T, repel=T)

DefaultAssay(M_Hep_T)<-"RNA"
FeaturePlot(M_Hep_T, features=c("GLS2", "HAL", "SDS", "ALDOB"), order=T, min.cutoff = 0, raster=F)
FeaturePlot(M_Hep_T, features=c("GLUL", "CYP1A2", "CYP3A4", "ADH1B"), order=T, min.cutoff = 0, raster=F)
DimPlot(M_Hep_T, group.by = "orig.ident")
FeaturePlot(M_Hep_T, features=c("nCount_RNA", "nFeature_RNA"), order=T, min.cutoff = 0, raster=F) # no clear effects from this
FeaturePlot(M_Hep_T, features=c("percent.mt"), order=T, min.cutoff = 0, raster=F) # some differences here
FeaturePlot(M_Hep_T, features=c("ALB"), order=T, min.cutoff = 0, raster=F) # some differences here

# I'm struggling to get any clear zonation for either background or tumor, will want to talk with Katie about how she wants to proceed.


#######################
#                     #
# Cycling Hepatocytes #
#                     #
#######################

# Subset to Cycling and plot with higher resolution
Idents(M2)<-"Resolution_1"

myCols=hue_pal()(3)
pdf(file = "2_Cycling_for_highlight_3clus_UMAP.pdf", width = 8.5, height = 7)
par(mar=c(2, 2, 2, 2))
  DimPlot(M2, raster=F, label=T, repel=T, cols=c("grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", myCols[1], "grey", "grey", "grey", myCols[2], myCols[3], "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey"))
dev.off()

myCols=hue_pal()(11)
pdf(file = "2_Cycling_for_highlight_1clus_UMAP.pdf", width = 8.5, height = 7)
par(mar=c(2, 2, 2, 2))
  DimPlot(M2, group.by="subcluster", raster=F, label=T, repel=T, cols=c("grey", "grey", myCols[3], "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey"))
dev.off()

M_Cyc<-subset(M2, idents=c(8,12,13))
DimPlot(M_Cyc, raster=F, label=T, repel=T)

DefaultAssay(M_Cyc)<-"integrated"
M_Cyc <- RunPCA(M_Cyc, npcs = 30, verbose = FALSE)
M_Cyc <- RunUMAP(M_Cyc, reduction = "pca", dims = 1:30)

pdf(file = "2_Cycling_81213_UMAP.pdf", width = 7.5, height = 7)
par(mar=c(2, 2, 2, 2))
  DimPlot(M_Cyc, group.by="Resolution_1", raster=F, label=T, repel=T)
dev.off()

# Cell cycle scoring
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
DefaultAssay(M_Cyc)<-"RNA"
M_Cyc <- CellCycleScoring(M_Cyc, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

pdf(file = "2_CC_Scoring_altView_UMAP.pdf", width = 7.5, height = 7)
par(mar=c(2, 2, 2, 2))
  DimPlot(M_Cyc, raster=F)
dev.off()

pdf(file = "2_CC_Scoring_UMAP.pdf", width = 14.5, height = 7)
par(mar=c(2, 2, 2, 2))
  FeaturePlot(M_Cyc, features=c("S.Score", "G2M.Score"), order=T, min.cutoff = 0, raster=F)
dev.off()

Idents(M_Cyc)<-"Resolution_1"
results=FindAllMarkers(M_Cyc, only.pos = T)
write.table(results, "81213_Cluster_DE_Results.txt", sep="\t", quote=F)
# These are largely going to be cell cycle phase associated genes (as expected)

FeaturePlot(M_Cyc, features=c("EZH2"), order=T, min.cutoff = 0, raster=F)

pdf(file = "2_Dot_Plot_BandT_SampleType.pdf", width = 9, height = 5)
par(mar=c(2, 2, 2, 2))
  DotPlot(M_Cyc, assay="RNA", group.by="SampleType", features=c("EZH2", "SUZ12", "EED", "JARID2", "MKI67", "AURKB"), dot.scale=12, scale=F) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()

pdf(file = "2_Dot_Plot_BandT_Cluster.pdf", width = 9, height = 5)
par(mar=c(2, 2, 2, 2))
  DotPlot(M_Cyc, assay="RNA", group.by="Resolution_1", features=c("EZH2", "SUZ12", "EED", "JARID2", "MKI67", "AURKB"), dot.scale=12, scale=F)+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()

table(M_Cyc$SampleType, M_Cyc$Resolution_1)
#               8     12    13
# Background    168   137   162
# Tumor         4431  2687  2640

Idents(M_Cyc)<-"SampleType"
M_Cyc_B<-subset(M_Cyc, idents="Background")
M_Cyc_T<-subset(M_Cyc, idents="Tumor")

pdf(file = "2_Dot_Plot_B_Cluster.pdf", width = 9, height = 5)
par(mar=c(2, 2, 2, 2))
  DotPlot(M_Cyc_B, assay="RNA", group.by="Resolution_1", features=c("EZH2", "SUZ12", "EED", "JARID2", "MKI67", "AURKB"), dot.scale=12, scale=F)+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()

pdf(file = "2_Dot_Plot_T_Cluster.pdf", width = 9, height = 5)
par(mar=c(2, 2, 2, 2))
  DotPlot(M_Cyc_T, assay="RNA", group.by="Resolution_1", features=c("EZH2", "SUZ12", "EED", "JARID2", "MKI67", "AURKB"), dot.scale=12, scale=F)+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()



##################
# 13 December 2024

# Load in object
# M2 <- readRDS("/data/GI-Informatics/DePasquale/Projects/Bondoc_HB/Resolution_0.5noLung/Seurat_Liver_25_backsub_backGroups_multiRes_BTonly_subcluster.rds")

setwd("/data/GI-Informatics/DePasquale/Projects/Bondoc_HB/Resolution_0.5noLung/Prelim_Figs/")

Idents(M2)<-"orig.ident"
A=unique(M2@active.ident)
M2=subset(M2, ident=setdiff(A, c("96UES")))

myCols=hue_pal()(11)
pdf(file = "2_AllHep_for_highlight_1clus_UMAP.pdf", width = 8.5, height = 7)
par(mar=c(2, 2, 2, 2))
  DimPlot(M2, group.by="subcluster", raster=F, cols=c("grey", myCols[1], myCols[1], myCols[1], "grey", "grey", "grey", "grey", "grey", "grey", "grey"))
dev.off()

pdf(file = "2_AllHep_for_highlight_3clus_UMAP.pdf", width = 8.5, height = 7)
par(mar=c(2, 2, 2, 2))
  DimPlot(M2, group.by="subcluster", raster=F, cols=c("grey", myCols[2], myCols[3], myCols[4], "grey", "grey", "grey", "grey", "grey", "grey", "grey"))
dev.off()

# Subset to Hepatocyte and plot with higher resolution
Idents(M2)<-"subcluster"
M_Hep<-subset(M2, idents=c("Cycling Hepatocyte", "Hepatocyte", "Cholangiocyte-like"))

DefaultAssay(M_Hep)<-"integrated"
M_Hep <- RunPCA(M_Hep, npcs = 30, verbose = FALSE)
M_Hep <- RunUMAP(M_Hep, reduction = "pca", dims = 1:30)

pdf(file = "2_AllHep_new_3clus_UMAP.pdf", width = 8.5, height = 7)
par(mar=c(2, 2, 2, 2))
  DimPlot(M_Hep, raster=F, label=T, repel=T, cols=c(myCols[4], myCols[2], myCols[3])) #Panel A
dev.off()

data_long=as.data.frame(cbind(Sample=M_Hep@meta.data[["SampleType"]], Cluster=as.character(M_Hep@meta.data[["subcluster"]])))
data=as.data.frame.matrix(table(data_long))
data=t(data) # Read and wrangle data
data <- sweep(data, 2, colSums(data), "/")*100
myColors=c(myCols[2], myCols[3], myCols[4])

pdf("2_OPT1_Sub_Frequencies.pdf", width = 6, height = 7)
par(mar = c(8,4,4,16), xpd = T)
  barplot(as.matrix(data[nrow(data):1,]), col = rev(myColors), xaxt = "n", ylab = "Population frequency (%)", border = NA)
  axis(side = 1, at = seq(1,ncol(data))*1.2-0.5, labels = colnames(data), las = 2)
  legend(x = ncol(data)*1.2+0.5, y = 100, legend = row.names(data), fill = myColors, bty = "n", border = NA)
dev.off()

pdf("2_OPT2_Sub_Split.pdf", width = 14, height = 7)
par(mar = c(8,4,4,16), xpd = T)
  DimPlot(M_Hep, raster=F, label=T, repel=T, split.by="SampleType", cols=c(myCols[4], myCols[2], myCols[3])) #Panel A
dev.off()

Idents(M_Hep)<-"SampleType"
M_Hep_T<-subset(M_Hep, idents="Tumor")
Idents(M_Hep_T)<-"subcluster"

DefaultAssay(M_Hep_T)<-"RNA"

pdf("2_Dotplot_from_Katie_25Jan2025.pdf", width = 9, height = 4)
par(mar = c(8,4,4,16), xpd = T)
  DotPlot(M_Hep_T, features=c("ALB","AFP", "STAT3", "KRT19", "JARID2", "EZH2", "SUZ12", "EED", "MKI67", "AURKB", "ASPM", "TOP2A", "HELLS"), dot.scale=12) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()
# DotPlot(M_Hep_T, features=c("AFP", "KRT19", "JARID2", "EZH2", "SUZ12", "EED", "MKI67", "AURKB"), scale=F, dot.scale=12)

pdf("2_Dotplot_from_Katie_for_grant_2May2025.pdf", width = 9, height = 4)
par(mar = c(8,4,4,16), xpd = T)
  DotPlot(M_Hep_T, features=c("ALB","AFP", "GPC3", "STAT3", "CDH1", "MLXIPL", "KRT19", "JARID2", "PTCH1", "VCAN", "EZH2", "SUZ12", "EED", "MKI67", "AURKB", "ASPM"), dot.scale=12) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()

pdf(file = "2_SubHep_for_highlight_1clus_UMAP.pdf", width = 8.5, height = 7)
par(mar=c(2, 2, 2, 2))
  DimPlot(M_Hep_T, group.by="subcluster", raster=F,  cols=c("grey", myCols[3], "grey"))
dev.off()

pdf(file = "2_SubHep_KatieGenes_UMAP_25Jan2025.pdf", width = 35, height = 21)
par(mar=c(2, 2, 2, 2))
  FeaturePlot(M_Hep_T, features=c("ALB","AFP", "STAT3", "KRT19", "JARID2", "EZH2", "SUZ12", "EED", "MKI67", "AURKB", "ASPM", "TOP2A", "HELLS"), ncol=5, order=T, min.cutoff = 0, raster=F)
dev.off()

# Subset to Cycling and plot with higher resolution
Idents(M_Hep_T)<-"Resolution_1"

M_Cyc<-subset(M_Hep_T, idents=c(8,12,13))
DimPlot(M_Cyc, raster=F, label=T, repel=T)

DefaultAssay(M_Cyc)<-"integrated"
M_Cyc <- RunPCA(M_Cyc, npcs = 30, verbose = FALSE)
DimHeatmap(M_Cyc, dims = 1:15, cells = 500, balanced = TRUE)
ElbowPlot(M_Cyc) # UMAP with only 1:5 makes no difference

M_Cyc <- RunUMAP(M_Cyc, reduction = "pca", dims = 1:30)

pdf(file = "2_Cycling_81213_UMAP.pdf", width = 7.5, height = 7)
par(mar=c(2, 2, 2, 2))
  DimPlot(M_Cyc, group.by="Resolution_1", raster=F, label=T, repel=T)
dev.off()

# Cell cycle scoring
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
DefaultAssay(M_Cyc)<-"RNA"
M_Cyc <- CellCycleScoring(M_Cyc, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

pdf(file = "2_CC_Scoring_altView_UMAP.pdf", width = 7.5, height = 7)
par(mar=c(2, 2, 2, 2))
  DimPlot(M_Cyc, raster=F) 
dev.off()

pdf(file = "2_CC_Scoring_UMAP.pdf", width = 14.5, height = 7)
par(mar=c(2, 2, 2, 2))
  FeaturePlot(M_Cyc, features=c("S.Score", "G2M.Score"), order=T, min.cutoff = 0, raster=F, keep.scale="all") 
dev.off()

pdf(file = "2_CC_Scoring_UMAP_small.pdf", width = 8, height = 3.5)
par(mar=c(2, 2, 2, 2))
  FeaturePlot(M_Cyc, features=c("S.Score", "G2M.Score"), order=T, min.cutoff = 0, raster=F, keep.scale="all") 
dev.off()

