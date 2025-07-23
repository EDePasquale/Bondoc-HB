#=========================================================================================
#Set resolution
resolution=0.3

dir.create(paste0("/data/GI-Informatics/DePasquale/Projects/Bondoc_HB/Resolution_", resolution, "background_no17"))
setwd(paste0("/data/GI-Informatics/DePasquale/Projects/Bondoc_HB/Resolution_", resolution, "background_no17"))

# Load libraries
print("Loading Libraries..")
suppressMessages(library(reticulate))
suppressMessages(library(Seurat))
suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(Matrix))
suppressMessages(library(gplots))
suppressMessages(library(ggplot2))
suppressMessages(library(gdata))
print("DONE! Now starting Analysis..")
#=========================================================================================

getFiles <- function(fileName,SampleName,grp){
seuratseuratObject <- Read10X(data.dir = fileName, gene.column = 2, unique.features = TRUE)
sample       <- CreateSeuratObject(counts = seuratseuratObject, project = SampleName, min.cells = 5)
sample$stim  <- SampleName
sample$grp  <- grp
sample[["percent.mt"]] <- PercentageFeatureSet(sample, pattern = "^MT-")
sample 	   <- subset(sample, subset =  percent.mt <= 20 & nFeature_RNA >= 500 & nCount_RNA >= 1000 )
return(sample)
}
#=========================================================================================

#sample1<-getFiles("/data/GI-Informatics/DePasquale/Projects/Bondoc_HB/Cell_Ranger/10X-Bondoc-81T/81T/strainedCounts", "81T", "Tumor")
#sample2<-getFiles("/data/GI-Informatics/DePasquale/Projects/Bondoc_HB/Cell_Ranger/10X-Bondoc-96UES/96UES/strainedCounts", "96UES", "Tumor")
#sample3<-getFiles("/data/GI-Informatics/DePasquale/Projects/Bondoc_HB/Cell_Ranger/10X-Bondoc-HB17-Background-20191001-3v3hg/HB17B/strainedCounts", "HB17B", "Background") #MT
#sample4<-getFiles("/data/GI-Informatics/DePasquale/Projects/Bondoc_HB/Cell_Ranger/10X-Bondoc-HB17-PDX-F2-20191002-3v3hg/HB17P/strainedCounts", "HB17P", "PDX")
#sample5<-getFiles("/data/GI-Informatics/DePasquale/Projects/Bondoc_HB/Cell_Ranger/10X-Bondoc-HB17-Tumor-8-16-19-20190821-3v3hg/HB17T/strainedCounts", "HB17T", "Tumor")
#sample6<-getFiles("/data/GI-Informatics/DePasquale/Projects/Bondoc_HB/Cell_Ranger/10X-Bondoc-HB30B-20210105-3v3-1hg/HB30B/strainedCounts", "HB30B", "Background") # Too few cells
#sample7<-getFiles("/data/GI-Informatics/DePasquale/Projects/Bondoc_HB/Cell_Ranger/10X-Bondoc-HB30-PDX-20201123-3v3-1hg/HB30P/strainedCounts", "HB30P", "PDX")
#sample8<-getFiles("/data/GI-Informatics/DePasquale/Projects/Bondoc_HB/Cell_Ranger/10X-Bondoc-HB30T-20201109-3v3-1hg/HB30T/strainedCounts", "HB30T", "Tumor")
sample9<-getFiles("/data/GI-Informatics/DePasquale/Projects/Bondoc_HB/Cell_Ranger/10X-Bondoc-HB31B-Background-Liver-20210621-3v3-1hg/HB31B/strainedCounts", "HB31B", "Background")
#sample10<-getFiles("/data/GI-Informatics/DePasquale/Projects/Bondoc_HB/Cell_Ranger/10X-Bondoc-HB31T-20210622-3v3-1hg/HB31T/strainedCounts", "HB31T", "Tumor")

#sample11<-getFiles("/data/GI-Informatics/DePasquale/Projects/Bondoc_HB/Cell_Ranger/10X-Bondoc-HB35PDX-20210623-3v3-1hg/HB35P/strainedCounts", "HB35P", "PDX")
sample12<-getFiles("/data/GI-Informatics/DePasquale/Projects/Bondoc_HB/Cell_Ranger/10X-Bondoc-HB53-Background-Liver-10-6-2020-20201030-3v3-1hg/HB53B/strainedCounts", "HB53B", "Background")
#sample13<-getFiles("/data/GI-Informatics/DePasquale/Projects/Bondoc_HB/Cell_Ranger/10X-Bondoc-HB53-Tumor-53T-10-6-2020-20201030-3v3-1hg/HB53T/strainedCounts", "HB53T", "Tumor")
#sample14<-getFiles("/data/GI-Informatics/DePasquale/Projects/Bondoc_HB/Cell_Ranger/10X-Bondoc-HB60B-20210312-3v3-1hg/HB60B/strainedCounts", "HB60B", "Background") # Lung background
#sample15<-getFiles("/data/GI-Informatics/DePasquale/Projects/Bondoc_HB/Cell_Ranger/10X-Bondoc-HB60-PDX-20210319-3v3-1hg/HB60P/strainedCounts", "HB60P", "PDX")
#sample16<-getFiles("/data/GI-Informatics/DePasquale/Projects/Bondoc_HB/Cell_Ranger/10X-Bondoc-HB60T-20210409-3v3-1hg/HB60T/strainedCounts", "HB60T", "Tumor")
#sample17<-getFiles("/data/GI-Informatics/DePasquale/Projects/Bondoc_HB/Cell_Ranger/10X-Bondoc-HB70T-20210624-3v3-1hg/HB70T/strainedCounts", "HB70T", "Tumor")
#sample18<-getFiles("/data/GI-Informatics/DePasquale/Projects/Bondoc_HB/Cell_Ranger/10X-Bondoc-HB82B-20220401-3v3-1DI-hg/HB82B/strainedCounts", "HB82B", "Background") # Lung background
#sample19<-getFiles("/data/GI-Informatics/DePasquale/Projects/Bondoc_HB/Cell_Ranger/10X-Bondoc-HB82T-20220404-3v3-1DI-hg/HB82T/strainedCounts", "HB82T", "Tumor")
sample20<-getFiles("/data/GI-Informatics/DePasquale/Projects/Bondoc_HB/Cell_Ranger/10X-Bondoc-HB83B-20220328-3v3-1DI-hg/HB83B/strainedCounts", "HB83B", "Background")

#sample21<-getFiles("/data/GI-Informatics/DePasquale/Projects/Bondoc_HB/Cell_Ranger/10X-Bondoc-HB83T-20220330-3v3-1DI-mm/HB83T/strainedCounts", "HB83T", "Tumor")
sample22<-getFiles("/data/GI-Informatics/DePasquale/Projects/Bondoc_HB/Cell_Ranger/10X-Bondoc-HB85B-20220408-3v3-1DI-hg/HB85B/strainedCounts", "HB85B", "Background")
#sample23<-getFiles("/data/GI-Informatics/DePasquale/Projects/Bondoc_HB/Cell_Ranger/10X-Bondoc-HB85-T-20220414-3v3-1DI-mm/HB85T/strainedCounts", "HB85T", "Tumor")
#sample24<-getFiles("/data/GI-Informatics/DePasquale/Projects/Bondoc_HB/Cell_Ranger/10X-Bondoc-HLO-Control-20230106-3v3-1DI-hg/HLOC/strainedCounts", "HLOC", "Tumor") # Condition not included
#sample25<-getFiles("/data/GI-Informatics/DePasquale/Projects/Bondoc_HB/Cell_Ranger/10X-Bondoc-HLO-HAR/HLOHAR/strainedCounts", "HLOHAR", "HLO") # Condition not included
#sample26<-getFiles("/data/GI-Informatics/DePasquale/Projects/Bondoc_HB/Cell_Ranger/10X-Bondoc-HLO-PCO-20220603-3v3-1DI-hg/HLOPCO/strainedCounts", "HLOPCO", "HLO") # Condition not included
#sample27<-getFiles("/data/GI-Informatics/DePasquale/Projects/Bondoc_HB/Cell_Ranger/10X-Bondoc-HLO-PPO-20220602-3v3-1DI-hg/HLOPPO/strainedCounts", "HLOPPO", "HLO") # SoupX didn't work,  # Condition not included
#sample28<-getFiles("/data/GI-Informatics/DePasquale/Projects/Bondoc_HB/Cell_Ranger/10X-Bondoc-HLO-Wnt-20230109-3v3-1DI-hg/HLOW/strainedCounts", "HLOW", "HLO") # Condition not included
sample29<-getFiles("/data/GI-Informatics/DePasquale/Projects/Bondoc_HB/Cell_Ranger/10X-Bondoc-Liver-80B-20220113-3v3-1DI-hg/80B/strainedCounts", "80B", "Background")
#sample30<-getFiles("/data/GI-Informatics/DePasquale/Projects/Bondoc_HB/Cell_Ranger/10X-Bondoc-Liver-80T-20220126-3v3-1DI-hg/80T/strainedCounts", "80T", "Tumor")


sample.all <- merge(sample9,
y = c(sample12, sample20, sample22, sample29),
add.cell.ids = c("HB31B", "HB53B", "HB83B", "HB85B", "80B"),
project = "sampleAll")


pdf(file="QC_features.pdf", width = 15, height = 8)
par(mar=c(4, 4, 4, 4))
  print(VlnPlot(sample.all, features = c("nFeature_RNA"), ncol = 1, group.by = "stim"))
dev.off()

pdf(file="QC_count.pdf", width = 15, height = 8)
par(mar=c(4, 4, 4, 4))
  print(VlnPlot(sample.all, features = c("nCount_RNA"), ncol = 1, group.by = "stim"))
dev.off()

pdf(file="QC_MT.pdf", width = 15, height = 8)
par(mar=c(4, 4, 4, 4))
  print(VlnPlot(sample.all, features = c("percent.mt"), ncol = 1, group.by = "stim"))
dev.off()

#=========================================================================================

#Integration

scrna.list <- SplitObject(sample.all, split.by = "stim")
scrna.list <- lapply(X = scrna.list, FUN = NormalizeData, normalization.method = "LogNormalize", scale.factor = 10000)
intersect_all_genes=NULL
for(i in 1:length(scrna.list)){
  scrna.list[[i]] <- NormalizeData(scrna.list[[i]])
  scrna.list[[i]] <- FindVariableFeatures(scrna.list[[i]], selection.method = "vst", nfeatures = 2000)
}

features <- SelectIntegrationFeatures(object.list = scrna.list)

save.image("environ_1_seurat.RData")
# resolution=0.5
# setwd(paste0("/data/GI-Informatics/DePasquale/Projects/Bondoc_HB/Resolution_", resolution, "noLung"))
# load("environ_1_seurat.RData")

message("before FindIntegrationAnchors")
immune.anchors <- FindIntegrationAnchors(object.list = scrna.list, anchor.features = features) # normally CCA
message("after FindIntegrationAnchors")

save.image("environ_2_seurat.RData")

message("before IntegrateData")
immune.combined2 <- IntegrateData(anchorset = immune.anchors)
message("after IntegrateData")

save.image("environ_3_seurat.RData") 

message("before ScaleData")
immune.combined2 <- ScaleData(object = immune.combined2)
message("after ScaleData")

save.image("environ_4_seurat.RData")

message("before RunPCA")
immune.combined2 <- RunPCA(immune.combined2, verbose = FALSE)
message("after RunPCA")

save.image("environ_5_seurat.RData")

immune.combined2 <- RunUMAP(immune.combined2, reduction = "pca", dims = 1:40)
immune.combined2 <- RunTSNE(immune.combined2, reduction = "pca", dims = 1:40)
immune.combined2 <- FindNeighbors(immune.combined2, reduction = "pca", dims = 1:40)
immune.combined2 <- FindClusters(immune.combined2, resolution = as.numeric(resolution))

#=========================================================================================
#Save plots and RDS

### --- save RDS file for further analysis ---###
message("saving the rds file for downstream analysis ...")
saveRDS(immune.combined2, file = "Log-Integrated.rds")

pdf("Integrated-Batch_umap.pdf", width = 11, height = 10)
DimPlot(immune.combined2, group.by = "grp", raster=FALSE, shuffle=T) + ggtitle("Integrated UMAP Plot")
dev.off()


pdf("Integrated-Clusters_umap.pdf", width = 10, height = 10)
DimPlot(immune.combined2, group.by = "seurat_clusters", label = T, raster=FALSE, shuffle=T) + NoLegend()
dev.off()

### ---save umap cordinate --- ###
message("saving umap cordinate to txt file...")
umap = as.data.frame(immune.combined2@reductions[["umap"]]@cell.embeddings)
setDT(umap, keep.rownames = "UID")[]
message("saving umap cordinates file ...")
write.table(umap,"umap_cordinate.txt", col.names = T, sep ="\t", row.names = F, quote = F)


### cellbarcode to cluster associatoion ###
seurat_clusters <- immune.combined2@meta.data[['seurat_clusters']]
Dataset <- immune.combined2@meta.data[['stim']]
UID<- rownames(immune.combined2@meta.data)
cell_clusters <- cbind(UID, Dataset, seurat_clusters)
message("saving cellbarcode_to_cluster.txt file ...")
write.table(cell_clusters, "cellbarcode_to_cluster.txt", col.names = T, sep ="\t", row.names = F, quote = F)

### --- findMarkers for each cluster --- ###
message("Finding markers for the integrated object ...")
DefaultAssay(immune.combined2) <- "RNA"
M.markers <- FindAllMarkers(immune.combined2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(M.markers, paste0("Integration_Markers_res_", resolution ,".txt"), col.names = T, sep = "\t", row.names = F, quote = F)


### -- save metadata information ---###
message("writting the meta data file ...")
Meta = as.data.frame(immune.combined2@meta.data)
write.table(Meta, "Meta_Data.txt", col.names = T, sep = "\t", row.names = F, quote = F)


#######################################
# <><><><> Erica's Additions <><><><> #
#######################################

# Visualization
p1 <- DimPlot(immune.combined2, reduction = "umap", group.by = "orig.ident", raster=FALSE, shuffle=T)
p2 <- DimPlot(immune.combined2, reduction = "umap", label = TRUE, repel = TRUE, raster=FALSE, shuffle=T)
p3 <- DimPlot(immune.combined2, reduction = "tsne", group.by = "orig.ident", raster=FALSE, shuffle=T)
p4 <- DimPlot(immune.combined2, reduction = "tsne", label = TRUE, repel = TRUE, raster=FALSE, shuffle=T)

pdf(file = "UMAP_sample_and_cluster.pdf", width = 18, height = 8.5)
par(mar=c(4, 4, 4, 4))
p1 + p2
dev.off()

pdf(file = "UMAP_split_by_sample.pdf", width = 80, height = 8.5)
par(mar=c(4, 4, 4, 4))
DimPlot(immune.combined2, reduction = "umap", split.by = "orig.ident", raster=FALSE, shuffle=T)
dev.off()

pdf(file = "TSNE_sample_and_cluster.pdf", width = 18, height = 8.5)
par(mar=c(4, 4, 4, 4))
p3 + p4
dev.off()

pdf(file = "TSNE_split_by_sample.pdf", width = 80, height = 8.5)
par(mar=c(4, 4, 4, 4))
DimPlot(immune.combined2, reduction = "tsne", split.by = "orig.ident", raster=FALSE, shuffle=T)
dev.off()

#############

source("/data/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/200116_FunctionsGeneral.R")

signatures <- read.xls("/data/GI-Informatics/DePasquale/Projects/Peters_5PrimeTCRBCR/220420_signatures_top150.xlsx", header = T, sheet = 1)
cexsize <- round(max(c(0.3, 0.5-nrow(immune.combined2@assays[["RNA"]]@data)/40000)),2)

# Color by signature score
message("\nCalculating signature scores")
# Signatures
signatures <- lapply(as.list(signatures), function(y) y[!is.na(y)])
signatures <- lapply(signatures, intersect, rownames(immune.combined2@assays[["RNA"]]@data))
# Average gene expression for scoreSignature
CM.mean <- rowMeans(immune.combined2@assays[["RNA"]]@data)
signScore <- lapply(names(signatures), function(g) {
  message(g)
  scoreSignature(CM = immune.combined2@assays[["RNA"]]@data, signatures = signatures[[g]], CM.mean = CM.mean, verbose = T)
})
names(signScore) <- names(signatures)

# Plot all signatures - UMAP
pdf(file = "UMAP_Liver_25_Seurat.pdf", width = 10, height = 10)
par(mar=c(4, 4, 4, 4))

for (n in names(signScore)) {
  mycol <- colItay(signScore[[n]])
  plot(immune.combined2@reductions[["umap"]]@cell.embeddings, pch = 16, cex = cexsize, col = mycol, main = n)
}

dev.off()

# Plot all signatures - TSNE
pdf(file = "TSNE_Liver_25_Seurat.pdf", width = 10, height = 10)
par(mar=c(4, 4, 4, 4))

for (n in names(signScore)) {
  mycol <- colItay(signScore[[n]])
  plot(immune.combined2@reductions[["tsne"]]@cell.embeddings, pch = 16, cex = cexsize, col = mycol, main = n)
}

dev.off()

# Plot Seurat UMAP
pdf(file = "UMAP_Liver_25_Seurat_Groups.pdf", width = 7, height = 6)
par(mar=c(2, 2, 2, 2))
DimPlot(object = immune.combined2, reduction = "umap",label= TRUE, raster=FALSE)
dev.off()

# Plot Seurat TSNE
pdf(file = "TSNE_Liver_25_Seurat_Groups.pdf", width = 7, height = 6)
par(mar=c(2, 2, 2, 2))
DimPlot(object = immune.combined2, reduction = "tsne",label= TRUE, raster=FALSE)
dev.off()

saveRDS(immune.combined2, file = paste0("Seurat_Log_res", resolution ,".rds"))

temp=as.data.frame(cbind(immune.combined2@meta.data[["orig.ident"]],immune.combined2@meta.data[[paste0("integrated_snn_res.", resolution)]]))
write.table(table(temp), paste0("Seurat_Log_res", resolution ,".txt"), sep="\t", quote=F)






