#################
#               #
#  Single Cell  #
#    Figure 1   #
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

# Load object
M<-readRDS("/data/GI-Informatics/DePasquale/Projects/Bondoc_HB/Resolution_0.5noLung/Seurat_Liver_25_backsub_backGroups_multiRes.rds")
DefaultAssay(M)<-"RNA"

# Remove the 96UES sample
Idents(M)<-"orig.ident"
M<-subset(M, idents=c(setdiff(unique(M$orig.ident), "96UES")))
unique(M$orig.ident) # good, it's removed

# Make healthy only subset
Idents(M)<-"SampleType"
M<-subset(M, idents=c("Background", "Tumor"))
M_B<-subset(M, idents="Background")
M_T<-subset(M, idents="Tumor")

# Make results folder
dir.create("/data/GI-Informatics/DePasquale/Projects/Bondoc_HB/Resolution_0.5noLung/Prelim_Figs/")
setwd("/data/GI-Informatics/DePasquale/Projects/Bondoc_HB/Resolution_0.5noLung/Prelim_Figs/")

# Panel A
#DimPlot(M, raster=F)
#DimPlot(M, group.by="Resolution_0.1", label=T, repel=T, raster=F)


##########
# Subcluster (12/11)

Idents(M)<-"Resolution_0.1"
DefaultAssay(M)<-"integrated"

### Cluster 4 (Immune)
M_4<-subset(M, idents=4)
M_4 <- FindClusters(M_4, resolution = 0.29) #need this to give 5 clusters
#DotPlot(M_4, assay="RNA", features=c("CD8A", "CD4", "PTPRC", "CD3G", "CD68", "SLC8A1", "AOAH", "VCAN", "S100A8", "CD14"))
temp=cbind(as.numeric(as.character(M@active.ident)),as.numeric(as.character(M@active.ident)))
row.names(temp)=M@assays[["RNA"]]@data@Dimnames[[2]]
new_clusters=mapvalues(as.numeric(as.character(M_4@active.ident)), from=c(0,1,2,3,4), to=c("Kupffer", "Kupffer", "Kupffer", "Monocyte", "Lymphocyte"))
temp[names(M_4@active.ident),2]<-new_clusters
M@meta.data[["subcluster"]] <- temp[,2]
table(M@meta.data[["subcluster"]])

### Cluster 5 (Endothelial)
M_5<-subset(M, idents=5)
M_5 <- FindClusters(M_5, resolution = 0.3)
DotPlot(M_5, assay="RNA", features=c("STAB1", "STAB2", "PDPN", "CLEC14A", "SELP"))
temp=cbind(as.character(M$subcluster),as.character(M$subcluster))
row.names(temp)=M@assays[["RNA"]]@data@Dimnames[[2]]
new_clusters=mapvalues(as.numeric(as.character(M_5@active.ident)), from=c(0,1,2), to=c("Sinusoid", "Venous", "Sinusoid"))
temp[names(M_5@active.ident),2]<-new_clusters
M@meta.data[["subcluster"]] <- temp[,2]
table(M@meta.data[["subcluster"]])

### Cluster 6 (Mesenchymal)
M_6<-subset(M, idents=6)
M_6 <- FindClusters(M_6, resolution = 0.3)
DotPlot(M_6, assay="RNA", features=c("PRKG1", "DPT", "FBLN1", "FBLN2", "VCAN", "MYH11", "MSLN"))
temp=cbind(as.character(M$subcluster),as.character(M$subcluster))
row.names(temp)=M@assays[["RNA"]]@data@Dimnames[[2]]
new_clusters=mapvalues(as.numeric(as.character(M_6@active.ident)), from=c(0,1,2), to=c("Portal Fibroblast", "Portal Fibroblast", "Stellate"))
temp[names(M_6@active.ident),2]<-new_clusters
M@meta.data[["subcluster"]] <- temp[,2]
table(M@meta.data[["subcluster"]])

### Rename other clusters
M$subcluster[c(which(M$subcluster==0),which(M$subcluster==1),which(M$subcluster==8))] <- "Hepatocyte"
M$subcluster[which(M$subcluster==2)] <- "Cycling Hepatocyte"
M$subcluster[which(M$subcluster==3)] <- "Cholangiocyte-like"
M$subcluster[which(M$subcluster==7)] <- "Cholangiocyte"
table(M@meta.data[["subcluster"]])

### Plot
DimPlot(M, group.by = "subcluster", label=T, repel=T, raster=F)

### Save new object
saveRDS(M, "/data/GI-Informatics/DePasquale/Projects/Bondoc_HB/Resolution_0.5noLung/Seurat_Liver_25_backsub_backGroups_multiRes_BTonly_subcluster.rds")

#
##########

pdf(file = "1A_Broad_Clusters_UMAP.pdf", width = 8.5, height = 7)
par(mar=c(2, 2, 2, 2))
  DimPlot(M, group.by="broad_cluster_names", label=T, repel=T, raster=F) + ggtitle("Broad Clusters")
dev.off()

pdf(file = "1A_Sub_Clusters_UMAP.pdf", width = 8.5, height = 7)
par(mar=c(2, 2, 2, 2))
  DimPlot(M, group.by="subcluster", label=T, repel=T, raster=F) + ggtitle("Broad Clusters")
dev.off()

# Panel B
pdf(file = "1B_Samples_UMAP.pdf", width = 7.5, height = 7)
par(mar=c(2, 2, 2, 2))
  DimPlot(M, group.by = "orig.ident", raster=F) + ggtitle("Samples")
dev.off()

# Panel C
pdf(file = "1C_Broad_Clusters_Split_UMAP.pdf", width = 14.5, height = 7)
par(mar=c(2, 2, 2, 2))
  DimPlot(M, group.by = "broad_cluster_names", split.by = "SampleType", raster=F, label=T, repel=T) + ggtitle(NULL)
dev.off()

pdf(file = "1C_Sub_Clusters_Split_UMAP.pdf", width = 14.5, height = 7)
par(mar=c(2, 2, 2, 2))
  DimPlot(M, group.by = "subcluster", split.by = "SampleType", raster=F, label=T, repel=T) + ggtitle(NULL)
dev.off()

# Panel D
data_long=as.data.frame(cbind(Sample=M@meta.data[["SampleType"]], Cluster=as.character(M@meta.data[["broad_cluster_names"]])))
data=as.data.frame.matrix(table(data_long))
data=t(data) # Read and wrangle data
data <- sweep(data, 2, colSums(data), "/")*100
data=data[c(5,3,2,6,4,7,1),] # Order clusters
myColors=hue_pal()(7)

pdf("1D_Frequencies.pdf", width = 6, height = 7)
par(mar = c(8,4,4,16), xpd = T)
  barplot(as.matrix(data[nrow(data):1,]), col = rev(myColors), xaxt = "n", ylab = "Population frequency (%)", border = NA)
  axis(side = 1, at = seq(1,ncol(data))*1.2-0.5, labels = colnames(data), las = 2)
  legend(x = ncol(data)*1.2+0.5, y = 100, legend = row.names(data), fill = myColors, bty = "n", border = NA)
dev.off()

data_long=as.data.frame(cbind(Sample=M@meta.data[["SampleType"]], Cluster=as.character(M@meta.data[["subcluster"]])))
data=as.data.frame.matrix(table(data_long))
data=t(data) # Read and wrangle data
data <- sweep(data, 2, colSums(data), "/")*100
myColors=hue_pal()(11)

pdf("1D_Sub_Frequencies.pdf", width = 6, height = 7)
par(mar = c(8,4,4,16), xpd = T)
  barplot(as.matrix(data[nrow(data):1,]), col = rev(myColors), xaxt = "n", ylab = "Population frequency (%)", border = NA)
  axis(side = 1, at = seq(1,ncol(data))*1.2-0.5, labels = colnames(data), las = 2)
  legend(x = ncol(data)*1.2+0.5, y = 100, legend = row.names(data), fill = myColors, bty = "n", border = NA)
dev.off()

# Panel E
data_long=as.data.frame(cbind(Sample=M@meta.data[["orig.ident"]], Cluster=as.character(M@meta.data[["broad_cluster_names"]])))
data=as.data.frame.matrix(table(data_long))
data=cbind(sample=row.names(data), data)
data=data[order(data$sample),]
data=t(data)

x=as.data.frame(cbind(M@meta.data$orig.ident, M@meta.data$SampleType))
y=distinct(x)
y=y[order(y$V1),]
sum(y$V1!=colnames(data)) # 0, meaning they match

data=data[-1,]
rnames=row.names(data)
rnames2=gsub(" ", "", rnames)
rnames3=gsub("-", "", rnames2)
rnames4=gsub("/", "", rnames3)
data=apply(data, 2, as.numeric)
data <- sweep(data, 2, colSums(data), "/")*100
row.names(data)<-rnames4

for(i in row.names(data)){
  z=cbind(y, data[i,])
  colnames(z)<-c("sample", "SampleType", i)
  z$SampleType=factor(z$SampleType, levels = c("Background", "Tumor"))
  formula <- paste(i, "~ SampleType")
  print(compare_means(as.formula(formula),  data = z))
}
# .y.               group1      group2    p       p.adj   p.format  p.signif  method  
#   1 Cholangiocyte Background  Tumor     0.276   0.28    0.28      ns        Wilcoxon

# .y.                   group1      group2    p          p.adj     p.format  p.signif  method  
#   1 Cholangiocytelike Background  Tumor     0.000431   0.00043   0.00043   ***       Wilcoxon

# .y.                   group1      group2    p        p.adj   p.format  p.signif  method  
#   1 CyclingHepatocyte Background  Tumor     0.0182   0.018   0.018     *         Wilcoxon

# .y.             group1      group2    p       p.adj   p.format  p.signif  method  
#   1 Endothelial Background  Tumor     0.250   0.25    0.25      ns        Wilcoxon

# .y.             group1      group2    p         p.adj     p.format  p.signif  method  
#   1 Hepatocyte  Background  Tumor     0.000215  0.00022   0.00022   ***       Wilcoxon

# .y.         group1      group2   p    p.adj   p.format  p.signif  method  
#   1 Immune  Background  Tumor    1    1       1         ns        Wilcoxon

# .y.             group1      group2    p       p.adj   p.format  p.signif  method  
#   1 Mesenchymal Background  Tumor     0.102   0.1     0.1       ns        Wilcoxon

my_comparisons <- list( c("Background", "Tumor"))

i="Hepatocyte"
z=cbind(y, data[i,])
colnames(z)<-c("sample", "SampleType", i)
z$SampleType=factor(z$SampleType, levels = c("Background", "Tumor"))

pdf(paste0("1E_Hepatocyte_Stats.pdf"), width = 3, height = 5)
par(mar = c(8,4,4,16), xpd = T)
ggboxplot(z, x = "SampleType", y = "Hepatocyte", color = "SampleType", add = "jitter", legend = "none")+ 
  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  labs(x="Sample Type", y = paste0("Proportion Hepatocyte (%)"))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),axis.title.x=element_blank())+
  ggtitle("Hepatocyte")
dev.off()

i="Cholangiocytelike"
z=cbind(y, data[i,])
colnames(z)<-c("sample", "SampleType", i)
z$SampleType=factor(z$SampleType, levels = c("Background", "Tumor"))

pdf(paste0("1E_CholangiocyteLike_Stats.pdf"), width = 3, height = 5)
par(mar = c(8,4,4,16), xpd = T)
ggboxplot(z, x = "SampleType", y = "Cholangiocytelike", color = "SampleType", add = "jitter", legend = "none")+ 
  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  labs(x="Sample Type", y = paste0("Proportion Cholangiocyte-like (%)"))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),axis.title.x=element_blank())+
  ggtitle("Cholangiocyte-like")
dev.off()

i="CyclingHepatocyte"
z=cbind(y, data[i,])
colnames(z)<-c("sample", "SampleType", i)
z$SampleType=factor(z$SampleType, levels = c("Background", "Tumor"))

pdf(paste0("1E_CyclingHepatocyte_Stats.pdf"), width = 3, height = 5)
par(mar = c(8,4,4,16), xpd = T)
ggboxplot(z, x = "SampleType", y = "CyclingHepatocyte", color = "SampleType", add = "jitter", legend = "none")+ 
  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  labs(x="Sample Type", y = paste0("Proportion Cycling Hepatocyte (%)"))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),axis.title.x=element_blank())+
  ggtitle("Cycling Hepatocyte")
dev.off()

data_long=as.data.frame(cbind(Sample=M@meta.data[["orig.ident"]], Cluster=as.character(M@meta.data[["subcluster"]])))
data=as.data.frame.matrix(table(data_long))
data=cbind(sample=row.names(data), data)
data=data[order(data$sample),]
data=t(data)

x=as.data.frame(cbind(M@meta.data$orig.ident, M@meta.data$SampleType))
y=distinct(x)
y=y[order(y$V1),]
sum(y$V1!=colnames(data)) # 0, meaning they match

data=data[-1,]
rnames=row.names(data)
rnames2=gsub(" ", "", rnames)
rnames3=gsub("-", "", rnames2)
rnames4=gsub("/", "", rnames3)
data=apply(data, 2, as.numeric)
data <- sweep(data, 2, colSums(data), "/")*100
row.names(data)<-rnames4

for(i in row.names(data)){
  z=cbind(y, data[i,])
  colnames(z)<-c("sample", "SampleType", i)
  z$SampleType=factor(z$SampleType, levels = c("Background", "Tumor"))
  formula <- paste(i, "~ SampleType")
  print(compare_means(as.formula(formula),  data = z))
}

# .y.           group1     group2     p p.adj p.format p.signif method  
#   1 Cholangiocyte Background Tumor  0.276  0.28 0.28     ns       Wilcoxon

# .y.               group1     group2        p   p.adj p.format p.signif method  
#   1 Cholangiocytelike Background Tumor  0.000431 0.00043 0.00043  ***      Wilcoxon

# .y.               group1     group2      p p.adj p.format p.signif method  
#   1 CyclingHepatocyte Background Tumor  0.0182 0.018 0.018    *        Wilcoxon

# .y.        group1     group2        p   p.adj p.format p.signif method  
#   1 Hepatocyte Background Tumor  0.000215 0.00022 0.00022  ***      Wilcoxon

# .y.     group1     group2     p p.adj p.format p.signif method  
#   1 Kupffer Background Tumor  0.820  0.82 0.82     ns       Wilcoxon

# .y.        group1     group2     p p.adj p.format p.signif method  
#   1 Lymphocyte Background Tumor  0.962  0.96 0.96     ns       Wilcoxon

# .y.      group1     group2     p p.adj p.format p.signif method  
#   1 Monocyte Background Tumor  0.276  0.28 0.28     ns       Wilcoxon

# .y.              group1     group2     p p.adj p.format p.signif method  
#   1 PortalFibroblast Background Tumor  0.494  0.49 0.49     ns       Wilcoxon

# .y.      group1     group2      p p.adj p.format p.signif method  
#   1 Sinusoid Background Tumor  0.0529 0.053 0.053    ns       Wilcoxon

# .y.      group1     group2       p  p.adj p.format p.signif method  
#   1 Stellate Background Tumor  0.00229 0.0023 0.0023   **       Wilcoxon

# .y.    group1     group2      p p.adj p.format p.signif method  
#   1 Venous Background Tumor  0.0320 0.032 0.032    *        Wilcoxon

my_comparisons <- list( c("Background", "Tumor"))

i="Stellate"
z=cbind(y, data[i,])
colnames(z)<-c("sample", "SampleType", i)
z$SampleType=factor(z$SampleType, levels = c("Background", "Tumor"))

pdf(paste0("1E_Stellate_Stats.pdf"), width = 3, height = 5)
par(mar = c(8,4,4,16), xpd = T)
ggboxplot(z, x = "SampleType", y = "Stellate", color = "SampleType", add = "jitter", legend = "none")+ 
  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  labs(x="Sample Type", y = paste0("Proportion Stellate (%)"))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),axis.title.x=element_blank())+
  ggtitle("Stellate")
dev.off()

i="Venous"
z=cbind(y, data[i,])
colnames(z)<-c("sample", "SampleType", i)
z$SampleType=factor(z$SampleType, levels = c("Background", "Tumor"))

pdf(paste0("1E_Venous_Stats.pdf"), width = 3, height = 5)
par(mar = c(8,4,4,16), xpd = T)
ggboxplot(z, x = "SampleType", y = "Venous", color = "SampleType", add = "jitter", legend = "none")+ 
  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  labs(x="Sample Type", y = paste0("Proportion Venous Hepatocyte (%)"))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),axis.title.x=element_blank())+
  ggtitle("Venous")
dev.off()


# Panel F
Idents(M)<-"broad_cluster_names"
levels(x = M) <- rev(c("Hepatocyte", "Cycling Hepatocyte", "Cholangiocyte", "Cholangiocyte-like", "Immune", "Endothelial", "Mesenchymal"))

pdf(file = "1F_Dot_Plot.pdf", width = 9, height = 7)
par(mar=c(2, 2, 2, 2))
  DotPlot(M, features=c("ARG1", "MKI67", "KRT7", "NKD1", "PTPRC", "STAB2", "PRKG1"), dot.scale=12)
dev.off()

results=FindMarkers(M, ident.1 = "Cholangiocyte-like", only.pos=T) # find good marker for Cholangiocyte-like
write.table(results, "1F_CholangiocyteLike_Markers.txt", sep="\t", quote=F)

Idents(M)<-"subcluster"
levels(x = M) <- rev(c("Hepatocyte", "Cycling Hepatocyte", "Cholangiocyte", "Cholangiocyte-like", "Lymphocyte", "Kupffer", "Monocyte", "Portal Fibroblast", "Stellate", "Sinusoid", "Venous"))

pdf(file = "1F_Sub_Dot_Plot.pdf", width = 9, height = 5)
par(mar=c(2, 2, 2, 2))
  DotPlot(M, assay="RNA", features=c("ARG1", "MKI67", "KRT7", "NKD1", "PTPRC", "CD3G", "CD14", "SLC8A1", "VCAN", "PRKG1", "STAB2", "CLEC14A"), dot.scale=12) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()

############
# Plots for Katie poster
pdf(file = "Poster_UMAP_EZH2.pdf", width = 5.25, height = 5)
par(mar=c(2, 2, 2, 2))
  FeaturePlot(M_T, features="EZH2", min.cutoff = 0, raster=F)
dev.off()

pdf(file = "Poster_UMAP_AURKB.pdf", width = 5.25, height = 5)
par(mar=c(2, 2, 2, 2))
  FeaturePlot(M_T, features="AURKB", min.cutoff = 0, raster=F)
dev.off()

pdf(file = "Poster_UMAP_MKI67.pdf", width = 5.25, height = 5)
par(mar=c(2, 2, 2, 2))
  FeaturePlot(M_T, features="MKI67", min.cutoff = 0, raster=F)
dev.off()
