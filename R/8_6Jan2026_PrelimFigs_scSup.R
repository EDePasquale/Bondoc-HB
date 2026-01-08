################################################
#                                              #
# Additional figure based on reviewer comments #
#                Bondoc - HB                   #
#               6 January 2026                 #
#                                              #
################################################

# Load libraries
library(Seurat)
library(plyr)
library(ggplot2)
library(scales)
library(tidyr)
library(dplyr)
library(ggpubr)
library(speckle)

# Load in object
M2 <- readRDS("/data/GI-Informatics/DePasquale/Projects/Bondoc_HB/Resolution_0.5noLung/Seurat_Liver_25_backsub_backGroups_multiRes_BTonly_subcluster.rds")
setwd("/data/GI-Informatics/DePasquale/Projects/Bondoc_HB/Resolution_0.5noLung/Prelim_Figs/")

# Make new supplemental figure
pdf("Supplemental_split_by_sample.pdf", width = 18, height = 9.5)
par(mar = c(2,2,2,2))
  DimPlot(M2, group.by = "subcluster", split.by = "orig.ident", ncol = 6, raster = F, pt.size = 0.1) + ggtitle("Cluster by Sample")
dev.off()

# New method for cell type proportion differences (propeller)
prop_results <- propeller(x = M2, 
                          clusters = M2$subcluster, 
                          sample = M2$orig.ident, 
                          group = M2$SampleType)
knitr::kable(print(prop_results))
#   |                   |BaselineProp.clusters | BaselineProp.Freq| PropMean.Background| PropMean.Tumor| PropRatio| Tstatistic|   P.Value|       FDR|
#   |:------------------|:---------------------|-----------------:|-------------------:|--------------:|---------:|----------:|---------:|---------:|
#   |Cholangiocyte-like |Cholangiocyte-like    |         0.0845661|           0.0039880|      0.1198565| 0.0332733| -6.1942589| 0.0000051| 0.0000563|
#   |Stellate           |Stellate              |         0.0040542|           0.0001690|      0.0050708| 0.0333317| -4.6283896| 0.0001689| 0.0009291|
#   |Hepatocyte         |Hepatocyte            |         0.7142313|           0.8616270|      0.6364671| 1.3537652|  3.8484410| 0.0010284| 0.0037707|
#   |Cycling Hepatocyte |Cycling Hepatocyte    |         0.0819540|           0.0227360|      0.1305887| 0.1741040| -2.9940156| 0.0072652| 0.0199792|
#   |Venous             |Venous                |         0.0143029|           0.0086139|      0.0162687| 0.5294755| -2.1687639| 0.0425491| 0.0936080|
#   |Sinusoid           |Sinusoid              |         0.0300208|           0.0423711|      0.0230387| 1.8391311|  1.9808436| 0.0617709| 0.1132467|
#   |Cholangiocyte      |Cholangiocyte         |         0.0042990|           0.0001755|      0.0063850| 0.0274879| -0.9163945| 0.3708963| 0.5828371|
#   |Monocyte           |Monocyte              |         0.0062762|           0.0013632|      0.0060157| 0.2266142| -0.5843491| 0.5658303| 0.7780166|
#   |Lymphocyte         |Lymphocyte            |         0.0033286|           0.0041409|      0.0021288| 1.9451554|  0.3732409| 0.7130831| 0.8702943|
#   |Portal Fibroblast  |Portal Fibroblast     |         0.0174773|           0.0130900|      0.0199957| 0.6546390| -0.1660898| 0.8697823| 0.8702943|
#   |Kupffer            |Kupffer               |         0.0394896|           0.0417253|      0.0341844| 1.2205953|  0.1654303| 0.8702943| 0.8702943|
  
  