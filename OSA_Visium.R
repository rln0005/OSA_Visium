#### Spatial Transcriptome (Visium) Analysis of Canine OSA ####
# Written by: Rebecca Nance-Richey, 07/17/2024 
# Using R version 4.4.0 ("Puppy Cup")

library(SpotClean)
library(Seurat)
library(tidyverse)
library(SpaCET)
library(ggplot2)

#### 1. SpotClean (~2 hr run time) ####
basePath <- "U:/GitHub/OSA_Visium/raw_data/"
sampleNames <- list.files(path = basePath, pattern = NULL, all.files = FALSE, full.names = F)
files <- setNames(paste0(basePath, sampleNames), sampleNames)
#SpotClean function for each sample
lapply(1:length(files), function(i){
  sampleName <- names(files[i])
  matrix <- read10xRaw(paste0(files[i], "/raw_feature_bc_matrix"))
  slide <- read10xSlide(tissue_csv_file=paste0(files[i], "/spatial/tissue_positions.csv"),
                        tissue_img_file=paste0(files[i], "/spatial/tissue_lowres_image.png"),
                        scale_factor_file=paste0(files[i], "/spatial/scalefactors_json.json"))
  imagedir <- paste0(files[i], "/spatial")
  obj <- createSlide(count_mat=matrix, slide_info=slide)
  decont_obj <- spotclean(obj)
  seurat_obj <- convertToSeurat(decont_obj, slice=sampleName, image_dir = imagedir )
  seurat_obj <- AddMetaData(seurat_obj, metadata = sampleName, col.name = "sample")
  
  #save the RDS file for subsequent processing
  saveRDS(seurat_obj, paste0("U:/GitHub/OSA_Visium/processed_data/", sampleName, ".rds"))
})


#### 2. Individual Seurat Processing (~30 min run time) ####
basePath <- "U:/GitHub/OSA_Visium/processed_data/"
files <- list.files(path = basePath, pattern = ".rds", all.files = FALSE, full.names = F)
sampleNames <- sub('\\.rds$', '', files) 
files <- setNames(paste0(basePath, sampleNames, ".rds") , sampleNames)

lapply(1:length(files), function(i){
  inPath <- unname(files[i])
  sampleName <- names(files[i])
  
  seurat_obj <- readRDS(inPath)
  seurat_obj <- subset(x = seurat_obj, subset = nFeature_Spatial > 1500 & nCount_Spatial > 1000)
  seurat_obj <- seurat_obj %>% SCTransform(., assay="Spatial") %>% RunPCA() %>% FindNeighbors(., dims = 1:30, reduction = "pca") %>% FindClusters(., resolution = 0.4) %>% RunUMAP(., reduction="pca", dims=1:30)
  
  DimPlot(seurat_obj,  pt.size = 2.2, reduction = "umap", label = TRUE)
  ggsave(paste0("U:/GitHub/OSA_Visium/results/", sampleName, "_umap.png"), height = 6, width = 8)
  p <- SpatialDimPlot(seurat_obj, repel=TRUE, label=TRUE, pt.size.factor = 2.2, label.size=20) + NoLegend()
  p + scale_fill_manual(values=c("#F8766D", "#B79F00", "#00BA38", "#00BFC4", "#619CFF", "#F564E3"))
  ggsave(paste0("U:/GitHub/OSA_Visium/results/", sampleName, "_image.png"), height = 12, width = 14)
  saveRDS(seurat_obj, paste0("U:/GitHub/OSA_Visium/processed_data/", sampleName, ".rds"))
})


#### 3. Integrated Merged (Seurat) Analysis (~40 min run time) ####
basePath <- "U:/GitHub/OSA_Visium/processed_data/"
rds_files <- list.files(basePath, pattern = "\\.rds$", full.names = TRUE)
# Load each RDS file and assign it to a variable named after the file (without extension)
for (file in rds_files) {
  # Extract the base name without extension
  var_name <- tools::file_path_sans_ext(basename(file))
  # Read the RDS file and assign it to the variable in the environment
  assign(var_name, readRDS(file), envir = .GlobalEnv)
}
#merge samples (raw counts) into one seurat object
merge <- merge(x=OSA1, y=c(OSA1b,OSA2,OSA2b,OSA3,OSA4,OSA5,OSA6,OSA7,OSA8,OSA9,OSA9b,OSA10,OSA11), 
               add.cell.ids=c("OSA1","OSA1b","OSA2","OSA2b","OSA3","OSA4","OSA5","OSA6","OSA7","OSA8","OSA9","OSA9b","OSA10","OSA11"))
rm(list=setdiff(ls(), "merge"))
dim(merge)
merge <- subset(x = merge, subset = nFeature_Spatial > 1500 & nCount_Spatial > 1000)
#quality control of merged dataset
VlnPlot(merge, features=c("nFeature_Spatial", "nCount_Spatial"))
VlnPlot(merge, group.by="sample", features=c("nCount_Spatial"))
FeatureScatter(merge, feature1="nCount_Spatial", feature2="nFeature_Spatial")
#transform data
merge <- SCTransform(merge, assay="Spatial")
top20 <- head(VariableFeatures(merge), 20)
merge <- RunPCA(merge, assay="SCT")
#QC of normalized data--examine PCA heatmaps and elbow plot to evaluate dimensionality of dataset
DimHeatmap(merge, dims=1:9)
ElbowPlot(merge)
DimPlot(merge)

#identify integration anchors between samples, cluster, find neighbors, UMAP visualization (RPCA integration)
merge <- IntegrateLayers(object=merge, method=RPCAIntegration, normalization.method="SCT")
#join the data layers
merge[["Spatial"]] <- JoinLayers(merge[["Spatial"]])
merge <- FindNeighbors(merge,reduction="integrated.dr", dims=1:30)
merge <- FindClusters(merge, reduction="integrated.dr", resolution=0.4)
merge <- RunUMAP(merge, dims=1:30, reduction="integrated.dr")
DimPlot(merge, reduction="umap", group.by=c("sample", "seurat_clusters"))
p1 <- DimPlot(merge, reduction="umap", group.by="sample") + ggtitle("Integrated Datasets")
p2 <- DimPlot(merge, reduction="umap", group.by="seurat_clusters", label=TRUE, label.size=4) + ggtitle("Clusters")
p1 + p2
saveRDS(merge, paste0("U:/GitHub/OSA_Visium/processed_data/merge.rds"))
#get cluster totals
table(Idents(merge))
prop.table(table(Idents(merge)))
c <- table(Idents(merge), merge$sample)
write.csv(c, "U:/GitHub/OSA_Visium/results/cluster_totals.csv")


#### 4. Identify Spatially Variable Features {Moran's I} ####
basePath <- "U:/GitHub/OSA_Visium/processed_data/"
sampleNames <- list.files(path = basePath, pattern = ".rds", all.files = FALSE, full.names = F)
sampleNames <- sub('\\.rds$', '', sampleNames) 
sampleNames <- setdiff(sampleNames, "merge")
files <- setNames(paste0(basePath, sampleNames, ".rds") , sampleNames)
lapply(1:length(files), function(i){
  inPath <- unname(files[i])
  sampleName <- names(files[i])
  
  seurat_obj <- readRDS(inPath)
  seurat_obj <- FindSpatiallyVariableFeatures(seurat_obj, assay="SCT", features=VariableFeatures(seurat_obj)[1:1000], selection.method=c("moransi"))
  moransI <- seurat_obj@assays[["SCT"]]@meta.features
  moransI <- subset(x=moransI,subset = moransi.spatially.variable==TRUE) 
  moransI <- subset(x=moransI,subset = MoransI_p.value < 0.05) 
  write.csv(moransI, paste0("U:/GitHub/OSA_Visium/results/", sampleName, "_moransI.csv"))
  
  moransI <- moransI[order(moransI$MoransI_observed),]
  top.up.features <- moransI[1:6,]
  top.up.features <- rownames(top.up.features)
  SpatialFeaturePlot(seurat_obj, features=top.up.features, ncol=3, pt.size.factor = 2.2)
  ggsave(paste0("U:/GitHub/OSA_Visium/results/", sampleName, "_moransI_dispersed.png"), height = 10, width = 14)
  
  moransI <- moransI[order(-moransI$MoransI_observed), ]
  top.down.features <- moransI[1:6,]
  top.down.features <- rownames(top.down.features)
  SpatialFeaturePlot(seurat_obj, features=top.down.features, ncol=3, pt.size.factor = 2.2)
  ggsave(paste0("U:/GitHub/OSA_Visium/results/", sampleName, "_moransI_clustered.png"), height = 10, width = 14)
})


######## Find the common spatially variable genes across all samples
### Find what's in common between the samples
basePath <- "U:/GitHub/OSA_Visium/results/"
files <- list.files(path = basePath, pattern = "_moransI.csv", all.files = F, full.names=T)
#function to read first column
read_first_column <- function(file_name) {
  data <- read.csv(file_name)
  first_column <- data[, 1]
  return(first_column)
}
# Read first columns of each dataset
first_columns <- lapply(files, read_first_column)
# Find common elements
common_elements <- as.data.frame(Reduce(intersect, first_columns))
colnames(common_elements) <- "gene"
common_elements

### Subset common elements and values for each sample for the 10 genes that are in common to all
basePath <- "U:/GitHub/OSA_Visium/results/"
files <- list.files(path = basePath, pattern = "_moransI.csv", all.files = F, full.names=T)
sampleNames <- sub("^(?:[^/]*/){4}([^_]*).*", "\\1", files)
files <- setNames(files, sampleNames)
lapply(1:length(files), function(i){
  inPath <- unname(files[i])
  sampleName <- names(files[i])
  data <- read.csv(inPath)
  data <- data %>%
    filter(str_detect(X, "IGFBP7|COL4A2|COL1A1|COL18A1|APOE|CD74|HLA-DRB1|C1QA|C1QC|ACTA2"))
  data <- data[, c("X", "MoransI_observed", "MoransI_p.value")] 
  colnames(data) <- c("gene", paste0(sampleName, "_MoransI"), paste0(sampleName, "_MoransI_pvalue"))
  write.csv(data, row.names=FALSE, paste0("U:/GitHub/OSA_Visium/results/", sampleName, "_common.csv"))
})

#combine them into 1 file while simultaneously extracting their values
basePath <- "U:/GitHub/OSA_Visium/results/"
files <- list.files(path = basePath, pattern = "_common.csv", all.files = F, full.names=T)
list <- lapply(files, function(file){
  read.csv(file)
})
MyMerge <- function(x, y) {
  merge(x, y, by = "gene", all = TRUE)
}
combine <- Reduce(MyMerge, list)
write.csv(combine, row.names = F, "U:/GitHub/OSA_Visium/results/common_combined.csv") 

### Make figure of the common elements overlaid on each image
basePath <- "U:/GitHub/OSA_Visium/processed_data/"
sampleNames <- list.files(path = basePath, pattern = ".rds", all.files = FALSE, full.names = F)
sampleNames <- sub('\\.rds$', '', sampleNames) 
sampleNames <- setdiff(sampleNames, "merge")
files <- setNames(paste0(basePath, sampleNames, ".rds") , sampleNames)
lapply(1:length(files), function(i){
  inPath <- unname(files[i])
  sampleName <- names(files[i])
  
  seurat_obj <- readRDS(inPath)
  SpatialFeaturePlot(seurat_obj, features = c("APOE", "C1QC", "C1QA", "CD74", "COL1A1", "HLA-DRB1", "IGFBP7", "COL4A2", "ACTA2", "COL18A1"), pt.size.factor = 2.2, image.alpha=0.8, ncol=5, crop = TRUE)
  ggsave(paste0("U:/GitHub/OSA_Visium/results/", sampleName, "_common.png"), height = 7, width = 12) 
})


### Make example figure of Moran's I clustered vs dispersed with OSA5
a <- readRDS("U:/GitHub/OSA_Visium/processed_data/OSA5.rds")
SpatialFeaturePlot(a, features = c("S100A6", "ATRX"), pt.size.factor = 2.2, image.alpha=0.8, crop = TRUE)
ggsave(paste0("U:/GitHub/OSA_Visium/results/OSA5_example_MoransI.png"), height = 6, width = 8)



#### 5. SpaCET ####
basePath <- "U:/GitHub/OSA_Visium/processed_data/"
sampleNames <- list.files(path = basePath, pattern = NULL, all.files = FALSE, full.names = F)
files <- setNames(paste0(basePath, sampleNames, "/outs/") , sampleNames)

lapply(1:length(files), function(i){
  #set sample specific variables
  inPath <- unname(files[i])
  sampleName <- names(files[i])
  SpaCET_obj <- create.SpaCET.object.10X(visiumPath = inPath)
  
  #calculate and plot the QC metrics
  SpaCET_obj <- SpaCET.quality.control(SpaCET_obj, min.genes=1500)
  SpaCET.visualize.spatialFeature(
    SpaCET_obj, 
    spatialType = "QualityControl", 
    spatialFeatures=c("UMI","Gene"),
    imageBg = TRUE
  )
  ggsave(paste0("U:/GitHub/OSA_Visium/results", sampleName, "_SPACET_qc.png"), height = 7, width = 14)
  
  #deconvolute using PANCAN dictionary 
  SpaCET_obj <- SpaCET.deconvolution(
    SpaCET_obj,
    cancerType="PANCAN",
    coreNo = 8
  )
  #visualize the cell type proportion 
  p <- SpaCET.visualize.spatialFeature(
    SpaCET_obj, 
    spatialType = "CellFraction", 
    spatialFeatures = "All", 
    sameScaleForFraction = TRUE,
    pointSize = 0.1, 
    nrow = 5
  )
  p
  ggsave(paste0("U:/GitHub/OSA_Visium/results", sampleName, "_SPACET_spatialFeature.png"), height = 7, width = 14)
  
  #identify & Visualize the Tumor-Stroma Interface
  SpaCET_obj <- SpaCET.identify.interface(SpaCET_obj)
  SpaCET.visualize.spatialFeature(SpaCET_obj, pointSize=3, spatialType = "Interface", spatialFeature = "Interface")
  ggsave(paste0("U:/GitHub/OSA_Visium/results", sampleName, "_SPACET_interface.png"), height = 12, width = 14)
  
  saveRDS(SpaCET_obj, paste0("U:/GitHub/OSA_Visium/processed_data/", sampleName, ".rds"))
})




