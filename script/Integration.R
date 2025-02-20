suppressPackageStartupMessages({
library(doParallel)
library(Seurat)
library(SeuratWrappers)
library(harmony)
library(Rcpp)
library(ggplot2)
library(gridExtra)
library(dplyr)
})
rm(list = ls())
dim = 30
resolution = 0.9

data.dir <- "Specify your data stored path"

# load data
#### load data
object <- readRDS(paste0(data.dir,"/10XAll.pc30.res0.9.rds"))

cat("construct data_list")
ifnb_list <- list()
unique_ids <- unique(object@meta.data$ID)
unique_ids <- as.character(unique_ids)
unique_ids <- unique_ids[!unique_ids %in% c("10XS3", "10XS4", "10XS5")]
for (id in unique_ids) {
  idx <- which(object@meta.data$ID == id)
  subset_object <- object[,idx]
  subset_object@meta.data$orig.ident <- id
  ifnb_list[[id]] <- subset_object
}

## load ARD 2023 OA data 
## Senescent cell population with ZEB1 transcription factor as its main regulator promotes osteoarthritis in cartilage and meniscus
## link: https://pubmed.ncbi.nlm.nih.gov/36564153/
# CartOA1
object <- readRDS(paste0(data.dir,"/CartOA1.rds"))
object@meta.data$orig.ident <- "CartOA1"
CartOA1 <-object

# CartOA2
object<- readRDS(paste0(data.dir,"/CartOA2.rds"))
object@meta.data$orig.ident <- "CartOA2"
CartOA2 <- object

# CartOA3
object <- readRDS("/data/public/OA/hannah/CartOA3.rds")
object@meta.data$orig.ident <- "CartOA3"
CartOA3 <- object

# CartOA4
object <- readRDS("/data/public/OA/hannah/CartOA4.rds")
object@meta.data$orig.ident <- "CartOA4"
CartOA4 <- object

# CartOA5
object <- readRDS("/data/public/OA/hannah/CartOA5.rds")
object@meta.data$orig.ident <- "CartOA5"
CartOA5 <- object

# CartOA6
object <- readRDS("/data/public/OA/hannah/CartOA6.rds")
object@meta.data$orig.ident <- "CartOA6"
CartOA6 <- object
rm(object)

ifnb_list$CartOA1 <- CartOA1
ifnb_list$CartOA2 <- CartOA2
ifnb_list$CartOA3 <- CartOA3
ifnb_list$CartOA4 <- CartOA4
ifnb_list$CartOA5 <- CartOA5
ifnb_list$CartOA6 <- CartOA6
combined <- merge(x=ifnb_list[[1]],y=ifnb_list[-1],project='CartOA',add.cell.ids = names(ifnb_list),merge.data = T)
DefaultAssay(combined) <- "RNA"
combined$group <- as.factor(combined$orig.ident)
levels(combined$group) <- c(rep("Own",2),rep("ARD",6),rep("Own",14))
combined@meta.data$seurat_clusters <- NULL


# run harmony
cat("running harmony")
features <- 2500
combined <- NormalizeData(combined)
combined <- FindVariableFeatures(combined, selection.method = "vst", nfeatures = features) 
combined <- ScaleData(combined)
combined <- RunPCA(combined,reduction.name='pca',npcs = dim) 
combined <- RunHarmony(combined,reduction = 'pca',group.by.vars ='orig.ident', reduction.save = 'harmony')
combined <- RunUMAP(combined,reduction='harmony',dims=1:dim) 
combined <- FindNeighbors(combined, reduction = "harmony", dims = 1:dim)
combined <- FindClusters(combined, resolution = resolution)
saveRDS(combined, file = "/data/public/OA/hannah/all.dim_30_res_0.9.harmony.rds")

# Findmarkers
library(openxlsx)
count <- 1
markers_wilcox.list = list()
cat("running All FindMarkers")
for (i in levels(combined@meta.data$seurat_clusters)) {
  res <- FindMarkers(combined,slot="data",
                     ident.1=i,
                     ident.2=NULL,
                     logfc.threshold=0.5,
                     test.use='wilcox',
                     min.pct=0.25, 
                     min.diff.pct=-Inf, 
                     group.by='seurat_clusters')
  markers_wilcox.list[[count]] <- res[order(res$avg_log2FC, decreasing=TRUE), ]
  count <- count + 1
}## end for
