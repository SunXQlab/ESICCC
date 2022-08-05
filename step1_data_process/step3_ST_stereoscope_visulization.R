setwd("./step1_data_process/")

rm(list = ls())
library(Seurat)
library(tibble)
library(stringr)
set.seed(123)

# Create Seurat Object for ST data
expr.data <- Read10X("./dataset/ST.Data/filtered_count_matrices/CID4465_filtered_count_matrix/", 
                             gene.column = 1)
se.st <- CreateSeuratObject(counts = expr.data, assay = 'Spatial', min.cells = 10, min.features = 500)

# Load the image data
img <- Read10X_Image(image.dir = './dataset/ST.Data/spatial/CID4465_spatial/')
DefaultAssay(object = img) <- 'Spatial'
img <- img[colnames(se.st)]
se.st[['image']] <- img

# Load the orignal metadata
meta.tmp <- read.csv("./dataset/ST.Data/metadata/CID4465_metadata.csv")
meta.tmp <- column_to_rownames(meta.tmp, var = "X")
meta.tmp <- meta.tmp[colnames(se.st),3:5]
identical(rownames(meta.tmp), colnames(se.st))
se.st <- AddMetaData(se.st, meta.tmp)

# Load the celltype metadata from stereoscope
meta.tmp <- read.table("./result/stereoscope_result_CID4465/CID4465_spatial_cnt/W.2022-05-15135505.473576.tsv", 
                       header = TRUE, sep = "\t")
meta.tmp <- column_to_rownames(meta.tmp, var = "X")
colnames(meta.tmp) <- str_replace_all(colnames(meta.tmp), "\\.", "_")
identical(rownames(meta.tmp), colnames(se.st))

meta.tmp$celltype <- apply(meta.tmp, 1, function(x){
  colnames(meta.tmp)[which.max(x)]
})
se.st <- AddMetaData(se.st, meta.tmp)
save(se.st, file = './se_data/ST.se.RData')

p1 <- SpatialDimPlot(se.st, group.by = "celltype") 
p2 <- SpatialDimPlot(se.st, group.by = "Classification") 

p1|p2
