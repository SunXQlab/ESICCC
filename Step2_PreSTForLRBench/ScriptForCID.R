rm(list = ls())
suppressMessages(library(Seurat))
source('./Script/Step2_PreSTForLRBench/function.R')
set.seed(123)

# CID datasets
img.path <- list.dirs('~/1-Datasets/NatureGenetics_2021_BC/ST_Data/spatial')[-1]
matrix.path <- list.dirs('~/1-Datasets/NatureGenetics_2021_BC/ST_Data/filtered_count_matrices')[-1]
meta.path <- list.dirs('~/1-Datasets/NatureGenetics_2021_BC/ScriptForCCC/ST_deconvoulution', recursive = FALSE)
samples <- substring(meta.path, 76)

for (i in seq(img.path)) {
  output.path <- paste0('./Data/Step2_PreSTForLRBench/', samples[i])
  if(!dir.exists(output.path)){
    dir.create(output.path)
  }
  
  # Save SeuratObject of ST dataset
  if(T){
    count.st <- Read10X(matrix.path[i], gene.column = 1)
    img <- Read10X_Image(image.dir = img.path[i])
    meta.st <- read.csv(paste0(meta.path[i], '/c2l_means_result.csv'))
    meta.st <- tibble::column_to_rownames(meta.st, var = "X")
    colnames(meta.st) <- gsub('meanscell_abundance_w_sf_', '', colnames(meta.st))
    meta.st$celltype <- apply(meta.st, 1, function(x){
      colnames(meta.st)[which.max(x)]
    })
    
    shared.barcode <- intersect(colnames(count.st), rownames(meta.st))
    shared.barcode <- intersect(shared.barcode, rownames(img@coordinates))
    
    meta.st <- meta.st[shared.barcode, ]
    count.st <- count.st[, rownames(meta.st)]
    img <- img[colnames(count.st)]
    DefaultAssay(object = img) <- 'Spatial'
    
    identical(rownames(meta.st), colnames(count.st))
    
    ser <- CreateSeuratObject(count.st, meta.data = meta.st,
                              min.cells = 1, min.features = 1, 
                              assay = 'Spatial')
    ser <- SCTransform(ser, assay = "Spatial")
    ser[['image']] <- img
    
    saveRDS(ser, file = paste0(output.path, '/STser.rds'))
  }
  
  close_distant <- CloDistCP(ser)
  saveRDS(close_distant, file = paste0(output.path, '/CloDistCP.rds'))
}