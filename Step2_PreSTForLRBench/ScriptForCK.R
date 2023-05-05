rm(list = ls())
suppressWarnings(library(Seurat))
source('./Script/Step2_PreSTForLRBench/function.R')
set.seed(123)

samples <- readRDS('~/1-Datasets/Nature_2022_MI/dataset_id.rds')
samples <- samples$ST

for(sample in samples){
  output.path <- paste0('./Data/Step2_PreSTForLRBench/', sample)
  if(!dir.exists(output.path)){
    dir.create(output.path)
  }
  
  ser2 <- readRDS(paste0('~/1-Datasets/Nature_2022_MI/ST_Data/Visium-', sample, '.rds'))
  meta <- ser2@meta.data
  meta <- meta[, 15, drop = FALSE]; colnames(meta) <- 'celltype'
  meta$celltype <- gsub('Cycling.cells', 'Cycling', meta$celltype)
  rm(ser2);gc()
  
  st.path <- paste0('~/1-Datasets/Nature_2022_MI/ST_Data/', sample)
  counts <- Read10X_h5(paste0(st.path, '/filtered_feature_bc_matrix.h5'))
  img <- Read10X_Image(paste0(st.path, '/spatial'))
 
  share.barcode <- intersect(rownames(meta), colnames(counts))
  share.barcode <- intersect(share.barcode, rownames(img@coordinates))
  
  meta <- meta[share.barcode, , drop = FALSE]
  counts <- counts[, rownames(meta)]
  img <- img[colnames(counts)]
  DefaultAssay(object = img) <- 'Spatial'

  identical(rownames(meta), colnames(counts))
  
  ser <- CreateSeuratObject(counts, meta.data = meta,
                            min.cells = 1, min.features = 1, 
                            assay = 'Spatial')
  ser <- SCTransform(ser, assay = "Spatial")
  ser[['image']] <- img
  saveRDS(ser, file = paste0(output.path, '/STser.rds'))
  
  close_distant <- CloDistCP(ser)
  saveRDS(close_distant, file = paste0(output.path, '/CloDistCP.rds'))
}
