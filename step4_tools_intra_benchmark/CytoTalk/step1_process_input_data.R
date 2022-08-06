rm(list = ls())
library(Seurat)
set.seed(123)

# load ST data
if(T){
  load('./step1_data_process/se_data/ST.se.RData')
  
  # load gene expression matrix of ST data
  se.st <- SCTransform(se.st, assay = "Spatial")
  
  # load ST metadata
  meta <- se.st@meta.data
  
  meta$celltype <- stringr::str_replace_all(meta$celltype, "B_cells", "B-cells")
  meta$celltype <- stringr::str_replace_all(meta$celltype, "T_cells", "T-cells")
  meta$celltype <- stringr::str_replace_all(meta$celltype, "Cancer_Epithelial", "Cancer Epithelial")
  
  se.st <- AddMetaData(se.st, metadata = meta)
  Idents(se.st) <- "celltype"
  
  rm(meta)
}

setwd("./step4_tools_intra_benchmark/CytoTalk/")
norm.matrix <- GetAssayData(se.st, "data", "SCT")
norm.matrix <- as.matrix(norm.matrix)
write.table(norm.matrix, paste0("input/","CID4465", '_norm_count.txt'), sep='\t', quote=F)

meta.data <- cbind(rownames(se.st@meta.data), se.st@meta.data$celltype)  
meta.data <- as.matrix(meta.data)
write.table(meta.data, paste0("input/","CID4465", '_meta.txt'), sep='\t', quote=F, row.names=F)



