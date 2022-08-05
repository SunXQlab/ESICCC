rm(list = ls())
library(Seurat)
set.seed(123)

load("./step1_data_process/se_data/step2_scRNA_10X.RData")
setwd('./step2_sc_tools_benchmark/12-scConnect/')

norm.matrix <- GetAssayData(se.sc, "data", "RNA")
norm.matrix <- as.matrix(norm.matrix)
write.csv(norm.matrix, paste0("input/",sample, '_count.csv'), quote=F)
  
cell.meta <- data.frame(Cell = rownames(se.sc@meta.data), Annotation = se.sc$celltype_major)
write.csv(cell.meta, paste0("input/",sample, '_meta.csv'), quote=F, row.names = FALSE)