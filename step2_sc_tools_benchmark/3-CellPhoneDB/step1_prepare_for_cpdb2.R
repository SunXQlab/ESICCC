rm(list = ls())

library(Seurat)
set.seed(123)

load("./step1_data_process/se_data/step2_scRNA_10X.RData")
setwd("./step2_sc_tools_benchmark/3-CellPhoneDB/")

sample <- 'CID4465'
norm.matrix <- GetAssayData(se.sc, "data", "RNA")
norm.matrix <- as.matrix(norm.matrix)
write.table(norm.matrix, paste0("input/",sample, '_cellphonedb_count.txt'), sep='\t', quote=F)

meta.data <- cbind(rownames(se.sc@meta.data), se.sc@meta.data$celltype_major)  
meta.data <- as.matrix(meta.data)
write.table(meta.data, paste0("input/",sample, '_cellphonedb_meta.txt'), sep='\t', quote=F, row.names=F)


