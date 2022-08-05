rm(list = ls())
library(Seurat)
set.seed(123)

load("./step1_data_process/se_data/step2_scRNA_10X.RData")

setwd("./step2_sc_tools_benchmark/19-NATMI/")
id <- 'CID4465'
fpath.mat <- paste0("./input/", id, "_em.csv")
write.csv(100 * (exp(as.matrix(GetAssayData(object = se.sc, assay = "RNA", slot = "data"))) - 1), 
          fpath.mat, row.names = T)
meta <- data.frame(Cell = colnames(se.sc), Annotation = se.sc$celltype_major)
fpath.meta <- paste0("./input/", id, "_metadata.csv")
write.csv(meta,fpath.meta)