rm(list = ls())

library(Seurat)
set.seed(123)

load('./step1_data_process/se_data/step2_scRNA_10X.RData')
raw.count <- GetAssayData(se.sc, 'counts', 'RNA')
raw.count <- as.matrix(raw.count)
write.table(raw.count, file = './step2_sc_tools_benchmark/17-Cellinker/count.txt', quote = FALSE, sep = '\t')

meta <- data.frame(rownames(se.sc@meta.data), se.sc$celltype_major)
write.table(meta, file = './step2_sc_tools_benchmark/17-Cellinker/meta.txt', quote = FALSE, sep = '\t', 
            row.names = FALSE, col.names = FALSE)
