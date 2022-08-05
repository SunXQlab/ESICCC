rm(list = ls())

library(Seurat)
set.seed(123)

# gene expresseion matrix and cell type meta were the same as the input of CellPhoneDB v2

load("./step1_data_process/se_data/step2_scRNA_10X.RData")
setwd("./step2_sc_tools_benchmark/18-CellPhoneDB3/")

DEG.result <- list()
Idents(se.sc) <- se.sc$celltype_major
DEGs <- FindAllMarkers(se.sc, 
                       test.use = 't',
                       verbose = F, 
                       only.pos = T, 
                       random.seed = 123, 
                       logfc.threshold = 0.15, 
                       min.pct = 0.05, 
                       return.thresh = 0.05)
fDEGs = subset(DEGs, p_val_adj < 0.05 & avg_log2FC > 0.15)
fDEGs = fDEGs[, c('cluster', 'gene', 'p_val_adj', 'p_val', 'avg_log2FC', 'pct.1', 'pct.2')]
write.table(fDEGs, file = './input/degs.txt', sep = '\t', quote = F, row.names = F)
