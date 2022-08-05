rm(list = ls())

library(Seurat)
set.seed(123)

load("./step1_data_process/se_data/step2_scRNA_10X.RData")

write.csv(t(as.matrix(se.sc@assays$RNA@counts)),
          file = "./step2_sc_tools_benchmark/15-Domino/input_for_SCENIC/count.csv")
