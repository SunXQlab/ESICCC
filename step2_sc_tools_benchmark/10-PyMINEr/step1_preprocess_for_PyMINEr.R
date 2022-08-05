rm(list = ls())
library(Seurat)
set.seed(123)

load("./step1_data_process/se_data/step2_scRNA_10X.RData")

sample <- 'CID4465'
meta.sc <- data.frame(barcode = colnames(se.sc), celltype = se.sc$celltype_major)
i <- 0
for(ct in unique(meta.sc$celltype)){
  meta.sc[which(meta.sc$celltype == ct), "ct_num"] <- i
  i <- i+1
}
meta.sc <- meta.sc[,-2]
fpath <- paste0("./step2_sc_tools_benchmark/10-PyMINEr/input_and_output/",sample, "_meta.txt")
write.table(meta.sc, file = fpath,row.names = F, col.names = F, sep = "\t")

matrix.sc <- GetAssayData(se.sc, "data", "RNA")
matrix.sc <- as.data.frame(matrix.sc)
matrix.sc <- tibble::rownames_to_column(matrix.sc, var = "gene")
fpath <- paste0("./step2_sc_tools_benchmark/10-PyMINEr/input_and_output/",sample, "_count.txt")
write.table(matrix.sc, file = fpath, row.names = FALSE, sep = "\t")

