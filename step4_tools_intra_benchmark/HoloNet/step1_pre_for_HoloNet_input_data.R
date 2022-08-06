rm(list = ls())

library(Seurat)
library(tidyverse)
setwd(123)

# load ST data
# load gene expression matrix
load('./step1_data_process/se_data/ST.se.RData')
raw.count <- GetAssayData(se.st, "count", "Spatial")
raw.count <- as.matrix(raw.count)
write.csv(raw.count, "./step4_tools_intra_benchmark/HoloNet/input/CID4465_raw_count.csv", quote=F)

# load cell type metadata
if(T){
  meta <- se.st@meta.data
  meta$celltype <- stringr::str_replace_all(meta$celltype,"B_cells", "B-cells")
  meta$celltype <- stringr::str_replace_all(meta$celltype, "T_cells", "T-cells")
  meta$celltype <- stringr::str_replace_all(meta$celltype, "Cancer_Epithelial", "Cancer Epithelial")
}

meta <- data.frame(Barcodes = rownames(meta), 
                   celltype = meta$celltype, row.names = rownames(meta))
write.csv(meta, "./step4_tools_intra_benchmark/HoloNet/input/CID4465_meta.csv",quote=F, row.names = FALSE)


# load genes of interest
icgs <- readRDS("./step1_data_process/result/Giotto_result/CID4465_bc_ICGs.rds")
icgs <- icgs[["Cancer Epithelial"]]
goi <- c()
for (icg in icgs) {
  goi <- c(goi, as.character(icg))
}

goi <- unique(goi)
write.csv(goi, file = "./step4_tools_intra_benchmark/HoloNet/input/goi.csv", quote = FALSE, row.names = FALSE)
