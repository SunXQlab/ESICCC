rm(list = ls())
library(SingleCellSignalR)
library(Seurat)
set.seed(123)

load("./step1_data_process/se_data/step2_scRNA_10X.RData")

matrix.sc <- GetAssayData(se.sc, "data", "RNA")
meta.sc <- data.frame(celltype = se.sc$celltype_major, row.names = colnames(se.sc))
i <- 1
for(ct in unique(meta.sc$celltype)){
  meta.sc[which(meta.sc$celltype == ct), "ct_num"] <- i
  i <- i+1
}
c.names <- unique(meta.sc$celltype)
rm(i, ct)

signal <- cell_signaling(data = matrix.sc, genes = rownames(matrix.sc), int.type = "paracrine",
                         species = "homo sapiens", cluster = meta.sc$ct_num, c.names = c.names, write = FALSE)

inter.net <- inter_network(data = matrix.sc, signal = signal, genes = rownames(matrix.sc), 
                           cluster = meta.sc$ct_num, c.names = c.names, write = FALSE)
result.nor.SCSR.para <- inter.net

save(result.nor.SCSR.para, file = "./step2_sc_tools_benchmark/2-SingleCellSignalR/result.SingleCellSignal.para.RData")

