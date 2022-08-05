rm(list = ls())

library(iTALK)
library(Seurat)
set.seed(123)

load("./step1_data_process/se_data/step2_scRNA_10X.RData")

result.iTALK <- lapply(c(10, 25, 50), function(top_genes){
  matrix.sc <- GetAssayData(se.sc, "data", "RNA")
  matrix.sc <- as.matrix(matrix.sc)
  matrix.sc <- as.data.frame(t(matrix.sc))
  matrix.sc$cell_type <- se.sc$celltype_major
  
  # find top 50 percent highly expressed genes
  highly_exprs_genes<-rawParse(matrix.sc,top_genes=top_genes,stats='mean')
  
  # find the ligand-receptor pairs from highly expressed genes
  comm.list<-c('growth factor','other','cytokine','checkpoint')
  
  res<-NULL
  for(comm.type in comm.list){
    res.tmp<-FindLR(highly_exprs_genes,datatype='mean count',comm_type=comm.type)
    res.tmp<-res.tmp[order(res.tmp$cell_from_mean_exprs*res.tmp$cell_to_mean_exprs,decreasing=T),]
    res<-rbind(res,res.tmp)
  }
  res
})

names(result.iTALK) <- c(10, 25, 50)

save(result.iTALK, file = './step2_sc_tools_benchmark/5-iTALK/result.iTALK.RData')
