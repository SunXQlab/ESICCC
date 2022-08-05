rm(list = ls())
library(Seurat)
library(cellcall)
library(stringr)
set.seed(123)

# run CellCall
if(T){
  load("./step1_data_process/se_data/step2_scRNA_10X.RData")
  
  result <- list()
  matrix.sc <- GetAssayData(se.sc, "counts", "RNA")
  matrix.sc <- as.matrix(matrix.sc)
  colnames(matrix.sc) <- paste(colnames(matrix.sc), str_replace_all(se.sc$celltype_major, "-", " "), sep = "_")
  matrix.sc <- as.data.frame(matrix.sc)
  
  mt <- CreateNichConObject(data=matrix.sc,
                            names.field = 3,
                            names.delim = "_",
                            source = "UMI",
                            scale.factor = 10^6,
                            Org = "Homo sapiens",
                            project = "Microenvironment")
  
  mt <- TransCommuProfile(object = mt,
                          pValueCor = 0.05,
                          CorValue = 0.1,
                          topTargetCor=1,
                          p.adjust = 0.05,
                          use.type="mean",
                          method="weighted",
                          IS_core = TRUE,
                          Org = 'Homo sapiens')
  
  n <- mt@data$expr_l_r_log2_scale
  
  pathway.hyper.list <- lapply(colnames(n), function(i){
    print(i)
    tmp <- getHyperPathway(data = n, object = mt, cella_cellb = i, Org="Homo sapiens")
    return(tmp)
  })
  names(pathway.hyper.list) <- colnames(n)
  
  result[["lr"]] <- n
  result[["pathway"]] <- pathway.hyper.list
  result[["all"]] <- mt
  
  save(result, file = './step2_sc_tools_benchmark/9-CellCall/output/output.CellCall.RData')
}

# handle result
if(T){
  load("./step2_sc_tools_benchmark/9-CellCall/output/output.CellCall.RData")
  
  lr <- result$lr
  lr <- as.data.frame(lr)
  lr <- tibble::rownames_to_column(lr, var = "LR")
  lr <- tidyr::pivot_longer(lr,cols = -LR, names_to = "sr", values_to = "value")
  lr <- lr[which(lr$value > 0), ]
  lr <- tidyr::separate(data = lr, col = sr, into = c("sender", "reciever"), sep = "-")
  lr <- tidyr::separate(data = lr, col = LR, into = c("Ligand", "Receptor"), sep = "-")
  result.CellCall <- lr 
  
  save(result.CellCall, file = "./step2_sc_tools_benchmark/9-CellCall/result.CellCall.RData")
  
}

