rm(list = ls())

library(Seurat)
library(stMLnet)
library(Giotto)
library(tidyverse)
set.seed(123)

# load ST data
if(F){
  # load gene expression matrix
  load('./step1_data_process/se_data/ST.se.RData')
  se.st <- SCTransform(se.st, assay = "Spatial")
  
  # load cell type metadata
  if(T){
    meta <- se.st@meta.data
    meta$celltype <- stringr::str_replace_all(meta$celltype, "B_cells", "B-cells")
    meta$celltype <- stringr::str_replace_all(meta$celltype, "T_cells", "T-cells")
    meta$celltype <- stringr::str_replace_all(meta$celltype, "Cancer_Epithelial", "Cancer Epithelial")
    meta <- meta[colnames(se.st),]
  }
  identical(rownames(meta), colnames(se.st))
  se.st <- AddMetaData(se.st, metadata = meta)
  Idents(se.st) <- "celltype"
  
  rm(meta)
  save(se.st, file = "./step4_tools_intra_benchmark/stMLnet/input/CID4465_seruat.RData")
}

load("./step4_tools_intra_benchmark/stMLnet/input/CID4465_seruat.RData")
# define highly differentially expressed ligands
if(F){
  ## ligands
  Databases <- readRDS('./step4_tools_intra_benchmark/MistyR/prior_knowledge/Databases.rds')
  ligs_in_db <- Databases$LigRec.DB$source %>% unique()
  ligs_in_db <- intersect(ligs_in_db, rownames(se.st))
  
  clusters <- se.st@active.ident %>% as.character() %>% unique()
  df_markers_ligs <- lapply(clusters, function(cluster){
    
    df <- FindMarkers(se.st, ident.1 = cluster, features = ligs_in_db, only.pos = T, 
                      min.pct = 0.05)
    df$gene <- rownames(df)
    df$ident.1 <- cluster
    df
    
  }) %>% do.call('rbind',.)
  
  Ligs_up_list <- split(df_markers_ligs$gene,df_markers_ligs$ident.1)
  str(Ligs_up_list)
  
  saveRDS(Ligs_up_list, "./step4_tools_intra_benchmark/stMLnet/input/CID4465_Ligs_up.rds")
}

# define highly expressed receptors with mean greater than 0.1 and pct greater than 0.05
if(F){
  data <- se.st@assays$SCT@data
  BarCluTable <- data.frame(Barcode = colnames(se.st), Cluster = se.st$celltype)
  
  ## 参数
  expr.ct <- 0.1
  pct.ct <- 0.05
  
  ## 载入LigRec数据库
  Databases <- readRDS('./step4_tools_intra_benchmark/MistyR/prior_knowledge/Databases.rds')
  recs_in_db <- Databases$LigRec.DB$target %>% unique()
  
  ## 计算mean和pct
  clusters <- BarCluTable$Cluster %>% as.character() %>% unique()
  
  meanExpr_of_LR <- lapply(clusters, function(cluster){
    
    cluster.ids <- BarCluTable$Barcode[BarCluTable$Cluster == cluster]
    source_mean <- rowMeans(data[,cluster.ids])
    names(source_mean) <- rownames(data)
    source_mean
    
  }) %>% do.call('cbind',.) %>% as.data.frame()
  colnames(meanExpr_of_LR) <- clusters
  
  pct_of_LR <- lapply(clusters, function(cluster){
    
    cluster.ids <- BarCluTable$Barcode[BarCluTable$Cluster == cluster]
    dat <- data[,cluster.ids]
    pct <- rowSums(dat>0)/ncol(dat)
    names(pct) <- rownames(data)
    pct
    
  }) %>% do.call('cbind',.) %>% as.data.frame()
  colnames(pct_of_LR) <- clusters
  
  Recs_expr_list <- lapply(clusters, function(cluster){
    
    recs <- rownames(data)[meanExpr_of_LR[,cluster] >= expr.ct & pct_of_LR[,cluster] >= pct.ct]
    intersect(recs, recs_in_db)
    
  })
  names(Recs_expr_list) <- clusters
  str(Recs_expr_list)
  
  saveRDS(Recs_expr_list, "./step4_tools_intra_benchmark/stMLnet/input/CID4465_Recs_expr.rds")
}


