rm(list = ls())

setwd("./step1_data_process/")

library(Seurat)
library(stringr)
library(tibble)
set.seed(123)

## scRNA
if(F){
  file.path <- './dataset/scRNA.Data/CID4465/'
  meta.path <- './dataset/scRNA.Data/CID4465_metadata.csv'
  
  matrix.sc <- Read10X(data.dir = file.path, gene.column = 1)
  meta.sc <- read.csv(file = meta.path, header = TRUE)
  meta.sc <- column_to_rownames(meta.sc, var = 'X')
  meta.sc <- meta.sc[colnames(matrix.sc), ]
  identical(rownames(meta.sc), colnames(matrix.sc))
  
  se.sc <- CreateSeuratObject(counts = matrix.sc, meta.data = meta.sc, 
                              min.cells = 1, min.features = 0)
  se.sc[["percent.mt"]] <- PercentageFeatureSet(se.sc, pattern = "^MT-")
  se.sc <- subset(x = se.sc, subset = nFeature_RNA > 200)
  se.sc <- subset(x = se.sc, subset = nCount_RNA > 250)
  se.sc <- subset(x = se.sc, subset = percent.mt < 20)
  
  se.sc <- NormalizeData(object = se.sc, display.progress = F)
  se.sc <- FindVariableFeatures(object = se.sc, do.plot = F)
  se.sc <- ScaleData(object = se.sc, vars.to.regress = c("nFeature_RNA", "nCount_RNA"), 
                     display.progress = F)
  
  rm(matrix.sc, meta.sc, file.path, meta.path)
  
  save(se.sc, file = "./se_data/step0_scRNA_10X.RData")
}

load("./se_data/step1_se_sc.RData")

## for stereoscope #scRNA #sampling
if(T) {
  id <- 'CID4465'
  count <- GetAssayData(se.sc, slot = "counts", assay = "RNA")
  
  meta <- data.frame(bio_celltype = se.sc$celltype_major, row.names = colnames(se.sc))
  meta <- rownames_to_column(meta, var="cell")
  
  ct.meta <- as.data.frame(table(meta$bio_celltype))
  temp.meta <- NULL
  meta.filter <- NULL
  for (ct in ct.meta$Var1) {
    if(ct.meta[ct.meta$Var1 == ct, 2] < 25){
      temp.meta <- NULL
    }else if (ct.meta[ct.meta$Var1 == ct, 2] > 500){
      cell <- meta[meta$bio_celltype == ct, 1]
      set.seed(123)
      id.rand <- sample(x = cell, size = 500, replace = FALSE)
      temp.meta <- meta[which(meta$cell %in% id.rand), ]
    }else{
      temp.meta <- meta[which(meta$bio_celltype == ct),]
    }
    meta.filter <- rbind(meta.filter, temp.meta)
  }
  write.table(meta.filter,file=paste0(paste0("./input_data_stereoscope/scRNA_sampling/sample_",id), "_singlecell_mta.tsv"), 
              row.names = FALSE,col.names = TRUE,sep="\t")
  
  count <- as.matrix(count)
  count.filter <- count[, meta.filter$cell]
  count.filter <- t(count.filter)
  count.filter <- rownames_to_column(as.data.frame(count.filter), var = "cell")
  
  count.filter <- as.data.frame(count.filter, stringsAsFactors = FALSE)
  
  write.table(count.filter,file=paste0(paste0("./input_data_stereoscope/scRNA_sampling/sample_",id), "_singlecell_cnt.tsv"), 
              row.names = FALSE,col.names = TRUE,sep="\t")
}

## for stereoscope #ST
if(T){
  library(Seurat)
  set.seed(123)
  
  file.path <- paste0("./dataset/ST.Data/filtered_count_matrices/",paste0(id, "_filtered_count_matrix/"))
  matrix.st <- Read10X(data.dir = file.path, gene.column = 1)
  print(paste0("raw genes: ", dim(matrix.st)[1]))
  print(paste0("raw spots: ", dim(matrix.st)[2]))
  
  se.st <- CreateSeuratObject(matrix.st, min.cells = 10, min.features = 500)
  matrix.st.filter <- GetAssayData(se.st, "counts", "RNA")
  print(paste0("filter genes: ", dim(matrix.st.filter)[1]))
  print(paste0("filter spots: ", dim(matrix.st.filter)[2]))
  matrix.st.filter <- t(as.matrix(matrix.st.filter))
  
  write.table(matrix.st.filter,file = paste0(paste0("./input_data_stereoscope/ST/",id), "_spatial_cnt.tsv"),
              row.names = TRUE,col.names = TRUE,sep="\t")
}

