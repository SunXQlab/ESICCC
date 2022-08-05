setwd("./step1_data_process/")
rm(list = ls())
library(Seurat)
library(stringr)
library(dplyr)
library(rowr)
set.seed(123)

file.path <- list.dirs("./dataset/scRNA.Data")[-1]
meta.path <- paste("./dataset/scRNA.Data",list.files("./dataset/scRNA.Data")[2], sep = "/")
sample.id <- str_split(list.files("./dataset/scRNA.Data")[2], "_")[[1]][1]

# Normalize | FindVariableFeatures | ScaleData
if(T){
  matrix.sc <- Read10X(data.dir = file.path, gene.column = 1)
  meta.sc <- read.csv(file = meta.path, header = TRUE)
  meta.sc <- tibble::column_to_rownames(meta.sc, var = 'X')
  identical(colnames(matrix.sc), rownames(meta.sc))
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
  rm(matrix.sc, meta.sc)
}
save(se.sc, file = "./se_data/step0_scRNA_10X.RData")
load("se_data/step0_scRNA_10X.RData")

# RunPCA | JACKSTRAW ANALYSIS | DIMENSIONAL REDUCTION TSNE & UMAP
if(T){
  se <- RunPCA(object = se, npcs = 100)
  
  temp_PC_1 <- 0.0001
  temp_PC_2 <- 0.01
  temp_PC_3 <- 0.05
  
  # JACKSTRAW ANALYSIS
  if(T){
    se <- JackStraw(object = se, dims=100, verbose=T)
    se <- ScoreJackStraw(se, dims=1:100)
    temp_PC <- se@reductions$pca@jackstraw$overall.p.values
    temp_PC <- as.data.frame(temp_PC)
    temp_PC <- arrange(temp_PC, Score)
    
    temp_PCA_df <- NULL
    for(i in 1:3) {
      
      temp_jackstraw_cutoff <- 
        get(paste0("temp_PC_",i))
      
      print(paste0("jackstraw significance cutoff 1 = ",temp_jackstraw_cutoff))
      
      temp_PC_filtered <-
        temp_PC[temp_PC$Score < temp_jackstraw_cutoff,]
      
      temp_PC_filtered <- 
        temp_PC_filtered$PC
      
      print(paste0("Number of PCs used = ", length(temp_PC_filtered)))
      
      temp_PCA_df_subset <- 
        data.frame(PCA = temp_PC_filtered)
      names(temp_PCA_df_subset) <- paste0("PCA_cutoff_",temp_jackstraw_cutoff)
      
      temp_PCA_df <- 
        cbind.fill(temp_PCA_df,
                   temp_PCA_df_subset,
                   fill=NA)
      
      n <- paste0("temp_PC_",i)
      assign(n, temp_PC_filtered)
    }
    temp_PCA_df <- temp_PCA_df[,!colnames(temp_PCA_df) %in% "init", drop=F]
  }
  
  # TSNE
  for(i in 1:3) {
    temp_PC <- get(paste0("temp_PC_",i))
    se <- RunTSNE(object = se, dims = temp_PC,
                  reduction.key = paste0("TSNE",i,"_"),
                  reduction.name = paste0("TSNE",i))
  }
  
  # UMAP
  for(i in 1:3) {
    temp_PC <- get(paste0("temp_PC_",i))
    se <- RunUMAP(se, dims = temp_PC,
                  reduction.key = paste0("UMAP",i,"_"),
                  reduction.name = paste0("UMAP",i), verbose = F)
  }
  
  se.sc <- se
}
save(se.sc, file = "./se_data/step1_scRNA_10X.RData")
load("./se_data/step1_scRNA_10X.RData")

#filter the number of cells in celltype less than 25
ct.stat <- as.data.frame(table(se.sc$celltype_major))
ct.remove <- as.character(ct.stat[which(ct.stat$Freq < 25), 1])
Idents(se.sc) <- "celltype_major"
se.sc <- subset(se.sc, idents = ct.remove, invert = TRUE)

save(se.sc, file = "./se_data/step2_scRNA_10X.RData")
