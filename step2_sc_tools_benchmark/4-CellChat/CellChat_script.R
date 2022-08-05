
rm(list = ls())

library(Seurat)
library(CellChat)
library(tidyverse)
set.seed(123)

load("./step1_data_process/se_data/step2_scRNA_10X.RData")

matrix.sc <- GetAssayData(se.sc, "data", "RNA")
matrix.sc <- as.matrix(matrix.sc)
meta.sc <- data.frame(labels = se.sc$celltype_major, row.names = colnames(se.sc))
identical(colnames(matrix.sc), rownames(meta.sc))

cellchat <- createCellChat(object = matrix.sc, meta = meta.sc, group.by = "labels")
cellchat@DB <- CellChatDB.human

cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multiprocess", workers = 10) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

cellchat.trim <- computeCommunProb(cellchat) %>% filterCommunication(., min.cells = 10) %>% 
  subsetCommunication(.)
cellchat.trun05 <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.05) %>% 
  filterCommunication(., min.cells = 10) %>% subsetCommunication(.)
cellchat.trun10 <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.10) %>%
  filterCommunication(., min.cells = 10) %>% subsetCommunication(.)
cellchat.trun15 <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.15) %>%
  filterCommunication(., min.cells = 10) %>% subsetCommunication(.)
cellchat.trun20 <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.20) %>%
  filterCommunication(., min.cells = 10) %>% subsetCommunication(.)
result.CellChat <- list(trim = cellchat.trim, trun_05 = cellchat.trun05, 
                        trun_10 = cellchat.trun10, trun_15 = cellchat.trun15, 
                        trun_20 = cellchat.trun20)

save(result.CellChat, file = './step2_sc_tools_benchmark/4-CellChat/result.CellChat.RData')
