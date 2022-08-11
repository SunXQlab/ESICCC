rm(list = ls())
set.seed(123)
library(dplyr)
library(stringr)

# handle result
if(F){
  setwd('./step3_run_lr_benchmark/input/tools_result_for_benchmark/')
  # SingleCellSignalR
  load('./SingleCellSignalR_result.RData')
  result.SingleCellSignalR <- do.call(rbind, result.SingleCellSignalR)
  result.SingleCellSignalR$all <- paste(result.SingleCellSignalR$ligand, result.SingleCellSignalR$receptor,
                                        result.SingleCellSignalR$sr, sep = '_')
  result.SingleCellSignalR$rank <- rank(-result.SingleCellSignalR$LRscore)
  save(result.SingleCellSignalR, file = '../tools_result_for_rsi/SingleCellSignalR.result.RData')
  rm(result.SingleCellSignalR)
  
  # scSeqComm
  load("./scSeqComm_result.RData")
  result.scSeqComm <- result.scSeqComm$`0.8`
  result.scSeqComm <- do.call(rbind, result.scSeqComm)
  result.scSeqComm$all <- paste(result.scSeqComm$ligand, result.scSeqComm$receptor,
                                result.scSeqComm$sr, sep = "_")
  result.scSeqComm$rank <- rank(-result.scSeqComm$S_inter)
  save(result.scSeqComm, file = '../tools_result_for_rsi/scSeqComm.result.RData')
  rm(result.scSeqComm)
  
  # scConnect
  load('./scConnect_result.RData')
  result.scConnect <- do.call(rbind, result.scConnect)
  result.scConnect$all <- paste(result.scConnect$LR, result.scConnect$sr, sep = '_')
  result.scConnect$rank <- rank(-result.scConnect$score)
  save(result.scConnect, file = '../tools_result_for_rsi/scConnect.result.RData')
  rm(result.scConnect)
  
  # NicheNet
  load('./NicheNet_result.RData')
  result.NicheNet <- do.call(rbind, result.NicheNet)
  result.NicheNet$sr <- rownames(result.NicheNet)
  result.NicheNet$sr <- str_replace_all(result.NicheNet$sr, '\\.[0-9]+', '')
  result.NicheNet$all <- paste(result.NicheNet$ligand, result.NicheNet$receptor,
                               result.NicheNet$sr, sep = '_')
  result.NicheNet$rank <- rank(-result.NicheNet$weight)
  save(result.NicheNet, file = '../tools_result_for_rsi/NicheNet.result.RData')
  rm(result.NicheNet)
  
  # NATMI
  load('./NATMI_result.RData')
  result.NATMI <- do.call(rbind, result.NATMI)
  result.NATMI$all <- paste(result.NATMI$ligand, result.NATMI$receptor, 
                            result.NATMI$sr, sep = "_")
  result.NATMI$rank <- rank(-result.NATMI$Edge.average.expression.weight)
  save(result.NATMI, file = '../tools_result_for_rsi/NATMI.result.RData')
  rm(result.NATMI)
  
  # iTALK
  load('./iTALK_result.RData')
  result.iTALK <- result.iTALK$`10`
  result.iTALK <- do.call(rbind, result.iTALK)
  result.iTALK$all <- paste(result.iTALK$ligand, result.iTALK$receptor, 
                            result.iTALK$sr, sep = "_")
  result.iTALK$rank <- rank(-result.iTALK$LRscore)
  save(result.iTALK, file = '../tools_result_for_rsi/iTALK.result.RData')
  rm(result.iTALK)
  
  # ICELLNET
  load('./ICELLNET_result.RData')
  result.ICELLNET <- do.call(rbind, result.ICELLNET)
  result.ICELLNET$all <- paste(result.ICELLNET$ligand, result.ICELLNET$receptor, 
                               result.ICELLNET$sr, sep = '_')
  result.ICELLNET$rank <- rank(-result.ICELLNET$product_value)
  save(result.ICELLNET, file = '../tools_result_for_rsi/ICELLNET.result.RData')
  rm(result.ICELLNET)
  
  # CytoTalk
  load("~/CCC-benchmark/step2_sc_tools_benchmark/7-CytoTalk/result.CytoTalk.RData")
  result_CytoTalk[["CAFs-CAFs"]] <- NULL
  result_CytoTalk <- lapply(result_CytoTalk, function(res){res$pcst$final_network})
  result_CytoTalk <- do.call(rbind, result_CytoTalk)
  result_CytoTalk <- result_CytoTalk[which(result_CytoTalk$node1_type!=result_CytoTalk$node2_type),]
  result_CytoTalk$node1 <- toupper(result_CytoTalk$node1)
  result_CytoTalk$node2 <- toupper(result_CytoTalk$node2)
  result_CytoTalk$node1 <- str_replace_all(result_CytoTalk$node1, "C14ORF1", "C14orf1")
  result_CytoTalk$node2 <- str_replace_all(result_CytoTalk$node2, "C14ORF1", "C14orf1")
  result_CytoTalk$all <- paste(result_CytoTalk$node1, result_CytoTalk$node2, 
                               result_CytoTalk$node1_type, result_CytoTalk$node2_type, sep = '_')
  result_CytoTalk <- dplyr::distinct(result_CytoTalk, all, .keep_all = TRUE)
  
  load('./CytoTalk_result.RData')
  result.CytoTalk <- do.call(rbind, result.CytoTalk)
  result.CytoTalk$all <- paste(result.CytoTalk$ligand, result.CytoTalk$receptor, 
                               result.CytoTalk$sr, sep = '_')
  result.CytoTalk <- merge(result.CytoTalk, result_CytoTalk, by = 'all')
  rm(result_CytoTalk)
  result.CytoTalk$rank <- rank(result.CytoTalk$cost)
  save(result.CytoTalk, file = '../tools_result_for_rsi/CytoTalk.result.RData')
  rm(result.CytoTalk)
  
  # Connectome
  load('./Connectome_result.RData')
  result.Connectome <- do.call(rbind, result.Connectome)
  result.Connectome$all <- paste(result.Connectome$ligand, result.Connectome$receptor,
                                 result.Connectome$sr, sep = '_')
  result.Connectome$rank <- rank(-result.Connectome$weight_norm)
  save(result.Connectome, file = '../tools_result_for_rsi/Connectome.result.RData')
  rm(result.Connectome)
  
  # CellTalker
  load('./CellTalker_result.RData')
  result.CellTalker <- result.CellTalker$`0`
  result.CellTalker <- do.call(rbind, result.CellTalker)
  result.CellTalker$all <- paste(result.CellTalker$ligand, result.CellTalker$receptor,
                                 result.CellTalker$sr, sep = '_')
  result.CellTalker$rank <- rank(-result.CellTalker$value)
  save(result.CellTalker, file = '../tools_result_for_rsi/CellTalker.result.RData')
  rm(result.CellTalker)
  
  # CellPhoneDB3
  load('./CellPhoneDB3_result.RData')
  result.cpdb3 <- do.call(rbind, result.cpdb3)
  result.cpdb3$all <- paste(result.cpdb3$ligand, result.cpdb3$receptor,
                            result.cpdb3$sr, sep = '_')
  result.cpdb3$rank <- rank(-result.cpdb3$mean)
  save(result.cpdb3, file = '../tools_result_for_rsi/CellPhoneDB3.result.RData')
  rm(result.cpdb3)
  
  # CellPhoneDB2
  load('./CellPhoneDB2_result.RData')
  result.cpdb2 <- do.call(rbind, result.cpdb2)
  result.cpdb2$all <- paste(result.cpdb2$ligand, result.cpdb2$receptor, 
                            result.cpdb2$sr, sep = '_')
  result.cpdb2$rank <- rank(-result.cpdb2$mean)
  save(result.cpdb2, file = '../tools_result_for_rsi/CellPhoneDB2.result.RData')
  rm(result.cpdb2)
  
  # Cellinker
  load('./Cellinker_result.RData')
  result.Cellinker <- do.call(rbind, result.Cellinker)
  result.Cellinker$all <- paste(result.Cellinker$ligand, result.Cellinker$receptor,
                                result.Cellinker$sr, sep = '_')
  result.Cellinker$rank <- rank(-result.Cellinker$value)
  save(result.Cellinker, file = '../tools_result_for_rsi/Cellinker.result.RData')
  rm(result.Cellinker)
  
  # CellChat
  load('./CellChat_result.RData')
  result.CellChat <- result.CellChat$trim
  result.CellChat <- do.call(rbind, result.CellChat)
  result.CellChat$all <- paste(result.CellChat$ligand, result.CellChat$receptor, 
                               result.CellChat$sr, sep = '_')
  result.CellChat$rank <- rank(-result.CellChat$prob)
  save(result.CellChat, file = '../tools_result_for_rsi/CellChat.result.RData')
  rm(result.CellChat)
  
  # CellChat
  load('./CellCall_result.RData')
  result.CellCall <- do.call(rbind, result.CellCall)
  result.CellCall$all <- paste(result.CellCall$ligand, result.CellCall$receptor,
                               result.CellCall$sr, sep = '_')
  result.CellCall$rank <- rank(-result.CellCall$value)
  save(result.CellCall, file = '../tools_result_for_rsi/CellCall.result.RData')
  rm(result.CellCall)
}

# RSI
if(T){
  rm(list = ls())
  setwd('./step3_run_lr_benchmark/input/tools_result_for_rsi/')
  files <- list.files('.')
  for (file in files) { load(file) }
  rm(file, files); gc()
  results <- objects(pattern = '^result')
  rsi <- c()
  for (result.1 in results) {
    result_1 <- get(result.1)
    rsi.1 <- c()
    for(result.2 in results){
      result_2 <- get(result.2)
      overlap.lr <- intersect(result_1$all, result_2$all)
      if(length(overlap.lr)!=0){
        rank.1 <- result_1[which(result_1$all %in% overlap.lr),]$rank/dim(result_1)[1]
        rank.2 <- result_2[which(result_2$all %in% overlap.lr), ]$rank/dim(result_2)[1]
        mean.rank <- mean(abs(rank.1-rank.2))
        temp.rsi <- 1-mean.rank
      }else{temp.rsi <- 0}
      rsi.1 <- c(rsi.1, temp.rsi)
    }
    rsi <- cbind(rsi, rsi.1)
  }
  rsi <- as.data.frame(rsi)
  colnames(rsi) <- substring(results, 8)
  rownames(rsi) <- substring(results, 8)
  rsi <- rsi[, c('CytoTalk', 'NATMI', 'Connectome', 'iTALK', 
                 'CellCall', 'ICELLNET', 'scSeqComm', 
                 'SingleCellSignalR', 'CellTalker', 
                 'scConnect', 'cpdb2','cpdb3', 'Cellinker', 
                 'CellChat', 'NicheNet')]
  rsi <- rsi[c('CytoTalk', 'NATMI', 'Connectome', 'iTALK', 
               'CellCall', 'ICELLNET', 'scSeqComm', 
               'SingleCellSignalR', 'CellTalker', 
               'scConnect', 'cpdb2','cpdb3', 'Cellinker', 
               'CellChat', 'NicheNet'),]
  rsi <- as.matrix(rsi)
  
  pheatmap::pheatmap(rsi,#cluster_rows = FALSE, cluster_cols = FALSE,
                     treeheight_row = 0,
                     treeheight_col = 0,
                     display_numbers = TRUE,
                     fontsize=10,
                     #annotation_col = anno,
                     number_format = "%.2f",
                     border_color=NA,
                     angle_col = "45",
                     colorRampPalette(c("#438197","#FFFFFF", "#CF7660"))(50))
}
