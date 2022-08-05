rm(list = ls())
set.seed(123)
setwd('./step2_sc_tools_benchmark/')

source("../step3_run_lr_benchmark/code/bench.ST.cal.R")
load('../step3_run_lr_benchmark/input/norm_cd.RData')
close_distant[["0.5"]] <- NULL

# 2-SingleCellSignalR
if(T){
  # load and handle the result of SingleCellSignalR in CID4465
  if(F){
    load("./2-SingleCellSignalR/result.SingleCellSignal.para.RData")
    result.SingleCellSignalR <- result.nor.SCSR.para$`individual-networks`
    rm(result.nor.SCSR.para)
    
    result.SingleCellSignalR <- lapply(result.SingleCellSignalR, function(result){
      result <- result[,c(3:6,8)]
      colnames(result) <- c("ligand", "receptor", "sender", "reciever", "LRscore")
      result$sr <- paste(result$sender, result$reciever, sep = '_')
      result <- dplyr::distinct(result, ligand, receptor, .keep_all = T)
      result
    })
   
    names(result.SingleCellSignalR) <- stringr::str_replace_all(names(result.SingleCellSignalR), "-","_")
    names(result.SingleCellSignalR) <- stringr::str_replace_all(names(result.SingleCellSignalR), "T_cells","T-cells")
    names(result.SingleCellSignalR) <- stringr::str_replace_all(names(result.SingleCellSignalR), "B_cells","B-cells")
    
    save(result.SingleCellSignalR, file = "../step3_run_lr_benchmark/input/tools_result_for_benchmark/SingleCellSignalR_result.RData")
  }
  load("../step3_run_lr_benchmark/input/tools_result_for_benchmark/SingleCellSignalR_result.RData")
  
  # benchmark——mutual information | pcc
  bench.result.2 <- lapply(close_distant, function(cd){
    bench.result <- bench.ST(cd, norm.data, result.SingleCellSignalR)
    bench.result
  })
  
  save(bench.result.2, file = "../step3_run_lr_benchmark/result/2_SingleCellSignalR.RData")
  rm(result.SingleCellSignalR)
}

# 3-CellPhoneDB
if(T){
  # load and handle the result of cellphonedb in CID4465
  if(F){
    # load the result
    load("./3-CellPhoneDB/result.cpdb2.RData")
    result.cpdb2 <- tidyr::separate(result.cpdb2, "interacting_pair", 
                                   into = c("ligand", "receptor"), sep = "_")
    
    # handle the complex information
    if(F){
      cpdb.complex <-read.csv("./15-Domino/input_for_Domino/cpdb/complexes.csv")
      cpdb.complex <- cpdb.complex[,1:5]
      cpdb.gene <- read.csv("./15-Domino/input_for_Domino/cpdb/genes.csv")
      cpdb.gene <- dplyr::distinct(cpdb.gene, gene_name, uniprot, hgnc_symbol, .keep_all = TRUE)
      cpdb.gene <- cpdb.gene[,1:2]
      
      tmp.cpdb.complex <- merge(cpdb.gene, cpdb.complex, by.y = "uniprot_1", by.x = "uniprot")
      tmp.cpdb.complex <- tmp.cpdb.complex[,-1]
      colnames(tmp.cpdb.complex)[1] <- "gene_1"
      tmp.cpdb.complex <- merge(cpdb.gene, tmp.cpdb.complex, by.y = "uniprot_2", by.x = "uniprot")
      tmp.cpdb.complex <- tmp.cpdb.complex[,-1]
      colnames(tmp.cpdb.complex)[1] <- "gene_2"
      tmp.cpdb.complex <- tmp.cpdb.complex[,-5]
      tmp.cpdb.complex <- merge(cpdb.gene, tmp.cpdb.complex, by.y = "uniprot_3", by.x = "uniprot", all.y = TRUE)
      tmp.cpdb.complex <- tmp.cpdb.complex[,-1]
      colnames(tmp.cpdb.complex)[1] <- "gene_3"
      cpdb.complex <- tmp.cpdb.complex[,4:1]
      rm(cpdb.gene, tmp.cpdb.complex)
      
      cpdb.complex$gene <- paste(cpdb.complex$gene_1, cpdb.complex$gene_2, cpdb.complex$gene_3, sep = "&")
      cpdb.complex$gene <- stringr::str_replace_all(cpdb.complex$gene, "&NA", "")
      cpdb.complex <- cpdb.complex[,c("complex_name", "gene")]
    }
    
    # combine the complex information and result of cellphonedb
    if(F){
      tmp.result <- merge(result.cpdb2, cpdb.complex, by.x = "ligand", 
                          by.y = "complex_name", all.x = TRUE)
      tmp.result$ligand[which(grepl("complex",tmp.result$ligand))] <- 
        tmp.result$gene[which(grepl("complex",tmp.result$ligand))]
      tmp.result <- tmp.result[,-6]
      tmp.result <- merge(tmp.result, cpdb.complex, by.x = "receptor", 
                          by.y = "complex_name", all.x = TRUE)
      tmp.result$receptor[which(grepl("complex",tmp.result$receptor))] <- 
        tmp.result$gene[which(grepl("complex",tmp.result$receptor))]
      tmp.result <- tmp.result[,-6]
      result.cpdb2 <- tmp.result; rm(tmp.result)

      result.cpdb2$sr <- paste(result.cpdb2$sender, result.cpdb2$reciever, sep = "_")
      result.cpdb2 <- result.cpdb2[-which(result.cpdb2$sender==result.cpdb2$reciever),]
      result.cpdb2 <- split(result.cpdb2, result.cpdb2$sr)
      
      rm(cpdb.complex)
    }
    save(result.cpdb2, file = "../step3_run_lr_benchmark/input/tools_result_for_benchmark/CellPhoneDB2_result.RData")
  }
  load("../step3_run_lr_benchmark/input/tools_result_for_benchmark/CellPhoneDB2_result.RData")
  
  # benchmark——mutual information | pcc
  bench.result.3 <- lapply(close_distant, function(cd){
    bench.result <- bench.ST(cd, norm.data, result.cpdb2)
    bench.result
  })
  
  save(bench.result.3, file = "../step3_run_lr_benchmark/result/3_CellPhoneDB2.RData")
}

# 4-CellChat
if(T){
  # handle the result
  if(F){
    load('./4-CellChat/result.CellChat.RData')
    result.CellChat <- lapply(result.CellChat, function(res){
      res <- res[,c(1:6)]
      colnames(res)[1:4] <- c("sender", "reciever", "ligand", "receptor")
      res$receptor <- stringr::str_replace_all(res$receptor, "_", "&")
      res$sr <- paste(res$sender, res$reciever, sep = "_")
      res <- res[which(res$sender != res$reciever), ]
      res <- split(res, res$sr)
      res
    })
    save(result.CellChat, file = '../step3_run_lr_benchmark/input/tools_result_for_benchmark/CellChat_result.RData')
  }
  load("../step3_run_lr_benchmark/input/tools_result_for_benchmark/CellChat_result.RData")
  
  bench.result.4 <- lapply(close_distant, function(cd){
    bench.result <- lapply(result.CellChat, function(res){
      bench.result.para <- bench.ST(cd, norm.data, res)
      bench.result.para
    })
    bench.result
  })
  
  save(bench.result.4, file = "../step3_run_lr_benchmark/result/4_CellChat.RData")
}

# 5-iTALK
if(T){
  # load and handle the result of iTALK in CID4465
  if(F){
    load("./5-iTALK/result.iTALK.RData")
    result.iTALK <- lapply(result.iTALK, function(res){
      colnames(res)[c(4,6)] <- c("sender", "reciever")
      res$sr <- paste(res$sender, res$reciever, sep = "_")
      res <- res[which(res$sender != res$reciever), ]
      
      res$LRscore <- res$cell_from_mean_exprs*res$cell_to_mean_exprs
      res <- res[which(res$LRscore!=0),]
      
      res <- split(res, res$sr)
      res
    })
    
    save(result.iTALK,
         file = "../step3_run_lr_benchmark/input/tools_result_for_benchmark/iTALK_result.RData")
  }
  load("../step3_run_lr_benchmark/input/tools_result_for_benchmark/iTALK_result.RData")
  
  # benchmark——mutual information | pcc
  bench.result.5 <- lapply(close_distant, function(cd){
    bench.result <- lapply(result.iTALK, function(res){
      bench.result.tmp <- bench.ST(cd, norm.data, res)
      bench.result.tmp
    })
    bench.result
  })
  
  save(bench.result.5,
       file = "../step3_run_lr_benchmark/result/5_iTALK.RData")
}

# 7-CytoTalk
if(T){
  # load and handle the result of CytoTalk in CID4465
  if(F){
    load("./7-CytoTalk/result.CytoTalk.RData")
    names(result_CytoTalk) <- stringr::str_replace_all(names(result_CytoTalk), "-", "_")
    names(result_CytoTalk) <- stringr::str_replace_all(names(result_CytoTalk), "T_cells", "T-cells")
    names(result_CytoTalk) <- stringr::str_replace_all(names(result_CytoTalk), "B_cells", "B-cells")
    result_CytoTalk[["CAFs_CAFs"]] <- NULL
    
    tmp.result.CytoTalk <- lapply(result_CytoTalk, function(res){
      if(!is.null(res$pathways)){
        res <- res$pathways$df_pval
        res$ligand <- toupper(res$ligand)
        res$receptor <- toupper(res$receptor)
        res$ligand[which(res$ligand == "C5ORF55")] <- "C5orf55"
        res$ligand[which(res$ligand == "C5ORF55")] <- "C5orf55"
        
        res$receptor[which(res$receptor == "C14ORF1")] <- "C14orf1"
        res$receptor[which(res$receptor == "C14ORF1")] <- "C14orf1"
        colnames(res)[3:4] <- c("sender", "reciever")
      }else{
        res <- NULL
      }
      res
    })
    result.CytoTalk <- tmp.result.CytoTalk[!sapply(tmp.result.CytoTalk,is.null)]
    rm(result_CytoTalk, tmp.result.CytoTalk)
    
    result.CytoTalk <- do.call(rbind, result.CytoTalk)
    result.CytoTalk <- dplyr::distinct(result.CytoTalk, ligand, receptor, sender, 
                                       reciever, .keep_all = T)
    result.CytoTalk$sr <- paste(result.CytoTalk$sender, result.CytoTalk$reciever, sep = "_")
    result.CytoTalk <- result.CytoTalk[which(result.CytoTalk$sender != result.CytoTalk$reciever), ]
    result.CytoTalk <- split(result.CytoTalk, result.CytoTalk$sr)
    save(result.CytoTalk, file = "../step3_run_lr_benchmark/input/tools_result_for_benchmark/CytoTalk_result.RData")
  }
  load("../step3_run_lr_benchmark/input/tools_result_for_benchmark/CytoTalk_result.RData")
  
  bench.result.7 <- lapply(close_distant, function(cd){
    bench.result <- bench.ST(cd, norm.data, result.CytoTalk)
    bench.result
  })
  save(bench.result.7, file = "../step3_run_lr_benchmark/result/7_CytoTalk.RData")
}

# 8-scMLnet
if(T){
  # load and handle the result of scMLnet in CID4465
  if(F){
    load("./8-scMLnet/result.scMLnet.RData")
    result.scMLnet.2[unlist(lapply(result.scMLnet.2, function(res){all(is.na(res))})) ] <- NULL
    result.scMLnet <- lapply(result.scMLnet.2, function(res){
      res <- data.frame(lr = res$LigRec)
      colnames(res)[1:2] <- c("ligand", "receptor")
      res
    })
    
    save(result.scMLnet, file = "../step3_run_lr_benchmark/input/tools_result_for_benchmark/scMLnet_result.RData")
    
  }
  load("../step3_run_lr_benchmark/input/tools_result_for_benchmark/scMLnet_result.RData")
  
  bench.result.8 <- lapply(close_distant, function(cd){
    bench.result <- bench.ST(cd, norm.data, result.scMLnet)
    bench.result
  })

  save(bench.result.8,
       file = "../step3_run_lr_benchmark/result/8_scMLnet.RData")
}

# 9-CellCall
if(T){
  # load and handle the result of CellCall in CID4465
  if(F){
    load("./9-CellCall/result.CellCall.RData")
    result.CellCall$sender <- stringr::str_replace_all(result.CellCall$sender, "T cells", "T-cells")
    result.CellCall$sender <- stringr::str_replace_all(result.CellCall$sender, "B cells", "B-cells")
    result.CellCall$reciever <- stringr::str_replace_all(result.CellCall$reciever, "T cells", "T-cells")
    result.CellCall$reciever <- stringr::str_replace_all(result.CellCall$reciever, "B cells", "B-cells")
    result.CellCall$sr <- paste(result.CellCall$sender, result.CellCall$reciever, sep = "_")
    colnames(result.CellCall)[1:2] <- c("ligand", "receptor")
    result.CellCall$ligand <- stringr::str_replace_all(result.CellCall$ligand, ",", "&")
    result.CellCall$receptor <- stringr::str_replace_all(result.CellCall$receptor, ",", "&")
    result.CellCall <- result.CellCall[which(result.CellCall$sender!=result.CellCall$reciever),]
    
    result.CellCall <- split(result.CellCall, result.CellCall$sr)
    save(result.CellCall, file = "../step3_run_lr_benchmark/input/tools_result_for_benchmark/CellCall_result.RData")
  }
  load("../step3_run_lr_benchmark/input/tools_result_for_benchmark/CellCall_result.RData")
  
  bench.result.9 <- lapply(close_distant, function(cd){
    bench.result <- bench.ST(cd, norm.data, result.CellCall)
    bench.result
  })
  save(bench.result.9, file = "../step3_run_lr_benchmark/result/9_CellCall.RData")
}

# 10-PyMINEr
if(T){
  # load and handle the result of PyMINEr in CID4465
  if(F){
    load("./10-PyMINEr/result.PyMINEr.RData")
    result.PyMINEr <- result.PyMINEr[,c(2,4,5,6,8,9)]
    colnames(result.PyMINEr)[c(1:2,4:5)] <- c("sender", "ligand", "reciever", "receptor")
    result.PyMINEr$sr <- paste(result.PyMINEr$sender, result.PyMINEr$reciever, sep = "_")
    result.PyMINEr <- result.PyMINEr[which(result.PyMINEr$sender!=result.PyMINEr$reciever),]
    result.PyMINEr <- split(result.PyMINEr, result.PyMINEr$sr)
    save(result.PyMINEr, 
         file = "../step3_run_lr_benchmark/input/tools_result_for_benchmark/PyMINEr_result.RData")
  }
  load("../step3_run_lr_benchmark/input/tools_result_for_benchmark/PyMINEr_result.RData")
  
  bench.result.10 <- lapply(close_distant, function(cd){
    bench.result <- bench.ST(cd, norm.data, result.PyMINEr)
    bench.result
  })
  save(bench.result.10, file = "../step3_run_lr_benchmark/result/10_PyMINEr.RData")
}

# 13-scSeqComm
if(T){
  # load and handle the result of scSeqComm in CID4465
  if(F){
    load("./13-scSeqComm/result.scSeqComm.RData") 
    result.scSeqComm <- lapply(result.scSeqComm, function(res){
      res$receptor <- stringr::str_replace_all(res$receptor, ",", "&")
      res$ligand <- stringr::str_replace_all(res$ligand, ",", "&")
      res <- dplyr::distinct(res, ligand, receptor, cluster_L, cluster_R, .keep_all = T)
      res$sr <- paste(res$cluster_L, res$cluster_R, sep = "_")
      res <- res[which(res$cluster_L != res$cluster_R), ]
      res <- split(res, res$sr)
    })
    
    save(result.scSeqComm, 
         file = "../step3_run_lr_benchmark/input/tools_result_for_benchmark/scSeqComm_result.RData")
  }
  load("../step3_run_lr_benchmark/input/tools_result_for_benchmark/scSeqComm_result.RData")
  
  bench.result.13 <- lapply(close_distant, function(cd){
    bench.result <- lapply(result.scSeqComm, function(res){
      bench.result.tmp <- bench.ST(cd, norm.data, res)
      bench.result.tmp
    })
    bench.result
  })
  save(bench.result.13, 
       file = "../step3_run_lr_benchmark/result/13_scSeqComm.RData")
}

# 12-scConnect
if(T){
  # load and handle the result of scConnect in CID4465
  if(F){
    load("./12-scConnect/result/result.scConnect.RData")
    result.scConnect <- result.scConnect$filter
    colnames(result.scConnect)[c(8,11)] <- c("ligand", "receptor")
    result.scConnect$sr <- paste(result.scConnect$sender, result.scConnect$reciever, sep = "_")
    result.scConnect <- dplyr::distinct(result.scConnect, ligand, receptor, sr, .keep_all = T)
    result.scConnect <- result.scConnect[which(result.scConnect$sender != result.scConnect$reciever),]
    result.scConnect <- split(result.scConnect, result.scConnect$sr)
    
    save(result.scConnect, 
         file = "../step3_run_lr_benchmark/input/tools_result_for_benchmark/scConnect_result.RData")
  }
  load("../step3_run_lr_benchmark/input/tools_result_for_benchmark/scConnect_result.RData")
  
  bench.result.12 <- lapply(close_distant, function(cd){
    bench.result <- bench.ST(cd, norm.data, result.scConnect)
    bench.result
  })
  save(bench.result.12, file = "../step3_run_lr_benchmark/result/12_scConnect.RData")
}

# 14-Connectome
if(T){
  # load and handle the result of Connectome in CID4465
  if(F){
    load('./14-Connectome/result.Connectome.RData')
    result.Connectome <- result.Connectome[,c(1:4,9:18, 22:23)]
    colnames(result.Connectome)[1:4] <- c("sender", "reciever", "ligand", "receptor")
    result.Connectome$sr <- paste(result.Connectome$sender, result.Connectome$reciever, sep = "_")
    result.Connectome <- result.Connectome[which(result.Connectome$sender!=result.Connectome$reciever),]
    
    result.Connectome <- split(result.Connectome, result.Connectome$sr)
    save(result.Connectome, 
         file = "../step3_run_lr_benchmark/input/tools_result_for_benchmark/Connectome_result.RData")
  }
  load("../step3_run_lr_benchmark/input/tools_result_for_benchmark/Connectome_result.RData")
  
  bench.result.14 <- lapply(close_distant, function(cd){
    bench.result <- bench.ST(cd, norm.data, result.Connectome)
    bench.result
  })
  save(bench.result.14, 
       file = "../step3_run_lr_benchmark/result/14_Connectome.RData")
}

# 15-Domino
if(T){
  # load and handle the result of Domino in CID4465
  if(F){
    load("./15-Domino/result.Domino.RData")
    result.Domino <- result.Domino[which(result.Domino$sender!= result.Domino$reciever), ]
    result.Domino <- dplyr::distinct(result.Domino, ligand, receptor, sr, .keep_all = TRUE)
    result.Domino <- split(result.Domino, result.Domino$sr)
    
    save(result.Domino, 
         file = "../step3_run_lr_benchmark/input/tools_result_for_benchmark/Domino_result.RData")
  }
  load("../step3_run_lr_benchmark/input/tools_result_for_benchmark/Domino_result.RData")
  
  bench.result.15 <- lapply(close_distant, function(cd){
    bench.result <- bench.ST(cd, norm.data, result.Domino)
    bench.result
  })
  save(bench.result.15, 
       file = "../step3_run_lr_benchmark/result/15_Domino.RData")
}

# 16-CellTalker
if(T){
  # load and handle the result of CellTalker
  if(F){
    library(tidyr)
    load("./16-CellTalker/result.CellTalker.RData")
    result.CellTalker <- lapply(result.celltalker, function(res){
      res <- res$filter.fdr
      res <- tidyr::separate(res, col = interaction, into = c("ligand", "receptor"),sep = "_")
      res <- tidyr::separate(res, col = interaction_pairs, into = c("sender", "reciever"),sep = "_")
      res <- res[which(res$sender != res$reciever), ]
      res$sr <- paste(res$sender, res$reciever, sep = '_')
      res <- split(res, res$sr)
      res
      })
    rm(result.celltalker)
    save(result.CellTalker, 
         file = "../step3_run_lr_benchmark/input/tools_result_for_benchmark/CellTalker_result.RData")
  }
  load("../step3_run_lr_benchmark/input/tools_result_for_benchmark/CellTalker_result.RData")
  
  bench.result.16 <- lapply(close_distant, function(cd){
    bench.result <- lapply(result.CellTalker, function(res){
      bench.result.tmp <- bench.ST(cd, norm.data, res)
      bench.result.tmp
    })
    bench.result
  })
  save(bench.result.16, 
       file = "../step3_run_lr_benchmark/result/16_CellTalker.RData")
}

# 17-Cellinker
if(T){
  # load and handle the result of Cellinker in CID4465
  if(F){
    load("./17-Cellinker/result.Cellinker.RData")
    result.Cellinker$sr <- stringr::str_replace_all(result.Cellinker$sr, "T.cells", "T-cells")
    result.Cellinker$sr <- stringr::str_replace_all(result.Cellinker$sr, "B.cells", "B-cells")
    result.Cellinker$sr <- stringr::str_replace_all(result.Cellinker$sr, "Cancer.Epithelial", "Cancer Epithelial")
    result.Cellinker$sr <- stringr::str_replace_all(result.Cellinker$sr, "\\.", "_")
    result.Cellinker$LR <- stringr::str_replace_all(result.Cellinker$LR, "-", "_")
    result.Cellinker$LR <- stringr::str_replace_all(result.Cellinker$LR, "HLA_A", "HLA-A")
    result.Cellinker$LR <- stringr::str_replace_all(result.Cellinker$LR, "HLA_C", "HLA-C")
    result.Cellinker$LR <- stringr::str_replace_all(result.Cellinker$LR, "HLA_DRA", "HLA-DRA")
    result.Cellinker$LR <- stringr::str_replace_all(result.Cellinker$LR, "HLA_DRB1", "HLA-DRB1")
    result.Cellinker <- tidyr::separate(result.Cellinker, LR, sep = "_", into = c("ligand", "receptor"))
    result.Cellinker <- tidyr::separate(result.Cellinker, sr, sep = "_", into = c("sender", "reciever"),remove = FALSE)
    result.Cellinker <- result.Cellinker[which(result.Cellinker$sender!=result.Cellinker$reciever),]
    
    result.Cellinker <- split(result.Cellinker, result.Cellinker$sr)
    save(result.Cellinker, 
         file = "../step3_run_lr_benchmark/input/tools_result_for_benchmark/Cellinker_result.RData")
  }
  load("../step3_run_lr_benchmark/input/tools_result_for_benchmark/Cellinker_result.RData")
  
  bench.result.17 <- lapply(close_distant, function(cd){
    bench.result <- bench.ST(cd, norm.data, result.Cellinker)
    bench.result
  })
  save(bench.result.17, file = "../step3_run_lr_benchmark/result/17_Cellinker.RData")
}

# 18-CellPhoneDB3
if(T){
  # load and handle the result of CellPhoneDB3 in CID4465 | using the DB3
  if(F){
    load('./18-CellPhoneDB3/result.cpdb3.RData')
    result.cpdb <- separate(result.cpdb, interacting_pair, c('ligand', 'receptor'), '_')
    result.cpdb3 <- result.cpdb; rm(result.cpdb)
    # handle the complex information
    if(F){
      cpdb.complex <-read.csv("./15-Domino/input_for_Domino/cpdb/complexes.csv")
      cpdb.complex <- cpdb.complex[,1:5]
      cpdb.gene <- read.csv("./15-Domino/input_for_Domino/cpdb/genes.csv")
      cpdb.gene <- dplyr::distinct(cpdb.gene, gene_name, uniprot, .keep_all = FALSE)
      
      tmp.cpdb.complex <- merge(cpdb.gene, cpdb.complex, by.y = "uniprot_1", by.x = "uniprot")
      tmp.cpdb.complex <- tmp.cpdb.complex[,-1]
      colnames(tmp.cpdb.complex)[1] <- "gene_1"
      tmp.cpdb.complex <- merge(cpdb.gene, tmp.cpdb.complex, by.y = "uniprot_2", by.x = "uniprot")
      tmp.cpdb.complex <- tmp.cpdb.complex[,-1]
      colnames(tmp.cpdb.complex)[1] <- "gene_2"
      tmp.cpdb.complex <- merge(cpdb.gene, tmp.cpdb.complex, by.y = "uniprot_3", by.x = "uniprot", all.y = TRUE)
      tmp.cpdb.complex <- tmp.cpdb.complex[,-1]
      colnames(tmp.cpdb.complex)[1] <- "gene_3"
      tmp.cpdb.complex <- merge(cpdb.gene, tmp.cpdb.complex, by.y = "uniprot_4", by.x = "uniprot", all.y = TRUE)
      tmp.cpdb.complex <- tmp.cpdb.complex[,-1]
      colnames(tmp.cpdb.complex)[1] <- "gene_4"
      cpdb.complex <- tmp.cpdb.complex[,5:1]
      rm(cpdb.gene, tmp.cpdb.complex)
      
      cpdb.complex$gene <- paste(cpdb.complex$gene_1, cpdb.complex$gene_2, 
                                 cpdb.complex$gene_3, cpdb.complex$gene_4, sep = "&")
      cpdb.complex$gene <- stringr::str_replace_all(cpdb.complex$gene, "&NA", "")
      cpdb.complex <- cpdb.complex[,c("complex_name", "gene")]
    }
    
    # combine the complex information and result of cellphonedb
    if(F){
      tmp.result <- merge(result.cpdb3, cpdb.complex, by.x = "ligand", 
                          by.y = "complex_name", all.x = TRUE)
      tmp.result$ligand[which(grepl("complex",tmp.result$ligand))] <- 
        tmp.result$gene[which(grepl("complex",tmp.result$ligand))]
      tmp.result <- tmp.result[,-6]
      tmp.result <- merge(tmp.result, cpdb.complex, by.x = "receptor", 
                          by.y = "complex_name", all.x = TRUE)
      tmp.result$receptor[which(grepl("complex",tmp.result$receptor))] <- 
        tmp.result$gene[which(grepl("complex",tmp.result$receptor))]
      tmp.result <- tmp.result[,-6]
      result.cpdb3 <- tmp.result; rm(tmp.result)
      result.cpdb3 <- result.cpdb3[which(result.cpdb3$sender != result.cpdb3$reciever),]
      result.cpdb3$sr <- paste(result.cpdb3$sender, result.cpdb3$reciever, sep = "_")
      
      result.cpdb3 <- split(result.cpdb3, result.cpdb3$sr)
      rm(cpdb.complex)
    }
    save(result.cpdb3, 
         file = "../step3_run_lr_benchmark/input/tools_result_for_benchmark/CellPhoneDB3_result.RData")
  }
  load("../step3_run_lr_benchmark/input/tools_result_for_benchmark/CellPhoneDB3_result.RData")
  
  bench.result.18 <- lapply(close_distant, function(cd){
    bench.result <- bench.ST(cd, norm.data, result.cpdb3)
    bench.result
  })
  save(bench.result.18, file = "../step3_run_lr_benchmark/result/18_CellPhoneDB3.RData")
}

# 19-NATMI
if(T){
  # load and handle the result of NATMI in CID4465
  if(F){
    load('./19-NATMI/result.NATMI.RData')
    result.NATMI <- result.NATMI[which(result.NATMI$sender!= result.NATMI$reciever), ]
    result.NATMI <- split(result.NATMI, result.NATMI$sr)
    save(result.NATMI,
         file = "../step3_run_lr_benchmark/input/tools_result_for_benchmark/NATMI_result.RData")
  }
  load("../step3_run_lr_benchmark/input/tools_result_for_benchmark/NATMI_result.RData")
  
  bench.result.19 <- lapply(close_distant, function(cd){
    bench.result <- bench.ST(cd, norm.data, result.NATMI)
    bench.result
  })
  save(bench.result.19, file = "../step3_run_lr_benchmark/result/19_NATMI.RData")
}

# 20-iCELLNET
if(T){
  # load and handle the result of iCELLNET in CID4465
  if(F){
    load("./20-ICELLNET/result.ICELLNET.RData")
    
    result.ICELLNET$Ligand <- stringr::str_replace_all(result.ICELLNET$Ligand, " ", "")
    result.ICELLNET$Receptor <- stringr::str_replace_all(result.ICELLNET$Receptor, " ", "")
    result.ICELLNET$Ligand <- stringr::str_replace_all(result.ICELLNET$Ligand, "\\+", "&")
    result.ICELLNET$Receptor <- stringr::str_replace_all(result.ICELLNET$Receptor, "\\+", "&")
    colnames(result.ICELLNET)[1:2] <- c("ligand", "receptor")
    result.ICELLNET$sr <- paste(result.ICELLNET$sender, result.ICELLNET$reciever, sep = "_")
    result.ICELLNET <- result.ICELLNET[which(result.ICELLNET$sender!= result.ICELLNET$reciever),]
    
    result.ICELLNET <- split(result.ICELLNET, result.ICELLNET$sr)
    save(result.ICELLNET, file = "../step3_run_lr_benchmark/input/tools_result_for_benchmark/ICELLNET_result.RData")
  }
  load("../step3_run_lr_benchmark/input/tools_result_for_benchmark/ICELLNET_result.RData")
  bench.result.20 <- lapply(close_distant, function(cd){
    bench.result <- bench.ST(cd, norm.data, result.ICELLNET)
    bench.result
  })
  save(bench.result.20, 
       file = "../step3_run_lr_benchmark/result/20_ICELLNET.RData")
}

# 21-NicheNet
if(T){
  # load and handle the reuslt of NicheNet in CID4465
  if(F){
    load("./21-NicheNet/result.NicheNet.RData")
    result.NicheNet <- lapply(result.NicheNet, function(res){
      res <- res$lr
      colnames(res)[1:2] <- c("ligand", "receptor")
      res
    })
    save(result.NicheNet, 
         file = "../step3_run_lr_benchmark/input/tools_result_for_benchmark/NicheNet_result.RData")
  }
  load("../step3_run_lr_benchmark/input/tools_result_for_benchmark/NicheNet_result.RData")
  bench.result.21 <- lapply(close_distant, function(cd){
    bench.result <- bench.ST(cd, norm.data, result.NicheNet)
    bench.result
  })
  save(bench.result.21, 
       file = "../step3_run_lr_benchmark/result/21_NicheNet.RData")
}

rm(list = objects(pattern = "^result"))
rm(close_distant, norm.data, bench.ST, bench.ST.cal)
