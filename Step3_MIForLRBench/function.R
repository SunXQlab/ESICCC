EvaIndex1_1 <- function(CloDist, norm.data, result, cp, choose, logic){
  result.tmp <- result[which(result$sr == cp), ]
  if(!logic){
    cp.1 <- unlist(stringr::str_split(cp, "_"))[1]
    cp.2 <- unlist(stringr::str_split(cp, "_"))[2]
    cp <- paste(cp.2, cp.1, sep = "_")
  }
  
  if(choose == "close"){
    distant.tmp <- CloDist[[cp]]$close
  }else if(choose == 'distant'){
    distant.tmp <- CloDist[[cp]]$distant
  }
  
  tmp.result <- lapply(1:dim(result.tmp)[1], function(n){
    #print(n)
    ligands <- result.tmp[n,"Ligand"]
    ligands <- unlist(stringr::str_split(ligands, "&"))
    
    receptors <- result.tmp[n, "Receptor"]
    receptors <- unlist(stringr::str_split(receptors, "&"))
    
    tmp.lr <- expand.grid(ligands, receptors)
    colnames(tmp.lr) <- c('ligand', 'receptor')
    
    temp.result <- lapply(1:dim(tmp.lr)[1], function(m){
      ligand <- tmp.lr[m, "ligand"]
      receptor <- tmp.lr[m, "receptor"]
      
      if((ligand %in% rownames(norm.data)) & (receptor %in% rownames(norm.data))){
        if(logic){
          value.1 <- norm.data[ligand, distant.tmp$cell.1]
          value.2 <- norm.data[receptor, distant.tmp$cell.2]
        }else{
          value.1 <- norm.data[ligand, distant.tmp$cell.2]
          value.2 <- norm.data[receptor, distant.tmp$cell.1]
        }
        
        if(sum(value.1)!=0 & sum(value.2)!=0){
          tmp.sender.mat <- data.frame(value.1 = value.1, barcode.1 = names(value.1))
          tmp.reciever.mat <- data.frame(value.2 = value.2, barcode.2 = names(value.2))
          data.tmp = cbind(tmp.sender.mat, tmp.reciever.mat)
          data.tmp <- data.tmp[,c("value.1", "value.2")]
          rm(value.1, value.2, tmp.reciever.mat, tmp.sender.mat)
          
          ## calculate the mutual information
          data.disc <- infotheo::discretize(data.tmp)
          tmp.mi <- infotheo::mutinformation(data.disc[,1],data.disc[,2])
          
          ## calculate the Pearson Correlation Coefficient
          pcc.tmp <- cor.test(data.tmp[,1], data.tmp[,2], method = "pearson")
          tmp.pcc <- pcc.tmp$estimate[[1]]
          tmp.pcc.p <- pcc.tmp$p.value
          tmp <- c(tmp.mi, tmp.pcc, tmp.pcc.p)
          
          tmp
        }else{
          tmp.mi <- NA
          tmp.pcc <- NA
          tmp.pcc.p <- NA
          tmp <- c(tmp.mi, tmp.pcc, tmp.pcc.p)
          tmp <- as.data.frame(t(tmp))
          tmp
        }
      }else{
        tmp.mi <- NA
        tmp.pcc <- NA
        tmp.pcc.p <- NA
        tmp <- c(tmp.mi, tmp.pcc, tmp.pcc.p)
        tmp <- as.data.frame(t(tmp))
        tmp
      }

      tmp
      
    })
    names(temp.result) <- paste(tmp.lr$ligand, tmp.lr$receptor, sep = "_")
    #temp.result[which(is.na(temp.result))] <- NULL
    temp.result <- do.call(rbind, temp.result)
    temp.result <- as.data.frame(temp.result)
    temp.result <- tibble::rownames_to_column(temp.result, 'lr')
    
    temp.result
  })

  bench.result <- do.call(rbind, tmp.result)
  
  bench.result
}

EvaIndex1_2 <- function(CloDist, norm.data, result){
  result.bench <- lapply(unique(result$sr), function(cp){
    print(cp)
    if(cp %in% names(CloDist)){
      close.bench <- EvaIndex1_1(CloDist, norm.data, result, cp, "close", TRUE)
      distant.bench <- EvaIndex1_1(CloDist, norm.data, result, cp, "distant", TRUE)
      identical(distant.bench$lr, close.bench$lr)
      tmp.result <- cbind(close.bench, distant.bench)
      tmp.result <- tmp.result[,-5]
      colnames(tmp.result) <- c("lr", "closeMI", "closePCC", "closePCCPvalue",
                                "distantMI", "distantPCC","distantPCCPvalue")
    }else{
      close.bench <- EvaIndex1_1(CloDist, norm.data, result, cp, "close", FALSE)
      distant.bench <- EvaIndex1_1(CloDist, norm.data, result, cp,  "distant", FALSE)
      tmp.result <- cbind(close.bench, distant.bench)
      tmp.result <- tmp.result[,-5]
      colnames(tmp.result) <- c("lr","closeMI", "closePCC", "closePCCPvalue",
                                "distantMI", "distantPCC","distantPCCPvalue")
      
    }
    tmp.result <- as.data.frame(tmp.result)
    tmp.result
  })
  names(result.bench) <- unique(result$sr)
  result.bench
}


mm2hs <- function(result, genes){
  genes <- genes[, c('mouse', 'human')]
  
  result <- merge(result, genes, by.x = 'Ligand', by.y = 'mouse', all.x = TRUE)
  multi <- which(is.na(result$human) & grepl('&', result$Ligand))
  if(length(multi)>0){
    result$human[multi] <- unlist(lapply(stringr::str_split(result$Ligand[multi], '&'), function(x){
      tmp <- unlist(lapply(seq(x), function(y){genes$human[which(genes$mouse %in% x[y])]}))
      paste(tmp, collapse = '&')
    }))
  }
  result <- result[!(is.na(result$human) & !grepl('&', result$Ligand)), ]
  result$Ligand <- NULL
  colnames(result)[which(colnames(result) == 'human')] <- 'Ligand'
  
  result <- merge(result, genes, by.x = 'Receptor', by.y = 'mouse', all.x = TRUE)
  multi <- which(is.na(result$human) & grepl('&', result$Receptor))
  if(length(multi)>0){
    result$human[multi] <- unlist(lapply(stringr::str_split(result$Receptor[multi], '&'), function(x){
      tmp <- unlist(lapply(seq(x), function(y){genes$human[which(genes$mouse %in% x[y])]}))
      paste(tmp, collapse = '&')
    }))
  }
  result <- result[!(is.na(result$human) & !grepl('&', result$Receptor)), ]
  result$Receptor <- NULL
  colnames(result)[which(colnames(result) == 'human')] <- 'Receptor'
  result$all <- paste(result$Sender, result$Ligand, result$Receiver, result$Receptor, sep = '_')
  
  return(result)
}

hs2mm <- function(result, genes){
  genes <- genes[, c('mouse', 'human')]
  
  result <- merge(result, genes, by.x = 'Ligand', by.y = 'human', all.x = TRUE)
  multi <- which(is.na(result$mouse) & grepl('&', result$Ligand))
  if(length(multi)>0){
    result$mouse[multi] <- unlist(lapply(stringr::str_split(result$Ligand[multi], '&'), function(x){
      tmp <- unlist(lapply(seq(x), function(y){genes$mouse[which(genes$human %in% x[y])]}))
      paste(tmp, collapse = '&')
    }))
  }
  result <- result[!(is.na(result$mouse) & !grepl('&', result$Ligand)), ]
  result$Ligand <- NULL
  colnames(result)[which(colnames(result) == 'mouse')] <- 'Ligand'
  
  result <- merge(result, genes, by.x = 'Receptor', by.y = 'human', all.x = TRUE)
  multi <- which(is.na(result$mouse) & grepl('&', result$Receptor))
  if(length(multi)>0){
    result$mouse[multi] <- unlist(lapply(stringr::str_split(result$Receptor[multi], '&'), function(x){
      tmp <- unlist(lapply(seq(x), function(y){genes$mouse[which(genes$human %in% x[y])]}))
      paste(tmp, collapse = '&')
    }))
  }
  result <- result[!(is.na(result$mouse) & !grepl('&', result$Receptor)), ]
  result$Receptor <- NULL
  colnames(result)[which(colnames(result) == 'mouse')] <- 'Receptor'
  result$all <- paste(result$Sender, result$Ligand, result$Receiver, result$Receptor, sep = '_')
  
  return(result)
}


# further handle the results of MI
Eval1Process <- function(result.path){
  Eval1 <- readRDS(result.path)
  Eval1[which(is.na(Eval1))] <- NULL
  Eval1 <- lapply(Eval1, function(methods_eval){
    perc_temp <- lapply(methods_eval, function(perc_eval){
      perc_tmp <- do.call(rbind, perc_eval)
      perc_tmp <- tibble::rownames_to_column(perc_tmp, 'sr')
      perc_tmp$sr <- gsub('\\.[0-9]+', '', perc_tmp$sr)
      perc_tmp
    })
    names(perc_temp) <- c(10, 20, 30, 40)
    perc_temp <- do.call(rbind, perc_temp)
    perc_temp <- tibble::rownames_to_column(perc_temp, 'perc')
    perc_temp$perc <- gsub('\\.[0-9]+', '', perc_temp$perc)
    perc_temp
  })
  Eval1 <- do.call(rbind, Eval1)
  Eval1 <- tibble::rownames_to_column(Eval1, 'methods')
  Eval1$methods <- gsub('\\.[0-9]+', '', Eval1$methods)
  allna <- apply(Eval1[, c("closeMI", "closePCC", "distantMI", "distantPCC")], 
                 1, FUN = function(x){any(is.na(x))})
  Eval1 <- Eval1[which(!allna), ]
  Eval1$all <- paste(Eval1$methods, Eval1$perc, Eval1$sr, Eval1$lr, sep = '_')
  Eval1 <- dplyr::distinct(Eval1, all, .keep_all = TRUE)
  
  return(Eval1)
}

# Calculate the DLRC index of each method —— p.value(whether significant?) + median + weight
EvalIndex_DLRC <- function(Eval1Result){
  Eval1Result$class <- paste(Eval1Result$methods, Eval1Result$perc, sep = '_')
  
  temp <- lapply(unique(Eval1Result$class), function(class){
    result_sub <- Eval1Result[which(Eval1Result$class == class),]
    medianMI_close <- median(result_sub$closeMI)
    medianMI_distant <- median(result_sub$distantMI)
    medianPCC_close <- median(abs(result_sub$closePCC))
    medianPCC_distant <- median(abs(result_sub$distantPCC))
    PCC_pval <- wilcox.test (result_sub$closePCC, result_sub$distantPCC, 
                             alternative = "greater")$p.value
    MI_pval <- wilcox.test (result_sub$closeMI, result_sub$distantMI, 
                            alternative = "greater")$p.value
    return(data.frame(medianMI_close = medianMI_close, medianMI_distant = medianMI_distant, MI_pval = MI_pval,
                      medianPCC_close = medianPCC_close, medianPCC_distant = medianPCC_distant, PCC_pval = PCC_pval))
  })
  names(temp) <- unique(Eval1Result$class)
  temp <- do.call(rbind, temp)
  temp <- tibble::rownames_to_column(temp, 'class')
  temp <- tidyr::separate(temp, class, c('methods', 'perc'), '_')
  #temp <- temp[-which(temp$perc == 40),]
  
  tmp <- lapply(unique(temp$methods), function(method){
    result_sub <- temp[which(temp$methods == method), ]
    if(dim(result_sub)[1]==4){
      result_sub$MIIndex <- (0.5-as.numeric(result_sub$perc)*0.01)*
        ifelse(result_sub$MI_pval<0.05, 1, 0)*
        (result_sub$medianMI_close - result_sub$medianMI_distant)
      
      result_sub$PCCIndex <- (0.5-as.numeric(result_sub$perc)*0.01)*
        ifelse(result_sub$PCC_pval<0.05, 1, 0)*
        (result_sub$medianPCC_close - result_sub$medianPCC_distant)
      return(data.frame(MIIndex = sum(result_sub$MIIndex), PCCIndex = sum(result_sub$PCCIndex)))
    }
  })
  names(tmp) <- unique(temp$methods)
  tmp <- do.call(rbind, tmp)
  tmp <- tibble::rownames_to_column(tmp, 'methods')
  tmp$MIRank <- rank(-tmp$MIIndex)
  tmp$PCCRank <- rank(-tmp$PCCIndex)
  return(tmp)
}
