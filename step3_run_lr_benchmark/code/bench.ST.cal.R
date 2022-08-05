# the function of calculating the mutual information | pcc
bench.ST.cal <- function(close_distant, norm.data, result, cp, choose, logic){
  result.tmp <- result[[cp]]
  if(!logic){
    cp.1 <- unlist(stringr::str_split(cp, "_"))[1]
    cp.2 <- unlist(stringr::str_split(cp, "_"))[2]
    cp <- paste(cp.2, cp.1, sep = "_")
  }
  
  if(choose == "close"){
    distant.tmp <- close_distant[[cp]]$close
  }else{
    distant.tmp <- close_distant[[cp]]$distant
  }
  
  if(logic){
    sender.mat <- norm.data[,distant.tmp$cell.1]
    reciever.mat <- norm.data[, distant.tmp$cell.2]
  }else{
    sender.mat <- norm.data[,distant.tmp$cell.2]
    reciever.mat <- norm.data[, distant.tmp$cell.1]
  }
  
  sender.mat <- sender.mat[which(rowSums(sender.mat)!=0), ]
  reciever.mat <- reciever.mat[which(rowSums(reciever.mat)!=0), ]
  
  tmp.result <- lapply(1:dim(result.tmp)[1], function(n){
    ligands <- result.tmp[n,"ligand"]
    ligands <- unlist(stringr::str_split(ligands, "&"))
    
    receptors <- result.tmp[n, "receptor"]
    receptors <- unlist(stringr::str_split(receptors, "&"))

    if(all(ligands %in% rownames(sender.mat)) & all(receptors %in% rownames(reciever.mat))){
      tmp.lr <- data.frame()
      for (ligand in ligands) {
        temp <- data.frame(ligand = rep(ligand, length(receptors)), receptor = receptors )
        tmp.lr <- rbind(tmp.lr, temp)
      }
      rm(ligand, temp)
      temp.result <- lapply(1:dim(tmp.lr)[1], function(m){
        ligand <- tmp.lr[m, "ligand"]
        receptor <- tmp.lr[m, "receptor"]
        
        tmp.sender.mat <- data.frame(value.1 = sender.mat[ligand,], barcode.1 = colnames(sender.mat))
        tmp.reciever.mat <- data.frame(value.2 = reciever.mat[receptor,], barcode.2 = colnames(reciever.mat))
        tmp.mat = cbind(tmp.sender.mat, tmp.reciever.mat)
        
        data <- tmp.mat[,c("value.1", "value.2")]
        
        ## calculate the mutual information
        data.disc <- infotheo::discretize(data)
        tmp.mutinfo <- infotheo::mutinformation(data.disc[,1],data.disc[,2])
        
        ## calculate the Pearson Correlation Coefficient
        tmp.pcc <- cor.test(data[,1], data[,2], method = "pearson")$estimate[[1]]
        tmp.pcc.p <- cor.test(data[,1], data[,2], method = "pearson")$p.value
        tmp <- c(tmp.mutinfo, tmp.pcc, tmp.pcc.p)
        tmp
      })
      result.temp <- do.call(rbind, temp.result)
      result.temp <- as.data.frame(result.temp)
      result.temp$lr <- paste(tmp.lr$ligand, tmp.lr$receptor, sep = "_")
    }else{
      tmp.lr <- data.frame()
      for (ligand in ligands) {
        temp <- data.frame(ligand = rep(ligand, length(receptors)), receptor = receptors )
        tmp.lr <- rbind(tmp.lr, temp)
      }
      rm(ligand, temp)
      
      temp.result <- lapply(1:dim(tmp.lr)[1], function(m){
        tmp.mutinfo <- NA
        tmp.pcc <- NA
        tmp.pcc.p <- NA
        tmp <- c(tmp.mutinfo, tmp.pcc, tmp.pcc.p)
        tmp <- as.data.frame(t(tmp))
      })
      
      result.temp <- do.call(rbind, temp.result)
      result.temp$lr <- paste(tmp.lr$ligand, tmp.lr$receptor, sep = "_")
    }
    result.temp
  })
  bench.result <- do.call(rbind, tmp.result)
}

bench.ST <- function(close_distant, norm.data, result){
  result.bench <- lapply(names(result), function(cp){
    print(cp)
    if(cp %in% names(close_distant)){
      close.bench <- bench.ST.cal(close_distant, norm.data, result, 
                                  cp, "close", TRUE)
      distant.bench <- bench.ST.cal(close_distant, norm.data, result, 
                                     cp, "distant", TRUE)
      identical(distant.bench$lr, close.bench$lr)
      tmp.result <- cbind(close.bench, distant.bench)
      tmp.result <- tmp.result[,-4]
      colnames(tmp.result) <- c("close.mutinfo", "close.pcc", "close.pcc.p",
                                "distant.mutinfo", "distant.pcc","distant.pcc.p", "lr")
    }else{
      close.bench <- bench.ST.cal(close_distant, norm.data, result, 
                                  cp, "close", FALSE)
      distant.bench <- bench.ST.cal(close_distant, norm.data, result, 
                                     cp,  "distant", FALSE)
      tmp.result <- cbind(close.bench, distant.bench)
      tmp.result <- tmp.result[,-4]
      colnames(tmp.result) <- c("close.mutinfo", "close.pcc", "close.pcc.p",
                                "distant.mutinfo", "distant.pcc","distant.pcc.p", "lr")
      
    }
    tmp.result <- as.data.frame(tmp.result)
    tmp.result
  })
  names(result.bench) <- names(result)
  result.bench
}

