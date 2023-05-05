IDTransform <- function(data, method, result){
  if(dim(result)[1]!=0){
    if(grepl('CID', data) & method == 'RNAMagnet'){
      genes <- readRDS('~/1-Datasets/NatureGenetics_2021_BC/ScriptForCCC/hs2mm_genes.rds')
      result <- mm2hs(result, genes)
    }else if(grepl('CK', data) & method == 'RNAMagnet'){
      genes <- readRDS('~/1-Datasets/Nature_2022_MI/ScriptForCCC/hs2mm_genes.rds')
      result <- mm2hs(result, genes)
    }else if(grepl('Slide', data) & 
             method %in% c('CellPhoneDB2', 'CellPhoneDB3', 'CellTalker', 'NATMI', 
                           'ICELLNET', 'NicheNet', 'scMLnet', 'iTALK', 'cell2cell')){
      genes <- readRDS('~/1-Datasets/MouseEmbryo_GSE166692/sciSpaceData/mm2hs.rds')
      result <- hs2mm(result, genes)
    }
  }
  return(result)
}

CalSI_1 <- function(a, b) {
  intersection = length(intersect(unique(a), unique(b)))
  #union = length(a) + length(b) - intersection
  union <- min(length(unique(a)), length(unique(b)))
  return (intersection/union)
}

CalSI_2 <- function(data){
  methods <- c('CellPhoneDB2', 'CellPhoneDB3',	'CellTalker',	'Connectome',
               'NATMI',	'ICELLNET',	'scConnect', 'CellChat', 'SingleCellSignalR',
               'CellCall',	'scSeqComm',	'NicheNet', 'Domino',	'scMLnet', 'iTALK',	
               'cell2cell', 'RNAMagnet', 'PyMINEr', 'CytoTalk')
  
  # get result path
  if(grepl('CID', data)){
    result.path <- paste0('./Data/Step1_LRPredictionResult/NG_BC_', data)
  }else if(grepl('CK', data)){
    result.path <- paste0('./Data/Step1_LRPredictionResult/N_MI_', data)
  }else if(grepl('Slide', data)){
    result.path <- paste0('./Data/Step1_LRPredictionResult/MouseEmbryo_', data)
  }
  
  SI_res <- lapply(methods, function(method1){
    print(method1)
    if(method1 == 'scMLnet' | method1 == 'Domino'){
      result1 <- readRDS(paste0(result.path, '/', method1, '/result1.rds'))
    }else{
      result1 <- readRDS(paste0(result.path, '/', method1, '/result.rds'))
    }
    
    result1 <- result1$result
    result1 <- IDTransform(data, method1, result1)
    temp.si <- lapply(methods, function(method2){
      if(method2 == 'scMLnet' | method2 == 'Domino'){
        result2 <- readRDS(paste0(result.path, '/', method2, '/result1.rds'))
      }else{
        result2 <- readRDS(paste0(result.path, '/', method2, '/result.rds'))
      }
      
      result2 <- result2$result
      result2 <- IDTransform(data, method2, result2)
      if(dim(result1)[1]==0 | dim(result2)[1]==0){
        tmp.si <- NA
      }else{
        tmp.si <- CalSI_1(as.character(result1$all), as.character(result2$all))
      }
      
      tmp.si
    })
    temp.si <- data.frame(SI = unlist(temp.si), row.names = methods)
    temp.si
  })
  SI_res <- do.call(cbind, SI_res)
  colnames(SI_res) <- methods
  #SI_res <- SI_res[!apply(SI_res, 1, function(x){all(is.na(x))}), ]
  #SI_res <- SI_res[, !apply(SI_res, 2, function(x){all(is.na(x))})]
  return(SI_res)
}

CalRSI <- function(data){
  # remove cell2cell
  methods <- c('CellPhoneDB2', 'CellPhoneDB3',	'CellTalker',	'Connectome',
               'NATMI',	'ICELLNET',	'scConnect', 'CellChat', 'SingleCellSignalR',
               'CellCall',	'scSeqComm',	'NicheNet', 'Domino',	'scMLnet', 'iTALK',	
               'cell2cell', 'RNAMagnet', 'PyMINEr', 'CytoTalk')
  
  # get result path
  if(grepl('CID', data)){
    result.path <- paste0('./Data/Step1_LRPredictionResult/NG_BC_', data)
  }else if(grepl('CK', data)){
    result.path <- paste0('./Data/Step1_LRPredictionResult/N_MI_', data)
  }else if(grepl('Slide', data)){
    result.path <- paste0('./Data/Step1_LRPredictionResult/MouseEmbryo_', data)
  }
  
  RSI_res <- lapply(methods, function(method1){
    print(paste0('Method1: ', method1))
    if(method1 == 'scMLnet' | method1 == 'Domino'){
      result1 <- readRDS(paste0(result.path, '/', method1, '/result1.rds')) # recode the LR score
    }else{
      result1 <- readRDS(paste0(result.path, '/', method1, '/result.rds'))
    }
    
    result1 <- as.data.frame(result1$result)
    result1$LRscore <- as.numeric(result1$LRscore)
    result1 <- IDTransform(data, method1, result1)
    if(method1 == 'CytoTalk'){
      result1$rank <- rank(result1$LRscore)
    }else{
      result1$rank <- rank(-result1$LRscore)
    }
    
    temp.rsi <- lapply(methods, function(method2){
      print(paste0('Method2: ', method2))
      if(method2 == 'scMLnet' | method2 == 'Domino'){
        result2 <- readRDS(paste0(result.path, '/', method2, '/result1.rds'))
      }else{
        result2 <- readRDS(paste0(result.path, '/', method2, '/result.rds'))
      }
      
      result2 <- as.data.frame(result2$result)
      result2$LRscore <- as.numeric(result2$LRscore)
      result2 <- IDTransform(data, method2, result2)
      if(method2 == 'CytoTalk'){
        result2$rank <- rank(result2$LRscore)
      }else{
        result2$rank <- rank(-result2$LRscore)
      }
      
      if(dim(result1)[1]==0 | dim(result2)[1]==0){
        tmp.rsi <- NA
      }else{
        overlap.lr <- intersect(result1$all, result2$all)
        if(length(overlap.lr)!=0){
          rank1 <- result1[which(result1$all %in% overlap.lr),]$rank/dim(result1)[1]
          rank2 <- result2[which(result2$all %in% overlap.lr), ]$rank/dim(result2)[1]
          mean.rank <- mean(abs(rank1-rank2))
          tmp.rsi <- 1-mean.rank
        }else{tmp.rsi <- 0}
        }
      tmp.rsi
    })
    temp.rsi <- data.frame(RSI = unlist(temp.rsi), row.names = methods)
    temp.rsi
  })
  RSI_res <- do.call(cbind, RSI_res)
  colnames(RSI_res) <- methods
  #RSI_res <- RSI_res[!apply(RSI_res, 1, function(x){all(is.na(x))}), ]
  #RSI_res <- RSI_res[, !apply(RSI_res, 2, function(x){all(is.na(x))})]
  
  return(RSI_res)
}
