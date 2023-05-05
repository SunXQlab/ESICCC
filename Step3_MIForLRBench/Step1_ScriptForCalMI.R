rm(list = ls())
suppressMessages(library(Seurat))
source('./Script/Step3_MIForLRBench/function.R')
set.seed(123)

datasets <- c('CID4465', 'CID44971', 'CK357', 'CK358',
              'CK368', 'CK162', 'CK362', 'CK361',
              'CK161', 'CK165', 'Slide14')
methods <- c('CellPhoneDB2', 'CellPhoneDB3',	'CellTalker',	'Connectome',
             'NATMI',	'ICELLNET',	'scConnect', 'CellChat', 'SingleCellSignalR',
             'CellCall',	'scSeqComm',	'NicheNet', 'Domino',	'scMLnet', 'iTALK',	
             'cell2cell', 'RNAMagnet', 'PyMINEr', 'CytoTalk')
ck.dataset <- readRDS('~/1-Datasets/Nature_2022_MI/dataset_id.rds')
# loop: datasets
for(data in datasets){
  print(data)
  if(grepl('CK', data)){
    # get normalized data of ST
    STser <- readRDS(paste0('./Data/Step2_PreSTForLRBench/', ck.dataset$ST[which(ck.dataset$snRNA == data)], '/STser.rds'))
    # get close and distant cell pairs
    CloDistCP <- readRDS(paste0('./Data/Step2_PreSTForLRBench/', ck.dataset$ST[which(ck.dataset$snRNA == data)], '/CloDistCP.rds'))
  }else{
    # get normalized data of ST
    STser <- readRDS(paste0('./Data/Step2_PreSTForLRBench/', data, '/STser.rds'))
    # get close and distant cell pairs
    CloDistCP <- readRDS(paste0('./Data/Step2_PreSTForLRBench/', data, '/CloDistCP.rds'))
  }
  
  # get normalized data of ST
  norm.data <- GetAssayData(STser, 'data', 'SCT')
  #norm.data <- as.matrix(norm.data)
  # get close and distant cell pairs
  CloDistCP[['0.5']] <- NULL
  # get celltype of ST data
  celltype.st <- unique(STser$celltype)
  rm(STser); gc()
  
  # get result path
  if(grepl('CID', data)){
    result.path <- paste0('./Data/Step1_LRPredictionResult/NG_BC_', data)
  }else if(grepl('CK', data)){
    result.path <- paste0('./Data/Step1_LRPredictionResult/N_MI_', data)
  }else if(grepl('Slide', data)){
    result.path <- paste0('./Data/Step1_LRPredictionResult/MouseEmbryo_', data)
  }
  
  output.path <- paste0('./Data/Step3_MIPCCForLRBench/', data, '_result.rds')
  dataset_EvalIndex1 <- lapply(methods, function(method){
    print(method)
    
    # remove the LR pairs whose sender/receiver celltypes are not in celltypes of ST
    result <- readRDS(paste0(result.path, '/', method, '/result.rds'))
    result <- result$result
    result <- result[which(result$Sender %in% celltype.st), ]
    result <- result[which(result$Receiver %in% celltype.st), ]
    result$sr <- paste(result$Sender, result$Receiver, sep = '_')
    CellPairs <- names(CloDistCP$`0.1`)
    CellPairs1 <- as.data.frame(stringr::str_split(CellPairs, '_', simplify = TRUE))
    CellPairs1 <- paste(CellPairs1$V2, CellPairs1$V1, sep = '_')
    CellPairs <- c(CellPairs, CellPairs1)
    remove.sr <- setdiff(unique(result$sr), CellPairs)
    if(length(remove.sr) != 0){
      result <- result[-which(result$sr %in% remove.sr), ]
    }
    
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
      
      EvaIndex1_result <- lapply(CloDistCP, function(CloDist){
        result.tmp <- EvaIndex1_2(CloDist, norm.data, result)
        result.tmp
      })
      
    }else{
      EvaIndex1_result <- NA
    }
    
    return(EvaIndex1_result)
  })
  names(dataset_EvalIndex1) <- methods
  saveRDS(dataset_EvalIndex1, file = output.path)
}