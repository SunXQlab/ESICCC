rm(list = ls());gc()
suppressMessages(library(tidyverse))
suppressMessages(library(ggsci))
set.seed(123)

datasets <- list.files('./Data/Step6_LRBenchSamplingResult/', recursive = FALSE) %>%
  gsub('_[0-9]+', '',.) %>% unique()
methods <- c("CellPhoneDB2", "CellPhoneDB3", "CellTalker", 
             "Connectome", "ICELLNET",  "NATMI", "iTALK", 
             "scConnect", "SingleCellSignalR", "CellChat", 
             "RNAMagnet", "PyMINEr",
             "scSeqComm", "NicheNet",
             "scMLnet", "CellCall","CytoTalk", "Domino")

JaccardIndex <- function(a, b) {
  intersection = length(intersect(unique(a), unique(b)))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}

result <- lapply(datasets, function(data){
  print(data)
  temp <- lapply(methods, function(method){
    print(method)
    if(grepl('CK', data)){
      res100 <- readRDS(paste('./Data/Step1_LRPredictionResult',
                              paste0('N_MI_', data), method, 
                              ifelse(method=='iTALK','result1.rds','result.rds'), 
                              sep = '/')) %>% .$result
    }else if(grepl('CID', data)){
      res100 <- readRDS(paste('./Data/Step1_LRPredictionResult',
                              paste0('NG_BC_', data), method, 
                              ifelse(method=='iTALK','result1.rds','result.rds'), 
                              sep = '/')) %>% .$result
    }else if(grepl('Slide14', data)){
      res100 <- readRDS(paste('./Data/Step1_LRPredictionResult',
                              paste0('MouseEmbryo_', data), method, 
                              ifelse(method=='iTALK','result1.rds','result.rds'), 
                              sep = '/')) %>% .$result
    }else{
      res100 <- readRDS(paste('./Data/Step1_LRPredictionResult',
                              data, method, 
                              ifelse(method=='iTALK','result1.rds','result.rds'), 
                              sep = '/')) %>% .$result
    }
    
    if(dim(res100)[[1]]==0 | is.null(res100)){
      tmp <- rep(NA, 5)
    }else{
      res90 <- readRDS(paste('./Data/Step10_LRBenchSamplingResult',
                             paste0(data, '_90'), method, 'result.rds', sep = '/'))%>% .$result
      JI_90 <- JaccardIndex(res100$all, res90$all)
      res80 <- readRDS(paste('./Data/Step10_LRBenchSamplingResult',
                             paste0(data, '_80'), method, 'result.rds', sep = '/'))%>% .$result
      JI_80 <- JaccardIndex(res100$all, res80$all)
      res70 <- readRDS(paste('./Data/Step10_LRBenchSamplingResult',
                             paste0(data, '_70'), method, 'result.rds', sep = '/'))%>% .$result
      JI_70 <- JaccardIndex(res100$all, res70$all)
      res60 <- readRDS(paste('./Data/Step10_LRBenchSamplingResult',
                             paste0(data, '_60'), method, 'result.rds', sep = '/'))%>% .$result
      JI_60 <- JaccardIndex(res100$all, res60$all)
      res50 <- readRDS(paste('./Data/Step10_LRBenchSamplingResult',
                             paste0(data, '_50'), method, 'result.rds', sep = '/'))%>% .$result
      JI_50 <- JaccardIndex(res100$all, res50$all)
      tmp <- c(JI_50, JI_60, JI_70, JI_80, JI_90)
    }
    
    tmp
  })
  names(temp) <- methods
  temp <- do.call(rbind, temp)
  temp <- as.data.frame(temp)
  colnames(temp) <- c('S50', 'S60', 'S70', 'S80', 'S90')
  temp <- tibble::rownames_to_column(temp, 'methods')
  temp$datasets <- data
  temp
})
result <- do.call(rbind, result)
saveRDS(result, file = './Data/Step7_LRBenchSampling/JaccardIndex.rds')
