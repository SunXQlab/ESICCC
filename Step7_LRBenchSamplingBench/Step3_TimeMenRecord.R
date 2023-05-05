rm(list = ls());gc()
suppressMessages(library(tidyverse))
source('./Script/Step11_LRBenchSamplingBench/function.R')
set.seed(123)

DataRecord <- readRDS(file = './Data/Step7_LRBenchSampling/DatasetsRecord.rds')
DataRecord$ratio <- gsub('S', '', DataRecord$ratio)
datasets <- unique(DataRecord$datasets)
methods <- c('CellPhoneDB2', 'CellPhoneDB3',	'CellTalker',	'Connectome',
             'NATMI',	'ICELLNET',	'scConnect', 'CellChat', 'SingleCellSignalR',
             'CellCall',	'scSeqComm',	'NicheNet',	'scMLnet', 'iTALK',	
             'RNAMagnet', 'PyMINEr', 'CytoTalk', 'Domino')
ratios <- c(50, 60, 70, 80, 90, 100)

TimeMenRecord <- lapply(datasets, function(data){
  print(data)
  tmp <- lapply(ratios, function(ratio){
    print(ratio)
    temp <- RunTimeMemRecord(data, ratio, methods, DataRecord)
    temp <- do.call(rbind, temp)
    temp <- tibble::rownames_to_column(temp, 'methods')
    temp$ratios <- ratio
    temp
  })
  tmp <- do.call(rbind, tmp) %>% as.data.frame()
  tmp$datasets <- data
  tmp
})
TimeMenRecord <- do.call(rbind, TimeMenRecord) %>% as.data.frame()
TimeMenRecord$class <- paste(TimeMenRecord$datasets, TimeMenRecord$ratios, sep = '_')

DataRecord$class <- paste(DataRecord$datasets, DataRecord$ratio, sep = '_')
DataRecord <- DataRecord[, -c(1:2)]
TimeMenRecord <- merge(TimeMenRecord, DataRecord, by = 'class')

saveRDS(TimeMenRecord, file = paste0('./Data/Step7_LRBenchSampling/TimeMemRecord.rds'))
