rm(list = ls()); gc()
suppressMessages(library(ggplot2))
source('./Script/Step9_LRTBench/function.R')
set.seed(123)

# Record the number of cell types in each dataset
if(F){
  suppressMessages(library(Seurat))
  set.seed(123)
  ser1.path <- list.files('~/1-Datasets/NatureGenetics_2021_BC/ScriptForCCC/ST_Data', full.names = TRUE)
  ser2.path <- list.files('~/1-Datasets/Nature_2022_MI/ScriptForCCC/ST_Data', full.names = TRUE)
  ser3.path <- list.files('~/1-Datasets/CancerCell_2022_Glioma/ScriptForLRTBenchmark/ST_Data', full.names = TRUE)
  ser.path <- c(ser1.path, ser2.path, ser3.path)
  rm(ser1.path, ser2.path, ser3.path);gc()
  
  tmp <- c()
  for (i in seq(ser.path)) {
    ser <- readRDS(ser.path[i])
    remove.ct <- as.numeric(length(which(table(ser$celltype)<3)))
    temp <- length(unique(ser$celltype))-remove.ct
    tmp <- c(tmp, temp)
    rm(ser); gc()
  }
  ct_record <- data.frame(datasets = gsub('_ser.rds', '', 
                                          unlist(lapply(stringr::str_split(ser.path, '/'), 
                                                        function(x){x[8]}))),
                          num = as.numeric(tmp))
  
  saveRDS(ct_record, file = './Data/Step9_LRTBenchResult/DatasetCellTypeInfo.rds')
  rm(i, ser.path, temp, tmp, remove.ct, ct_record);gc()
}

# Record the running time of methods in each dataset
if(T){
  datasets <- list.dirs('./Data/Step8_LRTPredictionResult', recursive = FALSE, full.names = FALSE)
  methods <- c('CytoTalk', 'NicheNet', 'HoloNet', 'MISTy', 'stMLnet')
  DataRecord <- readRDS('./Data/Step9_LRTBenchResult/DatasetCellTypeInfo.rds')
  wd <- './Data/Step8_LRTPredictionResult/'
  
  TimeMemRecord <- lapply(datasets, function(data){
    print(data)
    temp <- RunTimeMemRecord(wd, data, methods)
    temp
  })
  names(TimeMemRecord) <- datasets
  TimeMemRecord <- lapply(TimeMemRecord, function(data){
    temp <- do.call(rbind, data)
    temp <- tibble::rownames_to_column(temp, 'methods')
    temp
  })
  TimeMemRecord <- do.call(rbind, TimeMemRecord)
  TimeMemRecord <- tibble::rownames_to_column(TimeMemRecord, 'datasets')
  TimeMemRecord$datasets <- gsub('\\.[0-9]+', '', TimeMemRecord$datasets)
  
  # Specific handle: the running time of CytoTalk in GBM datasets
  CytoTalk_TimeMemRecord <- TimeMemRecord[which(TimeMemRecord$methods == 'CytoTalk' & grepl('UKF', TimeMemRecord$datasets)), ]
  CytoTalk_TimeMemRecord <- merge(CytoTalk_TimeMemRecord, DataRecord, by = 'datasets')
  CytoTalk_TimeMemRecord$total_num <- 2*(CytoTalk_TimeMemRecord$num-1)
  CytoTalk_TimeMemRecord$actual_num <- 2*(CytoTalk_TimeMemRecord$num-1)-1
  CytoTalk_TimeMemRecord$clock_time <- (CytoTalk_TimeMemRecord$clock_time/CytoTalk_TimeMemRecord$total_num)*CytoTalk_TimeMemRecord$actual_num
  CytoTalk_TimeMemRecord$linux_time <- (CytoTalk_TimeMemRecord$linux_time/CytoTalk_TimeMemRecord$total_num)*CytoTalk_TimeMemRecord$actual_num
  CytoTalk_TimeMemRecord <- CytoTalk_TimeMemRecord[, -c(7:9)]
  TimeMemRecord <- TimeMemRecord[-which(TimeMemRecord$methods == 'CytoTalk' & grepl('UKF', TimeMemRecord$datasets)), ]
  TimeMemRecord <- rbind(TimeMemRecord, CytoTalk_TimeMemRecord)
  rm(CytoTalk_TimeMemRecord, DataRecord, datasets, methods, wd); gc()
}