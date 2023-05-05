rm(list = ls());gc()
suppressMessages(library(Seurat))
suppressMessages(library(tidyverse))
set.seed(123)

files <- list.dirs('./Data/Step1_LRPredictionResult', recursive = FALSE)
methods <- c('scMLnet', 'Domino')

for (file in files) {
  print(file)
  if(grepl('Slide', file)){
    ser <- readRDS('~/1-Datasets/MouseEmbryo_GSE166692/ScriptForCCC/DataForCCC/Slide14_hs_ser.rds')
  }else if(grepl('N_MI', file)){
    sample <- substring(file, 38)
    ser <- readRDS(paste0('~/1-Datasets/Nature_2022_MI/ScriptForCCC/scRNA_Data/', sample, '_ser.rds'))
  }else if(grepl('NG_BC', file)){
    sample <- substring(file, 39)
    ser <- readRDS(paste0('~/1-Datasets/NatureGenetics_2021_BC/ScriptForCCC/scRNA_Data/', sample, '_ser.rds'))
  }else if(grepl('pbmc', file)){
    sample <- substring(file, 33)
    ser <- readRDS(paste0('~/1-Datasets/10XGenomics_PBMC/ScriptForCellTypeAnn/DataForCCC/', sample, '_ser.rds'))
  }else if(grepl('Sper', file)){
    ser <- readRDS('~/1-Datasets/Spermatogenesis_GSE106487/ScriptForCCC/Sper_ser.rds')
  }
  
  data <- GetAssayData(ser, 'data', 'RNA') %>% 
    as.matrix(.) %>% t(.) %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column(., 'barcodes')
  meta <- ser@meta.data %>% .[, 'celltype', drop =FALSE] %>% tibble::rownames_to_column(., 'barcodes')
  rm(ser);gc()
  
  data <- merge(data, meta, by = 'barcodes')
  data$barcodes <- NULL
  data <- aggregate(.~celltype, data, mean)
  data <- pivot_longer(data, -celltype, 'genes')
  data$ct_genes <- paste(data$celltype, data$genes, sep = '_')
  data <- data[ -c(1:2)]
  rm(meta);gc()
  
  for (method in methods) {
    print(method)
    result <- readRDS(paste0(file, '/', method, '/result.rds'))
    result_tmp <- result$result
    result_tmp$sl <- paste(result_tmp$Sender, result_tmp$Ligand, sep = '_')
    result_tmp$rr <- paste(result_tmp$Receiver, result_tmp$Receptor, sep = '_')
    result_tmp <- merge(result_tmp, data, by.x = 'sl', by.y = 'ct_genes')
    colnames(result_tmp)[8] <- 'ligand_value'
    result_tmp <- merge(result_tmp, data, by.x = 'rr', by.y = 'ct_genes')
    colnames(result_tmp)[9] <- 'receptor_value'
    result_tmp$LRscore <- result_tmp$ligand_value*result_tmp$receptor_value
    result_tmp <- result_tmp[,c('Ligand', 'Receptor', 'Sender', 'Receiver', 'LRscore', 'all')]
    result$result <- result_tmp
    saveRDS(result, paste0(file, '/', method, '/result1.rds'))
  }
}
