rm(list = ls());gc()
setwd('~/2-CCC-benchmark-202204-202208/CCC-benchmark-202212/Data/Step8_LRTPredictionResult/')
library(dplyr)
library(ROCR)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(tidyverse)
source('../../Script/Step5_BenchBasedCAGEProteomic/function.R')
set.seed(123)


###########
## color ##
###########

scales::show_col(pal_aaas(palette = "default", alpha = 0.6)(10))
mycolors_aaas <- pal_aaas(palette = "default", alpha = 0.6)(10)

mycolor_software <- mycolors_aaas[c(1:5)]
names(mycolor_software) <- c('CytoTalk', "NicheNet", "MISTy", "HoloNet", "stMLnet")
scales::show_col(mycolor_software)

#####################
## load all result ##
#####################

if(F){
  data_dirs <- list.files('.', full.names = TRUE)
  all_result <- lapply(data_dirs, function(data_dir){
    print(data_dir)
    methods_dirs <- list.files(data_dir, full.names = TRUE)
    methods <- list.files(data_dir)
    used_methods <- c('CytoTalk', 'HoloNet', 'MISTy', 'NicheNet', 'stMLnet')
    index <- which(methods %in% used_methods)
    methods_dirs <- methods_dirs[index]
    methods <- methods[index]
    data_result <- lapply(seq(methods_dirs), function(i){
      print(i)
      res_files <- list.files(methods_dirs[i], '_result.rds', full.names = TRUE)
      if(length(res_files)>0){
        res <- lapply(res_files, function(res_file){
          res_temp <- readRDS(res_file)
          if(dim(res_temp)[[1]] != 0){
            res_temp$methods <- methods[i]
          }else{
            res_temp <- NA
          }
          res_temp
        })
        names(res) <- list.files(methods_dirs[i], '_result.rds') %>% gsub('_result.rds', '', .)
        res[which(is.na(res))] <- NULL
      }else{
        res <- NA
      }
      res
    })
    names(data_result) <- methods
    data_result[which(is.na(data_result))] <- NULL
    data_result
  })
  names(all_result) <- list.files('.')
  saveRDS(all_result, file = 'all_result.rds')
}

setwd('~/2-CCC-benchmark-202204-202208/CCC-benchmark-202212/Data')
source('../Script/Step9_LRTBench/function.R')
all_result <- readRDS('./Step8_LRTPredictionResult/all_result.rds')

###############
## Calculate ##
###############

result_index <- lapply(names(all_result), function(data){
  print(data)
  data_res <- all_result[[data]]
  
  if(grepl('CID', data)){
    
    ## groundtrue
    degs_ls <- readRDS('./Step0_SharedInfo/LRTBench/BC_celline/degs_ls_p.rds')
    
    methods_res <- lapply(data_res, function(x){x$CancerEpithelial})
    methods_res[unlist(lapply(methods_res, is.null))] <- NULL
    
    result <- getallmetrics(methods_res, degs_ls)
    result$celltype <- 'CancerEpithelial'
    
    result1 <- getallrecord(methods_res, degs_ls)
    saveRDS(result1, file = paste0('./Step8_LRTBenchResult/record/', data, '.rds'))
    
  }else if(grepl('control', data)){
    
    ## groundtrue
    degs_ls <- readRDS('~/1-Datasets/CellLinesDatasets/ProjectForLRTBench/result/degs_ls_p.rds')
    used <- which(grepl('GSE206947', names(degs_ls)) | 
                    grepl('GSE181575', names(degs_ls)) | 
                    grepl('GSE123018', names(degs_ls)))
    degs_ls_used <- degs_ls[used]
    
    methods_res <- lapply(data_res, function(x){x$Fibroblast})
    methods_res[unlist(lapply(methods_res, is.null))] <- NULL
    
    result <- getallmetrics(methods_res, degs_ls_used)
    result$celltype <- 'Fibroblast'
    
    result1 <- getallrecord(methods_res, degs_ls_used)
    saveRDS(result1, file = paste0('./Step8_LRTBenchResult/record/', data, '.rds'))
    
  }else if(grepl('UKF', data)){
    
    ## groundtrue
    degs_ls <- readRDS('~/1-Datasets/CellLinesDatasets/ProjectForLRTBench/result/degs_ls_p.rds')
    
    methods_res_macro <- lapply(data_res, function(x){x$macrophages})
    methods_res_macro[unlist(lapply(methods_res_macro, is.null))] <- NULL
    
    methods_res_mal <- lapply(data_res, function(x){x$Malignant})
    methods_res_mal[unlist(lapply(methods_res_mal, is.null))] <- NULL
    
    used <- which(grepl('GSE69104', names(degs_ls)) & grepl('TAM', names(degs_ls)))
    degs_ls_used <- degs_ls[used]
    result_macro <- getallmetrics(methods_res_macro, degs_ls_used)
    result_macro1 <- getallrecord(methods_res_macro, degs_ls_used)
    saveRDS(result_macro1, file = paste0('./Step8_LRTBenchResult/record/', data, '_macro.rds'))
    
    used <- which(grepl('GSE140145', names(degs_ls))|
                    grepl('GSE116414', names(degs_ls)))
    degs_ls_used <- degs_ls[used]
    result_mal <- getallmetrics(methods_res_mal, degs_ls_used)
    result_mal1 <- getallrecord(methods_res_mal, degs_ls_used)
    saveRDS(result_mal1, file = paste0('./Step8_LRTBenchResult/record/', data, '_mal.rds'))
    
    result_macro$celltype <- 'macrophages'
    result_mal$celltype <- 'malignant'
    result <- rbind(result_macro, result_mal)
    
  }
  result
})
names(result_index) <- names(all_result)
saveRDS(result_index,file = './Step8_LRTBenchResult/result_index.rds')
