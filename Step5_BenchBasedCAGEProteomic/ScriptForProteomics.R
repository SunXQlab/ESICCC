rm(list = ls());gc()
suppressMessages(library(data.table))
suppressMessages(library(tidyverse))
suppressMessages(library(ROCR))
source('./Script/Step3_MIForLRBench/function.R')
source('./Script/Step8_BenchBasedCAGEProteomic/function.R')
set.seed(123)

# process proteomics data
if(F){
  proteinGroup <- fread('~/1-Datasets/Proteomics_28/Data/proteinGroups.txt')
  pep <- proteinGroup[, c(7,549:816)]
  rm(proteinGroup); gc()
  
  pep <- pep[which(pep$Gene.names != ''), ]
  # pep
  if(F){
    colnames(pep) <- gsub('Peptides.', '', colnames(pep))
    pep <- as.data.frame(pep)
    pep <- pep[, c(1, which(grepl('steady-state', colnames(pep))))]
    pep <- pep[, which(!grepl('Library_single', colnames(pep)))]
    colnames(pep) <- gsub(' \\(steady-state\\)', '', colnames(pep))
    
    colnames(pep)[1] <- 'genes'
    pep <- tidyr::separate_rows(pep, genes, sep = ';')
    pep <- aggregate(.~genes, data = pep, FUN=mean) %>%
      tibble::column_to_rownames(., 'genes') %>% 
      as.matrix() %>% t() %>% as.data.frame()
    
    pep$celltype <- gsub('_[0-9]+', '', rownames(pep))
    pep <- aggregate(.~celltype, data = pep, FUN = mean) %>%
      tibble::column_to_rownames(., 'celltype') %>%
      as.matrix() %>% t() %>% as.data.frame()
    
    celltype <- c('MO.classical', 'MO.nonclassical', 'mDC', 'pDC',
                  'B.naive', 'B.memory',
                  'T8.naive', 'T8.CM', 'T8.EM', 'T8.EMRA',
                  'T4.CM', 'T4.EM', 'T4.naive',
                  'NK.dim', 'NK.bright', 'Thrombocyte'
    )
    
    pep <- pep[ , which(colnames(pep) %in% celltype)]
    
    pep <- as.matrix(pep) %>% t() %>% as.data.frame()
    pep$celltype <- c('B', 'B', 'mDC', 'CD14Mono', 'CD16Mono',
                      'NK', 'NK', 'pDC', 'memoryCD4T', 'memoryCD4T',
                      'naiveCD4T', 'CD8T', 'CD8T', 'CD8T', 'CD8T', 'Platelet')
    pep <- aggregate(.~celltype, data = pep, FUN = mean) %>%
      tibble::column_to_rownames(., 'celltype') %>%
      as.matrix() %>% t() %>% as.data.frame()
    
    pep_genes <- lapply(1:dim(pep)[[2]], function(i){
      x <- pep[,i]
      genes <- rownames(pep)[which(x>=2)]
      paste(genes, colnames(pep)[i], sep = '_')
    })
    pep_genes <- unlist(pep_genes)
    pep_genes_list <- list(all = rownames(pep), genes_ct = pep_genes)
    saveRDS(pep_genes_list, 
            file = '~/2-CCC-benchmark-202204-202208/CCC-benchmark-202212/Data/Step5_BenchBasedCAGEProteomic/pepdata.rds')
  }
}

rm(list = ls());gc()
#calculate index -- peptides
if(F){
  pep_genes <- readRDS('./Data/Step8_BenchBasedCAGEProteomic/pepdata.rds')
  
  datasets <- c('pbmc4k', 'pbmc6k', 'pbmc8k')
  
  for (data in datasets) {
    fpath <- paste0('./Data/Step1_LRPredictionResult/', data)
    result.files <- list.files(fpath, full.names = TRUE)
    result.files <- result.files[-which(grepl('cellinker', result.files))]
    result_index <- lapply(result.files, function(result.file){
      print(result.file)
      if(grepl('scMLnet', result.file) | grepl('Domino', result.file)){
        result <- readRDS(paste0(result.file, '/result1.rds'))
      }else{
        result <- readRDS(paste0(result.file, '/result.rds'))
      }
      result <- result$result
      
      if((dim(result)[[1]] != 0) & 'LRscore' %in% colnames(result)){
        if(grepl('RNAMagnet', result.file)){
          genes <- readRDS(paste0('~/1-Datasets/10XGenomics_PBMC/ScriptForCellTypeAnn/', data, '_hs2mm.rds'))
          result <- mm2hs(result, genes)
        }
        
        result$receptor_ct <- lapply(1:dim(result)[[1]], function(i){
          receptor <- result$Receptor[i]
          if(grepl('&', receptor)){
            receptor <- unlist(strsplit(receptor, '&'))
          }
          receiver <- result$Receiver[i]
          tmp <- paste(receptor, receiver, sep = '_')
          label <- ifelse(all(tmp %in% pep_genes$genes_ct), TRUE, FALSE)
          if(!label){
            label <- ifelse(all(receptor %in% pep_genes$all), FALSE, NA)
          }
          label
        }) %>% unlist()
        
        result$ligand_ct <- lapply(1:dim(result)[[1]], function(i){
          ligand <- result$Ligand[i]
          if(grepl('&', ligand)){
            ligand <- unlist(strsplit(ligand, '&'))
          }
          sender <- result$Sender[i]
          tmp <- paste(ligand, sender, sep = '_')
          label <- ifelse(all(tmp %in% pep_genes$genes_ct), TRUE, FALSE)
          if(!label){
            label <- ifelse(all(ligand %in% pep_genes$all), FALSE, NA)
          }
          label
        }) %>% unlist()
        
        result <- result[which(!is.na(result$receptor_ct)), ]
        result <- result[which(!is.na(result$ligand_ct)),]
        
        result$label <- lapply(1:dim(result)[[1]], function(i){
          all(result$ligand_ct[i] & result$receptor_ct[i])
        }) %>% unlist()
        
        result$sr <- paste(result$Sender, result$Receiver, sep = '_')
        result <- result[, c('LRscore', 'label', 'sr')]
        result <- split(result, result$sr)
        
        res_index <- lapply(result, function(res){
          index <- get_evaluate_metrics(as.numeric(res$LRscore), res$label)
          index <- index$perf_metrics
        })
        res_index <- do.call(rbind, res_index)
        res_index <- as.data.frame(res_index)
        
      }else{
        res_index <- NA
      }
      res_index
    })
    
    names(result_index) <- substring(result.files, 40)
    result_index[which(is.na(result_index))] <- NULL
    result_index <- do.call(rbind, result_index)
    result_index <- result_index[which(rowSums(result_index)!=0), ]
    result_index <- tibble::rownames_to_column(result_index, 'temp')
    result_index <- tidyr::separate(result_index, temp, c('methods', 'sr'), '\\.')
    result_index$dataset <- data
    result_index
    saveRDS(result_index, file = paste0('./Data/Step5_BenchBasedCAGEProteomic/', data, '_pep.rds'))
  }
}
