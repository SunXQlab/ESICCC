rm(list = ls());gc()
suppressMessages(library(ROCR))
suppressMessages(library(tidyverse))
source('./Script/Step5_BenchBasedCAGEProteomic/function.R')
source('./Script/Step3_MIForLRBench/function.R')
set.seed(123)

#process cage data 
if(F){
  celltype <- c('CD14+CD16- Monocytes', 'CD14-CD16+ Monocytes', 'CD19+ B cells', 
                'Dendritic Monocyte Immature derived', 'Dendritic Plasmacytoid',
                'CD4+CD25-CD45RA- memory conventional T cells', 
                'CD4+CD25-CD45RA+ naive conventional T cells',
                'CD8+ T cells', 'NK cells')
  cage_all <- data.table::fread('~/1-Datasets/Fantom5_CAGE/ExpressionGenes.txt', header = TRUE, sep = '\t')
  cage_all <- as.data.frame(cage_all)
  cage_all <- cage_all[ , c(1, which(colnames(cage_all) %in% celltype))]
  cage_all <- tibble::column_to_rownames(cage_all, 'ApprovedSymbol')
  colnames(cage_all) <- c('CD14Mono', 'CD16Mono', 'B', 
                          'naiveCD4T', 'memoryCD4T', 'CD8T',
                          'mDC', 'pDC', 'NK')
  cage_all <- as.matrix(cage_all)
  cage <- lapply(1:dim(cage_all)[[2]], function(i){
    x <- cage_all[,i]
    genes <- rownames(cage_all)[which(x>10)]
    paste(genes, colnames(cage_all)[i], sep = '_')
  })
  cage <- unlist(cage)
  cage_genes <- list(all = rownames(cage_all), genes_ct = cage)
  saveRDS(cage_genes, file = './Data/Step5_BenchBasedCAGEProteomic/cagedata.rds')
}

#calculate index
if(F){
  cage_genes <- readRDS('./Data/Step5_BenchBasedCAGEProteomic/cagedata.rds')
  
  datasets <- c('pbmc4k', 'pbmc6k', 'pbmc8k')
  
  for (data in datasets) {
    fpath <- paste0('./Data/Step1_LRPredictionResult/', data)
    result.files <- list.files(fpath, full.names = TRUE)
    result.files <- result.files[-which(grepl('cellinker', result.files))]
    result_index <- lapply(result.files, function(result.file){
      print(result.file)
      if(grepl('Domino', result.file) | grepl('scMLnet', result.file)){
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
          label <- ifelse(all(tmp %in% cage_genes$genes_ct), TRUE, FALSE)
          if(!label){
            label <- ifelse(all(receptor %in% cage_genes$all), FALSE, NA)
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
          label <- ifelse(all(tmp %in% cage_genes$genes_ct), TRUE, FALSE)
          if(!label){
            label <- ifelse(all(ligand %in% cage_genes$all), FALSE, NA)
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
    saveRDS(result_index, file = paste0('./Data/Step5_BenchBasedCAGEProteomic/', data, '_cage.rds'))
  }
}
