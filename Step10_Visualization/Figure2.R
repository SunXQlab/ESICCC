rm(list = ls()); gc()
suppressMessages(library(tidyverse))
suppressMessages(library(ROCR))
suppressMessages(library(ggplotify))
suppressMessages(library(pheatmap))
suppressMessages(library(reshape2))
suppressMessages(library(ggplot2))
suppressMessages(library(ggrepel))
suppressMessages(library(RColorBrewer))
set.seed(123)

##################################### 
##                                 ## each dataset: median AUPRC
## LR Benchmark based on CAGE data ## metric based on all datasets: median AUPRC of all datasets
##                                 ##
#####################################

if(F){
  
  source('./Script/Step3_MIForLRBench/function.R')
  source('./Script/Step5_BenchBasedCAGEProteomic/function.R')
  
  cage_genes <- readRDS('./Data/Step5_BenchBasedCAGEProteomic/cagedata.rds')
  methods <- c("CellPhoneDB2", "CellPhoneDB3", "CellTalker", 
               "Connectome", "ICELLNET",  "NATMI", "iTALK", 
               "scConnect", "SingleCellSignalR", "CellChat", 
               "RNAMagnet", "PyMINEr",
               "scSeqComm", "NicheNet","CytoTalk",
               "scMLnet", "Domino", "CellCall")
  
  datasets <- c('pbmc4k', 'pbmc6k', 'pbmc8k')
  
  cage_eval <- lapply(datasets, function(data){
    fpath <- paste0('./Data/Step1_LRPredictionResult/', data)
    result.files <- list.files(fpath, full.names = TRUE)
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
    result_index <- tibble::rownames_to_column(result_index, 'temp')
    result_index <- tidyr::separate(result_index, temp, c('methods', 'sr'), '\\.')
    result_index$dataset <- data
    result_index <- result_index[which(result_index$methods %in% methods), ]
    result_index
  })
  names(cage_eval) <- datasets
  saveRDS(cage_eval, file = './Data/Step10_Visualization/cage_eval_updata.rds')
  
  rm(list = ls()); gc()
  
}

cage_eval <- readRDS('./Data/Step10_Visualization/cage_eval_updata.rds')
cage_prc_summary <- lapply(cage_eval, function(eval){
  summary <- aggregate(AUPRC~methods, eval, median)
  summary
})
cage_eval <- do.call(rbind, cage_eval)
rownames(cage_eval) <- NULL
cage_prc_summary[['all']] <- aggregate(AUPRC~methods, cage_eval, median)
rm(cage_eval); gc()


########################################## 
##                                      ## each dataset: median AUPRC
## LR Benchmark based on Proteomic data ## metric based on all datasets: median AUPRC of all datasets
##                                      ##
##########################################

if(F){
  
  source('./Script/Step3_MIForLRBench/function.R')
  source('./Script/Step5_BenchBasedCAGEProteomic/function.R')
  
  pep_genes <- readRDS('./Data/Step5_BenchBasedCAGEProteomic/pepdata.rds')
  
  datasets <- c('pbmc4k', 'pbmc6k', 'pbmc8k')
  methods <- c("CellPhoneDB2", "CellPhoneDB3", "CellTalker", 
               "Connectome", "ICELLNET",  "NATMI", "iTALK", 
               "scConnect", "SingleCellSignalR", "CellChat", 
               "RNAMagnet", "PyMINEr",
               "scSeqComm", "NicheNet","CytoTalk",
               "scMLnet", "Domino", "CellCall")
  
  
  proteomic_eval <- lapply(datasets, function(data){
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
    result_index <- tibble::rownames_to_column(result_index, 'temp')
    result_index <- tidyr::separate(result_index, temp, c('methods', 'sr'), '\\.')
    result_index <- result_index[which(result_index$methods %in% methods), ]
    result_index$dataset <- data
    result_index
  })
  
  saveRDS(proteomic_eval, file = './Data/Step10_Visualization/proteomic_eval_updata.rds')
  
  rm(list = ls()[-which(ls() %in% 'cage_prc_summary')]); gc()
}

proteomic_eval <- readRDS('./Data/Step10_Visualization/proteomic_eval_updata.rds')
proteomic_prc_summary <- lapply(proteomic_eval, function(eval){
  summary <- aggregate(AUPRC~methods, eval, median)
  summary
})
names(proteomic_prc_summary) <- c('pbmc4k', 'pbmc6k', 'pbmc8k')
proteomic_eval <- do.call(rbind, proteomic_eval)
proteomic_prc_summary[['all']] <- aggregate(AUPRC~methods, proteomic_eval, median)
rm(proteomic_eval); gc()


############################## 
##                          ## each dataset: one dataset, one MIIndex
## LR Benchmark based on MI ## metric based on all datasets: 1, scaled from 0 to 1; 2, mean of all datasets
##                          ##
##############################

if(F){
  source('./Script/Step3_MIForLRBench/function.R')
  
  result.path <- list.files('./Data/Step3_MIForLRBench', pattern = '_result.rds', full.names = TRUE)
  samples <- gsub('_result.rds', '', substring(result.path, 30))
  
  tmp <- lapply(seq(result.path), function(i){
    print(i)
    Eval1Result <- Eval1Process(result.path[i])
    Eval1Result <- EvalIndex_DLRC(Eval1Result)
    Eval1Result
  })
  names(tmp) <- samples
  saveRDS(tmp, file = './Data/Step10_Visualization/MI_eval.rds')
  rm(list = ls()[-which(ls()%in%c('cage_prc_summary', 'proteomic_prc_summary'))]);gc()
}

mi_eval <- readRDS('./Data/Step10_Visualization/MI_eval.rds')
mi_summary <- lapply(mi_eval, function(eval){
  eval <- eval[, c('methods','MIIndex')]
  eval
})
mi_eval <- do.call(rbind, mi_summary)
mi_eval <- tibble::rownames_to_column(mi_eval, 'datasets')
mi_eval$datasets <- gsub('\\.[0-9]+', '', mi_eval$datasets)
mi_eval <- spread(mi_eval, 'datasets', 'MIIndex')
mi_eval <- tibble::column_to_rownames(mi_eval, 'methods')
mi_eval <- apply(mi_eval, 2, 
               function(x){(x-min(x, na.rm = TRUE))/(max(x, na.rm = TRUE)-min(x, na.rm = TRUE))})
mi_eval_summary <- data.frame(methods = rownames(mi_eval), 
                              MIIndex = rowMeans(mi_eval, na.rm = TRUE))
rownames(mi_eval_summary) <- NULL
mi_summary[['all']] <- mi_eval_summary
rm(mi_eval, mi_eval_summary); gc()


#################################### 
##                                ## each dataset: one dataset, one Stabilty
## LR Benchmark based on Sampling ## metric based on all datasets: 1, scaled from 0 to 1; 2, mean of all datasets
##                                ##
####################################

result <- readRDS('./Data/Step6_LRBenchSampling/JaccardIndex.rds')
result$stability <- lapply(seq(dim(result)[[1]]), function(x){
  tmp <- 1/sum(1 - result$S50[x], 1 - result$S60[x], 1 - result$S70[x], 
             1 - result$S80[x], 1 - result$S90[x])
  tmp
})%>% unlist

result <- spread(result[, c(1,7,8)], 'datasets', 'stability')
result <- tibble::column_to_rownames(result, 'methods')
result <- apply(result, 2, 
                 function(x){(x-min(x, na.rm = TRUE))/(max(x, na.rm = TRUE)-min(x, na.rm = TRUE))})
sampling_summary <- data.frame(methods = rownames(result), 
                              Stability = rowMeans(result, na.rm = TRUE))
rownames(sampling_summary) <- NULL

rm(result);gc()

#################################
##                             ## each datasets: score by defined score
## Running time and Max memory ## metric based on all datasets: 1, scaled from 0 to 1; 2, mean of all datasets
##                             ## Usability: mean(score_time, score_memory)
################################# 


TimeMenRecord <- readRDS('./Data/Step7_LRBenchSampling/TimeMemRecord.rds')
temp <- TimeMenRecord[, c(2,4,5,7,8)]
temp$linux_time <- temp$linux_time/60
temp$Score_time <- lapply(temp$linux_time, function(time){
  if(time <= 1){ score <- 1 
  }else if(time > 1 & time <= 5){ score <- 0.9 
  }else if(time > 5 & time <= 10){ score <- 0.8 
  }else if(time > 10 & time <= 30){ score <- 0.7 
  }else if(time > 30 & time <= 60){ score <- 0.6 
  }else if(time >60 & time <= 720){ score <- 0.5
  }else if(time>720 & time <= 1440){ score <- 0.4
  }else if(time >1440 & time <= 2880){ score <- 0.3 
  }else if(time > 2880 & time <=4320){ score <- 0.2
  }else if(time > 4320){ score <- 0.1 }
  score
})%>%unlist()

temp$Score_memory <- lapply(temp$max_memory, function(memory){
  if(memory <= 1){ score <- 1 
  }else if(memory > 1 & memory <= 5){ score <- 0.9 
  }else if(memory > 5 & memory <= 10){ score <- 0.8 
  }else if(memory > 10 & memory <= 15){ score <- 0.7 
  }else if(memory > 15 & memory <= 20){ score <- 0.6 
  }else if(memory > 20 & memory <= 25){ score <- 0.5
  }else if(memory > 25 & memory <= 30){ score <- 0.4
  }else if(memory > 30 & memory <= 35){ score <- 0.3 
  }else if(memory > 35 & memory <= 40){ score <- 0.2
  }else if(memory > 40){ score <- 0.1 }
  score
})%>%unlist()
temp$class <- paste(temp$datasets, temp$ratios, sep = '_')

Time_record <- spread(temp[, c(1,6,8)], 'class', 'Score_time', fill = NA)
Time_record <- tibble::column_to_rownames(Time_record, 'methods')
Time_record <- apply(Time_record, 2, 
                 function(x){(x-min(x, na.rm = TRUE))/(max(x, na.rm = TRUE)-min(x, na.rm = TRUE))})
Time_record <- data.frame(methods = rownames(Time_record), 
                          Score_time = rowMeans(Time_record,na.rm = TRUE))

Mem_record <- spread(temp[, c(1,7,8)], 'class', 'Score_memory', fill = NA)
Mem_record <- tibble::column_to_rownames(Mem_record, 'methods')
Mem_record <- apply(Mem_record, 2, 
                     function(x){(x-min(x, na.rm = TRUE))/(max(x, na.rm = TRUE)-min(x, na.rm = TRUE))})
Mem_record <- data.frame(methods = rownames(Mem_record), 
                          Score_memory = rowMeans(Mem_record,na.rm = TRUE))
TimeMen_summary <- merge(Time_record, Mem_record, by = 'methods')
TimeMen_summary$usability <- rowMeans(TimeMen_summary[,2:3])

rm(list = ls()[-which(ls() %in% objects(pattern = '_summary'))]);gc()

######################
## Merge all Result ##
######################

cage_tmp <- do.call(rbind, cage_prc_summary) 
cage_tmp <- tibble::rownames_to_column(cage_tmp, 'datasets')
cage_tmp$datasets <- gsub('\\.[0-9]+', '', cage_tmp$datasets)
cage_tmp <- spread(cage_tmp, 'datasets', 'AUPRC')
colnames(cage_tmp)[2:5] <- paste0(colnames(cage_tmp)[2:5], '_CAGE')
rm(cage_prc_summary);gc()

proteomic_tmp <- do.call(rbind, proteomic_prc_summary)
proteomic_tmp <- tibble::rownames_to_column(proteomic_tmp, 'datasets')
proteomic_tmp$datasets <- gsub('\\.[0-9]+', '', proteomic_tmp$datasets)
proteomic_tmp <- spread(proteomic_tmp, 'datasets', 'AUPRC')
colnames(proteomic_tmp)[2:5] <- paste0(colnames(proteomic_tmp)[2:5], '_Proteomic')
rm(proteomic_prc_summary);gc()

mi_tmp <- do.call(rbind, mi_summary)
mi_tmp <- tibble::rownames_to_column(mi_tmp, 'datasets')
mi_tmp$datasets <- gsub('\\.[0-9]+', '', mi_tmp$datasets)
mi_tmp <- spread(mi_tmp, 'datasets', 'MIIndex')
colnames(mi_tmp)[2:13] <- paste0(colnames(mi_tmp)[2:13], '_MI')
rm(mi_summary);gc()

sampling_tmp <- sampling_summary
rm(sampling_summary);gc()

TimeMen_tmp <- TimeMen_summary
rm(TimeMen_summary); gc()

all_summary <- merge(cage_tmp, proteomic_tmp, by = 'methods', all = TRUE)
all_summary <- merge(all_summary, mi_tmp, by = 'methods', all = TRUE)
all_summary <- merge(all_summary, sampling_tmp, by = 'methods', all = TRUE)
all_summary <- merge(all_summary, TimeMen_tmp, by = 'methods', all = TRUE)
rm(list = objects(pattern = '_tmp')); gc()
saveRDS(all_summary, file = './Data/Step10_Visualization/all_summary_updata.rds')

###################
## Visualization ##
###################

all_summary <- readRDS('./Data/Step10_Visualization/all_summary_updata.rds')
all_summary <- tibble::column_to_rownames(all_summary, 'methods')
all_summary <- as.matrix(all_summary)
all_summary <- all_summary[, c(1,5,9,21:24)]
all_summary <- apply(all_summary, 2, 
                     function(x){(x-min(x, na.rm = TRUE))/(max(x, na.rm = TRUE)-min(x, na.rm = TRUE))})
all_summary <- as.data.frame(all_summary)
all_summary$overallForBench <- rowMeans(all_summary[, 1:3])
all_summary$overallForBench <- (all_summary$overallForBench-min(all_summary$overallForBench, na.rm = TRUE))/(max(all_summary$overallForBench, na.rm = TRUE)-min(all_summary$overallForBench, na.rm = TRUE))
all_summary$overall <- rowMeans(all_summary[, c(4,7,8)], na.rm = TRUE)

  ###################
  ## based on rank ##
  ###################

temp1 <- apply(all_summary, 2, function(x){rank(-x, na.last = 'keep')})
temp1 <- as.data.frame(temp1)
temp1 <- as.matrix(temp1)
temp1 <- melt(temp1)
colnames(temp1) <- c('methods', 'Index', 'Rank')
temp1$class <- paste(temp1$methods, temp1$Index, sep = '-')
temp1 <- temp1[, -c(1:2)]

  ####################
  ## based on score ##
  ####################

temp2 <- all_summary
temp2 <- as.matrix(temp2)
temp2 <- melt(temp2)
colnames(temp2) <- c('methods', 'Index', 'Score')
temp2$class <- paste(temp2$methods, temp2$Index, sep = '-')
temp2 <- temp2[, -c(1:2)]

  ###########
  ## merge ##
  ###########
temp <- merge(temp1, temp2, by = 'class', all = TRUE)
temp <- separate(temp, class, c('methods', 'Index'), sep = '-')
temp$methods <- as.character(temp$methods) %>%
  factor(., levels = rev(c("CellPhoneDB2", "CellPhoneDB3", "CellTalker", 
                           "Connectome", "ICELLNET",  "NATMI", "iTALK", 
                           "scConnect", "SingleCellSignalR", "CellChat", 
                           "RNAMagnet", "PyMINEr",
                           "scSeqComm", "NicheNet","CytoTalk",
                           "scMLnet", "Domino", "CellCall")))


  ##########
  ## draw ##
  ##########

  #######################################
  ## all_CAGE | all_Proteomic | all_MI ##
  #######################################

tmp1 <- temp[which(temp$Index %in% c("all_CAGE", "all_Proteomic", "all_MI")), ]
tmp1$Index <- as.character(tmp1$Index) %>%
  factor(., levels = c("all_MI", "all_CAGE", "all_Proteomic"))

cols <- brewer.pal(9, "Blues")[c(1,2, 7,8,9)] #Blues Purples
pal <- colorRampPalette(cols)
mycolors <- pal(9) %>% rev(.)
scales::show_col(mycolors)
ggplot(tmp1, aes(x = Index, y = methods, size = Rank, fill=Score)) + geom_point(colour = "black", shape = 21)+
  scale_size(range = c(18, 3))+
  scale_fill_gradientn(colours  = mycolors)+
  geom_text_repel(data = tmp1, aes(label = Rank),
    size = 5, segment.color = "black", show.legend = FALSE )

  ###############
  ## Stability ##
  ###############

tmp2 <- temp[which(temp$Index %in% c('Stability')), ]
cols <- brewer.pal(9, "Oranges")[-c(8:9)] #Blues 
pal <- colorRampPalette(cols)
mycolors <- pal(9) %>% rev(.)
scales::show_col(mycolors)
ggplot(tmp2, aes(x = Index, y = methods, size = Rank, fill=Score)) + 
  geom_point(colour = "black", shape = 21)+
  scale_size(range = c(18, 3))+
  scale_fill_gradientn(colours  = mycolors)+
  geom_text_repel(data = tmp2, aes(label = Rank),
                  size = 5, segment.color = "black", show.legend = FALSE )

  ###############
  ## Usability ##
  ###############

tmp3 <- temp[which(temp$Index %in% c("Score_time", "Score_memory")), ]
cols <- brewer.pal(9, "YlOrBr")[c(1:4, 6, 7, 8)] #Blues 
pal <- colorRampPalette(cols)
mycolors <- pal(9) %>% rev(.)
scales::show_col(mycolors)
ggplot(tmp3, aes(x = Index, y = methods, size = Rank, fill=Score)) + 
  geom_point(colour = "black", shape = 21)+
  scale_size(range = c(18, 3))+
  scale_fill_gradientn(colours  = mycolors)+
  geom_text_repel(data = tmp3, aes(label = Rank),
                  size = 5, segment.color = "black", show.legend = FALSE )

  #####################
  ## overallForBench ##
  #####################
tmp4 <- temp[which(temp$Index %in% c('overallForBench')), ]
cols <- brewer.pal(9, "Blues")[c(1,2, 7,8,9)] #Blues Purples
pal <- colorRampPalette(cols)
mycolors <- pal(9) %>% rev(.)
scales::show_col(mycolors)
ggplot(tmp4, aes(x = Index, y = methods, size = Rank, fill=-Rank)) + 
  geom_point(colour = "black", shape = 22)+
  scale_size(range = c(18, 5))+
  scale_fill_gradientn(colours  = mycolors)+
  geom_text_repel(data = tmp4, aes(label = Rank),
                  size = 5, segment.color = "black", show.legend = FALSE )

  #########################
  ## Usability (Overall) ##
  #########################
tmp5 <- temp[which(temp$Index %in% c('usability')), ]
cols <- brewer.pal(9, "YlOrBr")[c(1:4, 6, 7, 8)]
pal <- colorRampPalette(cols)
mycolors <- pal(9) %>% rev(.)
scales::show_col(mycolors)
ggplot(tmp5, aes(x = Index, y = methods, size = Rank, fill=-Rank)) + 
  geom_point(colour = "black", shape = 22)+
  scale_size(range = c(18, 5))+
  scale_fill_gradientn(colours  = mycolors)+
  geom_text_repel(data = tmp5, aes(label = Rank),
                  size = 5, segment.color = "black", show.legend = FALSE )

  #############
  ## Overall ##
  #############
tmp6 <- temp[which(temp$Index %in% c('overall')), ]
cols <- brewer.pal(9, "Purples")
pal <- colorRampPalette(cols)
mycolors <- pal(9) %>% rev(.)
scales::show_col(mycolors)
ggplot(tmp6, aes(x = Index, y = methods, size = Rank, fill=-Rank)) + 
  geom_point(colour = "black", shape = 22)+
  scale_size(range = c(18, 5))+
  scale_fill_gradientn(colours  = mycolors)+
  geom_text_repel(data = tmp6, aes(label = Rank),
                  size = 5, segment.color = "black", show.legend = FALSE )