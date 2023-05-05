rm(list = ls());gc()
################
## TimeMemory ##
################

TimeMemRecord <- readRDS('./Data/Step9_LRBench/LRT_RunTimeMem.rds') 
temp <- TimeMemRecord
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


Time_record <- spread(temp[, c(1,2,7)], 'datasets', 'Score_time', fill = NA)
Time_record <- tibble::column_to_rownames(Time_record, 'methods')
Time_record <- apply(Time_record, 2, 
                     function(x){(x-min(x, na.rm = TRUE))/(max(x, na.rm = TRUE)-min(x, na.rm = TRUE))})
Time_record <- data.frame(methods = rownames(Time_record), 
                          Score_time = rowMeans(Time_record,na.rm = TRUE))

Mem_record <- spread(temp[, c(1,2,8)], 'datasets', 'Score_memory', fill = NA)
Mem_record <- tibble::column_to_rownames(Mem_record, 'methods')
Mem_record <- apply(Mem_record, 2, 
                    function(x){(x-min(x, na.rm = TRUE))/(max(x, na.rm = TRUE)-min(x, na.rm = TRUE))})
Mem_record <- data.frame(methods = rownames(Mem_record), 
                         Score_memory = rowMeans(Mem_record,na.rm = TRUE))
TimeMen_summary <- merge(Time_record, Mem_record, by = 'methods')
TimeMen_summary <- tibble::column_to_rownames(TimeMen_summary, 'methods')
TimeMen_summary <- apply(TimeMen_summary, 2, 
                    function(x){(x-min(x, na.rm = TRUE))/(max(x, na.rm = TRUE)-min(x, na.rm = TRUE))})
TimeMen_summary <- as.data.frame(TimeMen_summary)
TimeMen_summary$usability <- rowMeans(TimeMen_summary[,1:2])
TimeMen_summary <- tibble::rownames_to_column(TimeMen_summary, 'methods')

rm(list = ls()[-which(ls()%in%c('TimeMen_summary'))])

###########
## AUPRC ##
###########

result_index <- readRDS('./Data/Step9_LRTBenchResult/result_index.rds')
result_index <- do.call(rbind, result_index)
result_index <- tibble::rownames_to_column(result_index, 'datasets')
result_index$datasets <- gsub('\\.[0-9]+', '', result_index$datasets)
result_index <- result_index[which(result_index$methods %in% 
                                     c('CytoTalk', "NicheNet", "MISTy", "HoloNet", "stMLnet")), ]
tmp <- result_index[result_index$datasets %in% 
                      c("CID4465", "CID44971", "control_P7", "control_P8", 
                        'UKF260_T_ST', 'UKF266_T_ST', 'UKF334_T_ST', 'UKF243_T_ST'), ]
tmp <- tmp[which(!grepl('TAM_RB', tmp$celllines)), ]
tmp <- aggregate(AUPRC~methods, tmp, FUN = function(x){median(x, na.rm = TRUE)})
Evaluation_summary <- tmp
rm(list = ls()[-which(ls()%in%c('TimeMen_summary', 'Evaluation_summary'))])

###########
## Merge ##
###########
tmp <- merge(Evaluation_summary, TimeMen_summary, by = 'methods')
tmp <- tibble::column_to_rownames(tmp, 'methods')
tmp <- apply(tmp, 2, 
             function(x){(x-min(x, na.rm = TRUE))/(max(x, na.rm = TRUE)-min(x, na.rm = TRUE))})
tmp <- as.data.frame(tmp)
tmp$overall <- rowMeans(tmp[, c(1,4)], na.rm = TRUE)
tmp <- tibble::rownames_to_column(tmp, 'methods')
tmp <- melt(tmp)
colnames(tmp) <- c('methods', 'Index', 'Score')

cols <- brewer.pal(9, "Blues")[c(1,2, 7,8,9)] #Blues Purples
cols <- brewer.pal(9, "YlOrBr")[c(1:4, 6, 7, 8)] #Blues 
cols <- brewer.pal(9, "Purples") 
pal <- colorRampPalette(cols)
mycolors <- pal(9) %>% rev(.)
scales::show_col(mycolors)
ggplot(tmp, aes(x = Index, y = methods, size = 100, fill=Score)) + geom_point(colour = "black", shape = 21)+
  scale_fill_gradientn(colours  = mycolors)
 