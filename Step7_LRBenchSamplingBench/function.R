RunTimeMemRecord <- function(data, ratio, methods, DataRecord){
  # get result path
  if(ratio==100){
    if(grepl('CK', data)){
      result.path <- paste0('./Data/Step1_LRPredictionResult/N_MI_', data)
    }else if(grepl('CID', data)){
      result.path <- paste0('./Data/Step1_LRPredictionResult/NG_BC_', data)
    }else if(grepl('Slide14', data)){
      result.path <- paste0('./Data/Step1_LRPredictionResult/MouseEmbryo_', data)
    }else{
      result.path <- paste0('./Data/Step1_LRPredictionResult/', data)
    }
  }else{
    result.path <- paste0('./Data/Step10_LRBenchSamplingResult/', data, '_', ratio)
  }
  
  dataset_Record <- lapply(methods, function(method){
    print(method)
      
    if(method == 'CytoTalk'){
      record_dir <- paste0(result.path, '/', method)
      record_files <- list.files(record_dir, pattern = 'TimeMemRecord', full.names = TRUE)
      temp_record <- lapply(record_files, function(file){
        record <- tryCatch(read.table(file, header = FALSE, sep = '\t', fill = TRUE),
                        error=function(e){NA}
        )
        if(class(record)=='data.frame'){
          tmp_record <- RunTimeMemRecord_2(record)
        }else{
          tmp_record <- NA
        }
        tmp_record
      })
      temp_record[which(is.na(temp_record))] <- NULL
      temp_record <- do.call(rbind, temp_record) %>% as.data.frame()
      
      if(ratio == 100){
        number_record <- DataRecord$Celltypes[which((DataRecord$ratio == ratio) & (DataRecord$datasets==data))]
        actual_cp <- number_record * number_record - number_record - number_record
        combn_cp <- dim(combn(seq(number_record), 2))[2]
        temp_record$clock_time <- (temp_record$clock_time/actual_cp)*combn_cp
        temp_record$linux_time <- (temp_record$linux_time/actual_cp)*combn_cp
      }
      
      temp_record <- data.frame(clock_time = sum(temp_record$clock_time), 
                                linux_time = sum(temp_record$linux_time), 
                                max_memory = max(temp_record$max_memory), 
                                cpu_linux = paste0(round(mean(as.numeric(gsub('%', '', temp_record$cpu_linux))), 0), '%'))
     }else{
       record <- read.table(paste0(result.path, '/', method, '/', 'TimeMemRecord.txt'),
                             header = FALSE, sep = '\t', fill = TRUE)
       temp_record <- RunTimeMemRecord_2(record)
      }
    
    temp_record
  })
  names(dataset_Record) <- methods
  return(dataset_Record)
}


RunTimeMemRecord_2 <- function(record){
  
  record <- record[, 2, drop = FALSE]
  if(dim(record)[1]==24){
    record <- record[-1, , drop = FALSE]
  }
  record <- record[c(2,3,4,5,10),, drop = FALSE]
  
  
  # Linux_time: user+system time
  user_time <- as.numeric(gsub('User time \\(seconds\\): ', '', record[1,])) %>% as.numeric()
  sys_time <- as.numeric(gsub('System time \\(seconds\\): ', '', record[2,])) %>% as.numeric()
  linux_time <- user_time+sys_time
  
  # Clock_time
  clock_time <- gsub('Elapsed \\(wall clock\\) time \\(h:mm:ss or m:ss\\): ', '', record[4,])
  if(stringr::str_count(clock_time, ':') == 1){
    min_time <- as.numeric(unlist(stringr::str_split(clock_time, ':'))[1])
    min_time <- min_time*60
    sec_time <- as.numeric(unlist(stringr::str_split(clock_time, ':'))[2])
    clock_time <- min_time+sec_time; rm(sec_time, min_time)
  }else if(stringr::str_count(clock_time, ':') == 2){
    h_time <- as.numeric(unlist(stringr::str_split(clock_time, ':'))[1])
    h_time <- h_time*60*60
    min_time <- as.numeric(unlist(stringr::str_split(clock_time, ':'))[2])
    min_time <- min_time*60
    sec_time <- as.numeric(unlist(stringr::str_split(clock_time, ':'))[3])
    clock_time <- h_time+min_time+sec_time; 
    rm(sec_time, min_time, h_time);gc()
  }
  
  mem_linux <- gsub('Maximum resident set size \\(kbytes\\): ', '', record[5,]) 
  mem_linux <- round(as.numeric(mem_linux)/1024/1024, 2)
  cpu_linux <- gsub('Percent of CPU this job got: ', '', record[3,])
  
  return(data.frame(clock_time = clock_time, linux_time = linux_time, max_memory = mem_linux, cpu_linux = cpu_linux))
}
