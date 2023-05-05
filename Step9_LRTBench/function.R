getallmetrics <- function(methods_res, degs_ls){
  result <- lapply(methods_res, function(res){
    degs.tmp <- lapply(names(degs_ls), function(degs_names){
      
      key <- strsplit(degs_names,'_')[[1]][2]
      type <- ifelse(grepl('KO',degs_names),'Receptor','Ligand')
      
      res$regulon <- as.character(res$regulon)
      res$target <- as.character(res$target)
      res$value <- as.numeric(res$value)
      
      if(grepl('NicheNet', res$methods[1]) & type == 'Receptor'){
        if(key == 'AXL'){
          score <- res[which(res$regulon=='GAS6'),]
        }else if(key == 'CXCR4'){
          score <- res[which(res$regulon=='CXCL12'),]
        }else if(key == 'NRP1'){
          score <- res[which(res$regulon=='VEGFA'),]
        }else if(key == 'CSF1R'){
          score <- res[which(res$regulon=='CSF1'),]
        }else if(key == 'FGFR1'){
          score <- res[which(res$regulon=='FGF1'),]
        }else if(key == 'ALK'){
          score <- res[which(res$regulon %in% c('ALKAL1', 'ALKAL2')),]
          if(dim(score)[[1]]!=0){
            score <- aggregate(value~target, score, sum)
          }
        }
      }else{
        score <- res %>% filter(regulon == key & type == type)
      }     
      
      if(grepl('GSE181575',degs_names)){
        neg.genes <- setdiff(score$target, names(degs_ls[[degs_names]]))
        isDEGs <- rep(FALSE, length(neg.genes))
        names(isDEGs) <- neg.genes
        degs_ls_used <- append(degs_ls[[degs_names]], isDEGs)
      }else{
        degs_ls_used <- degs_ls[[degs_names]]
      }
      
      genes <- intersect(score$target, names(degs_ls_used))
      label <- degs_ls_used[genes]
      pred <- score$value[match(genes,score$target)]
      res.tmp <- get_evaluate_metrics(pred,label)
      res.tmp <- res.tmp$perf_metrics
      
      if(length(names(table(label))) == 2){
        AUPRCRatios <- res.tmp[2]/(as.numeric(table(label)['TRUE'])/length(label))
        names(AUPRCRatios) <- 'AUPRCRatios'
      }else{
        AUPRCRatios <- NA
        names(AUPRCRatios) <- 'AUPRCRatios'
      }
      res.tmp <- append(res.tmp, AUPRCRatios)
      res.tmp
      
    })
    names(degs.tmp) <- names(degs_ls)
    degs.tmp <- do.call(rbind, degs.tmp) %>% 
      as.data.frame() %>% 
      tibble::rownames_to_column(., 'celllines')
    degs.tmp
  })
  result <- do.call(rbind, result) %>% as.data.frame() %>% 
    tibble::rownames_to_column(., 'methods')
  result$methods <- gsub('\\.[0-9]+', '', result$methods)
  return(result)
}

getallrecord <- function(methods_res, degs_ls){
  result1 <- lapply(methods_res, function(res){
    degs.tmp <- lapply(names(degs_ls), function(degs_names){
      
      key <- strsplit(degs_names,'_')[[1]][2]
      type <- ifelse(grepl('KO',degs_names),'Receptor','Ligand')
      
      res$regulon <- as.character(res$regulon)
      res$target <- as.character(res$target)
      res$value <- as.numeric(res$value)
      
      if(grepl('NicheNet', res$methods[1]) & type == 'Receptor'){
        if(key == 'AXL'){
          score <- res[which(res$regulon=='GAS6'),]
        }else if(key == 'CXCR4'){
          score <- res[which(res$regulon=='CXCL12'),]
        }else if(key == 'NRP1'){
          score <- res[which(res$regulon=='VEGFA'),]
        }else if(key == 'CSF1R'){
          score <- res[which(res$regulon=='CSF1'),]
        }else if(key == 'FGFR1'){
          score <- res[which(res$regulon=='FGF1'),]
        }else if(key == 'ALK'){
          score <- res[which(res$regulon %in% c('ALKAL1', 'ALKAL2')),]
          if(dim(score)[[1]]!=0){
            score <- aggregate(value~target, score, sum)
          }
        }
      }else{
        score <- res %>% filter(regulon == key & type == type)
      }     
      
      if(grepl('GSE181575',degs_names)){
        neg.genes <- setdiff(score$target, names(degs_ls[[degs_names]]))
        isDEGs <- rep(FALSE, length(neg.genes))
        names(isDEGs) <- neg.genes
        degs_ls_used <- append(degs_ls[[degs_names]], isDEGs)
      }else{
        degs_ls_used <- degs_ls[[degs_names]]
      }
      
      genes <- intersect(score$target, names(degs_ls_used))
      label <- degs_ls_used[genes]
      pred <- score$value[match(genes,score$target)]
      
      res.tmp <- cbind(label, pred) %>% as.data.frame() %>% tibble::rownames_to_column(., 'genes')
    })
    names(degs.tmp) <- names(degs_ls)
    degs.tmp
  })
  return(result1)
}

RunTimeMemRecord <- function(wd, data, methods){
  # get result path
  result.path <- paste0(wd, data)
  
  dataset_Record <- lapply(methods, function(method){
    print(method)
    
    record <- read.table(paste0(result.path, '/', method, '/', 'TimeMemRecord.txt'),
                         header = FALSE, sep = '\t', fill = TRUE)
    record <- record[, 2, drop = FALSE]
    if(dim(record)[1]==24){
      record <- record[-1, , drop = FALSE]
    }
    record <- record[c(2,3,4,5,10),, drop = FALSE]
    
    
    # Linux_time: user+system time
    user_time <- as.numeric(gsub('User time \\(seconds\\): ', '', record[1,]))
    sys_time <- as.numeric(gsub('System time \\(seconds\\): ', '', record[2,]))
    time_linux <- user_time+sys_time
    
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
      clock_time <- h_time+min_time+sec_time; rm(sec_time, min_time, h_time)
    }
    
    # Max memory usage | CPU_used
    mem_linux <- gsub('Maximum resident set size \\(kbytes\\): ', '', record[5,]) 
    mem_linux <- round(as.numeric(mem_linux)/1024/1024, 2)
    cpu_linux <- gsub('Percent of CPU this job got: ', '', record[3,])
    
    return(data.frame(clock_time = clock_time, linux_time = time_linux, 
                      max_memory = mem_linux, cpu_linux = cpu_linux))
  })
  names(dataset_Record) <- methods
  return(dataset_Record)
}