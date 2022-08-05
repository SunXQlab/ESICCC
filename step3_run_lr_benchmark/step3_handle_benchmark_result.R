rm(list = ls())
library(tidyr)
set.seed(123)
setwd('./step3_run_lr_benchmark/result/')

## bench.st.handle.result.df.RData
if(T){
  source('../code/bench.handle.R')
  files <- list.files(".")
  for(file in files){
    file.path <- paste0("./", file)
    load(file.path)
    results.names <- objects(pattern = "^bench")
    if(substring(results.names, 14) %in% c(4, 5, 13, 16)){
      results <- get(results.names)
      results <- results$`0.4`
      result.names <- names(results)
      for (res.name in result.names){
        fin.result <- list()
        result <- results[[res.name]]
        result <- imputation.na(result) %>% remove.auto() %>% 
          remove.null()
        result <- do.call(rbind, result)
        result <- tibble::rownames_to_column(result, var = "sr")
        result$sr <- stringr::str_replace_all(result$sr, "\\.[0-9]+", "")
        result <- tidyr::separate(result, sr, c("sender", "reciever"), "_", remove = FALSE)
        tool <- paste(substring(file, gregexpr('_',file)[[1]][1]+1, nchar(file)-6), res.name, sep = "_")
        result$tool <- tool
        
        fin.result[["close"]] <- result[,c(1:5, 7, 9, 11)]
        colnames(fin.result[["close"]])[5:7] <- c("mutinfo", "pcc", "pcc.p")
        fin.result[["distant"]] <- result[,c(1:4, 6,8,10,11)]
        colnames(fin.result[["distant"]])[5:7] <- c("mutinfo", "pcc", "pcc.p")
        fin.result <- do.call(rbind, fin.result)
        fin.result <- tibble::rownames_to_column(fin.result, var = "groups")
        fin.result$groups <- stringr::str_replace_all(fin.result$groups, "\\.[0-9]+", "")
        assign(paste0("handle.result.", tool), fin.result)
      }
    }else{
      fin.result <- list()
      result <- get(results.names)
      result <- result$`0.4`
      result <- imputation.na(result) %>% remove.auto() %>% 
        remove.null()
      
      result <- do.call(rbind, result)
      result <- tibble::rownames_to_column(result, var = "sr")
      result$sr <- stringr::str_replace_all(result$sr, "\\.[0-9]+", "")
      result <- tidyr::separate(result, sr, c("sender", "reciever"), "_", remove = FALSE)
      tool <- substring(file, gregexpr('_',file)[[1]][1]+1, nchar(file)-6)
      result$tool <- tool
      
      fin.result[["close"]] <- result[,c(1:5, 7, 9, 11)]
      colnames(fin.result[["close"]])[5:7] <- c("mutinfo", "pcc", "pcc.p")
      fin.result[["distant"]] <- result[,c(1:4, 6,8,10,11)]
      colnames(fin.result[["distant"]])[5:7] <- c("mutinfo", "pcc", "pcc.p")
      fin.result <- do.call(rbind, fin.result)
      fin.result <- tibble::rownames_to_column(fin.result, var = "groups")
      fin.result$groups <- stringr::str_replace_all(fin.result$groups, "\\.[0-9]+", "")
      assign(paste0("handle.result.", tool), fin.result)
    }
    
    rm(list = results.names)
  }
  rm(list = ls.str(mode = 'character'))
  rm(list = ls.str(mode = "numeric"))
  rm(list = ls.str(mode = "function"))
  rm(result, results, fin.result)
  save.image("~/CCC-benchmark/step3_run_lr_benchmark/result.bench/result.benchmark.tools.4.RData")
}

tmp <- data.frame()
for (res in objects(pattern = "^handle")) {
  result <- get(res)
  tmp <- rbind(tmp, result)
}
rm(result, res)
rm(list = objects(pattern = "^handle"))
result.benchmark <- tmp
save(result.benchmark, file = '../result.bench/result.benchmark.4.RData')
