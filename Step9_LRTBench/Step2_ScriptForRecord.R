
files <- list.files('./Data/Step8_LRTBenchResult/record', full.names = TRUE)
record_data <- lapply(files, function(file){
  result <- readRDS(file)
  record_methods <- lapply(result, function(method){
    record_cellines <- lapply(method, function(cellines){
      if(dim(cellines)[[1]] == 0){
        text <- 'No rec/lig'
      }else{
        tmp <- as.data.frame(table(cellines$label))
        false <- as.numeric(tmp$Freq[which(tmp$Var1==0)])
        true <- as.numeric(tmp$Freq[which(tmp$Var1==1)])
        text <- paste0('F: ', false, '; T: ', true)
      }
      text
    })
    names(record_cellines) <- names(method)
    record_cellines
  })
  record_methods <- do.call(cbind, record_methods)
  record_methods <- as.data.frame(record_methods)
  record_methods
})
names(record_data) <- gsub('\\.rds', '', list.files('./Data/Step8_LRTBenchResult/record'))

record_data1 <- lapply(seq(5, 29, 2), function(i){
  print(i)
  marco <- record_data[[i]]
  mal <- record_data[[i+1]]
  marco.col <- setdiff(colnames(mal), colnames(marco))
  mal.col <- setdiff(colnames(marco), colnames(mal))
  
  if (length(marco.col)>0) {
    marco[, marco.col] <- '——'
    marco <- marco[, colnames(mal)]
  }else if(length(mal.col)>0){
    mal[, mal.col] <- '——'
    mal <- mal[, colnames(marco)]
  }
  tmp <- rbind(marco, mal)
  tmp
})

names(record_data1) <- gsub('_macro', '', names(record_data)[seq(5, 29, 2)])
