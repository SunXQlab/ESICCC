imputation.na <- function(result){
  tmp.result <- lapply(result, function(res){
    mutinfo <- res[,c("close.mutinfo", "distant.mutinfo", "lr")]
    mutinfo <- mutinfo[!apply(is.na(mutinfo[,c(1,2)]),1, all), ]
    mutinfo$close.mutinfo[is.na(mutinfo$close.mutinfo)] <- 0
    mutinfo$distant.mutinfo[is.na(mutinfo$distant.mutinfo)] <- 0
    
    pcc <- res[, c("close.pcc", "distant.pcc", "close.pcc.p", "distant.pcc.p","lr")]
    pcc <- pcc[!apply(is.na(pcc[,c(1:4)]),1, all), ]
    pcc$close.pcc[is.na(pcc$close.pcc)] <- 0
    pcc$distant.pcc[is.na(pcc$distant.pcc)] <- 0
    pcc$close.pcc.p[is.na(pcc$close.pcc.p)] <- 0
    pcc$distant.pcc.p[is.na(pcc$distant.pcc.p)] <- 0
    
    res <- merge(mutinfo, pcc, by = "lr")
    res
  })
  tmp.result
}
remove.auto <- function(result){
  remove <- lapply(names(result), function(cp){
    ct <- unique(unlist(stringr::str_split(cp, "_")))
    if(length(ct)==1){
      logic <- TRUE
    }else{
      logic <- FALSE
    }
    logic
  })
  remove <- unlist(remove)
  result[remove] <- NULL
  result
}
remove.null <- function(result){
  tmp.logic <- lapply(result, function(res){
    if(dim(res)[1]==0){
      logic <- TRUE
    }else{
      logic <- FALSE
    }
  })
  result[unlist(tmp.logic)] <- NULL
  result
}