rm(list = ls());gc()
source('./Script/Step3_MIForLRBench/function.R')
set.seed(123)

result.path <- list.files('./Data/Step3_MIPCCForLRBench', full.names = TRUE)
samples <- gsub('_result.rds', '', substring(result.path, 30))

tmp <- lapply(seq(result.path), function(i){
  Eval1Result <- Eval1Process(result.path[i])
  Eval1Result <- EvalIndex_DLRC(Eval1Result)
  Eval1Result
})
names(tmp) <- samples
tmp <- do.call(rbind, tmp)
tmp <- tibble::rownames_to_column(tmp, 'datasets')
tmp$datasets <- gsub('\\.[0-9]+', '', tmp$datasets)