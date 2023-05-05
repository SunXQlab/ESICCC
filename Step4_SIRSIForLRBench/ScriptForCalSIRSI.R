rm(list = ls());gc()
suppressMessages(library(Seurat))
source('./Script/Step3_MIForLRBench/function.R')
source('./Script/Step4_SIRSIForLRBench/function.R')
set.seed(123)

datasets <- c('CID4465', 'CID44971', 'CK357', 'CK358',
              'CK368', 'CK162', 'CK362', 'CK361',
              'CK161', 'CK165', 'Slide14') 
ck.dataset <- readRDS('~/1-Datasets/Nature_2022_MI/dataset_id.rds')

# Calculate SI
for(data in datasets){
  print(data)
  SIResult <- CalSI_2(data)
  saveRDS(SIResult, file = paste0('./Data/Step4_SIRSIForLRBench/', data, '_SIResult.rds'))
}
rm(SIResult, data);gc()

#Calculate RSI
for (data in datasets) {
  print(data)
  RSIResult <- CalRSI(data)
  saveRDS(RSIResult, file = paste0('./Data/Step4_SIRSIForLRBench/', data, '_RSIResult.rds'))
}
rm(RSIResult, data);gc()
