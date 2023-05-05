suppressMessages(library(Seurat))
source('../Script/Step8_LRTToolsFunction/Step8_LRTToolsFunction.R')
set.seed(123)

args <- commandArgs()
tools <- args[6]
print(tools)
sampleID <- args[7]
print(sampleID)

if(grepl('CID', sampleID)){
  fpath.dataset <- paste0('/home/ljx/1-Datasets/NatureGenetics_2021_BC/ScriptForCCC/ST_Data/', sampleID, '_ser.rds')
  receiver <- 'CancerEpithelial'
}else if(grepl('UKF', sampleID)){
  fpath.dataset <- paste0('/home/ljx/1-Datasets/CancerCell_2022_Glioma/ScriptForLRTBenchmark/ST_Data/', sampleID, '_ser.rds')
  receiver <- c('macrophages', 'Malignant')
}else{#MI
  fpath.dataset <- paste0('/home/ljx/1-Datasets/Nature_2022_MI/ScriptForCCC/ST_Data/', sampleID, '_ser.rds')
  receiver <- c('Fibroblast')
}

dir.output <- paste0('./Step8_LRTPredictionResult/', sampleID)
if(!dir.exists(dir.output)){
  dir.create(dir.output)
}

fpath.tools <- paste(dir.output, tools, sep = '/')
if(!dir.exists(fpath.tools)){
  dir.create(fpath.tools)
}

ser <- readRDS(fpath.dataset)

if(tools == 'CytoTalk'){
  CytoTalk_function(ser, fpath.tools, receiver = receiver)
}else if(tools == 'HoloNet'){
  HoloNet_function(ser, fpath.tools, sampleID, receiver = receiver)
}else if(tools == 'MISTy'){
  MISTy_function(ser, fpath.tools, sampleID, receiver = receiver)
}else if(tools == 'NicheNet'){
  NicheNet_function(ser, fpath.tools, sampleID, receiver = receiver)
}else if(tools == 'stMLnet'){
  stMLnet_function(ser, fpath.tools, sampleID, receiver = receiver)
}