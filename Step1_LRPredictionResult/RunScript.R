suppressMessages(library(Seurat))
source('../Script/Step0_LRToolsFunction/Step1_LRToolsFunction.R')
set.seed(123)

args <- commandArgs()
tools <- args[6]
print(tools)
sampleID <- args[7]
print(sampleID)
species <- args[8]
print(species)

if(grepl('CID', sampleID) & tools == 'RNAMagnet'){
  fpath.dataset <- paste0('~/1-Datasets/NatureGenetics_2021_BC/ScriptForCCC/scRNA_Data_Mm/',
                          sampleID, '_ser.rds')
  dir.output <- paste0('~/2-CCC-benchmark-202204-202208/CCC-benchmark-202212/Data/Step1_LRPredictionResult/NG_BC_', sampleID)
}else if(grepl('CID', sampleID) & !(tools == 'RNAMagnet')){
  fpath.dataset <- paste0('~/1-Datasets/NatureGenetics_2021_BC/ScriptForCCC/scRNA_Data/',
                          sampleID, '_ser.rds')
  dir.output <- paste0('~/2-CCC-benchmark-202204-202208/CCC-benchmark-202212/Data/Step1_LRPredictionResult/NG_BC_', sampleID)
}else if(grepl('CK', sampleID) & tools == 'RNAMagnet'){
  fpath.dataset <- paste0('~/1-Datasets/Nature_2022_MI/ScriptForCCC/scRNA_Data_Mm/',
                          sampleID, '_ser.rds')
  dir.output <- paste0('~/2-CCC-benchmark-202204-202208/CCC-benchmark-202212/Data/Step1_LRPredictionResult/N_MI_', sampleID)
}else if(grepl('CK', sampleID) & !(tools == 'RNAMagnet')){
  fpath.dataset <- paste0('~/1-Datasets/Nature_2022_MI/ScriptForCCC/scRNA_Data/',
                          sampleID, '_ser.rds')
  dir.output <- paste0('~/2-CCC-benchmark-202204-202208/CCC-benchmark-202212/Data/Step1_LRPredictionResult/N_MI_', sampleID)
}else if(grepl('Slide', sampleID)){
  if(tools %in% c('CellPhoneDB2', 'CellPhoneDB3', 'CellTalker', 'NicheNet', 'iTALK', 'ICELLNET', 'scMLnet','NATMI', 'cell2cell')){
    fpath.dataset <- paste0('~/1-Datasets/MouseEmbryo_GSE166692/ScriptForCCC/DataForCCC/',
                            sampleID, '_hs_ser.rds')
  }else{
    fpath.dataset <- paste0('~/1-Datasets/MouseEmbryo_GSE166692/ScriptForCCC/DataForCCC/',
                            sampleID, '_mm_ser.rds')
  }
  dir.output <- paste0('~/2-CCC-benchmark-202204-202208/CCC-benchmark-202212/Data/Step1_LRPredictionResult/MouseEmbryo_', sampleID)
}else if(grepl('pbmc', sampleID) & !(tools == 'RNAMagnet')){
  fpath.dataset <- paste0('~/1-Datasets/10XGenomics_PBMC/ScriptForCellTypeAnn/DataForCCC/', sampleID, '_ser.rds')
  dir.output <- paste0('~/2-CCC-benchmark-202204-202208/CCC-benchmark-202212/Data/Step1_LRPredictionResult/', sampleID)
}else if(grepl('pbmc', sampleID) & (tools == 'RNAMagnet')){
  fpath.dataset <- paste0('~/1-Datasets/10XGenomics_PBMC/ScriptForCellTypeAnn/DataForCCC/', sampleID, '_mm_ser.rds')
  dir.output <- paste0('~/2-CCC-benchmark-202204-202208/CCC-benchmark-202212/Data/Step1_LRPredictionResult/', sampleID)
}

ser <- readRDS(fpath.dataset)

if(!dir.exists(dir.output)){
  dir.create(dir.output)
}

fpath.tools <- paste(dir.output, tools, sep = '/')
if(!dir.exists(fpath.tools)){
  dir.create(fpath.tools)
}

if(tools == 'cell2cell'){
  result <- cell2cell_function(ser, fpath = fpath.tools, species = species)
}else if(tools == 'CellCall'){
  if(grepl('CK', sampleID) | grepl('pbmc', sampleID)){
    n <- 2
  }else if(grepl('CID', sampleID)){
    n <- 3
  }else if(grepl('Slide', sampleID)){
    n <- 5
  }
  result <- CellCall_function(ser, n, species = species)
}else if(tools == 'CellChat'){
  result <- CellChat_function(ser, species = species)
}else if(tools == 'CellPhoneDB2'){
  result <- CellPhoneDB2_function(ser, fpath.tools)
}else if(tools == 'CellPhoneDB3'){
  result <- CellPhoneDB3_function(ser, fpath.tools)
}else if(tools == 'CellTalker'){
  result <- CellTalker_function(ser)
}else if(tools == 'Connectome'){
  result <- Connectome_function(ser, species = species)
}else if(tools == 'CytoTalk'){
  result <- CytoTalk_function(ser, fpath.tools, species = species)
}else if(tools == 'Domino'){
  if(sampleID %in% c('pbmc4k', 'pbmc8k')){
    species <- 'hg38'
  }else if(sampleID == 'pbmc6k'){
    species <- 'hg19'
  }else if(grepl('Slide' ,sampleID)){
    species <- 'mm10'
  }else{
    species <- 'hg38'
  }
  result <- Domino_function(ser, fpath.tools, species = species)
}else if(tools == 'ICELLNET'){
  result <- ICELLNET_function(ser)
}else if(tools == 'iTALK'){
  result <- iTALK_function(ser)
}else if(tools == 'NATMI'){
  result <- NATMI_function(ser, fpath.tools, species = species)
}else if(tools == 'NicheNet'){
  result <- NicheNet_function(ser, lr = TRUE)
}else if(tools == 'PyMINEr'){
  result <- PyMINEr_function(ser, fpath.tools, species)
}else if(tools == 'RNAMagnet'){
  result <- RNAMagnet_function(ser)
}else if(tools == 'scConnect'){
  result <- scConnect_function(ser, fpath.tools, species)
}else if(tools == 'scMLnet'){
  result <- scMLnet_function(ser)
}else if(tools == 'scSeqComm'){
  result <- scSeqComm_function(ser, species)
}else if(tools == 'SingleCellSignalR'){
  result <- SCSR_function(ser, species)
}

saveRDS(result, file = paste0(fpath.tools, '/result.rds'))