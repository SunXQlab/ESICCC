rm(list = ls()); gc()
suppressMessages(library(tidyverse))
set.seed(123)

datasets <- list.files('./Data/Step6_LRBenchSamplingResult/', recursive = FALSE) %>%
  gsub('_[0-9]+', '',.) %>% unique()
ratios <- c(50, 60, 70, 80, 90, 100)

record_cells <- lapply(datasets, function(dataset){
  tmp <- lapply(ratios, function(ratio){
    if(grepl('CID', dataset) & ratio == 100){
      fpath.dataset <- paste0('~/1-Datasets/NatureGenetics_2021_BC/ScriptForCCC/scRNA_Data/', dataset,  '_ser.rds')
    }else if(grepl('CID', dataset) & ratio != 100){
      fpath.dataset <- paste0('~/1-Datasets/NatureGenetics_2021_BC/ScriptForCCC/scRNA_Data_sampling/', dataset,'_', ratio,  '_ser.rds')
    }else if(grepl('CK', dataset) & ratio == 100){
      fpath.dataset <- paste0('~/1-Datasets/Nature_2022_MI/ScriptForCCC/scRNA_Data/', dataset, '_ser.rds')
    }else if(grepl('CK', dataset) & ratio != 100){
      fpath.dataset <- paste0('~/1-Datasets/Nature_2022_MI/ScriptForCCC/scRNA_Data_sampling/', dataset, '_', ratio , '_ser.rds')
    }else if(grepl('Slide', dataset) & ratio == 100){
      fpath.dataset <- paste0('~/1-Datasets/MouseEmbryo_GSE166692/ScriptForCCC/DataForCCC/', dataset, '_mm_ser.rds')
    }else if(grepl('Slide', dataset) & ratio != 100){
      fpath.dataset <- paste0('~/1-Datasets/MouseEmbryo_GSE166692/ScriptForCCC/scRNA_Data_Sampling/', dataset, '_', ratio , '_mm_ser.rds')
    }else if(grepl('pbmc', dataset) & ratio == 100){
      fpath.dataset <- paste0('~/1-Datasets/10XGenomics_PBMC/ScriptForCellTypeAnn/DataForCCC/', dataset, '_ser.rds')
    }else if(grepl('pbmc', dataset) & ratio != 100){
      fpath.dataset <- paste0('~/1-Datasets/10XGenomics_PBMC/ScriptForCellTypeAnn/scRNA_Data_sampling/', dataset, '_', ratio , '_ser.rds')
    }
    ser <- readRDS(fpath.dataset)
    
    cells <- dim(ser)[[2]]
  })%>% unlist()
  
  tmp
})
names(record_cells) <- datasets
record_cells <- do.call(rbind, record_cells)
colnames(record_cells) <- paste0('S', ratios)
record_cells <- as.data.frame(record_cells)
record_cells <- tibble::rownames_to_column(record_cells, 'datasets')
record_cells <- pivot_longer(record_cells, -datasets, 'SampleRate')
colnames(record_cells)[3] <- 'Cells'
record_cells$class <- paste(record_cells$datasets, record_cells$SampleRate, sep = '_')
record_cells <- record_cells[,-c(1:2)]


record_celltypes <- lapply(datasets, function(dataset){
  tmp <- lapply(ratios, function(ratio){
    if(grepl('CID', dataset) & ratio == 100){
      fpath.dataset <- paste0('~/1-Datasets/NatureGenetics_2021_BC/ScriptForCCC/scRNA_Data/', dataset,  '_ser.rds')
    }else if(grepl('CID', dataset) & ratio != 100){
      fpath.dataset <- paste0('~/1-Datasets/NatureGenetics_2021_BC/ScriptForCCC/scRNA_Data_sampling/', dataset,'_', ratio,  '_ser.rds')
    }else if(grepl('CK', dataset) & ratio == 100){
      fpath.dataset <- paste0('~/1-Datasets/Nature_2022_MI/ScriptForCCC/scRNA_Data/', dataset, '_ser.rds')
    }else if(grepl('CK', dataset) & ratio != 100){
      fpath.dataset <- paste0('~/1-Datasets/Nature_2022_MI/ScriptForCCC/scRNA_Data_sampling/', dataset, '_', ratio , '_ser.rds')
    }else if(grepl('Slide', dataset) & ratio == 100){
      fpath.dataset <- paste0('~/1-Datasets/MouseEmbryo_GSE166692/ScriptForCCC/DataForCCC/', dataset, '_mm_ser.rds')
    }else if(grepl('Slide', dataset) & ratio != 100){
      fpath.dataset <- paste0('~/1-Datasets/MouseEmbryo_GSE166692/ScriptForCCC/scRNA_Data_Sampling/', dataset, '_', ratio , '_mm_ser.rds')
    }else if(grepl('pbmc', dataset) & ratio == 100){
      fpath.dataset <- paste0('~/1-Datasets/10XGenomics_PBMC/ScriptForCellTypeAnn/DataForCCC/', dataset, '_ser.rds')
    }else if(grepl('pbmc', dataset) & ratio != 100){
      fpath.dataset <- paste0('~/1-Datasets/10XGenomics_PBMC/ScriptForCellTypeAnn/scRNA_Data_sampling/', dataset, '_', ratio , '_ser.rds')
    }
    ser <- readRDS(fpath.dataset)
    
    celltype <- length(unique(ser$celltype))
  })%>% unlist()
  
  tmp
})
names(record_celltypes) <- datasets
record_celltypes <- do.call(rbind, record_celltypes)
colnames(record_celltypes) <- paste0('S', ratios)
record_celltypes <- as.data.frame(record_celltypes)
record_celltypes <- tibble::rownames_to_column(record_celltypes, 'datasets')
record_celltypes <- pivot_longer(record_celltypes, -datasets, 'SampleRate')
colnames(record_celltypes)[3] <- 'Celltypes'
record_celltypes$class <- paste(record_celltypes$datasets, record_celltypes$SampleRate, sep = '_')
record_celltypes <- record_celltypes[,-c(1:2)]

record <- merge(record_cells, record_celltypes, by = 'class')
record <- separate(record, class, c('datasets', 'ratio'), sep = '_')
saveRDS(record, file = './Data/Step6_LRBenchSampling/DatasetsRecord.rds')
