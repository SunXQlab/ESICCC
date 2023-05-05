rm(list = ls())
suppressMessages(library(Seurat))
source('./Script/Step2_PreSTForLRBench/function.R')
set.seed(123)

counts.st <- Read10X('~/1-Datasets/MouseEmbryo_GSE166692/sciSpaceData/data/', 
                     gene.column = 1)

cellmeta <- read.csv('~/1-Datasets/MouseEmbryo_GSE166692/sciSpaceData/cellmeta.tsv', 
                     sep = '\t')
genes <- read.csv('~/1-Datasets/MouseEmbryo_GSE166692/sciSpaceData/genemeta.tsv', 
                  sep = '\t')

for (slide in unique(cellmeta$slide_id)) {
  output.path <- paste0('./Data/Step2_PreSTForLRBench/', gsub(' ', '', slide))
  if(!dir.exists(output.path)){
    dir.create(output.path)
  }
  
  if(T){
    cells_sel <- rownames(cellmeta)[which(cellmeta$slide_id==slide)]
    
    # counts: solve the problem of duplicated genes
    if(T){
      counts_sub <- counts.st[, cells_sel]
      counts_sub <- as.matrix(counts_sub)
      counts_sub <- as.data.frame(counts_sub)
      counts_sub$gene <- genes$gene_short_name; rownames(counts_sub) <- NULL
      dup_genes <- counts_sub$gene[which(duplicated(counts_sub$gene))]
      counts_sub1 <- counts_sub[-which(counts_sub$gene %in% dup_genes),]
      counts_sub2 <- counts_sub[which(counts_sub$gene %in% dup_genes),]
      counts_sub2 <- aggregate(.~gene,mean,data=counts_sub2)
      counts_sub <- rbind(counts_sub1, counts_sub2); rownames(counts_sub) <- NULL
      rm(counts_sub1, counts_sub2, dup_genes, cells_sel); gc()
      counts_sub <- tibble::column_to_rownames(counts_sub, 'gene')
      counts_sub <- as.matrix(counts_sub)
    }
    
    meta_sub <- cellmeta[which(cellmeta$slide_id==slide),]
    colnames(meta_sub)[c(9,10,19)] <- c('imagerow', 'imagecol', 'celltype')
    meta_sub$celltype <- gsub(' ', '', meta_sub$celltype)
    meta_sub$celltype <- gsub('Cardiacmusclelineages', 'CardiacMuscleLineages', meta_sub$celltype)
    
    shared.barcode <- intersect(colnames(counts_sub), rownames(meta_sub))
    
    meta_sub <- meta_sub[shared.barcode, ]
    counts_sub <- counts_sub[, rownames(meta_sub)]
    identical(rownames(meta_sub), colnames(counts_sub))
    
    ser <- CreateSeuratObject(counts = counts_sub, meta.data = meta_sub, 
                              assay = 'Spatial', min.cells = 1, min.features = 1)
    ser <- SCTransform(ser, assay = "Spatial")
    
    saveRDS(ser, file = paste0(output.path, '/STser.rds'))
    rm(counts_sub, meta_sub);gc()
  }
  
  ser <- readRDS(paste0(output.path, '/STser.rds'))
  close_distant <- CloDistCP(ser, image = FALSE)
  saveRDS(close_distant, file = paste0(output.path, '/CloDistCP.rds'))
}


