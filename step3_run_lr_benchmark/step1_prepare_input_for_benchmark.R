## load ST expression matrix and close_distance groups ##
rm(list = ls())
library(Seurat)
set.seed(123)

load('./step1_data_process/se_data/ST.se.RData')

setwd('./step3_run_lr_benchmark/')

# load gene expression matrix of ST data
se.st <- SCTransform(se.st, assay = "Spatial")
norm.data <- GetAssayData(se.st, "data", "SCT")
norm.data <- as.matrix(norm.data)

# load ST metadata
meta <- se.st@meta.data
coordinate <- se.st@images$image@coordinates
identical(rownames(meta), rownames(coordinate))
meta <- cbind(meta, coordinate)
meta <- meta[, c('celltype', 'imagerow', 'imagecol')]
rm(coordinate, se.st)

meta$celltype <- stringr::str_replace_all(meta$celltype, "B_cells", "B-cells")
meta$celltype <- stringr::str_replace_all(meta$celltype, "T_cells", "T-cells")
meta$celltype <- stringr::str_replace_all(meta$celltype, "Cancer_Epithelial", "Cancer Epithelial")
meta <- meta[order(meta$celltype),]

# calculate the euclidean distance between cell pairs
if(F){
  loc <- meta[, c("imagerow", "imagecol")]
  dist.cp <- as.matrix(dist(loc, method = "euclidean"))
  diag(dist.cp) <- 1
  dist.cp[!lower.tri(dist.cp, diag = TRUE)] <- 0
  dist.cp <- as.data.frame(dist.cp)
  dist.cp <- tibble::rownames_to_column(dist.cp, var = "cell.1")
  dist.cp <- tidyr::pivot_longer(data = dist.cp, cols = -cell.1, 
                                 names_to = "cell.2", values_to = "distance")
  dist.cp <- dist.cp[-which(dist.cp$distance == 0), ]
  
  ct <- data.frame(barcode = rownames(meta), celltype = meta$celltype)
  dist.cp <- merge(dist.cp, ct, by.x = "cell.1", by.y = "barcode")
  colnames(dist.cp)[4] <- "celltype.1"
  dist.cp <- merge(dist.cp, ct, by.x = "cell.2", by.y = "barcode")
  colnames(dist.cp)[5] <- "celltype.2"
  dist.cp$cp <- paste(dist.cp$celltype.1, dist.cp$celltype.2, sep = "_")
  dist.cp <- dist.cp[-which(dist.cp$celltype.1==dist.cp$celltype.2), ]
  dist.cp <- dist.cp[,c("cell.1", "cell.2", "distance", "cp")]
  
  perc <- c(0.1, 0.2, 0.3, 0.4, 0.5)
  close_distant <- lapply(perc, function(per){
    tmp.close_distant <- lapply(unique(dist.cp$cp), function(cp){
      dist <- list()
      tmp.dist <- dist.cp[which(dist.cp$cp == cp), ]
      tmp.dist <- tmp.dist[order(tmp.dist$distance, decreasing = T), ]
      pct.close <- floor(dim(tmp.dist)[1]*per)
      pct.distant <- floor(dim(tmp.dist)[1]*per)
      close.dist <- tail(tmp.dist, n = pct.close)
      distant.dist <- head(tmp.dist, n = pct.distant)
      dist[["close"]] <- close.dist
      dist[["distant"]] <- distant.dist
      dist
    })
    names(tmp.close_distant) <- unique(dist.cp$cp)
    tmp.close_distant
  })
  names(close_distant) <- perc
  
  rm(ct, dist.cp, loc, meta, perc)
  
  save(close_distant, norm.data, file = "./Input_for_benchmark/norm_cd.RData")
}


