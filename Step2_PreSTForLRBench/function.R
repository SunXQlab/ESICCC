# Calculate the euclidean distance between cell pairs
CloDistCP <- function(ser, image = TRUE){
  if(image){
    coord <- ser@images$image@coordinates[, c('imagerow', 'imagecol')]
  }else{
    coord <- ser@meta.data[, c('imagerow', 'imagecol')]
  }
  
  meta <- ser@meta.data[, 'celltype', drop = FALSE]
  identical(rownames(meta), rownames(coord))
  
  meta <- cbind(meta, coord)
  meta <- meta[order(meta$celltype),]
  
  dist.cp <- as.matrix(dist(meta[, 2:3], method = "euclidean"))
  dist.cp[dist.cp==0] <- 1
  dist.cp[!lower.tri(dist.cp, diag = TRUE)] <- 0
  dist.cp <- as.data.frame(dist.cp)
  dist.cp <- tibble::rownames_to_column(dist.cp, var = "cell.1")
  dist.cp <- tidyr::pivot_longer(data = dist.cp, cols = -cell.1, 
                                 names_to = "cell.2", values_to = "distance")
  dist.cp <- dist.cp[which(dist.cp$distance != 0), ]
  
  ct <- data.frame(barcode = rownames(meta), celltype = meta$celltype)
  dist.cp <- dplyr::inner_join(dist.cp, ct, by = c('cell.1' = 'barcode'))
  colnames(dist.cp)[4] <- "celltype.1"
  dist.cp <- dplyr::inner_join(dist.cp, ct, by = c('cell.2' = 'barcode'))
  colnames(dist.cp)[5] <- "celltype.2"
  dist.cp$cp <- paste(dist.cp$celltype.1, dist.cp$celltype.2, sep = "_")
  dist.cp <- dist.cp[which(dist.cp$celltype.1!=dist.cp$celltype.2), ]
  dist.cp <- dist.cp[,c("cell.1", "cell.2", "distance", "cp")]
  
  # remove cell pairs with the number less than 30
  remove.cp <- names(which(table(dist.cp$cp)<30))
  if(length(remove.cp)!=0){
    dist.cp <- dist.cp[-which(dist.cp$cp %in% remove.cp), ]
  }
  
  perc <- c(0.1, 0.2, 0.3, 0.4, 0.5)
  cellpairs <- unique(dist.cp$cp)
  library(doParallel)
  cl <- makeCluster(5)
  registerDoParallel(cl)
  close_distant <- foreach(per = perc) %dopar% {
    tmp.close_distant <- lapply(cellpairs, function(cp){
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
    return(tmp.close_distant)
  }
  stopCluster(cl)
  names(close_distant) <- perc
  return(close_distant)
}
