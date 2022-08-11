rm(list = ls())
set.seed(123)
setwd('./step3_run_lr_benchmark/input/tools_result_for_benchmark/')

# ST-benchmark——Jaccard index

file.names <- list.files(".")

for (file.name in file.names) {
  res.path <- paste0("./", file.name)
  load(file = res.path)
  res <- objects(pattern = "^result")
  result <- get(res)
  if(grepl('CellChat', file.name)){
    result <- result$trim
  }else if(grepl('scSeqComm', file.name)){
    result <- result$`0.8`
  }else if(grepl('iTALK',file.name)){
    result <- result$`10`
  }else if(grepl('CellTalker', file.name)){
    result <- result$`0`
  }
  result <- do.call(rbind, result)
  
  if(any(colnames(result) %in% c("sender", "reciever"))){
    result$all <- paste(result$sender, result$ligand, result$receptor, result$reciever, sep = "_")
  }else{
    result <- tibble::rownames_to_column(result, var = "senrec")
    result$senrec <- stringr::str_replace_all(result$senrec, "\\.[0-9]+", "")
    result <- tidyr::separate(result, senrec, c("sender", "reciever"), "_")
    result$all <- paste(result$sender, result$ligand, result$receptor, result$reciever, sep = "_")
  }
  result.names <- paste0("JC.", res)
  assign(result.names, result)
  rm(list = objects(pattern = "^result"))
}
rm(list = ls.str(mode = 'character'))

result <- matrix(nrow=18)
for(res.1 in objects(pattern = "^JC")){
  result.1 <- get(res.1)
  
  tmp.overlap <- c()
  for (res.2 in objects(pattern = "^JC")) {
    print(paste(res.1, res.2, sep = "_"))
    result.2 <- get(res.2)
    tmp <- intersect(result.1$all, result.2$all)
    tmp.overlap <- rbind(tmp.overlap, length(tmp))
  }
  result <- cbind(result, tmp.overlap)
}
rm(res.1, res.2, tmp.overlap, tmp, result.1, result.2)
result <- result[,-1]
result <- as.data.frame(result)
colnames(result) <- substring(objects(pattern = "^JC"), 11)
rownames(result) <- substring(objects(pattern = "^JC"), 11)
result <- result[, c("cpdb2", "cpdb3", "CellTalker", # expression value; mean
                     "Connectome", "NATMI", "ICELLNET", # expression value; product
                     "scConnect", "Cellinker", "CellChat", # expression value; others
                     "SingleCellSignalR", "CytoTalk", "CellCall", "scSeqComm", "NicheNet", # expression value; pathway-based
                     "Domino","scMLnet", "PyMINEr", "iTALK")]
result <- result[c("cpdb2", "cpdb3", "CellTalker", # expression value; mean
                   "Connectome", "NATMI", "ICELLNET", # expression value; product
                   "scConnect", "Cellinker", "CellChat", # expression value; others
                   "SingleCellSignalR", "CytoTalk", "CellCall", "scSeqComm", "NicheNet", # expression value; pathway-based
                   "Domino","scMLnet", "PyMINEr", "iTALK"),]
result <- as.matrix(result)

pheatmap::pheatmap(result,cluster_rows = FALSE, cluster_cols = FALSE,
                   treeheight_row = 0,
                   treeheight_col = 0,
                   display_numbers = TRUE,
                   fontsize=10,
                   #annotation_col = anno,
                   number_format = "%.0f",
                   #border_color=NA,
                   angle_col = "45",
                   colorRampPalette(c("#FFFFFF"))(50)
)
