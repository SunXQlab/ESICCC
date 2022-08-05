rm(list = ls())
set.seed(123)
setwd('./step3_run_lr_benchmark/input/tools_result_for_benchmark/')

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

counts <- c()
for(res in objects(pattern = "^JC")){
  result.tmp <- get(res)
  tmp <- dim(result.tmp)[1]
  counts <- c(counts, tmp)
}
result <- data.frame(tools = substring(objects(pattern = "^JC"), 11), 
                     LR_pairs = counts)
rownames(result) <- result$tools
rm(list = objects(pattern = '^JC'))
rm(result.tmp, counts, res, tmp)
result$tools <- factor(result$tools, levels = c("cpdb2", "cpdb3", "CellTalker", # expression value; mean
                                                "Connectome", "NATMI", "ICELLNET", # expression value; product
                                                "scConnect", "Cellinker", "CellChat", # expression value; others
                                                "SingleCellSignalR", "CytoTalk", "CellCall", "scSeqComm", "NicheNet", # expression value; pathway-based
                                                "Domino","scMLnet", "PyMINEr", "iTALK"))

library(ggplot2)

ggplot(data = result, mapping = aes(x = factor(tools), y = LR_pairs, fill = LR_pairs)) + 
  geom_bar(stat = 'identity', position = 'dodge', color=colorRampPalette(c("#CF7660"))(18))+
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))+
  xlab("") +
  scale_x_discrete(labels = c("cpdb2" = "CellPhoneDB2", "cpdb3" = "CellPhoneDB3"))+
  geom_text(mapping = aes(label = LR_pairs))
