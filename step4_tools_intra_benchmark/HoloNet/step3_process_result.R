rm(list = ls())

library(tidyr)
library(stringr)
set.seed(123)

setwd("./step4_tools_intra_benchmark/HoloNet/")

files <- list.files("./output/")
result <- lapply(files, function(file){
  lr <- str_split(file, "_")[[1]][1]
  target <- str_split(file, "_")[[1]][2] %>% str_replace(".csv", "")
  
  res <- read.csv(paste0("./output/", file)) %>%
    pivot_longer(cols = -X, names_to = "receiver", values_to = "value")
  res$receiver <- str_replace_all(res$receiver, "[.]", " ") %>%
    str_replace_all("B cells", "B-cells") %>%
    str_replace_all("T cells", "T-cells")
  colnames(res)[1] <- "sender"
  
  res <- res[which(res$receiver == "Cancer Epithelial"), ]
  res$lr <- rep(lr, n = nrow(res))
  res <- separate(res, lr, c("ligand", "receptor"), ":")
  res$target <- rep(target, n = nrow(res))
  res
})
result <- do.call(rbind, result)
result <- result[-which(result$sender == "Cancer Epithelial"), ]

result_holonet <- result; rm(result)

result_holonet$LRT <- paste(result_holonet$ligand, result_holonet$receptor, result_holonet$target, sep = "_")
result_holonet <- result_holonet[,c(3,7)]

result_holonet <- aggregate(value ~ LRT, data = result_holonet, sum)
result_holonet <- tidyr::separate(result_holonet, LRT, c('ligand', 'receptor', 'target'), "_")
result_holonet_ligand <- result_holonet[, c(1,3,4)]
result_holonet_receptor <- result_holonet[, c(2,3,4)]
result_holonet_ligand <- aggregate(value ~ ligand+target, data = result_holonet_ligand, sum)
result_holonet_receptor <- aggregate(value ~ receptor+target, data = result_holonet_receptor, sum)
colnames(result_holonet_ligand) <- c("regulon", 'target',"value")
colnames(result_holonet_receptor) <- c("regulon", 'target',"value")
result_holonet_ligand$type <- 'ligand'
result_holonet_receptor$type <- 'receptor'
result_holonet <- rbind(result_holonet_ligand, result_holonet_receptor)
save(result_holonet, file = './result_holonet.RData')
