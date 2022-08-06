rm(list = ls())
library(igraph)
library(tidyverse)

set.seed(123)

setwd("./step4_tools_intra_benchmark/CytoTalk/")

load("./result/result_CytoTalk.RData")
result_CytoTalk[["Cancer Epithelial_Cancer Epithelial"]] <- NULL

result <- lapply(names(result_CytoTalk), function(cp){
  browser()
  sender <-  gsub("_.*", "",cp)
  receiver <- "Cancer Epithelial"
  
  tmp <- list()
  res <- result_CytoTalk[[cp]]
  if(!is.null(res$pathways)){
    network <- res$pcst$final_network
    network <- network[!(network$node1_type == sender & network$node2_type == sender), ]
    network <- network[!(network$node1_type == receiver & network$node2_type == sender), ]
    network$node1 <- toupper(network$node1)
    network$node2 <- toupper(network$node2)
    LR <- network[which(network$node1_type == sender & network$node2_type == receiver), ]
    LR <- paste(LR$node1, LR$node2, sep = "_")
    Ligand <- network[which(network$node1_type == sender & network$node2_type == receiver), "node1"] %>%
      unique()
    Receptor <- network[which(network$node1_type == sender & network$node2_type == receiver), "node2"] %>%
      unique()
    Target <- c(network[which(network$node1_type == receiver & network$node2_type == receiver), "node1"],
                network[which(network$node1_type == receiver & network$node2_type == receiver), "node2"]) %>%
      unique()
    Target <- setdiff(Target, c(Ligand, Receptor))
    tmp <- list(Edge = network,
                LR = LR,
                Ligand = Ligand, 
                Receptor = Receptor, 
                Target = Target)
  }else{
    tmp <- NULL
  }
  tmp
})

names(result) <- names(result_CytoTalk); rm(result_CytoTalk)
result[[which(unlist(lapply(result, is.null)))]] <- NULL

edge_list <- lapply(result, function(res){
  res <- res$Edge[,c(1,2,10)]
  colnames(res) <- c('from', 'to', 'weight')
  res
})
edge_df <- do.call(rbind, edge_list)
edge_df <- distinct(edge_df, from, to)
intranetwork <- graph_from_edgelist(as.matrix(edge_df), directed = FALSE)

Ligand <- lapply(result, function(res){res$Ligand}) %>% unlist() %>% unique()
Receptor <- lapply(result, function(res){res$Receptor}) %>% unlist() %>% unique()
Target <- lapply(result, function(res){res$Target}) %>% unlist() %>% unique()

distance <- distances(intranetwork, 
                      v = c(Ligand,Receptor), 
                      to = Target)
distance <- reshape2::melt(distance)
distance <- distance[!is.infinite(distance$value),]
colnames(distance) <- c("from","to","distance")
distance$from_type <- ifelse(distance$from %in% Ligand,"Ligand","Receptor")

saveRDS(distance,file = './result/LRT_distance.rds')
