rm(list = ls())

setwd('./step2_sc_tools_benchmark/19-NATMI')

result.NATMI <- read.csv("./output/Edges_lrc2p.csv")
result.NATMI <- result.NATMI[,c(1:4, 6:20)]
result.NATMI <- dplyr::filter(result.NATMI, Ligand.detection.rate > 0.05, 
                              Receptor.detection.rate > 0.05)
colnames(result.NATMI)[1:4] <- c("sender", "ligand", "receptor", "reciever")
result.NATMI$sr <- paste(result.NATMI$sender, result.NATMI$reciever, sep = "_")

save(result.NATMI, file = "./result.NATMI.RData")
