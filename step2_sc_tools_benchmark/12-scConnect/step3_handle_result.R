rm(list = ls())
library(stringr)

setwd('./step2_sc_tools_benchmark/12-scConnect/')
LR <- read.csv("./result/output.scConnect.csv", header = TRUE)
LR <- LR[,-1]

ligands <- read.csv("./input/database/ligands.csv", header = T)
ligands <- ligands[,-c(1, 3, 9)]
length(unique(ligands$ligand)) # confirm non-duplicate
LR_updata <- merge(ligands, LR,by = "ligand")
colnames(LR_updata)[2] <- "ligand_gene"
LR_updata <- LR_updata[,-c(3:6)]

receptors <- read.csv('./input/database/receptors.csv')
receptors <- dplyr::distinct(receptors, receptor, type, gene, .keep_all = FALSE)
LR_updata <- merge(receptors, LR_updata, by = 'receptor')
colnames(LR_updata)[3] <- 'receptor_gene'

LR_updata$LR <- paste(LR_updata$ligand_gene, LR_updata$receptor_gene,sep = "_")
#LR_updata <- dplyr::distinct(LR_updata, LR, sender, reciever, .keep_all = TRUE)
LR_updata <- LR_updata[,c("sender", "reciever", "LR", "score", "log_score", "specificity", 
                  "importance", "ligand_gene", "ligand_zscore", "ligand_pval", "receptor_gene", 
                  "receptor_zscore", "receptor_pval")]

result.scConnect <- list()
tmp.result <- LR_updata[which(LR_updata$ligand_pval<0.05), ]
tmp.result <- tmp.result[which(tmp.result$receptor_pval<0.05),]
result.scConnect[["raw"]] <- LR_updata
result.scConnect[["filter"]] <- tmp.result
save(result.scConnect, file = "./result/result.scConnect.RData")
