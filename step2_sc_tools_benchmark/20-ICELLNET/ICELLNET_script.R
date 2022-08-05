rm(list = ls())

library(icellnet)
library(Seurat)
library(BiocGenerics)
library("org.Hs.eg.db")
library("hgu133plus2.db")
library(jetset)
library(ggplot2)
library(dplyr)
library(icellnet)
library(gridExtra)
set.seed(123)

db=as.data.frame(read.csv(curl::curl(url="https://raw.githubusercontent.com/soumelis-lab/ICELLNET/master/data/ICELLNETdb.tsv"), 
                          sep="\t",header = T, check.names=FALSE, stringsAsFactors = FALSE, na.strings = ""))

load("./step1_data_process/se_data/step2_scRNA_10X.RData")

## Retrieve gene expression matrix
# Taking into account the total nb of cells in each cluster
# filter out the percent of genes less than 10%
filter.perc=5
average.clean= sc.data.cleaning(object = se.sc, db = db, 
                                filter.perc = filter.perc, 
                                save_file = F)

## Apply icellnet pipeline on cluster of interest
data.icell=as.data.frame(gene.scaling(as.data.frame(average.clean), n=1, db=db))

PC.data=as.data.frame(data.icell[, colnames(data.icell)], row.names = rownames(data.icell))


PC.target=data.frame(Class = colnames(PC.data)[-dim(data.icell)[2]], 
                     ID = colnames(PC.data)[-dim(data.icell)[2]], 
                     Cell_type = colnames(PC.data)[-dim(data.icell)[2]])
rownames(PC.target) = colnames(PC.data)[-dim(data.icell)[2]]

my.selection = colnames(PC.data)[-dim(data.icell)[2]]

result.all <- lapply(colnames(PC.data)[-dim(data.icell)[2]], function(ct){
  result <- list()
  ## Compute intercellular communication scores
  score.computation = icellnet.score(direction = "out", PC.data = PC.data, 
                                     CC.data = as.data.frame(data.icell[,ct], row.names = rownames(data.icell)),  
                                     PC.target = PC.target, PC = my.selection, CC.type = "RNAseq", 
                                     PC.type = "RNAseq",  db = db)
  lr <- score.computation[[2]][apply(score.computation[[2]], 1, function(y) any(!is.na(y))),]
  lr <- lr[which(rowSums(lr) > 0),]
  result[["lr"]] = lr
  result[["sum"]] <- as.data.frame(score.computation[[1]])
  result
})

names(result.all) <- colnames(PC.data)[-dim(data.icell)[2]]


result.ICELLNET <- lapply(names(result.all), function(ct){
  result <- result.all[[ct]]
  lr <- result$lr
  lr <- as.data.frame(lr)
  lr <- tibble::rownames_to_column(lr, "LR")
  colnames(lr) <- c("LR", paste(ct, colnames(lr)[2:dim(lr)[2]], sep = "_"))
  result.lr <- lr %>% tidyr::pivot_longer(cols = -LR, names_to = "sr", values_to = "product_value")
  result.lr <- tidyr::separate(data = result.lr, col = sr, into = c("sender", "reciever"), sep = "_")
  result.lr <- tidyr::separate(data = result.lr, col = LR, into = c("Ligand", "Receptor"), sep = "/")
  result.lr <- result.lr[which(result.lr$product_value>0), ]
})
names(result.ICELLNET) <- names(result.all)
result.ICELLNET <- do.call(rbind, result.ICELLNET)
save(result.ICELLNET, file = './step2_sc_tools_benchmark/20-iCELLNET/result.ICELLNET.RData')
