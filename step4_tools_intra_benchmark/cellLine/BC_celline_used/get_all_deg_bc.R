rm(list = ls())

library(dplyr)
set.seed(123)

setwd('./step4_tools_intra_benchmark/cellLine/BC_celline_used/')

# logic
wd <- "./receptor/"
files <- list.files(wd)
degs_ls_receptor <- lapply(files, function(file){
  degs_df <- read.csv(paste0(wd,file), header = T)
  degs_df$isDEGs = degs_df$P.Value < 0.05
  deg_ls <- degs_df$isDEGs
  names(deg_ls) <- degs_df$X
  deg_ls
})
name <- gsub("DE_", "", files)
name <- gsub(".csv", "", name)
name <- paste0(name, "_KO")
names(degs_ls_receptor) <- name
rm(name, files, wd); gc()

wd <- "./ligand/"
files <- list.files(wd)
degs_ls_ligand <- lapply(files, function(file){
  print(file)
  degs_df <- read.csv(paste0(wd,file), header = T)
  degs_df$isDEGs = degs_df$P.Value < 0.05 #adj.P.Val < 0.05 & abs(degs_df$logFC)>1
  deg_ls <- degs_df$isDEGs
  names(deg_ls) <- degs_df$X
  deg_ls
})
name <- gsub("DE_", "", files)
name <- gsub(".csv", "", name)
name <- paste0(name, "_Treatment")
names(degs_ls_ligand) <- name
rm(name, files, wd); gc()

lapply(degs_ls_ligand, sum)
lapply(degs_ls_receptor, sum)

degs_ls <- append(degs_ls_receptor, degs_ls_ligand)
rm(list = objects(pattern = "^degs_ls_"))

degs_ls <- degs_ls[unlist(lapply(degs_ls, sum))>50]

save(degs_ls, file = "./degs_ls.RData")
