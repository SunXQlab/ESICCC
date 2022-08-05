setwd("./step2_sc_tools_benchmark/3-CellPhoneDB/")

rm(list = ls())
library(tidyr)

file.path <- "./output/significant_means.txt"

res <- read.table(file.path, header = TRUE, sep = "\t")
res <- res[, c(2, 13:dim(res)[2])]

result.df <- res %>% pivot_longer(cols = -interacting_pair, names_to = "sr", values_to = "mean")
result.df <- dplyr::filter(result.df,  !is.na(mean))
result.df$sr <- stringr::str_replace_all(result.df$sr,
                                         "T\\.cells", "T-cells")
result.df$sr <- stringr::str_replace_all(result.df$sr,
                                         "B\\.cells", "B-cells")
result.df$sr <- stringr::str_replace_all(result.df$sr, 
                                         "Cancer\\.Epithelial", "Cancer Epithelial")
result.df <- tidyr::separate(data = result.df, col = sr, into = c("sender", "reciever"), sep = "\\.")

result.cpdb2 <- result.df
save(result.cpdb2, file = './result.cpdb2.RData')
