rm(list = ls())
set.seed(123)

res.p <- read.table("./step2_sc_tools_benchmark/17-Cellinker/output/res_p.txt", 
                    header = TRUE, sep = "\t")
res.p <- tibble::rownames_to_column(res.p, var = "interacting_pair")
res.p <-tidyr::pivot_longer(res.p, cols = -interacting_pair, names_to = "sr", values_to = "pvalue")
res.p <- res.p[which(res.p$pvalue<0.05), ]
res.p$comb <- paste(res.p$interacting_pair, res.p$sr, sep = "_")

res.value <- read.table("./step2_sc_tools_benchmark/17-Cellinker/output/res_value.txt",
                        header = TRUE, sep = "\t")
res.value <- tibble::rownames_to_column(res.value, var = "interacting_pair")
res.value <-tidyr::pivot_longer(res.value, cols = -interacting_pair, names_to = "sr", values_to = "value")
res.value <- res.value[which(res.value$value>0), ]
res.value$comb <- paste(res.value$interacting_pair, res.value$sr, sep = "_")

result <- merge(res.p, res.value, by = "comb")
result <- result[,-c(1:3)]
colnames(result) <- c("pvalue","LR", "sr", "value")
result <- result[, c("LR", "sr", "value", "pvalue")]
result.Cellinker <- result
save(result.Cellinker, file = "./step2_sc_tools_benchmark/17-Cellinker/result.Cellinker.RData")
