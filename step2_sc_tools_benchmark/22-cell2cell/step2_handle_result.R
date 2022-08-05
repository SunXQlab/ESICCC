rm(list = ls())   
set.seed(123)

setwd('./step2_sc_tools_benchmark/22-cell2cell/')
ccc.pvalue <- read.csv('./output/ccc_pval.csv')
ccc.pvalue <- tidyr::pivot_longer(ccc.pvalue, -X, names_to = 'sr', values_to = 'pvalue')

table(ccc.pvalue$pvalue<0.05)
