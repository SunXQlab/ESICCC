rm(list = ls())
library(Seurat)
library(celltalker)
library(dplyr)
set.seed(123)

load("./step1_data_process/se_data/step2_scRNA_10X.RData")

min.expression <-c(0, 50, 100, 150, 200)

## Run celltalker
result.celltalker <- lapply(min.expression, function(min.expr){
  result_celltalker <- list()
  result <- celltalk(input_object=se.sc,
                     metadata_grouping="celltype_major",
                     ligand_receptor_pairs=ramilowski_pairs,
                     number_cells_required=1,
                     min_expression=min.expr,
                     max_expression=20000,
                     scramble_times=10)
  
  top_stats.1 <- result %>%
    mutate(fdr=p.adjust(p_val,method="fdr")) %>%
    #filter(fdr < 0.05) %>%
    filter(p_val < 0.05) %>%
    filter(interact_ratio > 0) %>%
    group_by(cell_type1) %>%
    #top_n(3,interact_ratio) %>%
    ungroup()
  top_stats.2 <- result %>%
    mutate(fdr=p.adjust(p_val,method="fdr")) %>%
    filter(fdr < 0.05) %>%
    filter(p_val < 0.05) %>%
    filter(interact_ratio > 0) %>%
    group_by(cell_type1) %>%
    #top_n(3,interact_ratio) %>%
    ungroup()
  result_celltalker[["unfilter"]] <- result
  result_celltalker[["filter.p"]] <- top_stats.1
  result_celltalker[["filter.fdr"]] <- top_stats.2
  result_celltalker
})

names(result.celltalker) <- min.expression
