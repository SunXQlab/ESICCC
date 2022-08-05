setwd("./step2_sc_tools_benchmark/10-PyMINEr/input_and_output/")
set.seed(123)

# handle result
if(T){
  result <- read.table("./autocrine_paracrine_signaling/extracellular_plasma_membrane_cell_type_specific_interactions_read.txt",
                       header = TRUE, sep = "\t")
  tmp.result <- result
  tmp.result$cell_type_1[which(tmp.result$cell_type_1 == "bool_sample_group_0")] <- "Endothelial"
  tmp.result$cell_type_1[which(tmp.result$cell_type_1 == "bool_sample_group_1")] <- "CAFs"
  tmp.result$cell_type_1[which(tmp.result$cell_type_1 == "bool_sample_group_2")] <- "PVL"
  tmp.result$cell_type_1[which(tmp.result$cell_type_1 == "bool_sample_group_3")] <- "B-cells"
  tmp.result$cell_type_1[which(tmp.result$cell_type_1 == "bool_sample_group_4")] <- "Plasmablasts"
  tmp.result$cell_type_1[which(tmp.result$cell_type_1 == "bool_sample_group_5")] <- "T-cells"
  tmp.result$cell_type_1[which(tmp.result$cell_type_1 == "bool_sample_group_6")] <- "Myeloid"
  tmp.result$cell_type_1[which(tmp.result$cell_type_1 == "bool_sample_group_7")] <- "Cancer Epithelial"
  
  tmp.result$cell_type_2[which(tmp.result$cell_type_2 == "bool_sample_group_0")] <- "Endothelial"
  tmp.result$cell_type_2[which(tmp.result$cell_type_2 == "bool_sample_group_1")] <- "CAFs"
  tmp.result$cell_type_2[which(tmp.result$cell_type_2 == "bool_sample_group_2")] <- "PVL"
  tmp.result$cell_type_2[which(tmp.result$cell_type_2 == "bool_sample_group_3")] <- "B-cells"
  tmp.result$cell_type_2[which(tmp.result$cell_type_2 == "bool_sample_group_4")] <- "Plasmablasts"
  tmp.result$cell_type_2[which(tmp.result$cell_type_2 == "bool_sample_group_5")] <- "T-cells"
  tmp.result$cell_type_2[which(tmp.result$cell_type_2 == "bool_sample_group_6")] <- "Myeloid"
  tmp.result$cell_type_2[which(tmp.result$cell_type_2 == "bool_sample_group_7")] <- "Cancer Epithelial"
}

temp.result <- dplyr::distinct(tmp.result, cell_type_1, gene_1_symbol, cell_type_2, gene_2_symbol,.keep_all = TRUE)
result.PyMINEr <- temp.result
save(result.PyMINEr, file = '../result.PyMINEr.RData')
