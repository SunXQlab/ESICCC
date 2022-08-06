rm(list = ls())

library(Seurat)
library(tidyverse)
library(nichenetr)
set.seed(123)

# load ST data
if(T){
  # load gene expression matrix
  expr.data <- Read10X("./step1_data_process/dataset/ST.Data/filtered_count_matrices/CID4465_filtered_count_matrix", 
                       gene.column = 1)
  se.st <- CreateSeuratObject(counts = expr.data, assay = 'Spatial', min.cells = 10, min.features = 500)
  se.st <- SCTransform(se.st, assay = "Spatial")
  
  # load cell type metadata
  if(T){
    metadata <- data.frame(se.st@meta.data, row.names = rownames(se.st@meta.data))
    load("./step3_benchmark/ST_benchmark_data/ST.CID4465.meta.RData")
    ST.CID4465.meta <- tibble::column_to_rownames(ST.CID4465.meta, var = "barcode")
    ST.CID4465.meta <- ST.CID4465.meta[order(ST.CID4465.meta$celltype),]
    ST.CID4465.meta$celltype <- stringr::str_replace_all(ST.CID4465.meta$celltype, 
                                                         "B_cells", "B-cells")
    ST.CID4465.meta$celltype <- stringr::str_replace_all(ST.CID4465.meta$celltype, 
                                                         "T_cells", "T-cells")
    ST.CID4465.meta$celltype <- stringr::str_replace_all(ST.CID4465.meta$celltype, 
                                                         "Cancer_Epithelial", "Cancer Epithelial")
  }
  
  se.st <- AddMetaData(se.st, metadata = ST.CID4465.meta)
  Idents(se.st) <- "celltype"
  
  rm(expr.data, ST.CID4465.meta)
}
load('./step4_tools_intra_benchmark/stMLnet/input/CID4465_seruat.RData')
# load prior data
if(T){
  ligand_target_matrix = readRDS("./step2_sc_tools_benchmark/21-NicheNet/input/ligand_target_matrix.rds")
  #ligand_target_matrix[1:5,1:5]
  
  lr_network = readRDS("./step2_sc_tools_benchmark/21-NicheNet/input/lr_network.rds")
  #head(lr_network)
  
  weighted_networks = readRDS("./step2_sc_tools_benchmark/21-NicheNet/input/weighted_networks.rds")
  weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network %>% distinct(from,to), by = c("from","to"))
  
}

celltype <- unique(se.st$celltype)

# define the sender genes
expressed_genes_sender <- lapply(celltype[-1], function(ct){
  expressed_genes_sender_ct <- get_expressed_genes(ct, se.st, pct = 0.05)
  expressed_genes_sender_ct
})
expressed_genes_sender = expressed_genes_sender %>% unlist() %>% unique()

# define the receiver genes
expressed_genes_receiver = get_expressed_genes(celltype[1], se.st, pct = 0.05)
background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

# define the genes of interest
geneset =  readRDS(file = "./step1_data_process/result/Giotto_result/CID4465_bc_ICGs.rds")
geneset <- geneset[['Cancer Epithelial']] %>% unlist() %>% unique()
geneset = geneset %>% .[. %in% rownames(ligand_target_matrix)]

# define the potential ligands
ligands = lr_network %>% pull(from) %>% unique()
receptors = lr_network %>% pull(to) %>% unique()

expressed_ligands = intersect(ligands,expressed_genes_sender)
expressed_receptors = intersect(receptors,expressed_genes_receiver)

potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()

# Perform NicheNet ligand activity analysis
ligand_activities = predict_ligand_activities(geneset = geneset, 
                                              background_expressed_genes = background_expressed_genes, 
                                              ligand_target_matrix = ligand_target_matrix, 
                                              potential_ligands = potential_ligands)
ligand_activities = ligand_activities %>% arrange(-pearson) %>% mutate(rank = rank(desc(pearson)))


# 选取所有配体作为下游分析目标方便验证
best_upstream_ligands = ligand_activities %>% arrange(-pearson) %>% pull(test_ligand)
head(best_upstream_ligands)

# Get the active ligand-target links
active_ligand_target_links_df = best_upstream_ligands %>% 
  lapply(get_weighted_ligand_target_links,
         geneset = geneset, 
         ligand_target_matrix = ligand_target_matrix, 
         n = nrow(ligand_target_matrix)) %>% 
  bind_rows()

active_ligand_target_links_df %>% group_by(ligand) %>% summarize(count=n())

NicheNet.lt.weights = active_ligand_target_links_df
saveRDS(NicheNet.lt.weights,"./step4_tools_intra_benchmark/NicheNet/result/tumor_LRpair_weight.rds")
