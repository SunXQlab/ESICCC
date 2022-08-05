rm(list = ls())

library(nichenetr)
library(Seurat)
library(tidyverse)
set.seed(123)

# load scRNA-seq data
load("./step1_data_process/se_data/step2_scRNA_10X.RData")
setwd('./step2_sc_tools_benchmark/21-NicheNet/')

# load data
if(T){
  ligand_target_matrix = readRDS("./input/ligand_target_matrix.rds")
  #ligand_target_matrix[1:5,1:5]
  
  lr_network = readRDS("./input/lr_network.rds")
  #head(lr_network)
  
  weighted_networks = readRDS("./input/weighted_networks.rds")
  weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network %>% distinct(from,to), by = c("from","to"))
  
}

# find DEGs in all celltypes
DE_table_receiver_all = FindAllMarkers(object = se.sc, min.pct = 0.05, 
                                       logfc.threshold = 0.15, test.use = "t", 
                                       return.thresh = 0.05)
celltype <- unique(se.sc$celltype_major)

result.NicheNet <- lapply(celltype, function(sender){
  tmp.result <- list()
  for (reciever in celltype[-which(celltype == sender)]) {
    # define the sender and reciever celltypes
    expressed_genes_receiver = get_expressed_genes(reciever, se.sc, pct = 0.05)
    background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
    
    expressed_genes_sender = get_expressed_genes(sender, se.sc, pct = 0.05)
    
    # find DEGs in reciever cells
    DE_table_receiver <- DE_table_receiver_all[which(DE_table_receiver_all$cluster == reciever), ]
    geneset_oi = DE_table_receiver %>% 
      filter(p_val_adj <= 0.05) %>% 
      pull(gene)
    geneset_oi = geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]
    
    # define the potential ligands
    ligands = lr_network %>% pull(from) %>% unique()
    receptors = lr_network %>% pull(to) %>% unique()
    
    expressed_ligands = intersect(ligands,expressed_genes_sender)
    expressed_receptors = intersect(receptors,expressed_genes_receiver)
    
    potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()
    
    # Perform NicheNet ligand activity analysis
    ligand_activities = predict_ligand_activities(geneset = geneset_oi, 
                                                  background_expressed_genes = background_expressed_genes, 
                                                  ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
    ligand_activities = ligand_activities %>% arrange(-pearson) %>% mutate(rank = rank(desc(pearson)))
    
    best_upstream_ligands = ligand_activities %>% top_n(20, pearson) %>% 
      arrange(-pearson) %>% pull(test_ligand) %>% unique()
    
    # Infer receptors and top-predicted target genes of top-ranked ligands
    ## Active target gene inference
    active_ligand_target_links_df = best_upstream_ligands %>% 
      lapply(get_weighted_ligand_target_links,geneset = geneset_oi, 
             ligand_target_matrix = ligand_target_matrix, n = 200) %>% 
      bind_rows() %>% drop_na()
    
    ## Receptors of top-ranked ligands
    lr_network_top = lr_network %>% 
      filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>%
      distinct(from,to)
    
    best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()
    
    lr_network_top_df_large = weighted_networks_lr %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)
    
    cp <- paste(sender, reciever, sep = "_")
    tmp.result[[cp]][["lr"]] <- lr_network_top_df_large
    tmp.result[[cp]][["lt"]] <- active_ligand_target_links_df
  }
  tmp.result
})

result <- list()
for (n in 1:length(result.NicheNet)) {
  tmp <- result.NicheNet[[n]]
  result <- c(result, tmp)
}
result.NicheNet <- result

save(result.NicheNet, file = "./result.NicheNet.RData")
