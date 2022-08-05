rm(list = ls())

library(domino)
library(Seurat)
library(dplyr)
set.seed(123)

# run domino
if(T){
  load("./step1_data_process/se_data/step2_scRNA_10X.RData")
  setwd('./step2_sc_tools_benchmark/15-Domino/')
  
  auc = t(read.table('./input_for_Domino/auc_mtx.csv', 
                     header = TRUE, row.names = 1, stringsAsFactors = FALSE, sep = ','))
  
  se.sc <- se.sc[, colnames(auc)]
  se.sc <- ScaleData(se.sc, features = rownames(se.sc))
  
  counts = se.sc@assays$RNA@counts
  z_scores = se.sc@assays$RNA@scale.data
  clusters = se.sc@active.ident
  rm(se.sc)
  
  source("./create_domino.R") 
  # line 200 of rscript of creat_domino function: change tf_genes = df[row, 10] to tf_genes = df[row, 9]
  tnbc_dom = create_domino(signaling_db = './input_for_Domino/cpdb/', 
                           features = auc, counts = counts, z_scores = z_scores, 
                           clusters = clusters, 
                           df = './input_for_Domino/reg.csv')
  rm(counts, z_scores, clusters, auc)
  rm(domino, add_rl_column, create_domino)
  
  dom = build_domino(tnbc_dom, min_tf_pval = 0.05, max_tf_per_clust = 10, max_rec_per_tf = 10)
  save(dom, tnbc_dom, file = "./output/domino.RData")
}
load("./output/domino.RData")
rm(tnbc_dom)

# bulid ligand-receptor dataframe
if(T){
  lig_to_rec <- data.frame()
  for (rec in names(dom@linkages$rec_lig)) {
    lig <- dom@linkages$rec_lig[[rec]]
    
    if(length(lig)>0){
      lig.rec <- data.frame(ligand = lig, receptor = rep(rec, length(lig)))
    }else{
      lig.rec <- data.frame(ligand = lig, receptor = rec)
    }
    
    lig_to_rec <- rbind(lig_to_rec, lig.rec)
  }
  rm(rec, lig, lig.rec)
  lig_to_rec <- lig_to_rec[-which(lig_to_rec$ligand==""),]
}

# handle result
if(T){
  domino.result <- lapply(levels(dom@clusters), function(reciever){
    print(reciever)
    # get ligands of sender cells which are communicated with reciever cells
    sender.ligands <- dom@cl_signaling_matrices[[reciever]] %>%
      as.data.frame(.) %>% tibble::rownames_to_column(var = "ligand") %>%
      tidyr::pivot_longer(cols = -ligand, names_to = "sender", values_to = "expression") %>%
      .[which(.$expression>0), ]
    sender.ligands$sender <- stringr::str_replace_all(sender.ligands$sender, "L_", "")
    sender.ligands <- sender.ligands[, -3]
    
    if(dim(sender.ligands)[1] != 0){
      # get tfs of reciever cells
      reciever.tfs <- dom@linkages$clust_tf[[reciever]]
      rec_tf <- data.frame()
      for(tf in reciever.tfs){
        rec <- dom@linkages$tf_rec[[tf]]
        
        if(length(rec)>0){
          rec.tf <- data.frame(receptor = rec, tf = rep(tf, length(rec)))
        }else{
          rec.tf <- data.frame()
        }
        rec_tf <- rbind(rec_tf, rec.tf)
      }
      rm(tf, rec, rec.tf)
      
      lig_rec_tf <- merge(lig_to_rec, rec_tf, by = "receptor")
      lig_rec_tf <- merge(sender.ligands, lig_rec_tf, by = "ligand")
      lig_rec_tf$reciever <- reciever
      lig_rec_tf$tf <- stringr::str_replace_all(lig_rec_tf$tf, "\\...", "")
    }else{
      lig_rec_tf <- NA
    }
    lig_rec_tf
  })
  
  names(domino.result) <- as.character(levels(dom@clusters))
  domino.result[["Cancer Epithelial"]] <- NULL
  
  tmp.result <- do.call(rbind, domino.result)
  tmp.result$sr <- paste(tmp.result$sender, tmp.result$reciever, sep = "_")
  
  result.Domino <- tmp.result
}
save(result.Domino, 
     file = "./result.Domino.RData")
