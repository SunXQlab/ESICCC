rm(list = ls())

library(scSeqComm)
library(Seurat)
library(dplyr)
set.seed(123)

load("./step1_data_process/se_data/step2_scRNA_10X.RData")

## scRNA-seq dataset
matrix <- GetAssayData(se.sc, "data", "RNA")
cell_cluster <- lapply(unique(se.sc$celltype_major), function(ct){
  meta <- se.sc@meta.data
  cells <- rownames(meta)[which(meta$celltype_major == ct)]
  cells
})
names(cell_cluster) <- unique(se.sc$celltype_major)

## Ligand - receptor pairs 
data("LR_pairs_Jin_2020")
LR_db <- LR_pairs_Jin_2020

## Transcriptional regulatory network——combined
if(F){
  TF_TG_db <- c(TF_TG_HTRIdb, TF_TG_RegNetwork_High, TF_TG_RegNetwork_Med_High, 
                TF_TG_TRRUSTv2, TF_TG_TRRUSTv2_HTRIdb_RegNetwork_High)
  tmp.db <- list()
  for (tf in unique(names(TF_TG_db))) {
    rep <- which(names(TF_TG_db) == tf)
    if(length(rep)>1){
      tmp.tg <- unique(unlist(TF_TG_db[rep]))
      tmp.db[[tf]] <- tmp.tg
    }else{
      tmp.db[[tf]] <- TF_TG_db[[tf]]
    }
  }
  TF_TG_db <- tmp.db
  rm(rep, tf, tmp.tg, tmp.db)
  save(TF_TG_db, file = "./step2_sc_tools_benchmark/13-scSeqComm/TF_TG_db.RData")
}
load("./step2_sc_tools_benchmark/13-scSeqComm/TF_TG_db.RData")

## Receptor-Transcription factor a-priori association from gene signaling networks
data("TF_PPR_KEGG_human")

## Intercellular and intracellular signaling analysis
scSeqComm_res <- scSeqComm_analyze(gene_expr = matrix,
                                   cell_group = cell_cluster,
                                   LR_pairs_DB = LR_db,
                                   TF_reg_DB = TF_TG_db, 
                                   R_TF_association = TF_PPR_KEGG_human,
                                   N_cores = 30)

result.scSeqComm <- scSeqComm_res$comm_results
#result.scSeqComm<- result.scSeqComm[-which(is.na(result.scSeqComm$S_intra)), ]

cutoff <- c(0.5, 0.6, 0.7, 0.8, 0.9)

result.scSeqComm <- lapply(cutoff, function(co){
  selected_comm <- dplyr::filter(result.scSeqComm, S_inter>co, S_intra>co)
})
names(result.scSeqComm) <- cutoff

save(result.scSeqComm, file = "./step2_sc_tools_benchmark/13-scSeqComm/result.scSeqComm.RData")
