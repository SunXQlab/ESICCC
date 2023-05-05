rm(list = ls())
# All the code run in the '/Data' directory

#######################
## SingleCellSignalR ##
#######################
SCSR_function <- function(ser, species = 'human'){
  suppressMessages(library(SingleCellSignalR))
  suppressMessages(library(Seurat))
  suppressMessages(library(tidyr))
  suppressMessages(library(dplyr))
  suppressMessages(library(lobstr))
  set.seed(123)
  
  start.time <- proc.time()
  
  if(species == 'human'){
    species <- "homo sapiens"
  }else{
    species <- "mus musculus"
  }
  
  matrix.sc <- GetAssayData(ser, "data", "RNA")
  matrix.sc <- as.matrix(matrix.sc)
  meta.sc <- data.frame(celltype = ser$celltype, row.names = colnames(ser))
  
  # Digitize Labels
  i <- 1
  for(ct in unique(meta.sc$celltype)){
    meta.sc[which(meta.sc$celltype == ct), "ct_num"] <- i
    i <- i+1
  }
  c.names <- as.character(unique(meta.sc$celltype))
  
  
  signal <- cell_signaling(data = matrix.sc, genes = rownames(matrix.sc), int.type = "paracrine",
                           species = species, cluster = meta.sc$ct_num, c.names = c.names, write = FALSE)
  inter.net <- inter_network(data = matrix.sc, signal = signal, genes = rownames(matrix.sc), 
                             cluster = meta.sc$ct_num, c.names = c.names, species = species,write = FALSE)
  
  used.time <- proc.time()-start.time
  used.memory <- mem_used()
  
  result <- inter.net$`full-network`
  result$interaction.type <- NULL
  result <- separate(result, 'ligand', c('Sender', 'Ligand'), sep = '\\.')
  result <- separate(result, 'receptor', c('Receiver','Receptor'), sep = '\\.')
  result$all <- paste(result$Sender, result$Ligand, result$Receiver, result$Receptor, sep = '_')
  result <- distinct(result, all, .keep_all = TRUE)
  rownames(result) <- NULL
  
  result_record <- list(result=result, 
                        used_time = paste0(round(used.time[3]/60,3), ' min'), 
                        used_memory =  round(used.memory/1024/1024/1024,3))
  return(result_record)
}

##############
## CellChat ##
##############
CellChat_function <- function(ser, species = 'human'){
  suppressMessages(library(Seurat))
  suppressMessages(library(CellChat))
  suppressMessages(library(tidyverse))
  suppressMessages(library(lobstr))
  set.seed(123)
  
  start.time <- proc.time()
  
  matrix.sc <- GetAssayData(ser, "data", "RNA")
  matrix.sc <- as.matrix(matrix.sc)
  meta.sc <- data.frame(celltype = ser$celltype, row.names = colnames(ser))
  meta.sc$celltype <- as.character(meta.sc$celltype)
  
  cellchat <- createCellChat(object = matrix.sc, meta = meta.sc, group.by = "celltype")
  if(species == 'human'){
    cellchat@DB <- CellChatDB.human
  }else{
    cellchat@DB <- CellChatDB.mouse
  }
  
  cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  result <- computeCommunProb(cellchat) %>% filterCommunication(.) %>% 
    subsetCommunication(.)
  
  used.time <- proc.time()-start.time
  used.memory <- mem_used()
  
  result <- result[,1:5]
  colnames(result) <- c('Sender','Receiver','Ligand', 'Receptor', 'LRscore')
  result <- result[which(result$Sender != result$Receiver),]
  result$Ligand <- gsub('_', '&', result$Ligand)
  result$Receptor <- gsub('_', '&', result$Receptor)
  result$all <- paste(result$Sender, result$Ligand, result$Receiver, result$Receptor, sep = '_')
  result <- dplyr::distinct(result, all, .keep_all = TRUE)
  
  
  result_record <- list(result=result, 
                        used_time = paste0(round(used.time[3]/60,3), ' min'), 
                        used_memory =  round(used.memory/1024/1024/1024,3))
  return(result_record)
}

########### 
## iTALK ## top_genes was set to 50 (default)
########### Only for human (LR prior databases); the LR prior databases can be designed by users (complex)
iTALK_function <- function(ser){
  suppressMessages(library(iTALK))
  suppressMessages(library(Seurat))
  suppressMessages(library(lobstr))
  set.seed(123)
  
  start.time <- proc.time()
  
  matrix.sc <- GetAssayData(ser, "data", "RNA")
  matrix.sc <- as.matrix(matrix.sc)
  matrix.sc <- as.data.frame(t(matrix.sc))
  matrix.sc$cell_type <- ser$celltype
  
  # find top 50 percent highly expressed genes
  highly_exprs_genes <- rawParse(matrix.sc,stats='mean')
  # find the ligand-receptor pairs from highly expressed genes
  comm.list<-c('growth factor','other','cytokine','checkpoint')
  
  result <- NULL
  for(comm.type in comm.list){
    res.tmp <- FindLR(highly_exprs_genes,datatype='mean count',comm_type=comm.type)
    res.tmp <- res.tmp[order(res.tmp$cell_from_mean_exprs*res.tmp$cell_to_mean_exprs,decreasing=T),]
    result<-rbind(result,res.tmp)
  }
  
  used.time <- proc.time()-start.time
  used.memory <- mem_used()
  
  result$LRscore <- result$cell_from_mean_exprs*result$cell_to_mean_exprs
  result[,c(3,5,7)] <- NULL
  colnames(result) <- c('Ligand', 'Receptor', 'Sender', 'Receiver', 'LRscore')
  result <- result[which(result$Receiver != result$Sender), ]
  result$all <- paste(result$Sender, result$Ligand, result$Receiver, result$Receptor, sep = '_')
  result <- dplyr::distinct(result, all, .keep_all = TRUE)
  
  result_record <- list(result=result, 
                        used_time = paste0(round(used.time[3]/60,3), ' min'), 
                        used_memory =  round(used.memory/1024/1024/1024,3))
  return(result_record)
}

############## Run in parallel with 10 cores 
## CytoTalk ## Run in R 4.1.0 version
############## 
CytoTalk_function <- function(ser, fpath, sender = NULL, receiver = NULL, species = 'human'){
  suppressMessages(library(CytoTalk))
  suppressMessages(library(Seurat))
  suppressMessages(library(lobstr))
  set.seed(123)
  
  # create conda environment
  if(F){
    library(reticulate)  # To install and call Python modules from R.
    conda_create(envname = "r_reticulate_CytoTalk", python_version = "3.7.3")  # Create a new Conda environment to facilitate the Python module installation.
    conda_install(envname = "r_reticulate_CytoTalk", "pybind11")  # Install two necessary Python modules for correctly compiling and using the "pcst_fast" Python module.
    conda_install(envname = "r_reticulate_CytoTalk", "numpy")
    conda_install(envname = "r_reticulate_CytoTalk", "git+https://github.com/fraenkel-lab/pcst_fast.git", pip = TRUE) # To install the "pcst_fast" module.
  }
  
  # input
  if(T){
    input_fpath <- paste0(fpath, '/input/')
    if(!dir.exists(input_fpath)){
      dir.create(input_fpath)
    }
    
    fpath.mat <- paste0(input_fpath, "counts.csv")
    fpath.meta <- paste0(input_fpath,"metadata.csv")
    
    if(T){
      norm.matrix <- as.matrix(GetAssayData(ser, "data", "RNA"))
      write.csv(norm.matrix, fpath.mat, quote=F)
      
      meta.data <- data.frame(rownames(ser@meta.data), ser@meta.data$celltype)
      colnames(meta.data) <- NULL
      meta.data <- as.matrix(meta.data)
      write.csv(meta.data, fpath.meta, quote=F, row.names = FALSE)
      rm(norm.matrix, meta.data, ser);gc()
    }
  }
  
  start.time <- proc.time()

  if(species == 'human'){
    pcg = CytoTalk::pcg_human
    lrp = CytoTalk::lrp_human
  }else{
    pcg = CytoTalk::pcg_mouse
    lrp = CytoTalk::lrp_mouse
  }
  
  lst.sc <- read_matrix_with_meta(fpath.mat, fpath.meta)
  celltype <- unique(lst.sc$cell_types)
  
  celltypes <- as.character(unique(ser$celltype))
  comb <- combn(celltypes, 2)
  comb <- t(comb)
  comb <- as.data.frame(comb)
  
  result <- list()
  for (i in 1:dim(comb)[[1]]) {
    tmp <- tryCatch(CytoTalk::run_cytotalk(lst.sc, comb[i, 1], comb[i, 2], cores = 10, 
                                           cutoff_a = 0.05, cutoff_b = 0.05,
                                           pcg = pcg, lrp = lrp),
                    error=function(e){NA}
    )
    cp <- paste(comb[i, 1], comb[i, 2], sep = "_")
    result[[cp]] <- tmp
  }
  
  used.time <- proc.time()-start.time
  used.memory <- mem_used()
  
  result[which(is.na(result))] <- NULL
  
  result <- lapply(result, function(res){res$pcst$final_network})
  result <- do.call(rbind, result)
  result <- result[which(result$node1_type!=result$node2_type),]
  if (species=='human') {
    result$node1 <- toupper(result$node1)
    result$node2 <- toupper(result$node2)
    result$node1 <- gsub("ORF", "orf", result$node1)
    result$node2 <- gsub("ORF", "orf", result$node2)
  }
  result <- result[,c(1:4,10)]
  colnames(result) <- c('Ligand', 'Receptor', 'Sender', 'Receiver', 'LRscore')
  result$all <- paste(result$Sender, result$Ligand, result$Receiver, result$Receptor, sep = '_')
  result <- dplyr::distinct(result, all, .keep_all = TRUE)
  rownames(result) <- NULL
  
  result_record <- list(result=result, 
                        used_time = paste0(round(used.time[3]/60,3), ' min'), 
                        used_memory =  round(used.memory/1024/1024/1024,3))
  return(result_record)
}

##############
## CellCall ##
##############
CellCall_function <- function(ser, n=3, species = 'human'){
  suppressMessages(library(Seurat))
  suppressMessages(library(cellcall))
  suppressMessages(library(lobstr))
  set.seed(123)
  
  start.time <- proc.time()
  
  if(species == 'human'){
    species <- "Homo sapiens"
  }else{
    species <- "Mus musculus"
  }
  
  matrix.sc <- GetAssayData(ser, "counts", "RNA")
  matrix.sc <- as.matrix(matrix.sc)
  colnames(matrix.sc) <- paste(colnames(matrix.sc), gsub("-", " ", ser$celltype), sep = "_")
  matrix.sc <- as.data.frame(matrix.sc)
  
  mt <- CreateNichConObject(data=matrix.sc,
                            names.field = n,
                            names.delim = "_",
                            source = "UMI",
                            scale.factor = 10^6,
                            Org = species,
                            project = "Microenvironment")
  
  mt <- TransCommuProfile(object = mt,
                          pValueCor = 0.05,
                          CorValue = 0.1,
                          topTargetCor=1,
                          p.adjust = 0.05,
                          use.type="mean",
                          method="weighted",
                          IS_core = TRUE,
                          Org = species)
  
  result <- mt@data$expr_l_r_log2_scale
  
  used.time <- proc.time()-start.time
  used.memory <- mem_used()
  
  result <- as.data.frame(result)
  result <- tibble::rownames_to_column(result, var = "LR")
  result <- tidyr::pivot_longer(result,cols = -LR, names_to = "sr", values_to = "value")
  result <- result[which(result$value > 0), ]
  result <- tidyr::separate(data = result, col = sr, into = c("Sender", "Receiver"), sep = "-")
  result <- tidyr::separate(data = result, col = LR, into = c("Ligand", "Receptor"), sep = "-")
  colnames(result)[5] <- 'LRscore'
  result <- result[which(result$Sender != result$Receiver), ]
  result$Ligand <- gsub(',', '&', result$Ligand)
  result$Receptor <- gsub(',', '&', result$Receptor)
  result$all <- paste(result$Sender, result$Ligand, result$Receiver, result$Receptor, sep = '_')
  result <- dplyr::distinct(result, all, .keep_all = TRUE)
  
  
  result_record <- list(result=result, 
                        used_time = paste0(round(used.time[3]/60,3), ' min'), 
                        used_memory =  round(used.memory/1024/1024/1024,3))
  return(result_record)
}

############### 
## scSeqComm ## S_inter and S_intra were set to 0.8
############### Run in parallel with 10 cores
scSeqComm_function <- function(ser, species = 'human'){
  suppressMessages(library(scSeqComm))
  suppressMessages(library(Seurat))
  suppressMessages(library(lobstr))
  set.seed(123)
  
  start.time <- proc.time()
  
  ## scRNA-seq dataset
  matrix.sc <- GetAssayData(ser, "data", "RNA")
  cell_cluster <- lapply(unique(ser$celltype), function(ct){
    meta <- ser@meta.data
    cells <- rownames(meta)[which(meta$celltype == ct)]
    cells
  })
  names(cell_cluster) <- unique(ser$celltype)
  
  if(species == 'human'){
    ## Ligand - receptor pairs
    LR_db <- LR_pairs_Jin_2020
    
    ## Transcriptional regulatory network——combined
    if(T){
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
      rm(rep, tf, tmp.tg, tmp.db);gc()
    }
    
    ## Receptor-Transcription factor a-priori association from gene signaling networks
    TF_PPR_db <- TF_PPR_KEGG_human
  }else{
    ## Ligand - receptor pairs
    LR_db <- LR_pairs_Jin_2020_mouse
    
    ## Transcriptional regulatory network——combined
    if(T){
      TF_TG_db <- c(TF_TG_RegNetwork_High_mouse, TF_TG_RegNetwork_Med_High_mouse, 
                    TF_TG_TRRUSTv2_mouse, TF_TG_TRRUSTv2_RegNetwork_High_mouse)
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
      rm(rep, tf, tmp.tg, tmp.db);gc()
    }
    
    ## Receptor-Transcription factor a-priori association from gene signaling networks
    TF_PPR_db <- TF_PPR_KEGG_mouse
    
  }
  
  ## Intercellular and intracellular signaling analysis
  scSeqComm_res <- scSeqComm_analyze(gene_expr = matrix.sc,
                                     cell_group = cell_cluster,
                                     LR_pairs_DB = LR_db,
                                     TF_reg_DB = TF_TG_db, 
                                     R_TF_association = TF_PPR_db,
                                     N_cores = 10)
  result <- scSeqComm_res$comm_results
  result <- dplyr::filter(result, S_inter>0.8, S_intra>0.8)
  
  used.time <- proc.time()-start.time
  used.memory <- mem_used()
  
  result <- result[which(result$cluster_L != result$cluster_R),]
  result <- result[,c(1:2,4:5,9)]
  colnames(result) <- c('Ligand', 'Receptor', 'Sender', 'Receiver', 'LRscore')
  result$Ligand <- gsub(',', '&', result$Ligand)
  result$Receptor <- gsub(',', '&', result$Receptor)
  result$all <- paste(result$Sender, result$Ligand, result$Receiver, result$Receptor, sep = '_')
  result <- dplyr::distinct(result, all, .keep_all = TRUE)
 
  
  result_record <- list(result=result, 
                        used_time = paste0(round(used.time[3]/60,3), ' min'), 
                        used_memory =  round(used.memory/1024/1024/1024,3))
  return(result_record)
}

################
## Connectome ##
################
Connectome_function <- function(ser, species = 'human'){
  suppressMessages(library(Seurat))
  suppressMessages(library(Connectome))
  suppressMessages(library(lobstr))
  set.seed(123)
  
  start.time <- proc.time()
  
  if(species == 'human'){
    connectome.genes <- union(Connectome::ncomms8866_human$Ligand.ApprovedSymbol,
                              Connectome::ncomms8866_human$Receptor.ApprovedSymbol)
  }else{
    connectome.genes <- union(Connectome::ncomms8866_mouse$Ligand.ApprovedSymbol,
                              Connectome::ncomms8866_mouse$Receptor.ApprovedSymbol)
    species <- 'mouse'
  }
  
  genes <- connectome.genes[connectome.genes %in% rownames(ser)]
  ser <- ScaleData(ser,features = genes)
  sc.con <- CreateConnectome(ser,species = species, calculate.DOR = TRUE)
  result <- FilterConnectome(sc.con, max.p = 0.05, min.pct = 0.05, remove.na = T)
  
  used.time <- proc.time()-start.time
  used.memory <- mem_used()
  
  result <- result[which(result$source != result$target),]
  result <- result[,c(1:4, 17)]
  colnames(result) <- c('Sender', 'Receiver', 'Ligand', 'Receptor', 'LRscore')
  result$all <- paste(result$Sender, result$Ligand, result$Receiver, result$Receptor, sep = '_')
  result <- dplyr::distinct(result, all, .keep_all = TRUE)
  
 
  result_record <- list(result=result, 
                        used_time = paste0(round(used.time[3]/60,3), ' min'), 
                        used_memory =  round(used.memory/1024/1024/1024,3))
  return(result_record)
}

################
## CellTalker ## min_expression was set to 0
################ Only for human (LR prior databases)
CellTalker_function <- function(ser){
  suppressMessages(library(Seurat))
  suppressMessages(library(celltalker))
  suppressMessages(library(dplyr))
  suppressMessages(library(lobstr))
  suppressMessages(library(tidyr))
  source('./Step0_SharedInfo/CellTalker/celltalk_code.R')
  set.seed(123)
  
  start.time <- proc.time()
  
  result <- celltalk(input_object=ser,
                     metadata_grouping="celltype",
                     ligand_receptor_pairs=ramilowski_pairs,
                     number_cells_required=1,
                     min_expression=0,
                     max_expression=20000,
                     scramble_times=10)
  result <- result %>%
    mutate(fdr=p.adjust(p_val,method="fdr")) %>%
    filter(fdr < 0.05) %>%
    filter(p_val < 0.05) %>%
    filter(interact_ratio > 0)
  
  used.time <- proc.time()-start.time
  used.memory <- mem_used()
  
  result <- tidyr::separate(result, 'interaction', c('Ligand', 'Receptor'), sep = '_')
  result <- tidyr::separate(result, 'interaction_pairs', c('Sender', 'Receiver'), sep = '_')
  result <- result[,c(1:5)]; colnames(result)[5] <- 'LRscore'
  result <- result[which(result$Sender != result$Receiver),]
  result$all <- paste(result$Sender, result$Ligand, result$Receiver, result$Receptor, sep = '_')
  result <- dplyr::distinct(result, all, .keep_all = TRUE)
  
  result_record <- list(result=result, 
                        used_time = paste0(round(used.time[3]/60,3), ' min'), 
                        used_memory =  round(used.memory/1024/1024/1024,3))
  return(result_record)
}

############## 
## ICELLNET ## filter.perc was set to 5
############## Only for human (LR prior databases)
ICELLNET_function <- function(ser){
  suppressMessages(library(icellnet))
  suppressMessages(library(Seurat))
  suppressMessages(library(BiocGenerics))
  suppressMessages(library("org.Hs.eg.db"))
  suppressMessages(library("hgu133plus2.db"))
  suppressMessages(library(jetset))
  suppressMessages(library(dplyr))
  suppressMessages(library(gridExtra))
  suppressMessages(library(lobstr))
  set.seed(123)
  
  start.time <- proc.time()
  
  db=as.data.frame(read.csv('./Step0_SharedInfo/ICELLNET/ICELLNETdb.tsv',  sep="\t",
                            header = T, check.names=FALSE, stringsAsFactors = FALSE, na.strings = ""))
  
  ## Retrieve gene expression matrix
  # Taking into account the total nb of cells in each cluster
  # filter out the percent of genes less than 5%
  filter.perc=5
  average.clean= sc.data.cleaning(object = ser, db = db, filter.perc = filter.perc, save_file = F)
  ## Apply icellnet pipeline on cluster of interest
  data.icell=as.data.frame(gene.scaling(as.data.frame(average.clean), n=1, db=db))
  
  PC.data=as.data.frame(data.icell[, colnames(data.icell)], row.names = rownames(data.icell))
  PC.target=data.frame(Class = colnames(PC.data)[-dim(data.icell)[2]], 
                       ID = colnames(PC.data)[-dim(data.icell)[2]], 
                       Cell_type = colnames(PC.data)[-dim(data.icell)[2]])
  rownames(PC.target) = colnames(PC.data)[-dim(data.icell)[2]]
  
  PC.ct <- colnames(PC.data)[-dim(data.icell)[2]]
  CC.ct <- colnames(PC.data)[-dim(data.icell)[2]]
  
  result.all <- lapply(CC.ct, function(ct){
    ## Compute intercellular communication scores
    score.computation = icellnet.score(direction = "out", PC.data = PC.data, 
                                       CC.data = as.data.frame(data.icell[,ct], row.names = rownames(data.icell)),  
                                       PC.target = PC.target, PC = PC.ct[which(PC.ct!=ct)], CC.type = "RNAseq", 
                                       PC.type = "RNAseq",  db = db)
    lr <- as.matrix(score.computation[[2]][apply(score.computation[[2]], 1, function(y) any(!is.na(y))),])
    lr <- as.matrix(lr[which(rowSums(lr) > 0),])
    lr
  })
  
  names(result.all) <- colnames(PC.data)[-dim(data.icell)[2]]
  
  used.time <- proc.time()-start.time
  used.memory <- mem_used()
  
  result <- lapply(names(result.all), function(ct){
    lr <- result.all[[ct]]
    lr <- as.data.frame(lr)
    lr <- tibble::rownames_to_column(lr, "LR")
    colnames(lr) <- c("LR", paste(ct, colnames(lr)[2:dim(lr)[2]], sep = "_"))
    result.lr <- lr %>% tidyr::pivot_longer(cols = -LR, names_to = "sr", values_to = "LRscore")
    result.lr <- tidyr::separate(data = result.lr, col = sr, into = c("Sender", "Receiver"), sep = "_")
    result.lr <- tidyr::separate(data = result.lr, col = LR, into = c("Ligand", "Receptor"), sep = " / ")
    result.lr <- result.lr[which(result.lr$LRscore>0), ]
    result.lr <- result.lr[which(result.lr$Sender != result.lr$Receiver),]
    result.lr$Ligand <- gsub(' \\+ ', '&', result.lr$Ligand)
    result.lr$Receptor <- gsub(' \\+ ', '&', result.lr$Receptor)
    result.lr$all <- paste(result.lr$Sender, result.lr$Ligand, result.lr$Receiver, result.lr$Receptor, sep = '_')
    result.lr <- distinct(result.lr, all, .keep_all = TRUE)
  })
  result <- do.call(rbind, result)
  
  result_record <- list(result=result, 
                        used_time = paste0(round(used.time[3]/60,3), ' min'), 
                        used_memory =  round(used.memory/1024/1024/1024,3))
  result_record
}

############## 
## NicheNet ## top 20 ligands, top 250 targets of ligand
############## Only for human (LR prior databases)
NicheNet_function <- function(ser, sender = NULL, receiver = NULL, lr = TRUE){
  suppressMessages(library(nichenetr))
  suppressMessages(library(Seurat))
  suppressMessages(library(tidyverse))
  suppressMessages(library(lobstr))
  set.seed(123)
  
  start.time <- proc.time()
  
  # load data
  if(T){
    ligand_target_matrix = readRDS("./Step0_SharedInfo/NicheNet/ligand_target_matrix.rds")
    
    lr_network = readRDS("./Step0_SharedInfo/NicheNet/lr_network.rds")
    
    weighted_networks = readRDS("./Step0_SharedInfo/NicheNet/weighted_networks.rds")
    weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network %>% distinct(from,to), by = c("from","to"))
  }
  
  # find DEGs in all celltypes
  DE_table_receiver_all = FindAllMarkers(object = ser, min.pct = 0.05, 
                                         logfc.threshold = 0.15, test.use = "t", 
                                         return.thresh = 0.05)
  if(is.null(sender) & is.null(receiver)){
    sender_ct <- unique(ser$celltype)
    receiver_ct <- DE_table_receiver_all[which(DE_table_receiver_all$p_val_adj<=0.05), ]$cluster %>%
      unique(.) %>% as.character(.)
  }else{
    sender_ct <- sender
    receiver_ct <- receiver
  }
  
  result <- lapply(sender_ct, function(sender){
    tmp.result <- list()
    for (receiver in receiver_ct[which(receiver_ct != sender)]) {
      # define the sender and receiver celltypes
      expressed_genes_receiver = get_expressed_genes(receiver, ser, pct = 0.05)
      background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
      
      expressed_genes_sender = get_expressed_genes(sender, ser, pct = 0.05)
      
      # find DEGs in receiver cells
      DE_table_receiver <- DE_table_receiver_all[which(DE_table_receiver_all$cluster == receiver), ]
      geneset_oi = DE_table_receiver %>% filter(p_val_adj <= 0.05) %>% pull(gene)
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
      colnames(active_ligand_target_links_df) <- c("Ligand", 'Target', 'Score')
      active_ligand_target_links_df$Sender <- sender
      active_ligand_target_links_df$Receiver <- receiver
      
      ## Receptors of top-ranked ligands
      lr_network_top = lr_network %>% 
        filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>%
        distinct(from,to)
      
      best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()
      
      lr_network_top_df_large = weighted_networks_lr %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)
      colnames(lr_network_top_df_large) <- c('Ligand', 'Receptor', 'LRscore')
      lr_network_top_df_large$Sender <- sender
      lr_network_top_df_large$Receiver <- receiver
      
      cp <- paste(sender, receiver, sep = "_")
      
      if(lr){
        tmp.result[[cp]] <- lr_network_top_df_large
      }else{
        tmp.result[[cp]] <- active_ligand_target_links_df
      }
    }
    tmp.result <- do.call(rbind, tmp.result)
    rownames(tmp.result) <- NULL
    tmp.result
  })
  
  used.time <- proc.time()-start.time
  used.memory <- mem_used()
  
  result <- do.call(rbind, result)
  result$all <- paste(result$Sender, result$Ligand, result$Receiver, result$Receptor, sep = '_')
  result <- dplyr::distinct(result, all, .keep_all = TRUE)
  
  result_record <- list(result=result, 
                        used_time = paste0(round(used.time[3]/60,3), ' min'), 
                        used_memory =  round(used.memory/1024/1024/1024,3))
  return(result_record)
}

################## Run in parallel with 10 cores
## CellPhoneDB2 ## Only for human (LR prior databases)
################## 
CellPhoneDB2_function <- function(ser, fpath){
  suppressMessages(library(tidyr))
  suppressMessages(library(Seurat))
  suppressMessages(library(lobstr))
  set.seed(123)
  
  # input
  if(T){
    input_fpath <- paste0(fpath, '/input/')
    if(!dir.exists(input_fpath)){
      dir.create(input_fpath)
    }
    
    fpath.mat <- paste0(input_fpath, "counts.txt")
    fpath.meta <- paste0(input_fpath,"metadata.txt")
    
    norm.matrix <- as.matrix(GetAssayData(ser, "data", "RNA"))
    write.table(norm.matrix, fpath.mat, sep='\t', quote=F)
    
    meta.data <- data.frame(rownames(ser@meta.data), ser@meta.data$celltype)
    colnames(meta.data) <- NULL
    meta.data <- as.matrix(meta.data)
    write.table(meta.data, fpath.meta, sep='\t', quote=F, row.names=F)
    rm(norm.matrix, meta.data, ser);gc()
    
    output_fpath <- paste0(fpath, '/output/')
    if(!dir.exists(output_fpath)){
      dir.create(output_fpath)
    }
  }
  
  start.time <- proc.time()
  
  command <- paste('sh ../Script/Step0_LRToolsFunction/CellPhoneDB2_shell.sh', 
                   fpath.meta, fpath.mat, output_fpath, sep = ' ')
  system(command)
  
  used.time <- proc.time()-start.time
  used.memory <- mem_used()
  
  file.path <- paste(output_fpath, 'significant_means.txt', sep = '/')
  result <- read.table(file.path, header = TRUE, sep = "\t")
  result <- result[, c(2, 13:dim(result)[2])]
  
  result <- result %>% pivot_longer(cols = -interacting_pair, names_to = "sr", values_to = "LRscore")
  result <- dplyr::filter(result,  !is.na(LRscore))
  result <- tidyr::separate(data = result, col = sr, into = c("Sender", "Receiver"), sep = "\\.")
  
  # handle the complex information
  if(T){
    cpdb.complex <-read.csv("./Step0_SharedInfo/CellPhoneDB/complexes.csv")
    rec.complex <- cpdb.complex$complex_name[which(cpdb.complex$receptor == TRUE)]
    lig.complex <- cpdb.complex$complex_name[which(cpdb.complex$receptor != TRUE)]
    cpdb.complex <- cpdb.complex[,1:5]
    cpdb.gene <- read.csv("./Step0_SharedInfo/CellPhoneDB/genes.csv")
    cpdb.gene <- dplyr::distinct(cpdb.gene, gene_name, uniprot, hgnc_symbol, .keep_all = TRUE)
    cpdb.gene <- cpdb.gene[,1:2]
    
    tmp.cpdb.complex <- merge(cpdb.gene, cpdb.complex, by.y = "uniprot_1", by.x = "uniprot")
    tmp.cpdb.complex <- tmp.cpdb.complex[,-1]
    colnames(tmp.cpdb.complex)[1] <- "gene_1"
    tmp.cpdb.complex <- merge(cpdb.gene, tmp.cpdb.complex, by.y = "uniprot_2", by.x = "uniprot")
    tmp.cpdb.complex <- tmp.cpdb.complex[,-1]
    colnames(tmp.cpdb.complex)[1] <- "gene_2"
    tmp.cpdb.complex <- tmp.cpdb.complex[,-5]
    tmp.cpdb.complex <- merge(cpdb.gene, tmp.cpdb.complex, by.y = "uniprot_3", by.x = "uniprot", all.y = TRUE)
    tmp.cpdb.complex <- tmp.cpdb.complex[-117,]
    tmp.cpdb.complex <- tmp.cpdb.complex[,-1]
    colnames(tmp.cpdb.complex)[1] <- "gene_3"
    cpdb.complex <- tmp.cpdb.complex[,4:1]
    rm(cpdb.gene, tmp.cpdb.complex)
    
    cpdb.complex$gene <- paste(cpdb.complex$gene_1, cpdb.complex$gene_2, cpdb.complex$gene_3, sep = "&")
    cpdb.complex$gene <- gsub("&NA", "", cpdb.complex$gene)
    cpdb.complex <- cpdb.complex[,c("complex_name", "gene")]
  }
  
  # combine the complex information and result of cellphonedb
  if(T){
    complexes <- cpdb.complex$complex_name[which(stringr::str_count(cpdb.complex$complex_name, pattern = '_')==1)]
    for (complex in complexes) {
      change.pair <- which(grepl(complex, result$interacting_pair))
      if(length(change.pair)>0){
        change.complex <- gsub('_', '*', complex)
        result$interacting_pair <- gsub(complex, change.complex, result$interacting_pair)
      }
    }
    
    result <- tidyr::separate(data = result, col = interacting_pair, into = c("Ligand", "Receptor"), sep = "_")
    result$Ligand <- gsub('\\*', '_', result$Ligand)
    result$Receptor <- gsub('\\*', '_', result$Receptor)
    result <- merge(result, cpdb.complex, by.x = "Ligand", by.y = "complex_name", all.x = TRUE)
    result$Ligand[!is.na(result$gene)] <- result$gene[!is.na(result$gene)]
    result <- result[,-6]
    result <- merge(result, cpdb.complex, by.x = "Receptor", by.y = "complex_name", all.x = TRUE)
    result$Receptor[!is.na(result$gene)] <- result$gene[!is.na(result$gene)]
    result <- result[,-6]
    
    result.rec <- which(result$Receptor %in% rec.complex)
    result.lig <- which(result$Ligand %in% lig.complex)
    
    if(length(result.lig)>0){
      result <- result[-result.lig,]
      print(paste0('Delate：', length(result.lig)))
    }else if(length(result.rec)>0){
      result <- result[-result.rec,]
      print(paste0('Delate：', length(result.rec)))
    }
    
    result <- result[which(result$Sender!=result$Receiver),]
    result$all <- paste(result$Sender, result$Ligand, result$Receiver, result$Receptor, sep = '_')
    result <- dplyr::distinct(result, all, .keep_all = TRUE)
  }
  
  result_record <- list(result=result, 
                        used_time = paste0(round(used.time[3]/60,3), ' min'), 
                        used_memory =  round(used.memory/1024/1024/1024,3))
  return(result_record)
}

################## 
## CellPhoneDB3 ## Only for human (LR prior databases)
################## Run in parallel with 10 cores
CellPhoneDB3_function <- function(ser, fpath){
  suppressMessages(library(tidyr))
  suppressMessages(library(Seurat))
  suppressMessages(library(lobstr))
  set.seed(123)
  
  # input
  if(T){
    input_fpath <- paste0(fpath, '/input/')
    if(!dir.exists(input_fpath)){
      dir.create(input_fpath)
    }
    
    fpath.mat <- paste0(input_fpath, "counts.txt")
    fpath.meta <- paste0(input_fpath,"metadata.txt")
    fpath.deg <- paste0(input_fpath, 'degs.txt')
    
    norm.matrix <- as.matrix(GetAssayData(ser, "data", "RNA"))
    write.table(norm.matrix, fpath.mat, sep='\t', quote=F)
    
    meta.data <- data.frame(rownames(ser@meta.data), ser@meta.data$celltype)
    colnames(meta.data) <- NULL
    meta.data <- as.matrix(meta.data)
    write.table(meta.data, fpath.meta, sep='\t', quote=F, row.names=F)
    
    Idents(ser) <- ser$celltype
    DEGs <- FindAllMarkers(ser, test.use = 't',
                           verbose = F, only.pos = T, 
                           random.seed = 123, logfc.threshold = 0.15, 
                           min.pct = 0.05, return.thresh = 0.05)
    fDEGs = subset(DEGs, p_val_adj < 0.05 & avg_log2FC > 0.15)
    fDEGs = fDEGs[, c('cluster', 'gene', 'p_val_adj', 'p_val', 'avg_log2FC', 'pct.1', 'pct.2')]
    write.table(fDEGs, file = fpath.deg, sep = '\t', quote = F, row.names = F)
    
    rm(norm.matrix, meta.data, DEGs, fDEGs, ser);gc()
    
    output_fpath <- paste0(fpath, '/output/')
    if(!dir.exists(output_fpath)){
      dir.create(output_fpath)
    }
  }
  
  start.time <- proc.time()
  
  command <- paste('sh ../Script/Step0_LRToolsFunction/CellPhoneDB3_shell.sh', fpath.meta, fpath.mat, fpath.deg, output_fpath, sep = ' ')
  system(command)
  
  used.time <- proc.time()-start.time
  used.memory <- mem_used()
  
  file.path <- paste(output_fpath, 'significant_means.txt', sep = '/')
  
  result <- read.table(file.path, header = TRUE, sep = "\t")
  result <- result[, c(2, 13:dim(result)[2])]
  
  result <- result %>% pivot_longer(cols = -interacting_pair, names_to = "sr", values_to = "LRscore")
  result <- dplyr::filter(result,  !is.na(LRscore))
  result <- tidyr::separate(data = result, col = sr, into = c("Sender", "Receiver"), sep = "\\.")
  
  # handle the complex information
  if(T){
    cpdb.complex <-read.csv("./Step0_SharedInfo/CellPhoneDB/complexes.csv")
    rec.complex <- cpdb.complex$complex_name[which(cpdb.complex$receptor == TRUE)]
    lig.complex <- cpdb.complex$complex_name[which(cpdb.complex$receptor != TRUE)]
    cpdb.complex <- cpdb.complex[,1:5]
    cpdb.gene <- read.csv("./Step0_SharedInfo/CellPhoneDB/genes.csv")
    cpdb.gene <- dplyr::distinct(cpdb.gene, gene_name, uniprot, hgnc_symbol, .keep_all = TRUE)
    cpdb.gene <- cpdb.gene[,1:2]
    
    tmp.cpdb.complex <- merge(cpdb.gene, cpdb.complex, by.y = "uniprot_1", by.x = "uniprot")
    tmp.cpdb.complex <- tmp.cpdb.complex[,-1]
    colnames(tmp.cpdb.complex)[1] <- "gene_1"
    tmp.cpdb.complex <- merge(cpdb.gene, tmp.cpdb.complex, by.y = "uniprot_2", by.x = "uniprot")
    tmp.cpdb.complex <- tmp.cpdb.complex[,-1]
    colnames(tmp.cpdb.complex)[1] <- "gene_2"
    tmp.cpdb.complex <- tmp.cpdb.complex[,-5]
    tmp.cpdb.complex <- merge(cpdb.gene, tmp.cpdb.complex, by.y = "uniprot_3", by.x = "uniprot", all.y = TRUE)
    tmp.cpdb.complex <- tmp.cpdb.complex[-117,]
    tmp.cpdb.complex <- tmp.cpdb.complex[,-1]
    colnames(tmp.cpdb.complex)[1] <- "gene_3"
    cpdb.complex <- tmp.cpdb.complex[,4:1]
    rm(cpdb.gene, tmp.cpdb.complex)
    
    cpdb.complex$gene <- paste(cpdb.complex$gene_1, cpdb.complex$gene_2, cpdb.complex$gene_3, sep = "&")
    cpdb.complex$gene <- gsub("&NA", "", cpdb.complex$gene)
    cpdb.complex <- cpdb.complex[,c("complex_name", "gene")]
  }
  
  # combine the complex information and result of cellphonedb
  if(T){
    complexes <- cpdb.complex$complex_name[which(stringr::str_count(cpdb.complex$complex_name, pattern = '_')==1)]
    for (complex in complexes) {
      change.pair <- which(grepl(complex, result$interacting_pair))
      if(length(change.pair)>0){
        change.complex <- gsub('_', '*', complex)
        result$interacting_pair <- gsub(complex, change.complex, result$interacting_pair)
      }
    }
    
    result <- tidyr::separate(data = result, col = interacting_pair, into = c("Ligand", "Receptor"), sep = "_")
    result$Ligand <- gsub('\\*', '_', result$Ligand)
    result$Receptor <- gsub('\\*', '_', result$Receptor)
    result <- merge(result, cpdb.complex, by.x = "Ligand", by.y = "complex_name", all.x = TRUE)
    result$Ligand[!is.na(result$gene)] <- result$gene[!is.na(result$gene)]
    result <- result[,-6]
    result <- merge(result, cpdb.complex, by.x = "Receptor", by.y = "complex_name", all.x = TRUE)
    result$Receptor[!is.na(result$gene)] <- result$gene[!is.na(result$gene)]
    result <- result[,-6]
    
    result.rec <- which(result$Receptor %in% rec.complex)
    result.lig <- which(result$Ligand %in% lig.complex)
    
    if(length(result.lig)>0){
      result <- result[-result.lig,]
      print(paste0('Delate：', length(result.lig)))
    }else if(length(result.rec)>0){
      result <- result[-result.rec,]
      print(paste0('Delate：', length(result.rec)))
    }
    
    result <- result[which(result$Sender!=result$Receiver),]
    result$all <- paste(result$Sender, result$Ligand, result$Receiver, result$Receptor, sep = '_')
    result <- dplyr::distinct(result, all, .keep_all = TRUE)
  }
  
  result_record <- list(result=result, 
                        used_time = paste0(round(used.time[3]/60,3), ' min'), 
                        used_memory =  round(used.memory/1024/1024/1024,3))
  return(result_record)
}

########### 
## NATMI ## Run in parallel with 10 cores
########### 
NATMI_function <- function(ser, fpath, species = 'human'){
  suppressMessages(library(Seurat))
  suppressMessages(library(lobstr))
  set.seed(123)
  
  input_fpath <- paste0(fpath, '/input/')
  if(!dir.exists(input_fpath)){
    dir.create(input_fpath)
  }
  
  output_fpath <- paste0(fpath, '/output/')
  if(!dir.exists(output_fpath)){
    dir.create(output_fpath)
  }
  
  fpath.mat <- paste0(input_fpath, "counts.csv")
  write.csv(100 * (exp(as.matrix(GetAssayData(object = ser, assay = "RNA", slot = "data"))) - 1), 
            fpath.mat, row.names = T)
  meta <- data.frame(Cell = colnames(ser), Annotation = ser$celltype)
  fpath.meta <- paste0(input_fpath, "metadata.csv")
  write.csv(meta,fpath.meta, row.names = FALSE)
  rm(ser);gc()
  
  start.time <- proc.time()
  
  if(species != 'human'){
    species = 'mouse'
  }
  
  command <- paste('sh ../Script/Step0_LRToolsFunction/NATMI_shell.sh', 
                   species, fpath.mat, fpath.meta, output_fpath, sep = ' ')
  system(command)
  
  used.time <- proc.time()-start.time
  used.memory <- mem_used()
  
  result_path <- paste0(output_fpath, 'Edges_lrc2p.csv')
  result <- read.csv(result_path)
  result <- result[,c(1:4, 6:20)]
  result <- dplyr::filter(result, Ligand.detection.rate > 0.05, Receptor.detection.rate > 0.05)
  result <- result[,c(1:4, 16)]
  colnames(result) <- c("Sender", "Ligand", "Receptor", "Receiver", "LRscore")
  result <- result[which(result$Sender!=result$Receiver),]
  result$all <- paste(result$Sender, result$Ligand, result$Receiver, result$Receptor, sep = '_')
  result <- dplyr::distinct(result, all, .keep_all = TRUE)
  
  result_record <- list(result=result, 
                        used_time = paste0(round(used.time[3]/60,3), ' min'), 
                        used_memory =  round(used.memory/1024/1024/1024,3))
  return(result_record)
}

############### 
## scConnect ## 
###############
scConnect_function <- function(ser, fpath, species = 'human'){
  suppressMessages(library(Seurat))
  suppressMessages(library(lobstr))
  set.seed(123)
  
  input_fpath <- paste0(fpath, '/input/')
  if(!dir.exists(input_fpath)){
    dir.create(input_fpath)
  }
  
  fpath.mat <- paste0(input_fpath, "counts.csv")
  norm.matrix <- as.matrix(GetAssayData(ser, "data", "RNA"))
  write.csv(norm.matrix, fpath.mat, quote=F)
  
  fpath.meta <- paste0(input_fpath, "metadata.csv")
  cell.meta <- data.frame(Cell = rownames(ser@meta.data), Annotation = ser$celltype)
  write.csv(cell.meta, fpath.meta, quote=F, row.names = FALSE)
  rm(norm.matrix, cell.meta, ser);gc()
  
  start.time <- proc.time()
  
  if(species == 'human'){
    species = 'hsapiens'
  }else{
    species = 'mmusculus'
  }
  
  output_fpath <- paste0(fpath, '/output/')
  if(!dir.exists(output_fpath)){
    dir.create(output_fpath)
  }
  output_fpath <- paste0(output_fpath, 'result.csv')
  
  command <- paste('sh ../Script/Step0_LRToolsFunction/scConnect_shell.sh', 
                   fpath.mat, fpath.meta, species, output_fpath, sep = ' ')
  system(command)
  
  used.time <- proc.time()-start.time
  used.memory <- mem_used()
  
  # handle result
  if(T){
    result <- read.csv(output_fpath)
    result <- result[,-1]
    
    ligands <- read.csv("./Step0_SharedInfo/scConnect/ligands.csv", header = T)
    ligands <- ligands[,-c(1, 3, 9)]
    result <- merge(ligands, result,by = "ligand")
    colnames(result)[2] <- "ligand_gene"
    result <- result[,-c(3:6)]
    
    
    receptors <- read.csv('./Step0_SharedInfo/scConnect/receptors.csv')
    receptors <- dplyr::distinct(receptors, receptor, gene, .keep_all = FALSE)
    result <- merge(receptors, result, by = 'receptor')
    colnames(result)[2] <- 'receptor_gene'
    result <- result[which(result$ligand_pval<0.05), ]
    result <- result[which(result$receptor_pval<0.05),]
    
    result <- result[,c("sender", "reciever", "ligand_gene", "receptor_gene", "score")]
    colnames(result) <- c("Sender", "Receiver", "Ligand", "Receptor", "LRscore")
    result <- result[which(result$Sender!=result$Receiver),]
    result$all <- paste(result$Sender, result$Ligand, result$Receiver, result$Receptor, sep = '_')
    result <- dplyr::distinct(result, all, .keep_all = TRUE)
  }
  
  
  result_record <- list(result=result, 
                        used_time = paste0(round(used.time[3]/60,3), ' min'), 
                        used_memory =  round(used.memory/1024/1024/1024,3))
  return(result_record)
}

###############
## cell2cell ## Only for human (LR prior databases); the LR prior databases can be designed by users (complex)
############### 
cell2cell_function <- function(ser, fpath, species = 'human'){
  suppressMessages(library(Seurat))
  suppressMessages(library(lobstr))
  set.seed(123)
  
  input_fpath <- paste0(fpath, '/input/')
  if(!dir.exists(input_fpath)){
    dir.create(input_fpath)
  }
  
  fpath.mat <- paste0(input_fpath, "counts.csv")
  norm.matrix <- as.matrix(GetAssayData(ser, "data", "RNA"))
  write.csv(norm.matrix, fpath.mat, quote=F)
  
  fpath.meta <- paste0(input_fpath, "metadata.csv")
  cell.meta <- data.frame(Cell = rownames(ser@meta.data), Annotation = ser$celltype)
  write.csv(cell.meta, fpath.meta, quote=F, row.names = FALSE)
  rm(norm.matrix, cell.meta, ser);gc()
  
  output_fpath <- paste0(fpath, '/output/')
  if(!dir.exists(output_fpath)){
    dir.create(output_fpath)
  }
  
  start.time <- proc.time()
  
  if(species == 'human'){
    fpath.lr <- './Step0_SharedInfo/cell2cell/lr_human.csv'
  }else{
    fpath.lr <- './Step0_SharedInfo/cell2cell/lr_human.csv'
  }
  
  command <- paste('sh ../Script/Step0_LRToolsFunction/cell2cell_shell.sh', fpath.mat, fpath.meta, 
                   fpath.lr, output_fpath, sep = ' ')
  system(command)
  used.time <- proc.time()-start.time
  used.memory <- mem_used()
  
  result.pval <- read.csv(paste0(output_fpath, 'ccc_pval.csv'))
  result.pval <- tidyr::pivot_longer(result.pval, -X, names_to = 'sr', values_to = 'pvalue')
  result.pval <- tidyr::separate(result.pval, X, c('Ligand', 'Receptor'), sep = ',')
  result.pval$Ligand <- gsub("\\(\\'", '', result.pval$Ligand)
  result.pval$Ligand <- gsub("\\'", '', result.pval$Ligand)
  result.pval$Receptor <- gsub("\\'\\)", '', result.pval$Receptor)
  result.pval$Receptor <- gsub("\\'", '', result.pval$Receptor)
  result.pval <- tidyr::separate(result.pval, sr, c('Sender', 'Receiver'), sep = '\\.')
  result.pval$all <- paste(result.pval$Sender, result.pval$Ligand, result.pval$Receiver, result.pval$Receptor, sep = '_')
  
  result.value <- read.csv(paste0(output_fpath, 'communication_matrix.csv'))
  result.value <- tidyr::pivot_longer(result.value, -X, names_to = 'sr', values_to = 'value')
  result.value <- tidyr::separate(result.value, X, c('Ligand', 'Receptor'), sep = ',')
  result.value$Ligand <- gsub("\\(\\'", '', result.value$Ligand)
  result.value$Ligand <- gsub("\\'", '', result.value$Ligand)
  result.value$Receptor <- gsub("\\'\\)", '', result.value$Receptor)
  result.value$Receptor <- gsub("\\'", '', result.value$Receptor)
  result.value <- tidyr::separate(result.value, sr, c('Sender', 'Receiver'), sep = '\\.')
  result.value$all <- paste(result.value$Sender, result.value$Ligand, result.value$Receiver, result.value$Receptor, sep = '_')
  result.value <- result.value[,c('all', 'value')]
  
  result <- merge(result.pval, result.value, by = 'all')
  result <- result[which(result$Sender != result$Receiver),]
  result <- dplyr::distinct(result, all, .keep_all = TRUE)
  result <- result[which(result$value>0), ]
  result <- result[which(result$pvalue<0.05),]
  result <- result[,1:5]
  
  result_record <- list(result=result, 
                        used_time = paste0(round(used.time[3]/60,3), ' min'), 
                        used_memory =  round(used.memory/1024/1024/1024,3))
  return(result_record) 
}


############
## Domino ## Parameters setting: max_tf_per_clust = 10, max_rec_per_tf = 10
############ Run in parallel with 30 cores (pyscenic)
Domino_function <- function(ser, fpath, species = 'hg38'){
  suppressMessages(library(Seurat))
  suppressMessages(library(domino))
  suppressMessages(library(dplyr))
  suppressMessages(library(lobstr))
  set.seed(123)
  
  input_fpath <- paste0(fpath, '/input/')
  if(!dir.exists(input_fpath)){
    dir.create(input_fpath)
  }
  
  fpath.mat <- paste0(input_fpath, 'counts.csv')
  write.csv(t(as.matrix(ser@assays$RNA@counts)), file = fpath.mat)
  
  start.time <- proc.time()
  
  if(species == 'hg38'){
    fpath.tfs <- './Step0_SharedInfo/Domino/hg38/hs_hgnc_curated_tfs.txt'
    fpath.feather.1 <- './Step0_SharedInfo/Domino/hg38/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather'
    fpath.feather.2 <- './Step0_SharedInfo/Domino/hg38/hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather'
    fpath.motif <- './Step0_SharedInfo/Domino/hg38/motifs-v9-nr.hgnc-m0.001-o0.0.tbl'
  }else if(species == 'hg38' & grepl('pbmc', fpath)){
    fpath.tfs <- './Step0_SharedInfo/Domino/hg38/hs_hgnc_curated_tfs.txt'
    fpath.feather.1 <- './Step0_SharedInfo/Domino/hg38/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather'
    fpath.feather.2 <- './Step0_SharedInfo/Domino/hg38/hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather'
    fpath.motif <- './Step0_SharedInfo/Domino/hg38/motifs-v9-nr.hgnc-m0.001-o0.0.tbl'
  }else if(species == 'hg19'){
    fpath.tfs <- './Step0_SharedInfo/Domino/hg38/hs_hgnc_curated_tfs.txt'
    fpath.feather.1 <- './Step0_SharedInfo/Domino/hg19/hg19-500bp-upstream-10species.mc9nr.genes_vs_motifs.rankings.feather'
    fpath.feather.2 <- './Step0_SharedInfo/Domino/hg19/hg19-500bp-upstream-7species.mc9nr.genes_vs_motifs.rankings.feather'
    fpath.feather.3 <- './Step0_SharedInfo/Domino/hg19/hg19-tss-centered-10kb-10species.mc9nr.genes_vs_motifs.rankings.feather'
    fpath.feather.4 <- './Step0_SharedInfo/Domino/hg19/hg19-tss-centered-10kb-7species.mc9nr.genes_vs_motifs.rankings.feather'
    fpath.feather.5 <- './Step0_SharedInfo/Domino/hg19/hg19-tss-centered-5kb-10species.mc9nr.genes_vs_motifs.rankings.feather'
    fpath.feather.6 <- './Step0_SharedInfo/Domino/hg19/hg19-tss-centered-5kb-7species.mc9nr.genes_vs_motifs.rankings.feather'
    fpath.motif <- './Step0_SharedInfo/Domino/hg38/motifs-v9-nr.hgnc-m0.001-o0.0.tbl'
  }else if(species == 'mm9'){
    fpath.tfs <- './Step0_SharedInfo/Domino/mm9/mm_mgi_tfs.txt'
    fpath.feather.1 <- './Step0_SharedInfo/Domino/mm9/mm9-500bp-upstream-10species.mc9nr.genes_vs_motifs.rankings.feather'
    fpath.feather.2 <- './Step0_SharedInfo/Domino/mm9/mm9-500bp-upstream-7species.mc9nr.genes_vs_motifs.rankings.feather'
    fpath.feather.3 <- './Step0_SharedInfo/Domino/mm9/mm9-tss-centered-10kb-10species.mc9nr.genes_vs_motifs.rankings.feather'
    fpath.feather.4 <- './Step0_SharedInfo/Domino/mm9/mm9-tss-centered-10kb-7species.mc9nr.genes_vs_motifs.rankings.feather'
    fpath.feather.5 <- './Step0_SharedInfo/Domino/mm9/mm9-tss-centered-5kb-10species.mc9nr.genes_vs_motifs.rankings.feather'
    fpath.feather.6 <- './Step0_SharedInfo/Domino/mm9/mm9-tss-centered-5kb-7species.mc9nr.genes_vs_motifs.rankings.feather'
    fpath.motif <- './Step0_SharedInfo/Domino/mm9/motifs-v9-nr.mgi-m0.001-o0.0.tbl'
  }else if(species == 'mm10'){
    fpath.tfs <- './Step0_SharedInfo/Domino/mm9/mm_mgi_tfs.txt'
    fpath.feather.1 <- './Step0_SharedInfo/Domino/mm10/mm10_refseq-r80_10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather'
    fpath.feather.2 <- './Step0_SharedInfo/Domino/mm10/mm10_refseq-r80_500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather'
    fpath.motif <- './Step0_SharedInfo/Domino/mm9/motifs-v9-nr.mgi-m0.001-o0.0.tbl'
  }

  fpath.adj <- paste0(input_fpath, 'adj.tsv')
  fpath.reg <- paste0(input_fpath, 'reg.csv')
  fpath.auc <- paste0(input_fpath, 'auc_mtx.csv')
  
  if(species == 'hg19' | species == 'mm9'){
    command <- paste('sh ../Script/Step0_LRToolsFunction/Domino_shell1.sh', 
                     fpath.adj, fpath.mat, fpath.reg, fpath.auc, 
                     fpath.tfs, fpath.feather.1, fpath.feather.2, fpath.feather.3,
                     fpath.feather.4, fpath.feather.5, fpath.feather.6,
                     fpath.motif, sep = ' ')
  }else if(species == 'hg38' | species == 'mm10'){
    command <- paste('sh ../Script/Step0_LRToolsFunction/Domino_shell2.sh', 
                     fpath.adj, fpath.mat, fpath.reg, fpath.auc, 
                     fpath.tfs, fpath.feather.1, fpath.feather.2,
                     fpath.motif, sep = ' ')
  }
  
  system(command)
  
  # run domino
  if(T){
    auc = t(read.table(fpath.auc, header = TRUE, row.names = 1, stringsAsFactors = FALSE, sep = ','))
    
    ser <- ser[, colnames(auc)]
    ser <- ScaleData(ser, features = rownames(ser))
    
    source("./Step0_SharedInfo/Domino/create_domino.R") 
    # line 200 of rscript of creat_domino function: change tf_genes = df[row, 10] to tf_genes = df[row, 9]
    dom = create_domino(signaling_db = './Step0_SharedInfo/CellPhoneDB', 
                        features = auc, counts = ser@assays$RNA@counts, z_scores = ser@assays$RNA@scale.data, 
                        clusters = ser@active.ident, df = fpath.reg)
    dom = build_domino(dom, min_tf_pval = 0.05, max_tf_per_clust = 10, max_rec_per_tf = 10)
  }
  
  used.time <- proc.time()-start.time
  used.memory <- mem_used()
  
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
    lig_to_rec <- lig_to_rec[which(lig_to_rec$ligand!=""),]
  }
  
  # handle result
  if(T){
    result <- lapply(levels(dom@clusters), function(reciever){
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
        
        lig_rec_tf <- merge(lig_to_rec, rec_tf, by = "receptor")
        lig_rec_tf <- merge(sender.ligands, lig_rec_tf, by = "ligand")
        lig_rec_tf$reciever <- reciever
        lig_rec_tf$tf <- stringr::str_replace_all(lig_rec_tf$tf, "\\...", "")
      }else{
        lig_rec_tf <- NA
      }
      lig_rec_tf
    })
    
    result[which(is.na(result))] <- NULL
    
    result <- do.call(rbind, result)
    result <- result[which(result$sender != result$reciever),]
    result <- result[,-4]
    colnames(result) <- c('Ligand', 'Sender', 'Receptor', 'Receiver')
    result$all <- paste(result$Sender, result$Ligand, result$Receiver, result$Receptor, sep = '_')
    result <- dplyr::distinct(result, all, .keep_all = TRUE)
  }
  
  result_record <- list(result=result, 
                        used_time = paste0(round(used.time[3]/60,3), ' min'), 
                        used_memory =  round(used.memory/1024/1024/1024,3))
  return(result_record)
}

#############
## scMLnet ## 
#############
scMLnet_function <- function(ser, sender = NULL, receiver = NULL){
  pacman::p_unload(Seurat)
  suppressMessages(library(Seurat, lib.loc = '/home/ljx/R/x86_64-redhat-linux-gnu-library/library2/'))
  suppressMessages(library(scMLnet))
  suppressMessages(library(lobstr))
  set.seed(123)
  
  start.time <- proc.time()
  
  # scMLnet v0.2.0
  # pacakge: https://github.com/SunXQlab/scMLnet2.0
  GCMat<- GetAssayData(ser, "counts", "RNA")
  BarCluTable <- data.frame(Barcode = colnames(ser), Cluster = ser$celltype)
  types <- unique(BarCluTable$Cluster)
  
  LigRecLib <- read.table("./Step0_SharedInfo/scMLnet/LigRec.txt", header = T)
  colnames(LigRecLib)[2:3] <- c("source", "target")
  TFTarLib <- read.table("./Step0_SharedInfo/scMLnet/TFTargetGene.txt", header = T)
  colnames(TFTarLib)[1:2] <- c("source", "target")
  RecTFLib <- read.table("./Step0_SharedInfo/scMLnet/RecTF.txt", header = T)
  colnames(RecTFLib)[1:2] <- c("source", "target")
  
  if(is.null(sender) & is.null(receiver)){
    sender_ct <- types
    receiver_ct <- types
  }else{
    sender_ct <- sender
    receiver_ct <- receiver
  }
  
  result <- list()
  for (LigClu in sender_ct) {
    for (RecClu in receiver_ct[which(receiver_ct != LigClu)]) {
      netList <- tryCatch(RunMLnet(data = GCMat, BarCluTable = BarCluTable, 
                                   RecClu = RecClu, LigClu = LigClu,
                                   LigRec.DB = LigRecLib, 
                                   TFTG.DB = TFTarLib,
                                   RecTF.DB = RecTFLib),
                          error=function(e){NA}
      )
      list.names <- paste(LigClu, RecClu, sep = "_")
      result[[list.names]] <- netList
    }
  }
  
  used.time <- proc.time()-start.time
  used.memory <- mem_used()
  
  result[which(is.na(result))] <- NULL
  
  result <- lapply(result, function(res){
    res <- res$LigRec
    res
  })
  result <- do.call(rbind, result)
  if(!is.null(result)){
    result <- tibble::rownames_to_column(result, 'sr')
    result$sr <- gsub('\\.[0-9]+', '', result$sr)
    result <- tidyr::separate(result, sr, c('Sender', 'Receiver'), sep = '_')
    result <- result[,-5]
    colnames(result)[3:4] <- c('Ligand', "Receptor")
    result$all <- paste(result$Sender, result$Ligand, result$Receiver, result$Receptor, sep = '_')
    result <- dplyr::distinct(result, all, .keep_all = TRUE)
    result <- result[which(result$Sender!=result$Receiver),]
  }
  
  result_record <- list(result=result, 
                        used_time = paste0(round(used.time[3]/60,3), ' min'), 
                        used_memory =  round(used.memory/1024/1024/1024,3))
  return(result_record)
}

#############
## PyMINEr ## 
#############
PyMINEr_function <- function(ser, fpath, species = 'human'){
  suppressMessages(library(Seurat))
  suppressMessages(library(lobstr))
  set.seed(123)
  
  # input
  if(T){
    input_fpath <- paste0(fpath, '/input/')
    if(!dir.exists(input_fpath)){
      dir.create(input_fpath)
    }
    
    fpath.mat <- paste0(input_fpath, 'counts.txt')
    matrix.sc <- GetAssayData(ser, "data", "RNA")
    matrix.sc <- as.data.frame(matrix.sc)
    matrix.sc <- tibble::rownames_to_column(matrix.sc, var = "gene")
    write.table(matrix.sc, file = fpath.mat, row.names = FALSE, sep = "\t")
    
    fpath.meta <- paste0(input_fpath, 'metadata.txt')
    meta.sc <- data.frame(barcode = colnames(ser), celltype = ser$celltype)
    i <- 0
    for(ct in unique(meta.sc$celltype)){
      meta.sc[which(meta.sc$celltype == ct), "ct_num"] <- i
      i <- i+1
    }
    meta <- meta.sc[,2:3]
    meta <- dplyr::distinct(meta)
    rownames(meta) <- NULL
    
    meta.sc <- meta.sc[,-2]
    write.table(meta.sc, file = fpath.meta,row.names = F, col.names = F, sep = "\t")
    rm(matrix.sc, meta.sc, ct, i, ser);gc()
  }
  
  start.time <- proc.time()
  
  if(species == 'human'){
    species <- 'hsapiens'
  }else{
    species <- 'mmusculus'
  }
  
  command <- paste('sh ../Script/Step0_LRToolsFunction/PyMINEr_shell.sh', 
                   fpath.mat, fpath.meta, species, sep = ' ')
  system(command)
  
  used.time <- proc.time()-start.time
  used.memory <- mem_used()
  
  result_path <- paste0(input_fpath, 
                        'autocrine_paracrine_signaling/extracellular_plasma_membrane_cell_type_specific_interactions.txt')
  result <- read.table(result_path, header = FALSE, sep = ',', fileEncoding = "utf-8")
  result <- tidyr::separate(result, V1, paste0('V', 2:19), sep = '\t')
  result <- result[,-c(12:17)]
  result[1, 12] <- 'LRscore'
  colnames(result) <- result[1,]
  result <- result[-1,]
  
  result <- result[which(result$cell_type_1!=result$cell_type_2),]
  result <- result[,c(2,4,6,8,12)]
  meta$ct_num <- paste0('bool_sample_group_', meta$ct_num)
  result <- merge(result, meta, by.x = 'cell_type_1', by.y = 'ct_num')
  colnames(result)[6] <- 'Sender'
  result <- merge(result, meta, by.x = 'cell_type_2', by.y = 'ct_num')
  colnames(result)[7] <- 'Receiver'
  result <- result[,-c(1:2)]
  colnames(result)[1:2] <- c('Ligand', 'Receptor')
  result$all <- paste(result$Sender, result$Ligand, result$Receiver, result$Receptor, sep = '_')
  result <- dplyr::distinct(result, all, .keep_all = TRUE)
  
  result_record <- list(result=result, 
                        used_time = paste0(round(used.time[3]/60,3), ' min'), 
                        used_memory =  round(used.memory/1024/1024/1024,3))
  return(result_record)
}

###############
## RNAMagnet ## 
############### Only for mouse (LR prior databases)
RNAMagnet_function <- function(ser){
  suppressMessages(library(RNAMagnet))
  suppressMessages(library(Seurat))
  suppressMessages(library(lobstr))
  set.seed(123)
  
  start.time <- proc.time()
  
  result <- RNAMagnetSignaling(ser, .version = "2.0.0", 
                               .cellularCompartment = c("Membrane", "ECM", "Both", "Secreted"))
  
  temp <- c()
  for (sender in unique(ser$celltype)) {
    for (receiver in unique(ser$celltype)) {
      tmp <- getRNAMagnetGenes(result, sender, receiver)
      if(dim(tmp)[[1]]!=0){
        tmp$Sender <- sender
        tmp$Receiver <- receiver
        temp <- rbind(temp,tmp)
      }
    }
  }
  
  used.time <- proc.time()-start.time
  used.memory <- mem_used()
  
  result <- temp
  result <- result[which(result$Sender!=result$Receiver),]
  result <- tidyr::separate(result, pair, c('Ligand', 'Receptor'), sep = '-')
  rownames(result) <- NULL
  need <- which(grepl('\\|', result$Ligand) | grepl('\\|', result$Receptor))
  if(length(need)>0){
    result_sub1 <- result[need, ]
    result_sub2 <- result[-need, ]
    
    result_tmp <- c()
    for (i in seq(need)) {
      ligands <- unlist(stringr::str_split(result_sub1[i, 'Ligand'], '\\|'))
      receptors <- unlist(stringr::str_split(result_sub1[i, 'Receptor'], '\\|'))
      tmp <- expand.grid(ligands, receptors)
      colnames(tmp) <- c('Ligand', 'Receptor')
      tmp$Sender <- rep(result_sub1[i, 'Sender'], n = dim(tmp)[1])
      tmp$Receiver <- rep(result_sub1[i, 'Receiver'], n = dim(tmp)[1])
      tmp$score <- rep(result_sub1[i, 'score'], n = dim(tmp)[1])
      tmp <- tmp[, c(5, 1:4)]
      result_tmp <- rbind(result_tmp, tmp)
    }
    result <- rbind(result_tmp, result_sub2)
  }
  
  result <- aggregate(score~Ligand+Receptor+Sender+Receiver, mean, data=result)
  result$all <- paste(result$Sender, result$Ligand, result$Receiver, result$Receptor, sep = '_')
  colnames(result)[5] <- 'LRscore'
  
  
  result_record <- list(result=result, 
                        used_time = paste0(round(used.time[3]/60,3), ' min'), 
                        used_memory =  round(used.memory/1024/1024/1024,3))
  return(result_record)
}
