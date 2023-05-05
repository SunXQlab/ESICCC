############## Run in parallel with 10 cores
## CytoTalk ## Run in R 4.1.0 version
############## Filter the spots with the number of cells less than 3
CytoTalk_function <- function(ser, fpath, sender = NULL, receiver = NULL, species = 'human'){
  suppressMessages(library(CytoTalk))
  suppressMessages(library(Seurat))
  suppressMessages(library(tidyverse))
  suppressMessages(library(igraph))
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
   
    norm.matrix <- as.matrix(GetAssayData(ser, "data", "SCT"))
    write.csv(norm.matrix, fpath.mat, quote=F)
    
    meta.data <- data.frame(rownames(ser@meta.data), ser@meta.data$celltype)
    colnames(meta.data) <- NULL
    meta.data <- as.matrix(meta.data)
    write.csv(meta.data, fpath.meta, quote=F, row.names = FALSE)
    rm(norm.matrix, meta.data);gc()
  }
  
  if(species == 'human'){
    pcg = CytoTalk::pcg_human
    lrp = CytoTalk::lrp_human
  }else{
    pcg = CytoTalk::pcg_mouse
    lrp = CytoTalk::lrp_mouse
  }
  
  lst.sc <- read_matrix_with_meta(fpath.mat, fpath.meta)
  
  sender_ct <- unique(ser$celltype)
  remove_ct <- as.character(names(which(table(ser$celltype)<3)))
  if (length(remove_ct)>0) {
    sender_ct <- sender_ct[-which(sender_ct %in% remove_ct)]
  }
  if(is.null(receiver)){
    receiver_ct <- unique(ser$celltype)
  }else if(!is.null(receiver)){
    receiver_ct <- receiver
  }
  rm(sender, receiver);gc()
  
  result <- list()
  for (sender in sender_ct) {
    for (receiver in receiver_ct[which(receiver_ct!=sender)]) {
      tmp <- tryCatch(CytoTalk::run_cytotalk(lst.sc, sender, receiver, cores = 10, 
                                             cutoff_a = 0.05, cutoff_b = 0.05,
                                             pcg = pcg, lrp = lrp),
                      error=function(e){NA}
      )
      cp <- paste(sender, receiver, sep = "_")
      result[[cp]] <- tmp
    }
  }
  
  result[which(is.na(result))] <- NULL
  saveRDS(result, file = paste0(fpath, '/', 'result.rds'))
  
  
  for (receiver in receiver_ct) {
    result.tmp <- lapply(names(result), function(cp){
      sender <-  gsub("_.*", "",cp)
      
      tmp <- list()
      res <- result[[cp]]
      if(!is.null(res$pathways)){
        network <- res$pcst$final_network
        network <- network[!(network$node1_type == sender & network$node2_type == sender), ]
        network <- network[!(network$node1_type == receiver & network$node2_type == sender), ]
        network$node1 <- toupper(network$node1)
        network$node2 <- toupper(network$node2)
        network$node1 <- gsub("ORF", "orf", network$node1)
        network$node2 <- gsub("ORF", "orf", network$node2)
        LR <- network[which(network$node1_type == sender & network$node2_type == receiver), ]
        LR <- paste(LR$node1, LR$node2, sep = "_")
        Ligand <- network[which(network$node1_type == sender & network$node2_type == receiver), "node1"] %>%
          unique()
        Receptor <- network[which(network$node1_type == sender & network$node2_type == receiver), "node2"] %>%
          unique()
        Target <- c(network[which(network$node1_type == receiver & network$node2_type == receiver), "node1"],
                    network[which(network$node1_type == receiver & network$node2_type == receiver), "node2"]) %>%
          unique()
        Target <- setdiff(Target, c(Ligand, Receptor))
        tmp <- list(Edge = network,
                    LR = LR,
                    Ligand = Ligand, 
                    Receptor = Receptor, 
                    Target = Target)
      }else{
        tmp <- NA
      }
      tmp
    })
    
    names(result.tmp) <- names(result);
    result.tmp[which(is.na(result.tmp))] <- NULL
    
    edge_list <- lapply(result.tmp, function(res){
      res <- res$Edge[,c(1,2)]
      colnames(res) <- c('from', 'to')
      res
    })
    edge_df <- do.call(rbind, edge_list)
    edge_df <- distinct(edge_df, from, to)
    intranetwork <- graph_from_edgelist(as.matrix(edge_df), directed = FALSE)
    
    Ligand <- lapply(result.tmp, function(res){res$Ligand}) %>% unlist() %>% unique()
    Receptor <- lapply(result.tmp, function(res){res$Receptor}) %>% unlist() %>% unique()
    Target <- lapply(result.tmp, function(res){res$Target}) %>% unlist() %>% unique()
    
    distance <- distances(intranetwork, 
                          v = c(Ligand,Receptor), 
                          to = Target)
    distance <- reshape2::melt(distance)
    distance <- distance[!is.infinite(distance$value),]
    colnames(distance) <- c("regulon", 'target',"value")
    distance$type <- ifelse(distance$regulon %in% Ligand,"ligand","receptor")
    
    saveRDS(distance, file = paste0(fpath, "/", receiver, '_result.rds'))
  }
}

##############
## NicheNet ## Filter the spots with the number of cells less than 3
############## 
NicheNet_function <- function(ser, fpath, sampleID, sender = NULL, receiver = NULL){
  suppressMessages(library(nichenetr))
  suppressMessages(library(Seurat))
  suppressMessages(library(tidyverse))
  set.seed(123)
 
  # load data
  if(T){
    ligand_target_matrix = readRDS("./Step0_SharedInfo/NicheNet/ligand_target_matrix.rds")
    
    lr_network = readRDS("./Step0_SharedInfo/NicheNet/lr_network.rds")
    
    weighted_networks = readRDS("./Step0_SharedInfo/NicheNet/weighted_networks.rds")
    weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network %>% distinct(from,to), by = c("from","to"))
  }
  
  sender_ct <- unique(ser$celltype)
  remove_ct <- as.character(names(which(table(ser$celltype)<3)))
  if (length(remove_ct)>0) {
    sender_ct <- sender_ct[-which(sender_ct %in% remove_ct)]
  }
  if(is.null(receiver)){
    receiver_ct <- unique(ser$celltype)
  }else if(!is.null(receiver)){
    receiver_ct <- receiver
  }
  rm(sender, receiver);gc()
  
  for (receiver in receiver_ct) {
    # define the sender genes
    expressed_genes_sender <- lapply(sender_ct[which(sender_ct != receiver)], function(ct){
      expressed_genes_sender_ct <- get_expressed_genes(ct, ser, pct = 0.05)
      expressed_genes_sender_ct
    })
    expressed_genes_sender = expressed_genes_sender %>% unlist() %>% unique()
    
    # define the receiver genes
    expressed_genes_receiver = get_expressed_genes(receiver, ser, pct = 0.05)
    background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
    
    # define the genes of interest
    if(grepl('CID', sampleID)){
      icgs_fpath <- paste0('~/1-Datasets/NatureGenetics_2021_BC/ScriptForCCC/Giotto_result/', sampleID, '_icgs.rds')
    }else if(grepl('UKF', sampleID)){
      icgs_fpath <- paste0('~/1-Datasets/CancerCell_2022_Glioma/ScriptForLRTBenchmark/Giotto_result/', sampleID, '_icgs.rds')
    }else{
      icgs_fpath <- paste0('~/1-Datasets/Nature_2022_MI/ScriptForCCC/Giotto_result/', sampleID,'_icgs.rds')
    }
    geneset =  readRDS(file = icgs_fpath)
    geneset <- geneset[[receiver]] %>% unlist() %>% unique()
    geneset = geneset %>% .[. %in% rownames(ligand_target_matrix)]
    
    # define the potential ligands
    ligands = lr_network %>% pull(from) %>% unique()
    receptors = lr_network %>% pull(to) %>% unique()
    
    expressed_ligands = intersect(ligands,expressed_genes_sender)
    expressed_receptors = intersect(receptors,expressed_genes_receiver)
    
    potential_ligands = lr_network %>% 
      filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% 
      pull(from) %>% 
      unique()
    
    # Perform NicheNet ligand activity analysis
    ligand_activities = predict_ligand_activities(geneset = geneset, 
                                                  background_expressed_genes = background_expressed_genes, 
                                                  ligand_target_matrix = ligand_target_matrix, 
                                                  potential_ligands = potential_ligands)
    ligand_activities = ligand_activities %>% arrange(-pearson) %>% mutate(rank = rank(desc(pearson)))
    
    # Get all the ligands for downstreams analysis
    best_upstream_ligands = ligand_activities %>% arrange(-pearson) %>% pull(test_ligand) %>% unique()
    head(best_upstream_ligands)
    
    # Get the active ligand-target links
    active_ligand_target_links_df = best_upstream_ligands %>% 
      lapply(get_weighted_ligand_target_links,
             geneset = geneset, 
             ligand_target_matrix = ligand_target_matrix, 
             n = nrow(ligand_target_matrix)) %>% 
      bind_rows() %>% drop_na()
    colnames(active_ligand_target_links_df) <- c("regulon", 'target',"value")
    active_ligand_target_links_df$type <- 'ligand'
    
    saveRDS(active_ligand_target_links_df, file = paste0(fpath, '/', receiver, '_result.rds'))
  }
}


############
## MISTy  ## Filter the spots with the number of cells less than 3
############ 
MISTy_function <- function(ser, fpath, sampleID, sender = NULL, receiver=NULL){
  # load package
  if(T){
    # MISTy
    suppressMessages(library(mistyR))
    suppressMessages(library(future))
    suppressMessages(library(future))
    
    # Seurat
    suppressMessages(library(Seurat))
    suppressMessages(library(Giotto))
    
    # data manipulation
    suppressMessages(library(Matrix))
    suppressMessages(library(tibble))
    suppressMessages(library(dplyr))
    suppressMessages(library(purrr))
    
    # normalization
    suppressMessages(library(sctransform))
    
    # resource
    suppressMessages(library(progeny))
    
    # setup parallel execution
    options(future.globals.maxSize = 1024^3)
    plan(multisession)
    set.seed(123)
  }
  
  source('./Step0_SharedInfo/MistyR/code.R')
  
  if(is.null(receiver)){
    receiver_ct <- unique(ser$celltype)
  }else if(!is.null(receiver)){
    receiver_ct <- receiver
  }
  rm(receiver);gc()
  
  # run normalization
  sct.data <- vst(GetAssayData(
    object = ser,
    slot = "counts",
    assay = "Spatial"
  ),
  verbosity = 0
  )
  
  ser[["SCT"]] <- CreateAssayObject(data = sct.data$y)
  
  gene.expression <- GetAssayData(ser, assay = "SCT")
  coverage <- rowSums(gene.expression > 0) / ncol(gene.expression)
  slide.markers <- names(which(coverage >= 0.05))
  
  Databases <- readRDS('./Step0_SharedInfo/stMLnet/Databases.rds')
  ligands <- Databases$LigRec.DB$source %>% unique() %>% .[. %in% slide.markers]
  receptors <- Databases$LigRec.DB$target %>% unique() %>% .[. %in% slide.markers]
  
  # load ICGs
  if(grepl('CID', sampleID)){
    icgs_fpath <- paste0('~/1-Datasets/NatureGenetics_2021_BC/ScriptForCCC/Giotto_result/', sampleID, '_icgs.rds')
  }else if(grepl('UKF', sampleID)){
    icgs_fpath <- paste0('~/1-Datasets/CancerCell_2022_Glioma/ScriptForLRTBenchmark/Giotto_result/', sampleID, '_icgs.rds')
  }else{
    icgs_fpath <- paste0('~/1-Datasets/Nature_2022_MI/ScriptForCCC/Giotto_result/', sampleID,'_icgs.rds')
  }
  icgs <- readRDS(icgs_fpath)
  
  celltypes <- ser$celltype %>% as.character()
  remove_ct <- as.character(names(which(table(ser$celltype)<3)))
  if (length(remove_ct)>0) {
    celltypes <- celltypes[-which(celltypes %in% remove_ct)]
  }
  
  geometry <- ser@images$image@coordinates %>% select(., row,col)
  
  output_fpath <- paste0(fpath, '/output')
  if(!dir.exists(output_fpath)){
    dir.create(output_fpath)
  }
  
  #########
  ## run ##
  #########
  
  for (receiver in receiver_ct) {
    ct <- receiver
    message(paste0('running jobs:',ct))
    targets <- icgs[[ct]] %>% unlist() %>% unique()
    
    ###############
    ## parameter ##
    ###############
    
    view.assays <- list(
      "main" = "SCT",
      "ligand" = "SCT"
    )
    
    # Define features for each view
    view.features <- list(
      "main" = c(targets,receptors),
      "ligand" = ligands
    )
    
    # Define spatial context for each view
    view.types <- list(
      "main" = "intra",
      "ligand" = "para"
    )
    
    # Define additional parameters (l in the case of paraview)
    view.params <- list(
      "main" = NULL,
      "ligand" = 10
    )
    
    # Define specific properties for each view
    view.properties <- list(
      "main" = ifelse(celltypes == ct,1,0),
      "ligand" = ifelse(celltypes != ct,1,0)
    )
    
    #########
    ## Run ##
    #########
    
    spot.ids = NULL
    out.alias = paste0(output_fpath, "/results_",ct,"_paraview_10")
    
    # Extracting data
    view.data <- map(view.assays,
                     extract_seurat_data,
                     geometry = geometry,
                     visium.slide = ser
    )
    str(view.data,max.level = 1)
    
    # Adding all spots ids in case they are not defined
    if (is.null(spot.ids)) {
      spot.ids <- rownames(view.data[[1]])
    }
    
    # First filter the features from the data
    view.data.filt <- map2(view.data, view.features, filter_data_features)
    str(view.data.filt,max.level = 1)
    view.data.filt[[1]][1:4,1:1]
    view.data.filt[[2]][1:4,1:1]
    
    # specific properties: celltype
    view.data.spec <- map2(view.data.filt, view.properties, add_specific_properties)
    str(view.data.spec,max.level = 1)
    view.data.spec[[1]][1:6,1:2]
    view.data.spec[[2]][1:6,1:2]
    
    # Create initial view
    views.main <- create_initial_view(view.data.spec[[1]] %>%
                                        rownames_to_column() %>%
                                        filter(rowname %in% spot.ids) %>%
                                        select(-rowname))
    str(views.main,max.level = 2)
    
    # Create other views
    view.names <- names(view.data.spec)
    
    all.views <- pmap(list(
      view.data.filt[-1],
      view.types[-1],
      view.params[-1],
      view.names[-1]
    ),
    create_default_views,
    spot.ids = spot.ids,
    geometry = geometry
    )
    str(all.views,max.level = 2)
    
    pline.views <- add_views(
      views.main,
      unlist(all.views, recursive = FALSE)
    )
    
    # Run MISTy
    run_misty(pline.views, out.alias)
    misty.results <- collect_results(out.alias)
    
    ############
    ## output ##
    ############
    
    misty_score <- misty.results$importances.aggregated %>% na.omit()
    misty_score <- misty_score %>% 
      select(Predictor,Target,Importance,view) %>%
      filter(Target %in% targets, Predictor %in% c(ligands,receptors)) %>%
      rename(regulon=Predictor,target=Target,value=Importance,type=view)
    misty_score$type <- gsub('_10','',misty_score$type)
    misty_score$type <- gsub('intra','receptor',misty_score$type)
    
    saveRDS(misty_score, paste0(fpath, '/', receiver, '_result.rds'))
  }
}


#############
## stMLnet ## Filter the spots with the number of cells less than 3
############# 
stMLnet_function <- function(ser, fpath, sampleID, sender = NULL, receiver = NULL){
  suppressMessages(library(Seurat))
  suppressMessages(library(stMLnet))
  suppressMessages(library(tidyverse))
  suppressMessages(library(doSNOW))
  set.seed(123)
  
  sender_ct <- unique(ser$celltype)
  if(is.null(receiver)){
    receiver_ct <- unique(ser$celltype)
  }else if(!is.null(receiver)){
    receiver_ct <- receiver
  }
  rm(receiver);gc()
  
  # load input data
  if(T){
    norm.matrix <- as.matrix(GetAssayData(ser, "data", "SCT"))
    metadata <- data.frame(Barcode = colnames(ser), Cluster = ser$celltype)
    
    loc <- ser@images$image@coordinates %>% select(., row,col)
    
    ex_databases <- readRDS("./Step0_SharedInfo/stMLnet/Databases.rds") 
    
    # load ICGs
    if(grepl('CID', sampleID)){
      icgs_fpath <- paste0('~/1-Datasets/NatureGenetics_2021_BC/ScriptForCCC/Giotto_result/', sampleID, '_icgs.rds')
    }else if(grepl('UKF', sampleID)){
      icgs_fpath <- paste0('~/1-Datasets/CancerCell_2022_Glioma/ScriptForLRTBenchmark/Giotto_result/', sampleID, '_icgs.rds')
    }else{
      icgs_fpath <- paste0('~/1-Datasets/Nature_2022_MI/ScriptForCCC/Giotto_result/', sampleID,'_icgs.rds')
    }
    ICG_list <- readRDS(icgs_fpath)
    
    # define highly differentially expressed ligands
    if(T){
      remove_ct <- as.character(names(which(table(ser$celltype)<3)))
      lig_clusters <- ser@active.ident %>% as.character() %>% unique()
      if (length(remove_ct)>0) {
        lig_clusters <- lig_clusters[-which(lig_clusters %in% remove_ct)]
      }
      
      if(F){
        ## ligands
        ligs_in_db <- ex_databases$LigRec.DB$source %>% unique()
        ligs_in_db <- intersect(ligs_in_db, rownames(ser))
        df_markers_ligs <- lapply(lig_clusters, function(cluster){
          print(cluster)
          df <- FindMarkers(ser, ident.1 = cluster, features = ligs_in_db, 
                            only.pos = T, min.pct = 0.05)
          if(dim(df)[[1]]!=0){
            df$gene <- rownames(df)
            df$ident.1 <- cluster
          }else{
            df <- NA
          }
          df
        })
        df_markers_ligs[which(is.na(df_markers_ligs))] <- NULL
        df_markers_ligs <- do.call(rbind, df_markers_ligs)
        Ligs_up_list <- split(df_markers_ligs$gene,df_markers_ligs$ident.1)
        str(Ligs_up_list)
      }
    }
    
    # define highly differentially expressed receptors
    if(T){
      BarCluTable <- data.frame(Barcode = rownames(ser@meta.data), 
                                Cluster = ser@meta.data$celltype)
      
      ## parameters
      expr.ct <- 0.1
      pct.ct <- 0.05
      
      ## receptors in prior database 
      recs_in_db <- ex_databases$LigRec.DB$target %>% unique()
      
      ## calculate mean and pct
      rec_clusters <- receiver_ct
      #clusters <- BarCluTable$Cluster %>% as.character() %>% unique()
      
      meanExpr_of_LR <- lapply(rec_clusters, function(cluster){
        
        cluster.ids <- BarCluTable$Barcode[BarCluTable$Cluster == cluster]
        source_mean <- rowMeans(norm.matrix[,cluster.ids])
        names(source_mean) <- rownames(norm.matrix)
        source_mean
        
      }) %>% do.call('cbind',.) %>% as.data.frame()
      colnames(meanExpr_of_LR) <- rec_clusters
      
      pct_of_LR <- lapply(rec_clusters, function(cluster){
        
        cluster.ids <- BarCluTable$Barcode[BarCluTable$Cluster == cluster]
        dat <- norm.matrix[,cluster.ids]
        pct <- rowSums(dat>0)/ncol(dat)
        names(pct) <- rownames(norm.matrix)
        pct
        
      }) %>% do.call('cbind',.) %>% as.data.frame()
      colnames(pct_of_LR) <- rec_clusters
      
      Recs_expr_list <- lapply(rec_clusters, function(cluster){
        
        recs <- rownames(norm.matrix)[meanExpr_of_LR[,cluster] >= expr.ct & pct_of_LR[,cluster] >= pct.ct]
        intersect(recs, recs_in_db)
        
      })
      names(Recs_expr_list) <- rec_clusters
      str(Recs_expr_list)
    }
    
  }
  
  output_fpath <- paste0(fpath, '/output/')
  if(!dir.exists(output_fpath)){
    dir.create(output_fpath)
  }
  setwd(output_fpath)
  
  # step1 create mulityayer network
  
  quan.cutoff <- 0.98
  Databases <- ex_databases
  Databases$RecTF.DB <- Databases$RecTF.DB %>% 
    .[.$score > quantile(.$score, quan.cutoff), ] %>% 
    dplyr::distinct(source, target)
  Databases$LigRec.DB <- Databases$LigRec.DB %>% 
    dplyr::distinct(source, target) %>% 
    dplyr::filter(target %in% Databases$RecTF.DB$source)
  Databases$TFTG.DB <- Databases$TFTG.DB %>% 
    dplyr::distinct(source, target) %>% 
    dplyr::filter(source %in% Databases$RecTF.DB$target)
  
  resMLnet <- runMLnet(ExprMat = norm.matrix, AnnoMat = metadata, Normalize = F, 
                       LigClus = lig_clusters, RecClus = rec_clusters,
                       logfc.ct = 0.1, pct.ct = 0.05, expr.ct = 0.1,
                       ProjectName = '2023-1', Databases = Databases,
                       TGList=ICG_list, RecList=Recs_expr_list) #LigList=Ligs_up_list, 
  rm(ICG_list, Ligs_up_list, Recs_expr_list, ser);gc()
  
  # step2 calculate Signal Activity
  if(T){
    ## calculate distant
    DistMat <- as.matrix(dist(loc))
    colnames(DistMat) <- rownames(loc)
    rownames(DistMat) <- rownames(loc)
    
    ## imputation
    exprMat.Impute <- runImputation(exprMat = norm.matrix)
    rm(norm.matrix);gc()
    
    # remove empty list
    ex_mulnetlist <- list()
    for (i in 1:length(resMLnet$mlnets)) {
      mlnets <- resMLnet$mlnets[[i]]
      for (j in 1:length(mlnets)) {
        mlnet <- mlnets[[j]]
        if(nrow(mlnet$LigRec)!=0) ex_mulnetlist[[names(mlnets)[j]]] = mlnet
      }
    }
    rm(i, j, mlnets, mlnet,loc);gc()
    
    Sender <- NULL
    for (receiver in receiver_ct) {
      resSigActList_tme <- getSiganlActivity(ExprMat = exprMat.Impute,
                                             DistMat = DistMat,
                                             AnnoMat = metadata,
                                             MulNetList = ex_mulnetlist,
                                             Receiver = receiver, Sender = Sender,
                                             Downsample = FALSE, ProjectName = '2023-1')
    }
  }
  
  # step3 getCPSiganlActivity
  if(T){
    inputDir <- paste0(getwd(),'/runModel/work_2023-1/')
    files <- list.files(inputDir)
    
    time_ls <- c()
    for(f in files){
      
      label <- paste(unlist(strsplit(f,'[_.]'))[3:4],collapse = '-')
      LRTG_allscore <- readRDS(paste0(getwd(),'/runModel/work_2023-1/',f))
      message(paste0('running jobs: ',label))
      
      t1 <- Sys.time()
      getSiganlImport(SiganlActivity = LRTG_allscore, Lable = label, ProjectName = '2023-1',
                      NCores = 1, AutoPara = TRUE, NTrees = 500, NTrys = 10,
                      TreeMethod = 'variance', NodeSize = 5,  NPert = 10)
      t2 <- Sys.time()
      time_ls <- c(time_ls,paste(signif(t2-t1,4),units(signif(t2-t1,4)),sep = ' '))
      
    }
  }
  
  #handle result
  files <- list.files('./getPIM/work_2023-1', full.names = TRUE, pattern = 'LRTG_pim_clean_')
  for (f in files) {
    result <- readRDS(f)
    result$pIM <- lapply(result$pIM,function(s){ifelse(s>=0,s,0)}) %>% unlist() %>% as.numeric()
    colnames(result) <- c("regulon", 'target',"value", 'type')
    saveRDS(result, file = paste0('../', gsub('.rds', '',substring(f, 41)), '_result.rds'))
  }
}


############# 
## HoloNet ##
############# 
HoloNet_function <- function(ser, fpath, sampleID, sender = NULL, receiver=NULL){
  suppressMessages(library(Seurat))
  suppressMessages(library(stringr))
  suppressMessages(library(tidyverse))
  suppressMessages(library(tidyr))
  set.seed(123)
  
  sender_ct <- unique(ser$celltype)
  if(is.null(receiver)){
    receiver_ct <- unique(ser$celltype)
  }else if(!is.null(receiver)){
    receiver_ct <- receiver
  }
  rm(receiver);gc()
  
  input_fpath <- paste0(fpath, '/input/')
  if(!dir.exists(input_fpath)){
    dir.create(input_fpath)
  }
  
  output_fpath <- paste0(fpath, '/output/')
  if (!dir.exists(output_fpath)) {
    dir.create(output_fpath)
  }
  
  ## load ICGs
  if(grepl('CID', sampleID)){
    icgs_fpath <- paste0('~/1-Datasets/NatureGenetics_2021_BC/ScriptForCCC/Giotto_result/', sampleID, '_icgs.rds')
  }else if(grepl('UKF', sampleID)){
    icgs_fpath <- paste0('~/1-Datasets/CancerCell_2022_Glioma/ScriptForLRTBenchmark/Giotto_result/', sampleID, '_icgs.rds')
  }else{
    icgs_fpath <- paste0('~/1-Datasets/Nature_2022_MI/ScriptForCCC/Giotto_result/', sampleID,'_icgs.rds')
  }
  icgs <- readRDS(icgs_fpath)
  goi <- lapply(receiver_ct, function(receiver){
    icgs.tmp <- icgs[[receiver]]
    goi.tmp <- unique(unlist(icgs.tmp))
    goi.tmp
  })
  goi <- unique(unlist(goi))
  goi_fpath <- paste0(input_fpath, 'goi.csv')
  write.csv(goi, file = goi_fpath, quote = FALSE, row.names = FALSE)
  
  ## run HoloNet
  # get gene expression matrix
  counts_fpath <- paste0(input_fpath, "raw_count.csv")
  raw.count <- as.matrix(GetAssayData(ser, "count", "Spatial"))
  write.csv(raw.count, file = counts_fpath, quote=F)
  
  # load cell type metadata
  meta_fpath <- paste0(input_fpath, "meta.csv")
  meta <- data.frame(Barcodes = rownames(ser@meta.data), 
                     celltype = ser@meta.data$celltype, 
                     row.names = rownames(ser@meta.data))
  write.csv(meta, file = meta_fpath,quote=F, row.names = FALSE)
  
  
  if(grepl('CID',sampleID)){
    
    img_fpath1 <- paste0('~/1-Datasets/NatureGenetics_2021_BC/ST_Data/spatial/', 
                        gsub('A', '', sampleID),'_spatial')
    img_fpath <- paste0(input_fpath, 'spatial')
    copy_commond <- paste0('cp -r ', img_fpath1, ' ', img_fpath)
    system(copy_commond)
    h5_file <- './Step0_SharedInfo/HoloNet/filtered_feature_bc_matrix.h5'
    copy_commond <- paste0('cp -r ', h5_file, ' ', input_fpath)
    system(copy_commond)
    command <- paste('sh ../Script/Step6_LRTToolsFunction/HoloNet_shell.sh', 
                     counts_fpath, meta_fpath, input_fpath, 
                     goi_fpath, output_fpath, sep = ' ')
    
  }else if(grepl('UKF', sampleID)){
    
    visumn_fpath <- paste0('/home/ljx/1-Datasets/CancerCell_2022_Glioma/10XVisium/',sampleID, '/outs/')
    command <- paste('sh ../Script/Step6_LRTToolsFunction/HoloNet_shell.sh', 
                     counts_fpath, meta_fpath, visumn_fpath, 
                     goi_fpath, output_fpath, sep = ' ')
  
  }else{ ## MI dataset
    
    visumn_fpath <- paste0('/home/ljx/1-Datasets/Nature_2022_MI/ST_Data/',sampleID)
    command <- paste('sh ../Script/Step6_LRTToolsFunction/HoloNet_shell.sh', 
                     counts_fpath, meta_fpath, visumn_fpath, 
                     goi_fpath, output_fpath, sep = ' ')
    
  }
  
  
  system(command)
  
  ## handle result
  setwd(fpath)
  if(T){
    files <- list.files("./output/")
    result <- lapply(files, function(file){
      lr <- str_split(file, "_")[[1]][1]
      target <- str_split(file, "_")[[1]][2] %>% str_replace(".csv", "")
      
      res <- read.csv(paste0("./output/", file)) %>%
        pivot_longer(cols = -X, names_to = "receiver", values_to = "value")
      colnames(res)[1] <- "sender"
      
      res <- res[which(res$receiver %in% receiver_ct), ]
      res$lr <- rep(lr, n = nrow(res))
      res <- separate(res, lr, c("ligand", "receptor"), ":")
      res$target <- rep(target, n = nrow(res))
      res
    })
    result <- do.call(rbind, result)
    result <- result[which(result$sender!=result$receiver), ]
    
    for (receiver in unique(result$receiver)) {
      icgs_tmp <- unique(unlist(icgs[[receiver]]))
      res <- result[which(result$receiver == receiver),]
      res <- res[which(res$target %in% icgs_tmp), ]
      res$LRT <- paste(res$ligand, res$receptor, res$target, sep = "_")
      res <- res[,c(3,7)]
      res <- aggregate(value ~ LRT, data = res, sum)
      res <- tidyr::separate(res, LRT, c('ligand', 'receptor', 'target'), "_")
      res_ligand <- res[, c(1,3,4)]
      res_receptor <- res[, c(2,3,4)]
      
      res_ligand <- aggregate(value ~ ligand+target, data = res_ligand, sum)
      res_receptor <- aggregate(value ~ receptor+target, data = res_receptor, sum)
      
      colnames(res_ligand) <- c("regulon", 'target',"value")
      colnames(res_receptor) <- c("regulon", 'target',"value")
      res_ligand$type <- 'ligand'
      res_receptor$type <- 'receptor'
      
      res <- rbind(res_ligand, res_receptor)
      saveRDS(res, file = paste0('./', receiver, '_result.rds'))
    }
  }
}


