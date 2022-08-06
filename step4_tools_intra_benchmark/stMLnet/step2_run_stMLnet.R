rm(list = ls())

library(Seurat)
library(stMLnet)
library(tidyverse)
library(doSNOW)
set.seed(123)

# load input data
if(T){
  load("./step4_tools_intra_benchmark/stMLnet/input/CID4465_seruat.RData")
  norm.matrix <- GetAssayData(se.st, "data", "SCT")
  norm.matrix <- as.matrix(norm.matrix)
  
  anno <- data.frame(Barcode = colnames(se.st), Cluster = se.st$celltype)
  
  loc <- se.st@meta.data
  loc <- loc[,7:8]
  
  ex_databases <- readRDS("./step4_tools_intra_benchmark/MistyR/prior_knowledge/Databases.rds")
  
  Ligs_up_list <- readRDS("step4_tools_intra_benchmark/stMLnet/input/CID4465_Ligs_up.rds")
  Recs_expr_list <- readRDS("./step4_tools_intra_benchmark/stMLnet/input/CID4465_Recs_expr.rds")
  ICG_list <- readRDS("./step1_data_process/result/Giotto_result/CID4465_bc_ICGs.rds")
}

setwd("./step4_tools_intra_benchmark/stMLnet/result/")

# step1 create mulityayer network
resMLnet <- runMLnet(ExprMat = norm.matrix, AnnoMat = anno, Normalize = F, 
                     logfc.ct = 0.1, pct.ct = 0.05, expr.ct = 0.1,
                     ProjectName = '0707-1', Databases = NULL,
                     TGList=ICG_list, LigList=Ligs_up_list, RecList=Recs_expr_list)
rm(ICG_list, Ligs_up_list, Recs_expr_list, ex_databases, se.st)
#str(resMLnet$mlnets$`Cancer Epithelial`, max.level = 2)

# step2 calculate Signal Activity
if(T){
  ## calculate distant
  DistMat <- as.matrix(dist(loc))
  colnames(DistMat) <- rownames(loc)
  rownames(DistMat) <- rownames(loc)
  
  ## imputation
  exprMat.Impute <- runImputation (exprMat = norm.matrix)
  rm(norm.matrix)
  
  # remove empty list
  ex_mulnetlist <- list()
  for (i in 1:length(resMLnet$mlnets)) {
    
    mlnets <- resMLnet$mlnets[[i]]
    for (j in 1:length(mlnets)) {
      mlnet <- mlnets[[j]]
      if(nrow(mlnet$LigRec)!=0) ex_mulnetlist[[names(mlnets)[j]]] = mlnet
    }
    
  }
  rm(i, j, mlnets, mlnet,loc)
  
  clusters <- anno$Cluster %>% unique() %>% as.character()
  
  Sender <- NULL
  resSigActList_tme <- getSiganlActivity(ExprMat = exprMat.Impute,
                                                    DistMat = DistMat,
                                                    AnnoMat = anno,
                                                    MulNetList = ex_mulnetlist,
                                                    Receiver = 'Cancer Epithelial', Sender = Sender,
                                                    Downsample = FALSE, ProjectName = '0707-1')
}

# step3 getCPSiganlActivity
if(T){
  inputDir <- paste0(getwd(),'/runModel/work_0707-1/')
  files <- list.files(inputDir)
  
  time_ls <- c()
  for(f in files){
    
    label <- paste(unlist(strsplit(f,'[_.]'))[3:4],collapse = '-')
    LRTG_allscore <- readRDS(paste0(getwd(),'/runModel/work_0707-1/',f))
    message(paste0('running jobs: ',label))
    
    t1 <- Sys.time()
    getSiganlImport(SiganlActivity = LRTG_allscore, Lable = label, ProjectName = '0707-1',
                    NCores = 50, AutoPara = TRUE, NTrees = 500, NTrys = 10,
                    TreeMethod = 'variance', NodeSize = 5,  NPert = 10)
    t2 <- Sys.time()
    time_ls <- c(time_ls,paste(signif(t2-t1,4),units(signif(t2-t1,4)),sep = ' '))
    
  }
}


