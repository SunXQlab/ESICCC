###################
## calculate_auc ##
###################

## get different metrices and ROC/PRC preformance object
get_evaluate_metrics <- function(pred,label)
{
  
  find_optimal_cutoff <- function(TPR, FPR, threshold){
    
    y = TPR - FPR
    Youden_index = which.max(y)
    # optimal_threshold = threshold[Youden_index]
    return(Youden_index)
    
  }
  
  if(length(which(label==TRUE))!=0 & length(which(label==FALSE))!=0){
    pred <- prediction(pred, label)
    
    perf_ROC <- performance(pred, measure = "tpr", x.measure = "fpr")
    ind_ROC <- find_optimal_cutoff(perf_ROC@y.values[[1]],perf_ROC@x.values[[1]],perf_ROC@alpha.values[[1]])
    cutoff_ROC <- perf_ROC@alpha.values[[1]][ind_ROC]
    
    perf_PRC <- performance(pred, measure = "prec", x.measure = "rec")
    ind_PRC <- find_optimal_cutoff(perf_PRC@y.values[[1]],perf_PRC@x.values[[1]],perf_PRC@alpha.values[[1]])
    cutoff_PRC <- perf_PRC@alpha.values[[1]][ind_PRC]
    
    ACC <- performance(pred, measure = "acc")@y.values[[1]] %>% .[ind_ROC] %>% signif(.,4)
    ERR <- performance(pred, measure = "err")@y.values[[1]] %>% .[ind_ROC]  %>% signif(.,4)
    PPV <- performance(pred, measure = "ppv")@y.values[[1]] %>% .[ind_ROC]  %>% signif(.,4)
    MCC <- performance(pred, measure = "mat")@y.values[[1]] %>% .[ind_ROC]  %>% signif(.,4)
    
    AUC <- performance(pred, measure = "auc")@y.values[[1]] %>% signif(.,4)
    AUCPR <- performance(pred, measure = "aucpr")@y.values[[1]] %>% signif(.,4)
    
    res = list(perf_ROC = perf_ROC,
               perf_PRC = perf_PRC,
               perf_metrics = c(ROC_AUC=AUC,PRC_AUC=AUCPR,
                                ACC=ACC,ERR=ERR,PPV=PPV,MCC=MCC))
    
  }else{
    
    res <- list(perf_ROC = NA, perf_PRC = NA, 
                perf_metrics = rep(0,6))
    names(res$perf_metrics) <- c('ROC_AUC','PRC_AUC','ACC','ERR','PPV','MCC')
    
  }
  
  
  return(res)
  
}