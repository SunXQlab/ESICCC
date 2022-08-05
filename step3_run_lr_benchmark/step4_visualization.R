rm(list = ls())
library(ggplot2)
set.seed(123)
setwd('./step3_run_lr_benchmark/')
load('./result.bench/result.benchmark.1.RData')
# Visualization —— CellChat/iTALK/scSeqComm/CellTalker —— only perc=0.1
if(F){
  # draw mutual information
  if(F){
    data <- result.benchmark[grepl("^CellChat",result.benchmark$tool), ]
    p1 <- ggplot(data, aes(tool, mutinfo, fill=factor(groups))) + 
      geom_boxplot(outlier.shape = NA, lwd=0.5)+
      theme_bw()+
      theme(plot.title = element_text(hjust = 0.5)) +
      ggpubr::stat_compare_means(aes(group=groups), label = "p.signif",
                                 method.args = list(alternative = "less"),
                                 label.y = 0.15)+
      coord_cartesian(ylim=c(0,0.15))+
      ylab("Mutual information")+
      ggtitle("CellChat")+
      scale_x_discrete("methods", labels = c("CellChat_trim" = "triM",
                                             "CellChat_trun_05" = "trunM_05",
                                             "CellChat_trun_10" = "trunM_10",
                                             "CellChat_trun_15" = "trunM_15",
                                             "CellChat_trun_20" = "trunM_20"))
    p1
    #comparisons = my_comparisons
    data <- result.benchmark[grepl("^iTALK",result.benchmark$tool), ]
    p2 <- ggplot(data, aes(tool, mutinfo, fill=factor(groups))) + 
      geom_boxplot(outlier.shape = NA, lwd = 0.5)+
      theme_bw()+
      theme(plot.title = element_text(hjust = 0.5)) +
      ggpubr::stat_compare_means(aes(group=groups), label = "p.signif", 
                                 method.args = list(alternative = "less"),
                                 label.y = 0.15)+
      coord_cartesian(ylim=c(0,0.15))+
      ylab("Mutual information")+
      ggtitle("iTALK")+
      scale_x_discrete("top_gene", labels = c("iTALK_10" = "10",
                                              "iTALK_25" = "25",
                                              "iTALK_50" = "50"))
    p2
    
    data <- result.benchmark[grepl("^scSeqComm",result.benchmark$tool), ]
    p3 <- ggplot(data, aes(tool, mutinfo, fill=factor(groups))) + 
      geom_boxplot(outlier.shape = NA, lwd = 0.5)+
      theme_bw()+
      theme(plot.title = element_text(hjust = 0.5)) +
      ggpubr::stat_compare_means(aes(group=groups), label = "p.signif", 
                                 method.args = list(alternative = "less"),
                                 label.y = 0.15)+
      coord_cartesian(ylim=c(0,0.15))+
      ylab("Mutual information")+
      ggtitle("scSeqComm")+
      scale_x_discrete("S_intra/S_inter cutoff", labels = c("scSeqComm_0.5" = "0.5",
                                                            "scSeqComm_0.6" = "0.6",
                                                            "scSeqComm_0.7" = "0.7",
                                                            "scSeqComm_0.8" = "0.8",
                                                            "scSeqComm_0.9" = "0.9"))
    p3
    data <- result.benchmark[grepl("^CellTalker",result.benchmark$tool), ]
    data$tool <- factor(data$tool,levels = c("CellTalker_0", "CellTalker_50", "CellTalker_100", "CellTalker_150",
                                             "CellTalker_200", "CellTalker_250", "CellTalker_500"))
    p4 <- ggplot(data, aes(tool, mutinfo, fill=factor(groups))) + 
      geom_boxplot(outlier.shape = NA, lwd = 0.5)+
      theme_bw()+
      theme(plot.title = element_text(hjust = 0.5)) +
      ggpubr::stat_compare_means(aes(group=groups), label = "p.signif", 
                                 method.args = list(alternative = "less"),
                                 label.y = 0.15)+
      coord_cartesian(ylim=c(0,0.15))+
      ylab("Mutual information")+
      ggtitle("CellTalker")+
      scale_x_discrete("min_expression", labels = c("CellTalker_0" = "0",
                                                    "CellTalker_50" = "50",
                                                    "CellTalker_100" = "100",
                                                    "CellTalker_150" = "150",
                                                    "CellTalker_200" = "200",
                                                    "CellTalker_250" = "250",
                                                    "CellTalker_500" = "500"))
    
    p4
    Rmisc::multiplot(p1, p3, p2, p4, cols = 2)
  }
  
  # draw pcc
  if(F){
    data <- result.benchmark[grepl("^CellChat",result.benchmark$tool), ]
    data$pcc <- abs(data$pcc)
    p1 <- ggplot(data, aes(tool, pcc, fill=factor(groups))) + 
      geom_boxplot(outlier.shape = NA, lwd = 0.5)+
      theme_bw()+
      theme(plot.title = element_text(hjust = 0.5)) +
      ggpubr::stat_compare_means(aes(group=groups), label = "p.signif", method.args = list(alternative = "less"),
                                 label.y = 0.3)+
      coord_cartesian(ylim=c(0,0.3))+
      ylab("Pearson correlation coefficient")+
      ggtitle("CellChat")+
      scale_x_discrete("methods", labels = c("CellChat_trim" = "triM",
                                             "CellChat_trun_05" = "trunM_05",
                                             "CellChat_trun_10" = "trunM_10",
                                             "CellChat_trun_15" = "trunM_15",
                                             "CellChat_trun_20" = "trunM_20"))
    p1
    
    data <- result.benchmark[grepl("^iTALK",result.benchmark$tool), ]
    data$pcc <- abs(data$pcc)
    p2 <- ggplot(data, aes(tool, pcc, fill=factor(groups))) + 
      geom_boxplot(outlier.shape = NA, lwd = 0.5)+
      theme_bw()+
      theme(plot.title = element_text(hjust = 0.5)) +
      ggpubr::stat_compare_means(aes(group=groups), label = "p.signif", method.args = list(alternative = "less"),
                                 label.y = 0.3)+
      coord_cartesian(ylim=c(0,0.3))+
      ylab("Pearson correlation coefficient")+
      ggtitle("iTALK")+
      scale_x_discrete("top_gene", labels = c("iTALK_10" = "10",
                                              "iTALK_25" = "25",
                                              "iTALK_50" = "50"))
    p2
    
    data <- result.benchmark[grepl("^scSeqComm",result.benchmark$tool), ]
    data$pcc <- abs(data$pcc)
    p3 <- ggplot(data, aes(tool, pcc, fill=factor(groups))) + 
      geom_boxplot(outlier.shape = NA, lwd = 0.5)+
      theme_bw()+
      theme(plot.title = element_text(hjust = 0.5)) +
      ggpubr::stat_compare_means(aes(group=groups), label = "p.signif", method.args = list(alternative = "less"),
                                 label.y = 0.3)+
      coord_cartesian(ylim=c(0,0.3))+
      ylab("Pearson correlation coefficient")+
      ggtitle("scSeqComm")+
      scale_x_discrete("S_intra/S_inter cutoff", labels = c("scSeqComm_0.5" = "0.5",
                                                            "scSeqComm_0.6" = "0.6",
                                                            "scSeqComm_0.7" = "0.7",
                                                            "scSeqComm_0.8" = "0.8",
                                                            "scSeqComm_0.9" = "0.9"))
    p3
    
    data <- result.benchmark[grepl("^CellTalker",result.benchmark$tool), ]
    data$pcc <- abs(data$pcc)
    data$tool <- factor(data$tool,levels = c("CellTalker_0", "CellTalker_50", "CellTalker_100", "CellTalker_150",
                                             "CellTalker_200", "CellTalker_250", "CellTalker_500"))
    p4 <- ggplot(data, aes(tool, pcc, fill=factor(groups))) + 
      geom_boxplot(outlier.shape = NA, lwd = 0.5)+
      theme_bw()+
      theme(plot.title = element_text(hjust = 0.5)) +
      ggpubr::stat_compare_means(aes(group=groups), label = "p.signif", method.args = list(alternative = "less"),
                                 label.y = 0.4)+
      coord_cartesian(ylim=c(0,0.4))+
      ylab("Pearson correlation coefficient")+
      ggtitle("CellTalker")+
      scale_x_discrete("min_expression", labels = c("CellTalker_0" = "0",
                                                    "CellTalker_50" = "50",
                                                    "CellTalker_100" = "100",
                                                    "CellTalker_150" = "150",
                                                    "CellTalker_200" = "200",
                                                    "CellTalker_250" = "250",
                                                    "CellTalker_500" = "500"))
    p4
    Rmisc::multiplot(p1, p3, p2, p4, cols = 2)
  }
}

#handle the result
if(T){
  temp <- result.benchmark[!grepl("^CellChat.trun", result.benchmark$tool),]
  temp <- temp[-which(temp$tool %in% c("iTALK_25", "iTALK_50")),]
  temp <- temp[-which(temp$tool %in% c("scSeqComm_0.7", "scSeqComm_0.6", "scSeqComm_0.9", "scSeqComm_0.5")),]
  temp <- temp[-which(temp$tool %in% unique(temp$tool)[7:12]),]
  temp$tool <- stringr::str_replace_all(temp$tool, "scSeqComm_0.8", "scSeqComm")
  temp$tool <- stringr::str_replace_all(temp$tool, "CellTalker_0", "CellTalker")
  temp$tool <- stringr::str_replace_all(temp$tool, "CellChat_trim", "CellChat")
  temp$tool <- stringr::str_replace_all(temp$tool, "iTALK_10", "iTALK")
}
rm(result.benchmark, data, p1,p2,p3,p4)

# Visualization —— all tools —— 0.1
if(T){
  
  data <- temp
  
  library(ggplot2)
  data$tool <- factor(data$tool,levels = c("CellPhoneDB2", "CellPhoneDB3", "CellTalker", # expression value; mean
                                           "Connectome", "NATMI", "ICELLNET", # expression value; product
                                           "scConnect", "Cellinker", "CellChat", # expression value; others
                                           "SingleCellSignalR", "CytoTalk", "CellCall", "scSeqComm", "NicheNet", # expression value; pathway-based
                                           "Domino","scMLnet", "PyMINEr", "iTALK"
  ))
  
  ggplot(data, aes(tool, mutinfo, fill=factor(groups))) + 
    geom_boxplot(outlier.shape = NA, lwd = 0.5)+
    theme_bw()+
    theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1), plot.title = element_text(hjust = 0.5)) +
    ggpubr::stat_compare_means(aes(group=groups), label = "p.signif", method.args = list(alternative = "less"),
                               label.y = 0.15)+
    coord_cartesian(ylim=c(0,0.15))+
    ylab("Mutual information")+
    ggtitle("top 10% close/distant cell pairs")+
    xlab("")
  
  data$pcc <- abs(data$pcc)
  ggplot(data, aes(tool, pcc, fill=factor(groups))) + 
    geom_boxplot(outlier.shape = NA, lwd = 0.5)+
    theme_bw()+
    theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1), plot.title = element_text(hjust = 0.5)) +
    ggpubr::stat_compare_means(aes(group=groups), label = "p.signif", method.args = list(alternative = "less"),
                               label.y = 0.4)+
    coord_cartesian(ylim=c(0,0.4))+
    ylab("Pearson correlation coefficient")+
    ggtitle("top 10% close/distant cell pairs")+
    xlab("")
  
}

# Visualization —— all tools —— 0.2
if(T){
  rm(list = ls())
  load('./result.bench/result.benchmark.2.RData')
  
  if(T){
    temp <- result.benchmark[!grepl("^CellChat.trun", result.benchmark$tool),]
    temp <- temp[-which(temp$tool %in% c("iTALK_25", "iTALK_50")),]
    temp <- temp[-which(temp$tool %in% c("scSeqComm_0.7", "scSeqComm_0.6", "scSeqComm_0.9", "scSeqComm_0.5")),]
    temp <- temp[-which(temp$tool %in% unique(temp$tool)[7:12]),]
    temp$tool <- stringr::str_replace_all(temp$tool, "scSeqComm_0.8", "scSeqComm")
    temp$tool <- stringr::str_replace_all(temp$tool, "CellTalker_0", "CellTalker")
    temp$tool <- stringr::str_replace_all(temp$tool, "CellChat_trim", "CellChat")
    temp$tool <- stringr::str_replace_all(temp$tool, "iTALK_10", "iTALK")
  }
  rm(result.benchmark)
  
  
  data <- temp
  
  library(ggplot2)
  data$tool <- factor(data$tool,levels = c("CellPhoneDB2", "CellPhoneDB3", "CellTalker", # expression value; mean
                                           "Connectome", "NATMI", "ICELLNET", # expression value; product
                                           "scConnect", "Cellinker", "CellChat", # expression value; others
                                           "SingleCellSignalR", "CytoTalk", "CellCall", "scSeqComm", "NicheNet", # expression value; pathway-based
                                           "Domino","scMLnet", "PyMINEr", "iTALK"
  ))
  
  ggplot(data, aes(tool, mutinfo, fill=factor(groups))) + 
    geom_boxplot(outlier.shape = NA, lwd = 0.5)+
    theme_bw()+
    theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1), plot.title = element_text(hjust = 0.5)) +
    ggpubr::stat_compare_means(aes(group=groups), label = "p.signif", method.args = list(alternative = "less"),
                               label.y = 0.1)+
    coord_cartesian(ylim=c(0,0.1))+
    ylab("Mutual information")+
    ggtitle("top 20% close/distant cell pairs")+
    xlab("")
  
  data$pcc <- abs(data$pcc)
  ggplot(data, aes(tool, pcc, fill=factor(groups))) + 
    geom_boxplot(outlier.shape = NA, lwd = 0.5)+
    theme_bw()+
    theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1), plot.title = element_text(hjust = 0.5)) +
    ggpubr::stat_compare_means(aes(group=groups), label = "p.signif", method.args = list(alternative = "less"),
                               label.y = 0.3)+
    coord_cartesian(ylim=c(0,0.3))+
    ylab("Pearson correlation coefficient")+
    ggtitle("top 20% close/distant cell pairs")+
    xlab("")
  
}

# Visualization —— all tools —— 0.3
if(T){
  rm(list = ls())
  load('./result.bench/result.benchmark.3.RData')
  
  if(T){
    temp <- result.benchmark[!grepl("^CellChat.trun", result.benchmark$tool),]
    temp <- temp[-which(temp$tool %in% c("iTALK_25", "iTALK_50")),]
    temp <- temp[-which(temp$tool %in% c("scSeqComm_0.7", "scSeqComm_0.6", "scSeqComm_0.9", "scSeqComm_0.5")),]
    temp <- temp[-which(temp$tool %in% unique(temp$tool)[7:12]),]
    temp$tool <- stringr::str_replace_all(temp$tool, "scSeqComm_0.8", "scSeqComm")
    temp$tool <- stringr::str_replace_all(temp$tool, "CellTalker_0", "CellTalker")
    temp$tool <- stringr::str_replace_all(temp$tool, "CellChat_trim", "CellChat")
    temp$tool <- stringr::str_replace_all(temp$tool, "iTALK_10", "iTALK")
  }
  rm(result.benchmark)
  
  
  data <- temp
  
  library(ggplot2)
  data$tool <- factor(data$tool,levels = c("CellPhoneDB2", "CellPhoneDB3", "CellTalker", # expression value; mean
                                           "Connectome", "NATMI", "ICELLNET", # expression value; product
                                           "scConnect", "Cellinker", "CellChat", # expression value; others
                                           "SingleCellSignalR", "CytoTalk", "CellCall", "scSeqComm", "NicheNet", # expression value; pathway-based
                                           "Domino","scMLnet", "PyMINEr", "iTALK"
  ))
  
  ggplot(data, aes(tool, mutinfo, fill=factor(groups))) + 
    geom_boxplot(outlier.shape = NA, lwd = 0.5)+
    theme_bw()+
    theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1), plot.title = element_text(hjust = 0.5)) +
    ggpubr::stat_compare_means(aes(group=groups), label = "p.signif", method.args = list(alternative = "less"),
                               label.y = 0.1)+
    coord_cartesian(ylim=c(0,0.1))+
    ylab("Mutual information")+
    ggtitle("top 30% close/distant cell pairs")+
    xlab("")
  
  data$pcc <- abs(data$pcc)
  ggplot(data, aes(tool, pcc, fill=factor(groups))) + 
    geom_boxplot(outlier.shape = NA, lwd = 0.5)+
    theme_bw()+
    theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1), plot.title = element_text(hjust = 0.5)) +
    ggpubr::stat_compare_means(aes(group=groups), label = "p.signif", method.args = list(alternative = "less"),
                               label.y = 0.3)+
    coord_cartesian(ylim=c(0,0.3))+
    ylab("Pearson correlation coefficient")+
    ggtitle("top 30% close/distant cell pairs")+
    xlab("")
  
}

# Visualization —— all tools —— 0.4
if(T){
  rm(list = ls())
  load('./result.bench/result.benchmark.4.RData')
  
  if(T){
    temp <- result.benchmark[!grepl("^CellChat.trun", result.benchmark$tool),]
    temp <- temp[-which(temp$tool %in% c("iTALK_25", "iTALK_50")),]
    temp <- temp[-which(temp$tool %in% c("scSeqComm_0.7", "scSeqComm_0.6", "scSeqComm_0.9", "scSeqComm_0.5")),]
    temp <- temp[-which(temp$tool %in% unique(temp$tool)[7:12]),]
    temp$tool <- stringr::str_replace_all(temp$tool, "scSeqComm_0.8", "scSeqComm")
    temp$tool <- stringr::str_replace_all(temp$tool, "CellTalker_0", "CellTalker")
    temp$tool <- stringr::str_replace_all(temp$tool, "CellChat_trim", "CellChat")
    temp$tool <- stringr::str_replace_all(temp$tool, "iTALK_10", "iTALK")
  }
  rm(result.benchmark)
  
  
  data <- temp
  
  library(ggplot2)
  data$tool <- factor(data$tool,levels = c("CellPhoneDB2", "CellPhoneDB3", "CellTalker", # expression value; mean
                                           "Connectome", "NATMI", "ICELLNET", # expression value; product
                                           "scConnect", "Cellinker", "CellChat", # expression value; others
                                           "SingleCellSignalR", "CytoTalk", "CellCall", "scSeqComm", "NicheNet", # expression value; pathway-based
                                           "Domino","scMLnet", "PyMINEr", "iTALK"
  ))
  
  ggplot(data, aes(tool, mutinfo, fill=factor(groups))) + 
    geom_boxplot(outlier.shape = NA, lwd = 0.5)+
    theme_bw()+
    theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1), plot.title = element_text(hjust = 0.5)) +
    ggpubr::stat_compare_means(aes(group=groups), label = "p.signif", method.args = list(alternative = "less"),
                               label.y = 0.06)+
    coord_cartesian(ylim=c(0,0.06))+
    ylab("Mutual information")+
    ggtitle("top 40% close/distant cell pairs")+
    xlab("")
  
  data$pcc <- abs(data$pcc)
  ggplot(data, aes(tool, pcc, fill=factor(groups))) + 
    geom_boxplot(outlier.shape = NA, lwd = 0.5)+
    theme_bw()+
    theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1), plot.title = element_text(hjust = 0.5)) +
    ggpubr::stat_compare_means(aes(group=groups), label = "p.signif", method.args = list(alternative = "less"),
                               label.y = 0.2)+
    coord_cartesian(ylim=c(0,0.2))+
    ylab("Pearson correlation coefficient")+
    ggtitle("top 40% close/distant cell pairs")+
    xlab("")
  
}