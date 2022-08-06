rm(list = ls())

#############
## library ##
#############

library(dplyr)
library(ROCR)
library(ggplot2)
library(ggsci)
library(tidyverse)
set.seed(123)
gc()

setwd('./step4_tools_intra_benchmark/')

source('./code.R')

###########
## color ##
###########

scales::show_col(pal_aaas(palette = "default", alpha = 0.6)(7))
mycolors_aaas <- pal_aaas(palette = "default", alpha = 0.6)(7)

mycolor_software <- mycolors_aaas[c(1:5)]
names(mycolor_software) <- c("CytoTalk", "NicheNet","stMLnet","MISTy", 'HoloNet')
scales::show_col(mycolor_software)

###################
## load-software ##
###################

## CytoTalk
cytotalk_score <- readRDS("./CytoTalk/result/LRT_distance.rds")
cytotalk_score$from <- as.vector(cytotalk_score$from)
head(cytotalk_score)

## NicheNet
nichenet_score <- readRDS("./NicheNet/result/tumor_LRpair_weight.rds")
nichenet_score$ligand <- as.vector(nichenet_score$ligand)
head(nichenet_score)

## stMLnet
stmlnet_score <- readRDS("./stMLnet/result/getPIM/work_0707-1/LRTG_pim_clean_TME-Cancer Epithelial.rds")
stmlnet_score$pIM <- lapply(stmlnet_score$pIM,function(s){ifelse(s>=0,s,0)}) %>% unlist() %>% as.numeric()
head(stmlnet_score)

## MISTy
misty_score <- readRDS("./MistyR/result/misty_score.rds")
misty_score <- as.data.frame(misty_score)
head(misty_score)

## HoloNet
load('./HoloNet/result_holonet.RData')
holonet_score <- result_holonet; rm(result_holonet)
head(holonet_score)

## groundtrue
load('./cellLine/BC_celline_used/degs_ls.RData')

#######################
## evaluate-software ##
#######################

# split

ress_ls <- lapply(1:length(degs_ls), function(i){
  
  print(i)
  deg_ls <- degs_ls[[i]]
  key <- strsplit(names(degs_ls)[i],'_')[[1]][2]
  type <- ifelse(grepl('KO',names(degs_ls)[i]),'Receptor','Ligand')
  
  score_holonet <- holonet_score %>% filter(regulon == key & type == type) 
  gene_holonet <- intersect(score_holonet$target,names(deg_ls))
  label_holonet <- deg_ls[gene_holonet]
  pred_holonet <- score_holonet$value[match(gene_holonet,score_holonet$target)]
  res_holonet <- get_evaluate_metrics(pred_holonet,label_holonet)
  
  score_stml <- stmlnet_score %>% filter(regulator == key & type == type) 
  gene_stml <- intersect(score_stml$Target,names(deg_ls))
  label_stml <- deg_ls[gene_stml]
  pred_stml <- score_stml$pIM[match(gene_stml,score_stml$Target)]
  res_stml <- get_evaluate_metrics(pred_stml,label_stml)
  
  if(type == 'Ligand'){
    score_nich <- nichenet_score[which(nichenet_score$ligand==key),]
  }else if(type == 'Receptor' & key == 'AXL'){
    score_nich <- nichenet_score[which(nichenet_score$ligand=='GAS6'),]
  }else if(type == 'Receptor' & key == 'CXCR4'){
    score_nich <- nichenet_score[which(nichenet_score$ligand=='CXCL12'),]
  }else if(type == 'Receptor' & key == 'NRP1'){
    score_nich <- nichenet_score[which(nichenet_score$ligand=='VEGFA'),]
  }
  gene_nich <- intersect(score_nich$target,names(deg_ls))
  label_nich <- deg_ls[gene_nich]
  pred_nich <- score_nich$weight[match(gene_nich,score_nich$target)]
  res_nich <- get_evaluate_metrics(pred_nich,label_nich)

  score_cyto <- cytotalk_score %>% filter(from == key & from_type == type)
  gene_cyto <- intersect(score_cyto$to,names(deg_ls))
  label_cyto <- deg_ls[gene_cyto]
  pred_cyto <- score_cyto$distance[match(gene_cyto,score_cyto$to)]
  res_cyto <- get_evaluate_metrics(pred_cyto,label_cyto)
  
  score_mist <- misty_score %>% filter(regulator == key & type == type)
  gene_mist <- intersect(score_mist$target,names(deg_ls))
  label_mist <- deg_ls[gene_mist]
  pred_mist <- score_mist$score[match(gene_mist,score_mist$target)]
  res_mist <- get_evaluate_metrics(pred_mist,label_mist)
  
  result <- list(res_cyto = res_cyto, res_nich=res_nich, 
                 res_stml=res_stml, res_mist=res_mist,res_holo = res_holonet )
  
  result
  
})

names(ress_ls) <- names(degs_ls)

# check

str(ress_ls,max.level = 1)

#####################
## output-software ##
#####################

## data 

evaluate_ls <- ress_ls
str(evaluate_ls,max.level = 1)
df_plot <- lapply(1:length(evaluate_ls), function(i){
  
  print(names(evaluate_ls)[i])
  res <- evaluate_ls[[i]]
  
  df <- do.call('rbind',list(res$res_cyto$perf_metrics,
                             res$res_nich$perf_metrics,
                             res$res_stml$perf_metrics,
                             res$res_mist$perf_metrics,
                             res$res_holo$perf_metrics)) %>% as.data.frame()
  df$method <- c('CytoTalk','NicheNet','stMLnet', 'MISTy', 'HoloNet')
  df$group <- names(evaluate_ls)[i]
  df
  
}) %>% do.call('rbind',.)

## barplot

df_plot$Dataset <- strsplit(df_plot$group,'_') %>% do.call('rbind',.) %>% .[,2]

rename_dataset <- function(x){
  switch(EXPR = x,
         "GSE120268_AXL_KO" = 'AXL', #'GSE120268',
         "GSE157680_NRP1_MDA-MB-231_KO" = 'NRP1', #'GSE157680',        
         "GSE15893_CXCR4_KO" = 'CXCR4',#'GSE15893',
         "GSE15893_CXCL12_Treatment" = 'CXCL12',
         "GSE160990_TGFB1_Treatment" = 'TGFB1',#'GSE160990',
         "GSE36051_DLL4_MCF7_Treatment" = 'DLL4(1)', #'GSE36051(1)',
         "GSE36051_DLL4_MDA-MB-231_Treatment" = 'DLL4(2)',#'GSE36051(2)',
         "GSE36051_JAG1_MDA-MB-231_Treatment" = 'JAG1',#'GSE36051(3)',
         "GSE65398_IGF1_Treatment" = 'IGF1(1)',#'GSE65398',       
         "GSE7561_IGF1_Treatment" = 'IGF1(2)' #'GSE7561', 
  )
}
df_plot$Dataset <- lapply(df_plot$group,rename_dataset) %>% unlist()
df_plot$ROC_AUC_label <- round(df_plot$ROC_AUC, digits = 2)
df_plot$PRC_AUC_label <- round(df_plot$PRC_AUC, digits = 2)

df_plot$Dataset <- factor(df_plot$Dataset,
                          levels = c('AXL', 'NRP1', 'CXCR4', 'CXCL12','TGFB1',
                                     'DLL4(1)', 'DLL4(2)', 'JAG1',
                                     'IGF1(1)', 'IGF1(2)'))

df_plot <- lapply(unique(df_plot$Dataset), function(ds){
  
  df <- df_plot[df_plot$Dataset == ds,]
  if(sum(df$ROC_AUC)==0){
    df = NULL
  }else{
    df = df
  }
  
})
delet_celllines <- which(lapply(df_plot, is.null) %>% unlist())
if(length(delet_celllines)>0) df_plot <- df_plot[-delet_celllines]
df_plot <- df_plot %>% do.call('rbind',.)
df_plot$method <-factor(df_plot$method,ordered=TRUE,levels=c('CytoTalk','NicheNet','stMLnet', 'MISTy', 'HoloNet'))

theme_bar <- function(..., bg='white'){
  require(grid)
  theme_classic(...) +
    theme(rect=element_rect(fill=bg),
          plot.margin=unit(rep(0.5,4), 'lines'),
          panel.background=element_rect(fill='transparent', color='black'),
          panel.border=element_rect(fill='transparent', color='transparent'),
          panel.grid=element_blank(),
          axis.title.x = element_blank(),
          axis.title.y=element_text(face = "bold",size = 14),
          axis.text = element_text(face = "bold",size = 14),
          legend.title=element_blank(),
          legend.position='top',
          legend.direction = "horizontal",
          legend.text = element_text(face = "bold",size = 12,margin = margin(r=20)),
          legend.background = element_rect( linetype="solid",colour ="black")
    )
  
}

p1 <- ggplot(df_plot, aes(x=Dataset,y=ROC_AUC,fill=method)) + 
  geom_bar(stat="identity",position=position_dodge(0.75),width=0.6) +
  scale_y_continuous(expand = c(0, 0.01)) + labs(y = 'AUCROC') +
  coord_cartesian(ylim=c(-0.02,1.1))+
  scale_fill_manual(values = mycolor_software, breaks=c('CytoTalk','NicheNet','stMLnet','MISTy', 'HoloNet')) + 
  theme_bar() + #coord_flip() +
  geom_text(aes(label=ROC_AUC_label),position=position_dodge(0.9),vjust = -0.5,size = 4) 
p1

p2 <- ggplot(df_plot, aes(x=Dataset,y=PRC_AUC,fill=method)) + 
  geom_bar(stat="identity",position=position_dodge(0.75),width=0.6) +
  scale_y_continuous(expand = c(0, 0.01)) + labs(y = 'AUCPRC') +
  coord_cartesian(ylim=c(-0.02,1.1))+
  scale_fill_manual(values = mycolor_software, breaks=c('CytoTalk','NicheNet','MISTy','stMLnet', 'HoloNet')) + 
  theme_bar() + #coord_flip() +
  geom_text(aes(label=PRC_AUC_label),position=position_dodge(0.9),vjust = -0.5,size = 4) 
p2

# violin-plot
if(T){
  try <- df_plot
  try$PRC_AUC[c(1,6, 11,16)] <- NA
  
  try <- try[-which(try$PRC_AUC==0),]
  ggpubr::ggviolin(try, "method", "PRC_AUC",  color = "method", palette  = mycolor_software,
                   add =  c( "jitter","boxplot"))
  
}
