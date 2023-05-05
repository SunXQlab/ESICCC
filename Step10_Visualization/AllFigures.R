rm(list = ls());gc()
suppressMessages(library(tidyverse))
suppressMessages(library(reshape2))
suppressMessages(library(ggrepel))
suppressMessages(library(RColorBrewer))
suppressMessages(library(ggplot2))
suppressMessages(library(gg.gap))
suppressMessages(library(patchwork))
suppressMessages(library(ggsci))
setwd('~/2-CCC-benchmark-202204-202208/CCC-benchmark-202212/')
set.seed(123)

#################
## Figure 3A-B ##
#################

if(T){
  # SI
  if(T){
    files <- list.files('./Data/Step4_SIRSIForLRBench', pattern = '_SIResult.rds', full.names = TRUE)
    samples <- list.files('./Data/Step4_SIRSIForLRBench', pattern = '_SIResult.rds') %>% gsub('_SIResult.rds', '', .)
    
    SI <- lapply(seq(samples), function(i){
      SI.tmp <- readRDS(files[i])
      SI.tmp[lower.tri(as.matrix(SI.tmp))] <- NA
      SI.tmp <- SI.tmp[-which(colnames(SI.tmp) %in% c('cellinker', 'cell2cell')),
                       -which(rownames(SI.tmp) %in% c('cellinker', 'cell2cell'))]
      SI.tmp <- tibble::rownames_to_column(SI.tmp, 'method1')
      SI.tmp <- tidyr::pivot_longer(SI.tmp, -method1, names_to = 'method2')
      SI.tmp <- SI.tmp[-which(is.na(SI.tmp$value)), ]
      SI.tmp$dataset <- samples[i]
      SI.tmp
    })
    
    SI <- do.call(rbind, SI)
    SI <- aggregate(value ~ method1+method2, data = SI, mean)
    SI.tmp <- SI
    colnames(SI.tmp)[1:2] <- c('method2', 'method1')
    SI.tmp <- SI.tmp[,c(2,1,3)]
    SI <- rbind(SI,SI.tmp)
    SI <- aggregate(value ~ method1+method2, data = SI, mean)
    SI <- tidyr::spread(SI, key = method1, value = value, fill = NA)
    SI <- tibble::column_to_rownames(SI, 'method2')
    SI <- as.matrix(SI)
    
    order.names <- c("CellPhoneDB2", "CellPhoneDB3", "CellTalker", 
                     "Connectome", "ICELLNET",  "NATMI", "iTALK", 
                     "scConnect", "SingleCellSignalR", "CellChat", 
                     "RNAMagnet", "PyMINEr",
                     "scSeqComm", "NicheNet","CytoTalk",
                     "scMLnet", "Domino", "CellCall")
    
    SI <- SI[order.names, order.names]
    
    bk <- c(seq(0,0.4,by=0.01),seq(0.5,1,by=0.01))
    p1 <- pheatmap::pheatmap(SI,cluster_rows = FALSE, cluster_cols = FALSE,
                             treeheight_row = 0,
                             treeheight_col = 0,
                             display_numbers = TRUE,
                             fontsize=13,
                             #annotation_col = anno,
                             number_format = "%.2f",
                             border_color=NA,
                             angle_col = "45",
                             color = c(colorRampPalette(colors = c("#438197","#FFFFFF"))(length(bk)/2),
                                       colorRampPalette(colors = c("#FFFFFF","#CF7660"))(length(bk)/2)),
                             #colorRampPalette(c("#438197","#FFFFFF", "#CF7660"))(50),
                             main = paste0('SI'),
                             #na_col = 'black', 
                             silent = TRUE,
                             legend_breaks=seq(0,1,0.5),
                             breaks=bk
    )
    
  }
  
  #RSI
  if(T){
    files <- list.files('./Data/Step4_SIRSIForLRBench', pattern = '_RSIResult.rds', full.names = TRUE)
    samples <- list.files('./Data/Step4_SIRSIForLRBench', pattern = '_RSIResult.rds')%>%gsub('_RSIResult.rds', '', .)
    
    RSI <- lapply(seq(samples), function(i){
      print(i)
      RSI.tmp <- readRDS(files[i])
      RSI.tmp[lower.tri(as.matrix(RSI.tmp))] <- NA
      if('cellinker' %in% colnames(RSI.tmp) | 'cell2cell' %in% colnames(RSI.tmp)){
        RSI.tmp <- RSI.tmp[-which(colnames(RSI.tmp) %in% c('cellinker', 'cell2cell')),
                           -which(rownames(RSI.tmp) %in% c('cellinker', 'cell2cell'))]
      }
      RSI.tmp <- tibble::rownames_to_column(RSI.tmp, 'method1')
      RSI.tmp <- tidyr::pivot_longer(RSI.tmp, -method1, names_to = 'method2')
      RSI.tmp <- RSI.tmp[-which(is.na(RSI.tmp$value)), ]
      RSI.tmp$dataset <- samples[i]
      RSI.tmp
    })
    
    RSI <- do.call(rbind, RSI)
    RSI <- aggregate(value ~ method1+method2, data = RSI, mean)
    RSI.tmp <- RSI
    colnames(RSI.tmp)[1:2] <- c('method2', 'method1')
    RSI.tmp <- RSI.tmp[,c(2,1,3)]
    RSI <- rbind(RSI,RSI.tmp)
    RSI <- aggregate(value ~ method1+method2, data = RSI, mean)
    RSI <- tidyr::spread(RSI, key = method1, value = value, fill = NA)
    RSI <- tibble::column_to_rownames(RSI, 'method2')
    RSI <- as.matrix(RSI)
    RSI <- RSI[order.names, order.names]
    
    bk <- c(seq(0,0.4,by=0.01),seq(0.5,1,by=0.01))
    p2 <- pheatmap::pheatmap(RSI, cluster_rows = FALSE, cluster_cols = FALSE,
                             treeheight_row = 0,
                             treeheight_col = 0,
                             display_numbers = TRUE,
                             fontsize = 13,
                             #annotation_col = anno,
                             number_format = "%.2f",
                             border_color=NA,
                             angle_col = "45",
                             color = c(colorRampPalette(colors = c("#438197","#FFFFFF"))(length(bk)/2),
                                       colorRampPalette(colors = c("#FFFFFF","#CF7660"))(length(bk)/2)),
                             #colorRampPalette(c("#438197","#FFFFFF", "#CF7660"))(50),
                             main = paste0('RSI'),
                             #na_col = 'black', 
                             silent = TRUE,
                             legend_breaks=seq(0,1,0.5),
                             breaks=bk
    )
    
  }
  
  cowplot::plot_grid(p1$gtable, p2$gtable, nrow = 1, labels = c('A', 'B'), label_size = 20)
}


###############
## Figure 3C ##
###############

if(T){
  all_summary <- readRDS('./Data/Step10_Visualization/all_summary_updata.rds')
  all_summary <- tibble::column_to_rownames(all_summary, 'methods')
  all_summary[, 9] <- all_summary[, 10:20] %>%
    apply(., 2, function(x){(x-min(x, na.rm = TRUE))/(max(x, na.rm = TRUE)-min(x, na.rm = TRUE))}) %>%
    rowMeans(., na.rm = TRUE)
  
  temp1 <- apply(all_summary, 2, function(x){rank(-x, na.last = 'keep')})
  temp1 <- melt(temp1)
  colnames(temp1) <- c('methods', 'Index', 'Rank')
  temp1$class <- paste(temp1$methods, temp1$Index, sep = '-')
  temp1 <- temp1[which(grepl('_MI',temp1$Index)),]
  temp1 <- temp1[, -c(1:2)]
  
  temp2 <- apply(all_summary, 2, 
                 function(x){(x-min(x, na.rm = TRUE))/(max(x, na.rm = TRUE)-min(x, na.rm = TRUE))})
  temp2 <- melt(temp2)
  colnames(temp2) <- c('methods', 'Index', 'Score')
  temp2$class <- paste(temp2$methods, temp2$Index, sep = '-')
  temp2 <- temp2[which(grepl('_MI',temp2$Index)),]
  temp2 <- temp2[, -c(1:2)]
  
  temp <- merge(temp1, temp2, by = 'class')
  temp <- separate(temp, class, c('methods', 'Index'), sep = '-')
  temp$methods <- as.character(temp$methods) %>%
    factor(., levels = c("CellCall", "Domino", "scMLnet", "CytoTalk", "NicheNet", "scSeqComm", # pathway-based
                         "PyMINEr", 'RNAMagnet',  # binary
                         "CellChat", "SingleCellSignalR", "scConnect",# expression value; others
                         "iTALK", "NATMI", "ICELLNET", "Connectome", # expression value; product
                         "CellTalker", "CellPhoneDB3", "CellPhoneDB2" # expression value; mean
    ))
  temp$Index <- as.character(temp$Index) %>%
    factor(., levels = c('all_MI','CID4465_MI', 'CID44971_MI', 'CK357_MI', 'CK358_MI', 
                         'CK161_MI', 'CK165_MI', 'CK361_MI', 'CK362_MI', 'CK162_MI', 'CK368_MI', 'Slide14_MI'))
  
  cols <- brewer.pal(9, "Blues")[c(1,2,3,7,8,9)] #Blues Purples
  pal <- colorRampPalette(cols)
  mycolors <- pal(9) %>% rev(.)
  scales::show_col(mycolors)
  p3 <- ggplot(temp[which(temp$Index != 'all_MI'),], aes(x = Index, y = methods, size = Rank, fill=Score)) + geom_point(colour = "black", shape = 21)+
    scale_size(range = c(10, 1))+
    scale_fill_gradientn(colours  = mycolors)+
    #geom_text_repel(data = temp, aes(label = Rank),
    #                size = 5, segment.color = "black", show.legend = FALSE )+
    theme(axis.text.y=element_text(size=13, color = 'black'), 
          axis.title.y=element_text(size=0),
          axis.text.x=element_text(size=13, angle = 45, vjust=1, hjust=1, color = 'black'), 
          axis.title.x =element_text(size=0))
  #rm(list = ls()[-which(ls()%in%paste0('p', 1:3))])
}



###############
## Figure 3D ##
###############

if(T){
  all_summary <- readRDS('./Data/Step10_Visualization/all_summary_updata.rds')
  all_summary <- tibble::column_to_rownames(all_summary, 'methods')
  temp1 <- apply(all_summary, 2, function(x){rank(-x, na.last = 'keep')})
  temp1 <- melt(temp1)
  colnames(temp1) <- c('methods', 'Index', 'Rank')
  temp1$class <- paste(temp1$methods, temp1$Index, sep = '-')
  temp1 <- temp1[which(grepl('_CAGE',temp1$Index)|grepl('_Proteomic',temp1$Index)),]
  temp1 <- temp1[, -c(1:2)]
  
  temp2 <- apply(all_summary, 2, 
                 function(x){(x-min(x, na.rm = TRUE))/(max(x, na.rm = TRUE)-min(x, na.rm = TRUE))})
  temp2 <- melt(temp2)
  colnames(temp2) <- c('methods', 'Index', 'Score')
  temp2$class <- paste(temp2$methods, temp2$Index, sep = '-')
  temp2 <- temp2[which(grepl('_CAGE',temp2$Index)|grepl('_Proteomic',temp2$Index)),]
  temp2 <- temp2[, -c(1:2)]
  
  temp <- merge(temp1, temp2, by = 'class')
  temp <- separate(temp, class, c('methods', 'Index'), sep = '-')
  temp$methods <- as.character(temp$methods) %>%
    factor(., levels = c("CellCall", "Domino", "scMLnet", "CytoTalk", "NicheNet", "scSeqComm", # pathway-based
                         "PyMINEr", 'RNAMagnet',  # binary
                         "CellChat", "SingleCellSignalR", "scConnect",# expression value; others
                         "iTALK", "NATMI", "ICELLNET", "Connectome", # expression value; product
                         "CellTalker", "CellPhoneDB3", "CellPhoneDB2" # expression value; mean
    ))
  temp$Index <- as.character(temp$Index) %>%
    factor(., levels = c("all_CAGE", "pbmc4k_CAGE", "pbmc6k_CAGE", "pbmc8k_CAGE", 
                         "all_Proteomic", "pbmc4k_Proteomic", "pbmc6k_Proteomic", "pbmc8k_Proteomic"))
  
  cols <- brewer.pal(9, "Blues")[c(1,2,3,7,8,9)] #Blues Purples
  pal <- colorRampPalette(cols)
  mycolors <- pal(9) %>% rev(.)
  scales::show_col(mycolors)
  p4 <- ggplot(temp[-which(temp$Index %in% c('all_CAGE', 'all_Proteomic')),], aes(x = Index, y = methods, size = Rank, fill=Score)) + geom_point(colour = "black", shape = 21)+
    scale_size(range = c(10, 1))+
    scale_fill_gradientn(colours  = mycolors)+
    #geom_text_repel(data = temp, aes(label = Rank),
    #                size = 5, segment.color = "black", show.legend = FALSE )+
    theme(axis.text.y=element_text(size=13, color = 'black'), 
          axis.title.y=element_text(size=0),
          axis.text.x=element_text(size=13, angle = 45, vjust=1, hjust=1, color = 'black'), 
          axis.title.x =element_text(size=0))
  #rm(list = ls()[-which(ls()%in%paste0('p', 1:4))])
}

cowplot::plot_grid(p3, p4, labels = LETTERS[3:4], label_size = 20)

#################
## Figure 3E-H ##
#################

if(F){
  ##########
  ## CAGE ##
  ##########
  setwd('~/2-CCC-benchmark-202204-202208/CCC-benchmark-202212/')
  result <- readRDS('./Data/Step10_Visualization/cage_eval_updata.rds')
  result <- do.call(rbind, result)
  result$methods <- factor(result$methods,
                           levels = c("CellPhoneDB2", "CellPhoneDB3", "CellTalker", 
                                      "Connectome", "ICELLNET",  "NATMI", "iTALK", 
                                      "scConnect", "SingleCellSignalR", "CellChat", 
                                      "RNAMagnet", "PyMINEr",
                                      "scSeqComm", "NicheNet","CytoTalk",
                                      "scMLnet", "Domino", "CellCall"))
  
  methods <- as.character(unique(result$methods))
  
  stat_result_cage <- lapply(methods, function(method1){
    result1 <- result[which(result$methods==method1), ]
    tmp <- lapply(methods, function(method2){
      result2 <- result[which(result$methods==method2), ]
      wilcox.test(result1$AUPRC, result2$AUPRC, alternative = "greater")$p.value
    })
    tmp <- do.call(cbind, tmp)
    colnames(tmp) <- methods
    tmp
  })
  stat_result_cage <- do.call(rbind, stat_result_cage)
  rownames(stat_result_cage) <- methods
  stat_binary_cage <- ifelse(stat_result_cage<0.05, 1, 0)
  
  order_names <- c("CellPhoneDB2", "CellPhoneDB3", "CellTalker", 
                   "Connectome", "ICELLNET",  "NATMI", "iTALK", 
                   "scConnect", "SingleCellSignalR", "CellChat", 
                   "RNAMagnet", "PyMINEr",
                   "scSeqComm", "NicheNet","CytoTalk",
                   "scMLnet", "Domino", "CellCall")
  stat_binary_cage <- stat_binary_cage[order_names, order_names]
  
  p5 <- ggplot(result, aes(methods, AUPRC, fill=methods)) + 
    geom_boxplot(lwd = 0.5, position=position_dodge(width=0.8))+ #outlier.shape = NA, 
    theme_bw()+
    ylab("AUPRC")+
    xlab("")+
    theme(text = element_text(size = 15),
          plot.title = element_text(size = 16, hjust = 0.5, face = 'bold'),
          axis.text.y=element_text(size=13, color = 'black'),
          axis.text.x=element_text(size=13, angle = 45, vjust=1, hjust=1, color = 'black'),
          legend.position="none")+
    ggtitle('CAGE')
  
  p6 <- pheatmap::pheatmap(stat_binary_cage,cluster_rows = FALSE, cluster_cols = FALSE,
                           treeheight_row = 0,
                           treeheight_col = 0,
                           display_numbers = TRUE,number_color = 'white',
                           fontsize=13,
                           #annotation_col = anno,
                           number_format = "%.0f",
                           border_color=NA,
                           angle_col = "45", 
                           color = colorRampPalette(c("#B4BFB9", "#6499A7"))(2),
                           main = 'CAGE',
                           #na_col = 'black', 
                           silent = FALSE, legend = FALSE
  )
 
  cowplot::plot_grid(p5, p6$gtable, ncol = 2, labels = LETTERS[5:6], label_size = 20)
  
  ################
  ## Proteomics ##
  ################
  
  setwd('~/2-CCC-benchmark-202204-202208/CCC-benchmark-202212/')
  result <- readRDS('./Data/Step10_Visualization/proteomic_eval_updata.rds')
  result <- do.call(rbind, result)
  result$methods <- factor(result$methods,
                           levels = c("CellPhoneDB2", "CellPhoneDB3", "CellTalker",
                                      "Connectome", "ICELLNET",  "NATMI", "iTALK",
                                      "scConnect", "SingleCellSignalR", "CellChat",
                                      "RNAMagnet", "PyMINEr",
                                      "scSeqComm", "NicheNet","CytoTalk",
                                      "scMLnet", "Domino", "CellCall"))
  
  methods <- as.character(unique(result$methods))
  
  stat_result_proteomic <- lapply(methods, function(method1){
    result1 <- result[which(result$methods==method1), ]
    tmp <- lapply(methods, function(method2){
      result2 <- result[which(result$methods==method2), ]
      wilcox.test(result1$AUPRC, result2$AUPRC, alternative = "greater")$p.value
    })
    tmp <- do.call(cbind, tmp)
    colnames(tmp) <- methods
    tmp
  })
  stat_result_proteomic <- do.call(rbind, stat_result_proteomic)
  rownames(stat_result_proteomic) <- methods
  stat_binary_proteomic <- ifelse(stat_result_proteomic<0.05, 1, 0)
  
  order_names <- c("CellPhoneDB2", "CellPhoneDB3", "CellTalker", 
                   "Connectome", "ICELLNET",  "NATMI", "iTALK", 
                   "scConnect", "SingleCellSignalR", "CellChat", 
                   "RNAMagnet", "PyMINEr",
                   "scSeqComm", "NicheNet","CytoTalk",
                   "scMLnet", "Domino", "CellCall")
  stat_binary_proteomic <- stat_binary_proteomic[order_names, order_names]
  
  p7 <- ggplot(result, aes(methods, AUPRC, fill=methods)) + 
    geom_boxplot(lwd = 0.5, position=position_dodge(width=0.8))+ #outlier.shape = NA, 
    theme_bw()+
    ylab("AUPRC")+
    xlab("")+
    theme(text = element_text(size = 15),
          plot.title = element_text(size = 16, hjust = 0.5, face = 'bold'),
          axis.text.y=element_text(size=13, color = 'black'),
          axis.text.x=element_text(size=13, angle = 45, vjust=1, hjust=1, color = 'black'),
          legend.position="none")+
    ggtitle('Proteomics')
  
  p8 <- pheatmap::pheatmap(stat_binary_proteomic,cluster_rows = FALSE, cluster_cols = FALSE,
                           treeheight_row = 0,
                           treeheight_col = 0,
                           display_numbers = TRUE,number_color = 'white',
                           fontsize=13,
                           #annotation_col = anno,
                           number_format = "%.0f",
                           border_color=NA,
                           angle_col = "45",
                           color = colorRampPalette(c("#B4BFB9", "#6499A7"))(2),
                           main = 'Proteomics',
                           #na_col = 'black', 
                           silent = FALSE, legend = FALSE
  )
  
  cowplot::plot_grid(p7, p8$gtable, ncol = 2, labels = LETTERS[7:8], label_size = 20)
}

#cowplot::plot_grid(p1$gtable, p2$gtable, p4, p5, nrow = 2, labels = c('A', 'B', 'C', 'D'))


#################
## Figure 4A-C ##
#################

result <- readRDS(file = './Data/Step7_LRBenchSampling/JaccardIndex.rds')
result$stability <- lapply(seq(dim(result)[[1]]), function(x){
  tmp <- 1/sum(1 - result$S50[x], 1 - result$S60[x], 1 - result$S70[x], 
             1 - result$S80[x], 1 - result$S90[x])
  tmp
})%>% unlist
result$methods <- factor(result$methods,levels = c("CellPhoneDB2", "CellPhoneDB3", "CellTalker", 
                                                   "Connectome", "ICELLNET",  "NATMI", "iTALK", 
                                                   "scConnect", "SingleCellSignalR", "CellChat", 
                                                   "RNAMagnet", "PyMINEr",
                                                   "scSeqComm", "NicheNet","CytoTalk",
                                                   "scMLnet", "Domino", "CellCall"))

temp <- reshape2::melt(result, id.vars = c('methods', 'datasets'), 
                       measure.vars = c('S50', 'S60', 'S70', 'S80', 'S90'),
                       variable.name = 'sampleRates',
                       value.name = 'jaccard')
temp$sampleRates <- factor(temp$sampleRates, levels = c('S90','S80', 'S70', 'S60', 'S50'))


scales::show_col(pal_aaas(palette = "default", alpha = 0.6)(10))
mycolors_aaas <- pal_aaas(palette = "default", alpha = 0.6)(10)

mycolor_software <- mycolors_aaas[c(1:5)]
names(mycolor_software) <- c('S90','S80', 'S70', 'S60', 'S50')
scales::show_col(mycolor_software)

p1 <- ggplot(temp, aes(methods, jaccard, fill=sampleRates)) +
  geom_boxplot(outlier.shape = NA, lwd = 0.5, position=position_dodge(width=0.8))+
  theme_bw()+
  scale_fill_manual(values = mycolor_software,
                    labels=c('S90','S80', 'S70', 'S60', 'S50'))+
  ylab("Jaccard Index")+
  xlab("")+
  theme(text = element_text(size = 15),
        plot.title = element_text(hjust = 0.5),
        axis.text.y=element_text(size=13, color = 'black'),
        axis.text.x=element_text(size=13, angle = 45, vjust=1, hjust=1, color = 'black'))

p2 <- ggplot(result, aes(methods, stability, fill=methods)) + 
  geom_boxplot(lwd = 0.5, position=position_dodge(width=0.8))+ #outlier.shape = NA, 
  theme_bw()+
  ylab("Stability")+
  xlab("")+
  theme(text = element_text(size = 15),
        plot.title = element_text(hjust = 0.5),
        axis.text.y=element_text(size=13, color = 'black'),
        axis.text.x=element_text(size=13, angle = 45, vjust=1, hjust=1, color = 'black'),
        legend.position="none")

methods <- unique(result$methods)%>%as.character()
stat_result <- lapply(methods, function(method1){
  result1 <- result[which(result$methods==method1), ]
  tmp <- lapply(methods, function(method2){
    result2 <- result[which(result$methods==method2), ]
    wilcox.test(result1$stability, result2$stability, alternative = "greater")$p.value
  })
  tmp <- do.call(cbind, tmp)
  colnames(tmp) <- methods
  tmp
})
stat_result <- do.call(rbind, stat_result)
rownames(stat_result) <- methods
stat_binary <- ifelse(stat_result<0.05, 1, 0)

order_names <- c("CellPhoneDB2", "CellPhoneDB3", "CellTalker", 
                 "Connectome", "ICELLNET",  "NATMI", "iTALK", 
                 "scConnect", "SingleCellSignalR", "CellChat", 
                 "RNAMagnet", "PyMINEr",
                 "scSeqComm", "NicheNet","CytoTalk",
                 "scMLnet", "Domino", "CellCall")
stat_binary <- stat_binary[order_names, order_names]

p3 <- pheatmap::pheatmap(stat_binary,cluster_rows = FALSE, cluster_cols = FALSE,
                         treeheight_row = 0,
                         treeheight_col = 0,
                         display_numbers = TRUE,number_color = 'white', 
                         fontsize=13,
                         #annotation_col = anno,
                         number_format = "%.0f",
                         border_color=NA,
                         angle_col = "45",
                         color = colorRampPalette(c("#B4BFB9", "#6499A7"))(2),
                         #main = 'Proteomics',
                         #na_col = 'black', 
                         silent = FALSE, legend = FALSE
)

cowplot::plot_grid(p1,labels = LETTERS[1], label_size = 20)
cowplot::plot_grid(p2,p3$gtable, labels = LETTERS[2:3], label_size = 20, ncol = 2)

#################
## Figure 5A-B ##
#################


if(T){
  record <- readRDS('./Data/Step7_LRBenchSampling/TimeMemRecord.rds')
  
  temp1 <- record
  temp1$linux_time <- temp1$linux_time/60/60/24
  temp1 <- temp1[which(temp1$methods %in% c('scMLnet','CytoTalk','Domino')),]
  temp1$linux_time <- round(temp1$linux_time, 2)
  p1 <- ggplot(temp1,aes(Cells,linux_time,group=methods,color=methods,shape=methods))+
    scale_shape_manual(values = 1:3)+
    geom_point(size=3.5)+
    #geom_text(aes(label = linux_time),vjust=-0.2)+
    geom_line(position = position_dodge(0.1),cex=1.3)+
    labs(x="Datasets",y = "User time + System time (day)")+
    theme(legend.position = "none")
  
  temp1 <- record
  temp1$linux_time <- temp1$linux_time/60/60
  temp1 <- temp1[which(temp1$methods %in% c('PyMINEr', 'NicheNet', 'CellCall', 'CellPhoneDB2')),]
  temp1$linux_time <- round(temp1$linux_time, 2)
  p2 <- ggplot(temp1,aes(Cells,linux_time,group=methods,color=methods,shape=methods))+
    scale_shape_manual(values = 4:7)+
    geom_point(size=3.5)+
    #geom_text(aes(label = linux_time),vjust=-0.2)+
    geom_line(position = position_dodge(0.1),cex=1.3)+
    labs(x="Cells",y = "User time + System time (hour)")+
    theme(legend.position = "none")
  
  temp1 <- record
  temp1$linux_time <- temp1$linux_time/60
  temp1 <- temp1[-which(temp1$methods %in% c('CytoTalk', 'scMLnet','Domino', 'PyMINEr', 'NicheNet', 'CellCall', 'CellPhoneDB2')),]
  p3 <- ggplot(temp1,aes(Cells,linux_time,group=methods,color=methods,shape=methods))+
    scale_shape_manual(values = 7:19)+
    geom_point(size=3.5)+
    geom_line(position = position_dodge(0.1),cex=1.3)+
    labs(x="Cells",y = "User time + System time (min)")+
    theme(legend.position = "none")
  p <- cowplot::plot_grid(p3, p2, p1, nrow = 3)
  cowplot::plot_grid(p, labels = LETTERS[1], label_size = 20)
  
  temp1 <- record
  p4 <- ggplot(temp1,aes(Cells,max_memory,group=methods,color=methods, shape = methods))+
    scale_shape_manual(values = 1:18)+
    geom_point(size=4)+
    geom_line(position = position_dodge(0.1),cex=1.3)+
    labs(x="Cells",y = "Maximum memory usage(Gb)")
  gg.gap(plot = p4,
         segments = c(c(41, 42)),
         tick_width = c(10,50),
         ylim = c(0, 300))
  
}


##############
## Figure 6 ##
##############

if(T){
  setwd('~/2-CCC-benchmark-202204-202208/CCC-benchmark-202212/')
  
  ###############
  ## Figure 6A ##
  ###############
  
  result_index <- readRDS('./Data/Step8_LRTBenchResult/result_index.rds')
  result_index <- do.call(rbind, result_index)
  result_index <- tibble::rownames_to_column(result_index, 'datasets')
  result_index$datasets <- gsub('\\.[0-9]+', '', result_index$datasets)
  result_index <- result_index[which(result_index$methods %in% 
                                       c('CytoTalk', "NicheNet", "MISTy", "HoloNet", "stMLnet")), ]
  result_index$methods <- factor(result_index$methods,
                                 levels = c('CytoTalk', "NicheNet", "MISTy", "HoloNet", "stMLnet"))
  
  tmp <- result_index[result_index$datasets %in% 
                        c("CID4465", "CID44971", "control_P7", "control_P8", 
                          'UKF260_T_ST', 'UKF266_T_ST', 'UKF334_T_ST', 'UKF243_T_ST'), ]
  tmp <- tmp[which(!grepl('TAM_RB', tmp$celllines)), ]
  
  scales::show_col(pal_aaas(palette = "default", alpha = 0.6)(10))
  mycolors_aaas <- pal_aaas(palette = "default", alpha = 0.6)(10)
  mycolor_software <- mycolors_aaas[c(1:5)]
  names(mycolor_software) <- c('CytoTalk', "NicheNet", "MISTy", "HoloNet", "stMLnet")
  scales::show_col(mycolor_software)
  
  p1 <- ggplot(tmp, aes(methods, AUPRC, fill=methods)) + 
    geom_boxplot(lwd = 0.5, position=position_dodge(width=0.8))+ #outlier.shape = NA, 
    theme_bw()+
    ylab("AUPRC")+
    xlab("")+
    theme(text = element_text(size = 15),
          plot.title = element_text(size = 16, hjust = 0.5, face = 'bold'),
          axis.text.y=element_text(size=13, color = 'black'),
          axis.text.x=element_text(size=13, angle = 45, vjust=1, hjust=1, color = 'black'),
          legend.position="none")+
    scale_fill_manual(values = mycolor_software, 
                     labels=c('CytoTalk', "NicheNet", "MISTy", "HoloNet", "stMLnet"))

  
  ###############
  ## Figure 6B ##
  ###############
  
  methods <- as.character(unique(result_index$methods))
  stat_result <- lapply(methods, function(method1){
    result1 <- result_index[which(result_index$methods==method1), ]
    tmp <- lapply(methods, function(method2){
      result2 <- result_index[which(result_index$methods==method2), ]
      wilcox.test(result1$AUPRC, result2$AUPRC, alternative = "greater")$p.value
    })
    tmp <- do.call(cbind, tmp)
    colnames(tmp) <- methods
    tmp
  })
  stat_result <- do.call(rbind, stat_result)
  rownames(stat_result) <- methods
  stat_binary <- ifelse(stat_result<0.05, 1, 0)
  order_names <- c('CytoTalk', "NicheNet", "MISTy", "HoloNet", "stMLnet")
  stat_binary <- stat_binary[order_names, order_names]
  
  p2 <- pheatmap::pheatmap(stat_binary,cluster_rows = FALSE, cluster_cols = FALSE,
                           treeheight_row = 0,
                           treeheight_col = 0,legend = FALSE,
                           display_numbers = TRUE,number_color = 'white',
                           fontsize=16,
                           #annotation_col = anno,
                           number_format = "%.0f",
                           border_color=NA,
                           angle_col = "45",
                           color = colorRampPalette(c("#B4BFB9", "#6499A7"))(2),
                           #main = paste0('SI in ', data),
                           #na_col = 'black', 
                           silent = FALSE
  )
  cowplot::plot_grid(p1, p2$gtable, labels = LETTERS[1:2], label_size = 20)
  
  
  #################
  ## Figure 6C-D ##
  ################# 
  
  if(F){
    rm(list = ls()[-which(ls()%in% paste0('p',1:3))]); gc()
    source('./Script/Step9_LRTBench/function.R')
    
    datasets <- list.dirs('./Data/Step8_LRTPredictionResult', recursive = FALSE, full.names = FALSE)
    methods <- c('CytoTalk', 'NicheNet', 'HoloNet', 'MISTy', 'stMLnet')
    DataRecord <- readRDS('./Data/Step8_LRTBenchResult/DatasetCellTypeInfo.rds')
    wd <- './Data/Step7_LRTPredictionResult/'
    
    TimeMemRecord <- lapply(datasets, function(data){
      print(data)
      temp <- RunTimeMemRecord(wd, data, methods)
      temp
    })
    names(TimeMemRecord) <- datasets
    TimeMemRecord <- lapply(TimeMemRecord, function(data){
      temp <- do.call(rbind, data)
      temp <- tibble::rownames_to_column(temp, 'methods')
      temp
    })
    TimeMemRecord <- do.call(rbind, TimeMemRecord)
    TimeMemRecord <- tibble::rownames_to_column(TimeMemRecord, 'datasets')
    TimeMemRecord$datasets <- gsub('\\.[0-9]+', '', TimeMemRecord$datasets)
    
    # Specific handle: the running time of CytoTalk in GBM datasets
    CytoTalk_TimeMemRecord <- TimeMemRecord[which(TimeMemRecord$methods == 'CytoTalk' & grepl('UKF', TimeMemRecord$datasets)), ]
    CytoTalk_TimeMemRecord <- merge(CytoTalk_TimeMemRecord, DataRecord, by = 'datasets')
    CytoTalk_TimeMemRecord$total_num <- 2*(CytoTalk_TimeMemRecord$num-1)
    CytoTalk_TimeMemRecord$actual_num <- 2*(CytoTalk_TimeMemRecord$num-1)-1
    CytoTalk_TimeMemRecord$clock_time <- (CytoTalk_TimeMemRecord$clock_time/CytoTalk_TimeMemRecord$total_num)*CytoTalk_TimeMemRecord$actual_num
    CytoTalk_TimeMemRecord$linux_time <- (CytoTalk_TimeMemRecord$linux_time/CytoTalk_TimeMemRecord$total_num)*CytoTalk_TimeMemRecord$actual_num
    CytoTalk_TimeMemRecord <- CytoTalk_TimeMemRecord[, -c(7:9)]
    TimeMemRecord <- TimeMemRecord[-which(TimeMemRecord$methods == 'CytoTalk' & grepl('UKF', TimeMemRecord$datasets)), ]
    TimeMemRecord <- rbind(TimeMemRecord, CytoTalk_TimeMemRecord)
    rm(CytoTalk_TimeMemRecord, DataRecord, datasets, methods, wd); gc()
    saveRDS(TimeMemRecord, file = './Data/Step10_Visualization/LRT_RunTimeMem.rds')
  }
  
  TimeMemRecord <- readRDS('./Data/Step10_Visualization/LRT_RunTimeMem.rds') 
  
  # visualization
  ## Running time
  if(T){
    scales::show_col(pal_aaas(palette = "default", alpha = 1)(10))
    mycolors_aaas <- pal_aaas(palette = "default", alpha = 1)(10)
    mycolor_software <- mycolors_aaas[c(1:5)]
    names(mycolor_software) <- c('CytoTalk', "NicheNet", "MISTy", "HoloNet", "stMLnet")
    scales::show_col(mycolor_software)
    
    temp1 <- TimeMemRecord;
    temp1$linux_time <- temp1$linux_time/60/60/24
    temp1 <- temp1[which(temp1$methods %in% c('CytoTalk', 'HoloNet','stMLnet')),]
    temp1$methods <- factor(temp1$methods, levels = c('CytoTalk', 'HoloNet','stMLnet'))
    temp1$linux_time <- round(temp1$linux_time, 3)
    p4 <- ggplot(temp1,aes(datasets,linux_time,group=methods,color=methods,shape=methods))+
      scale_shape_manual(values = 1:3)+
      geom_point(size=3.5)+
      geom_text(aes(label = linux_time),vjust=-0.2)+
      geom_line(position = position_dodge(0.1),cex=1.3)+
      theme(axis.title.x =element_text(size=16), 
            axis.title.y=element_text(size=16), 
            text = element_text(size=16), 
            axis.text.x = element_text(colour = 'black', angle = 45, vjust=1, hjust=1), 
            axis.text.y = element_text(colour = 'black'))+
      labs(x="Datasets",y = "User + System time (day)")
   
    p4
    
    temp2 <- TimeMemRecord;
    temp2$linux_time <- temp2$linux_time/60
    temp2 <- temp2[which(temp2$methods %in% c('NicheNet', 'MISTy')),]
    temp2$methods <- factor(temp2$methods, levels = c('NicheNet', 'MISTy'))
    temp2$linux_time <- round(temp2$linux_time, 2)
    p5 <- ggplot(temp2,aes(datasets,linux_time,group=methods,color=methods,shape=methods))+
      scale_shape_manual(values = 4:5)+
      geom_point(size=3.5)+
      geom_text(aes(label = linux_time),vjust=-0.2)+
      geom_line(position = position_dodge(0.1),cex=1.3)+
      theme(axis.title.x =element_text(size=16), 
            axis.title.y=element_text(size=16), 
            text = element_text(size=16), 
            axis.text.x = element_text(colour = 'black', angle = 45, vjust=1, hjust=1), 
            axis.text.y = element_text(colour = 'black'))+
      labs(x="Datasets",y = "User + System time (min)")
    p5
    
    p6 <- cowplot::plot_grid(p5, p4, nrow = 2)
  }
  
  # Maximum memory usage
  if(T){
    temp <- TimeMemRecord;
    temp$methods <- factor(temp$methods, levels = c('CytoTalk', 'HoloNet','stMLnet', 'NicheNet', 'MISTy'))
    p7 <- ggplot(temp,aes(datasets,max_memory,group=methods,color=methods, shape = methods))+
      scale_shape_manual(values = 1:5)+
      geom_point(size=4)+
      geom_line(position = position_dodge(0.1),cex=1.3)+
      geom_text(aes(label = max_memory),vjust=-0.2)+
      labs(x="Datasets",y = "Maximum memory usage(Gb)")+
      theme(axis.title.x =element_blank(), 
            axis.title.y=element_text(size=16), 
            axis.text.x = element_text(colour = 'black', angle = 45, vjust=1, hjust=1),
            axis.text.y = element_text(colour = 'black'),
            text = element_text(size=16))
    gg.gap(plot = p7,
           segments = c(55, 80),
           tick_width = c(25,25),
           ylim = c(0, 125),
    )
  }
  cowplot::plot_grid(p6, p7, nrow = 1, labels = LETTERS[3:4], label_size = 20)
}




############################
## Supplementary Figure 1 ##
############################

setwd('~/2-CCC-benchmark-202204-202208/CCC-benchmark-202212/Data/Step4_SIRSIForLRBench/')

datasets <- c('CID4465', 'CID44971', 'CK357', 'CK358',
              'CK368', 'CK162', 'CK362', 'CK361',
              'CK161', 'CK165', 'Slide14')

for (i in seq(datasets)) {
  data <- datasets[i]
  print(data)
  
  order.names <- c("CellPhoneDB2", "CellPhoneDB3", "CellTalker", 
                   "Connectome", "ICELLNET",  "NATMI", "iTALK", 
                   "scConnect", "SingleCellSignalR", "CellChat", 
                   "RNAMagnet", "PyMINEr",
                   "scSeqComm", "NicheNet","CytoTalk",
                   "scMLnet", "Domino", "CellCall")
  
  RSI <- readRDS(paste0(data, '_RSIResult.rds'))
  if('cellinker' %in% colnames(RSI) | 'cell2cell' %in% colnames(RSI)){
    RSI <- RSI[-which(colnames(RSI) %in% c('cellinker', 'cell2cell')),
               -which(rownames(RSI) %in% c('cellinker', 'cell2cell'))]
  }
  RSI <- RSI[!apply(RSI, 1, function(x){all(is.na(x))}), ]
  RSI <- RSI[, !apply(RSI, 2, function(x){all(is.na(x))})]
  order.names.rsi <- order.names[which(order.names %in% colnames(RSI))]
  RSI <- RSI[order.names.rsi, order.names.rsi]
  
  SI <- readRDS(paste0(data, '_SIResult.rds'))
  SI <- SI[!apply(SI, 1, function(x){all(is.na(x))}), ]
  SI <- SI[, !apply(SI, 2, function(x){all(is.na(x))})]
  if('cellinker' %in% colnames(SI) | 'cell2cell' %in% colnames(SI)){
    SI <- SI[-which(colnames(SI) %in% c('cellinker', 'cell2cell')),
             -which(rownames(SI) %in% c('cellinker', 'cell2cell'))]
  }
  order.names.si <- order.names[which(order.names %in% colnames(SI))]
  SI <- SI[order.names.si, order.names.si]
  
  bk <- c(seq(0,0.4,by=0.01),seq(0.5,1,by=0.01))
  p <- pheatmap::pheatmap(SI,cluster_rows = FALSE, cluster_cols = FALSE,
                           treeheight_row = 0,
                           treeheight_col = 0,
                           display_numbers = TRUE,
                           fontsize=13,
                           #annotation_col = anno,
                           number_format = "%.2f",
                           border_color=NA,
                           angle_col = "45",
                           color = c(colorRampPalette(colors = c("#438197","#FFFFFF"))(length(bk)/2),
                                     colorRampPalette(colors = c("#FFFFFF","#CF7660"))(length(bk)/2)),
                           #colorRampPalette(c("#438197","#FFFFFF", "#CF7660"))(50),
                           main = paste0('SI in ', data),
                           #na_col = 'black', 
                           silent = TRUE,
                           legend_breaks=seq(0,1,0.5),
                           breaks=bk
  )
  #assign(paste0('p', i), p)
  pp <- pheatmap::pheatmap(RSI,cluster_rows = FALSE, cluster_cols = FALSE,
                           treeheight_row = 0,
                           treeheight_col = 0,
                           display_numbers = TRUE,
                           fontsize=13,
                           #annotation_col = anno,
                           number_format = "%.2f",
                           border_color=NA,
                           angle_col = "45",
                           color = c(colorRampPalette(colors = c("#438197","#FFFFFF"))(length(bk)/2),
                                     colorRampPalette(colors = c("#FFFFFF","#CF7660"))(length(bk)/2)),
                           #colorRampPalette(c("#438197","#FFFFFF", "#CF7660"))(50),
                           main = paste0('RSI in ', data),
                           na_col = 'black', silent = TRUE,
                           legend_breaks=seq(0,1,0.5),
                           breaks=bk
  )
  ppp <- cowplot::plot_grid(p$gtable, pp$gtable, ncol= 2)
  assign(paste0('p', i), ppp)
}
cowplot::plot_grid(p1, p2, nrow= 2)
cowplot::plot_grid(p3, p4, nrow= 2)
cowplot::plot_grid(p5, p6, nrow= 2)
cowplot::plot_grid(p7, p8, nrow= 2)
cowplot::plot_grid(p9, p10, nrow= 2)


##############################
## Supplementary Figure 2-6 ##
##############################

if(F){
  setwd('~/2-CCC-benchmark-202204-202208/CCC-benchmark-202212/')
  source('./Script/Step3_MIForLRBench/function.R')
  set.seed(123)
  
  files <- list.files('./Data/Step3_MIForLRBench', pattern = '_result.rds', full.names = TRUE)
  Eval1Result <- lapply(files, function(file){
    temp <- Eval1Process(file)
    temp$datasets <- substring(file, 30) %>% gsub('_result.rds', '', .)
    temp
  })
  Eval1Result <- do.call(rbind, Eval1Result)
  data1 <- Eval1Result[,c(1:4, 10:12)]
  data1$group <- 'close'; colnames(data1)[4] <- 'MI'
  data2 <- Eval1Result[,c(1:3,7,10:12)]
  data2$group <- 'distant'; colnames(data2)[4] <- 'MI'
  data <- rbind(data1, data2)
  data$methods <- factor(data$methods,levels = c("CellPhoneDB2", "CellPhoneDB3", "CellTalker", 
                                                 "Connectome", "ICELLNET",  "NATMI", "iTALK", 
                                                 "scConnect", "SingleCellSignalR", "CellChat", 
                                                 "RNAMagnet", "PyMINEr",
                                                 "scSeqComm", "NicheNet","CytoTalk",
                                                 "scMLnet", "Domino", "CellCall"))
  rm(data1, data2);gc()
  
  #############
  ## CID4465 ##
  #############
  dataset <- 'CID4465'
  per <- 10
  data.tmp <- data[which(data$perc==per & data$datasets == dataset),]
  p1 <- ggplot(data.tmp, aes(methods, MI, fill=group)) + 
    geom_boxplot(outlier.shape = NA, lwd = 0.5)+
    theme_bw()+
    theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1), plot.title = element_text(hjust = 0.5)) +
    ggpubr::stat_compare_means(aes(group=group), label = "p.signif", method = 'wilcox.test',
                               method.args = list(alternative = "less"),
                               label.y = 0.05)+
    coord_cartesian(ylim=c(0,0.05))+
    ylab("Mutual information")+
    ggtitle(paste0("top", per, "% close/distant cell pairs in ", dataset))+
    xlab("")+
    theme(text = element_text(size = 13),
          axis.text.y=element_text(size=12, color = 'black'),
          axis.text.x=element_text(size=12, angle = 45, vjust=1, hjust=1, color = 'black'))
  
  per <- 20
  data.tmp <- data[which(data$perc==per & data$datasets == dataset),]
  p2 <- ggplot(data.tmp, aes(methods, MI, fill=group)) + 
    geom_boxplot(outlier.shape = NA, lwd = 0.5)+
    theme_bw()+
    theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1), plot.title = element_text(hjust = 0.5)) +
    ggpubr::stat_compare_means(aes(group=group), label = "p.signif", method = 'wilcox.test',
                               method.args = list(alternative = "less"),
                               label.y = 0.04)+
    coord_cartesian(ylim=c(0,0.04))+
    ylab("Mutual information")+
    ggtitle(paste0("top", per, "% close/distant cell pairs in ", dataset))+
    xlab("")+
    theme(text = element_text(size = 13),
          axis.text.y=element_text(size=12, color = 'black'),
          axis.text.x=element_text(size=12, angle = 45, vjust=1, hjust=1, color = 'black'))
  
  per <- 30
  data.tmp <- data[which(data$perc==per & data$datasets == dataset),]
  p3 <- ggplot(data.tmp, aes(methods, MI, fill=group)) + 
    geom_boxplot(outlier.shape = NA, lwd = 0.5)+
    theme_bw()+
    theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1), plot.title = element_text(hjust = 0.5)) +
    ggpubr::stat_compare_means(aes(group=group), label = "p.signif", method = 'wilcox.test',
                               method.args = list(alternative = "less"),
                               label.y = 0.03)+
    coord_cartesian(ylim=c(0,0.03))+
    ylab("Mutual information")+
    ggtitle(paste0("top", per, "% close/distant cell pairs in ", dataset))+
    xlab("")+
    theme(text = element_text(size = 13),
          axis.text.y=element_text(size=12, color = 'black'),
          axis.text.x=element_text(size=12, angle = 45, vjust=1, hjust=1, color = 'black'))
  
  per <- 40
  data.tmp <- data[which(data$perc==per & data$datasets == dataset),]
  p4 <- ggplot(data.tmp, aes(methods, MI, fill=group)) + 
    geom_boxplot(outlier.shape = NA, lwd = 0.5)+
    theme_bw()+
    theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1), plot.title = element_text(hjust = 0.5)) +
    ggpubr::stat_compare_means(aes(group=group), label = "p.signif", method = 'wilcox.test',
                               method.args = list(alternative = "less"),
                               label.y = 0.02)+
    coord_cartesian(ylim=c(0,0.02))+
    ylab("Mutual information")+
    ggtitle(paste0("top", per, "% close/distant cell pairs in ", dataset))+
    xlab("")+
    theme(text = element_text(size = 13),
          axis.text.y=element_text(size=12, color = 'black'),
          axis.text.x=element_text(size=12, angle = 45, vjust=1, hjust=1, color = 'black'))
  
  cowplot::plot_grid(p1, p2, p3, p4, nrow = 2, labels = LETTERS[1:4], label_size = 20)
  
  ##############
  ## CID44971 ##
  ##############
  
  dataset <- 'CID44971'
  per <- 10
  data.tmp <- data[which(data$perc==per & data$datasets == dataset),]
  p1 <- ggplot(data.tmp, aes(methods, MI, fill=group)) + 
    geom_boxplot(outlier.shape = NA, lwd = 0.5)+
    theme_bw()+
    theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1), plot.title = element_text(hjust = 0.5)) +
    ggpubr::stat_compare_means(aes(group=group), label = "p.signif", method = 'wilcox.test',
                               method.args = list(alternative = "less"),
                               label.y = 0.06)+
    coord_cartesian(ylim=c(0,0.06))+
    ylab("Mutual information")+
    ggtitle(paste0("top", per, "% close/distant cell pairs in ", dataset))+
    xlab("")+
    theme(text = element_text(size = 13),
          axis.text.y=element_text(size=12, color = 'black'),
          axis.text.x=element_text(size=12, angle = 45, vjust=1, hjust=1, color = 'black'))
  p1
  
  per <- 20
  data.tmp <- data[which(data$perc==per & data$datasets == dataset),]
  p2 <- ggplot(data.tmp, aes(methods, MI, fill=group)) + 
    geom_boxplot(outlier.shape = NA, lwd = 0.5)+
    theme_bw()+
    theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1), plot.title = element_text(hjust = 0.5)) +
    ggpubr::stat_compare_means(aes(group=group), label = "p.signif", method = 'wilcox.test',
                               method.args = list(alternative = "less"),
                               label.y = 0.05)+
    coord_cartesian(ylim=c(0,0.05))+
    ylab("Mutual information")+
    ggtitle(paste0("top", per, "% close/distant cell pairs in ", dataset))+
    xlab("")+
    theme(text = element_text(size = 13),
          axis.text.y=element_text(size=12, color = 'black'),
          axis.text.x=element_text(size=12, angle = 45, vjust=1, hjust=1, color = 'black'))
  p2
  
  per <- 30
  data.tmp <- data[which(data$perc==per & data$datasets == dataset),]
  p3 <- ggplot(data.tmp, aes(methods, MI, fill=group)) + 
    geom_boxplot(outlier.shape = NA, lwd = 0.5)+
    theme_bw()+
    theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1), plot.title = element_text(hjust = 0.5)) +
    ggpubr::stat_compare_means(aes(group=group), label = "p.signif", method = 'wilcox.test',
                               method.args = list(alternative = "less"),
                               label.y = 0.04)+
    coord_cartesian(ylim=c(0,0.04))+
    ylab("Mutual information")+
    ggtitle(paste0("top", per, "% close/distant cell pairs in ", dataset))+
    xlab("")+
    theme(text = element_text(size = 13),
          axis.text.y=element_text(size=12, color = 'black'),
          axis.text.x=element_text(size=12, angle = 45, vjust=1, hjust=1, color = 'black'))
  p3
  
  per <- 40
  data.tmp <- data[which(data$perc==per & data$datasets == dataset),]
  p4 <- ggplot(data.tmp, aes(methods, MI, fill=group)) + 
    geom_boxplot(outlier.shape = NA, lwd = 0.5)+
    theme_bw()+
    theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1), plot.title = element_text(hjust = 0.5)) +
    ggpubr::stat_compare_means(aes(group=group), label = "p.signif", method = 'wilcox.test',
                               method.args = list(alternative = "less"),
                               label.y = 0.03)+
    coord_cartesian(ylim=c(0,0.03))+
    ylab("Mutual information")+
    ggtitle(paste0("top", per, "% close/distant cell pairs in ", dataset))+
    xlab("")+
    theme(text = element_text(size = 13),
          axis.text.y=element_text(size=12, color = 'black'),
          axis.text.x=element_text(size=12, angle = 45, vjust=1, hjust=1, color = 'black'))
  p4
  
  cowplot::plot_grid(p1, p2, p3, p4, nrow = 2, labels = LETTERS[5:8], label_size = 20)
  
  ###########
  ## CK357 ##
  ###########
  
  dataset <- 'CK357'
  per <- 10
  data.tmp <- data[which(data$perc==per & data$datasets == dataset),]
  p1 <- ggplot(data.tmp, aes(methods, MI, fill=group)) + 
    geom_boxplot(outlier.shape = NA, lwd = 0.5)+
    theme_bw()+
    theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1), plot.title = element_text(hjust = 0.5)) +
    ggpubr::stat_compare_means(aes(group=group), label = "p.signif", method = 'wilcox.test',
                               method.args = list(alternative = "less"),
                               label.y = 0.03)+
    coord_cartesian(ylim=c(0,0.03))+
    ylab("Mutual information")+
    ggtitle(paste0("top", per, "% close/distant cell pairs in ", dataset))+
    xlab("")+
    theme(text = element_text(size = 13),
          axis.text.y=element_text(size=12, color = 'black'),
          axis.text.x=element_text(size=12, angle = 45, vjust=1, hjust=1, color = 'black'))
  p1
  
  per <- 20
  data.tmp <- data[which(data$perc==per & data$datasets == dataset),]
  p2 <- ggplot(data.tmp, aes(methods, MI, fill=group)) + 
    geom_boxplot(outlier.shape = NA, lwd = 0.5)+
    theme_bw()+
    theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1), plot.title = element_text(hjust = 0.5)) +
    ggpubr::stat_compare_means(aes(group=group), label = "p.signif", method = 'wilcox.test',
                               method.args = list(alternative = "less"),
                               label.y = 0.02)+
    coord_cartesian(ylim=c(0,0.02))+
    ylab("Mutual information")+
    ggtitle(paste0("top", per, "% close/distant cell pairs in ", dataset))+
    xlab("")+
    theme(text = element_text(size = 13),
          axis.text.y=element_text(size=12, color = 'black'),
          axis.text.x=element_text(size=12, angle = 45, vjust=1, hjust=1, color = 'black'))
  p2
  
  per <- 30
  data.tmp <- data[which(data$perc==per & data$datasets == dataset),]
  p3 <- ggplot(data.tmp, aes(methods, MI, fill=group)) + 
    geom_boxplot(outlier.shape = NA, lwd = 0.5)+
    theme_bw()+
    theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1), plot.title = element_text(hjust = 0.5)) +
    ggpubr::stat_compare_means(aes(group=group), label = "p.signif", method = 'wilcox.test',
                               method.args = list(alternative = "less"),
                               label.y = 0.015)+
    coord_cartesian(ylim=c(0,0.015))+
    ylab("Mutual information")+
    ggtitle(paste0("top", per, "% close/distant cell pairs in ", dataset))+
    xlab("")+
    theme(text = element_text(size = 13),
          axis.text.y=element_text(size=12, color = 'black'),
          axis.text.x=element_text(size=12, angle = 45, vjust=1, hjust=1, color = 'black'))
  p3
  
  per <- 40
  data.tmp <- data[which(data$perc==per & data$datasets == dataset),]
  p4 <- ggplot(data.tmp, aes(methods, MI, fill=group)) + 
    geom_boxplot(outlier.shape = NA, lwd = 0.5)+
    theme_bw()+
    theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1), plot.title = element_text(hjust = 0.5)) +
    ggpubr::stat_compare_means(aes(group=group), label = "p.signif", method = 'wilcox.test',
                               method.args = list(alternative = "less"),
                               label.y = 0.01)+
    coord_cartesian(ylim=c(0,0.01))+
    ylab("Mutual information")+
    ggtitle(paste0("top", per, "% close/distant cell pairs in ", dataset))+
    xlab("")+
    theme(text = element_text(size = 13),
          axis.text.y=element_text(size=12, color = 'black'),
          axis.text.x=element_text(size=12, angle = 45, vjust=1, hjust=1, color = 'black'))
  p4
  
  cowplot::plot_grid(p1, p2, p3, p4, nrow = 2, labels = LETTERS[1:4], label_size = 20)
  
  ###########
  ## CK358 ##
  ###########
  
  dataset <- 'CK358'
  per <- 10
  data.tmp <- data[which(data$perc==per & data$datasets == dataset),]
  p1 <- ggplot(data.tmp, aes(methods, MI, fill=group)) + 
    geom_boxplot(outlier.shape = NA, lwd = 0.5)+
    theme_bw()+
    theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1), plot.title = element_text(hjust = 0.5)) +
    ggpubr::stat_compare_means(aes(group=group), label = "p.signif", method = 'wilcox.test',
                               method.args = list(alternative = "less"),
                               label.y = 0.04)+
    coord_cartesian(ylim=c(0,0.04))+
    ylab("Mutual information")+
    ggtitle(paste0("top", per, "% close/distant cell pairs in ", dataset))+
    xlab("")+
    theme(text = element_text(size = 13),
          axis.text.y=element_text(size=12, color = 'black'),
          axis.text.x=element_text(size=12, angle = 45, vjust=1, hjust=1, color = 'black'))
  p1
  
  per <- 20
  data.tmp <- data[which(data$perc==per & data$datasets == dataset),]
  p2 <- ggplot(data.tmp, aes(methods, MI, fill=group)) + 
    geom_boxplot(outlier.shape = NA, lwd = 0.5)+
    theme_bw()+
    theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1), plot.title = element_text(hjust = 0.5)) +
    ggpubr::stat_compare_means(aes(group=group), label = "p.signif", method = 'wilcox.test',
                               method.args = list(alternative = "less"),
                               label.y = 0.03)+
    coord_cartesian(ylim=c(0,0.03))+
    ylab("Mutual information")+
    ggtitle(paste0("top", per, "% close/distant cell pairs in ", dataset))+
    xlab("")+
    theme(text = element_text(size = 13),
          axis.text.y=element_text(size=12, color = 'black'),
          axis.text.x=element_text(size=12, angle = 45, vjust=1, hjust=1, color = 'black'))
  p2
  
  per <- 30
  data.tmp <- data[which(data$perc==per & data$datasets == dataset),]
  p3 <- ggplot(data.tmp, aes(methods, MI, fill=group)) + 
    geom_boxplot(outlier.shape = NA, lwd = 0.5)+
    theme_bw()+
    theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1), plot.title = element_text(hjust = 0.5)) +
    ggpubr::stat_compare_means(aes(group=group), label = "p.signif", method = 'wilcox.test',
                               method.args = list(alternative = "less"),
                               label.y = 0.02)+
    coord_cartesian(ylim=c(0,0.02))+
    ylab("Mutual information")+
    ggtitle(paste0("top", per, "% close/distant cell pairs in ", dataset))+
    xlab("")+
    theme(text = element_text(size = 13),
          axis.text.y=element_text(size=12, color = 'black'),
          axis.text.x=element_text(size=12, angle = 45, vjust=1, hjust=1, color = 'black'))
  p3
  
  per <- 40
  data.tmp <- data[which(data$perc==per & data$datasets == dataset),]
  p4 <- ggplot(data.tmp, aes(methods, MI, fill=group)) + 
    geom_boxplot(outlier.shape = NA, lwd = 0.5)+
    theme_bw()+
    theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1), plot.title = element_text(hjust = 0.5)) +
    ggpubr::stat_compare_means(aes(group=group), label = "p.signif", method = 'wilcox.test',
                               method.args = list(alternative = "less"),
                               label.y = 0.015)+
    coord_cartesian(ylim=c(0,0.015))+
    ylab("Mutual information")+
    ggtitle(paste0("top", per, "% close/distant cell pairs in ", dataset))+
    xlab("")+
    theme(text = element_text(size = 13),
          axis.text.y=element_text(size=12, color = 'black'),
          axis.text.x=element_text(size=12, angle = 45, vjust=1, hjust=1, color = 'black'))
  p4
  
  cowplot::plot_grid(p1, p2, p3, p4, nrow = 2, labels = LETTERS[5:8], label_size = 20)
  
  ###########
  ## CK161 ##
  ###########
  
  dataset <- 'CK161'
  per <- 10
  data.tmp <- data[which(data$perc==per & data$datasets == dataset),]
  p1 <- ggplot(data.tmp, aes(methods, MI, fill=group)) + 
    geom_boxplot(outlier.shape = NA, lwd = 0.5)+
    theme_bw()+
    theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1), plot.title = element_text(hjust = 0.5)) +
    ggpubr::stat_compare_means(aes(group=group), label = "p.signif", method = 'wilcox.test',
                               method.args = list(alternative = "less"),
                               label.y = 0.0065)+
    coord_cartesian(ylim=c(0,0.007))+
    ylab("Mutual information")+
    ggtitle(paste0("top", per, "% close/distant cell pairs in ", dataset))+
    xlab("")+
    theme(text = element_text(size = 13),
          axis.text.y=element_text(size=10, color = 'black'),
          axis.text.x=element_text(size=10, angle = 45, vjust=1, hjust=1, color = 'black'))
  p1
  
  per <- 20
  data.tmp <- data[which(data$perc==per & data$datasets == dataset),]
  p2 <- ggplot(data.tmp, aes(methods, MI, fill=group)) + 
    geom_boxplot(outlier.shape = NA, lwd = 0.5)+
    theme_bw()+
    theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1), plot.title = element_text(hjust = 0.5)) +
    ggpubr::stat_compare_means(aes(group=group), label = "p.signif", method = 'wilcox.test',
                               method.args = list(alternative = "less"),
                               label.y = 0.0037)+
    coord_cartesian(ylim=c(0,0.004))+
    ylab("Mutual information")+
    ggtitle(paste0("top", per, "% close/distant cell pairs in ", dataset))+
    xlab("")+
    theme(text = element_text(size = 13),
          axis.text.y=element_text(size=10, color = 'black'),
          axis.text.x=element_text(size=10, angle = 45, vjust=1, hjust=1, color = 'black'))
  p2
  
  per <- 30
  data.tmp <- data[which(data$perc==per & data$datasets == dataset),]
  p3 <- ggplot(data.tmp, aes(methods, MI, fill=group)) + 
    geom_boxplot(outlier.shape = NA, lwd = 0.5)+
    theme_bw()+
    theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1), plot.title = element_text(hjust = 0.5)) +
    ggpubr::stat_compare_means(aes(group=group), label = "p.signif", method = 'wilcox.test',
                               method.args = list(alternative = "less"),
                               label.y = 0.0018)+
    coord_cartesian(ylim=c(0,0.002))+
    ylab("Mutual information")+
    ggtitle(paste0("top", per, "% close/distant cell pairs in ", dataset))+
    xlab("")+
    theme(text = element_text(size = 13),
          axis.text.y=element_text(size=10, color = 'black'),
          axis.text.x=element_text(size=10, angle = 45, vjust=1, hjust=1, color = 'black'))
  p3
  
  per <- 40
  data.tmp <- data[which(data$perc==per & data$datasets == dataset),]
  p4 <- ggplot(data.tmp, aes(methods, MI, fill=group)) + 
    geom_boxplot(outlier.shape = NA, lwd = 0.5)+
    theme_bw()+
    theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1), plot.title = element_text(hjust = 0.5)) +
    ggpubr::stat_compare_means(aes(group=group), label = "p.signif", method = 'wilcox.test',
                               method.args = list(alternative = "less"),
                               label.y = 0.0014)+
    coord_cartesian(ylim=c(0,0.0015))+
    ylab("Mutual information")+
    ggtitle(paste0("top", per, "% close/distant cell pairs in ", dataset))+
    xlab("")+
    theme(text = element_text(size = 13),
          axis.text.y=element_text(size=10, color = 'black'),
          axis.text.x=element_text(size=10, angle = 45, vjust=1, hjust=1, color = 'black'))
  p4
  
  cowplot::plot_grid(p1, p2, p3, p4, nrow = 2, labels = LETTERS[1:4], label_size = 20)
  
  ###########
  ## CK165 ##
  ###########
  
  dataset <- 'CK165'
  per <- 10
  data.tmp <- data[which(data$perc==per & data$datasets == dataset),]
  p1 <- ggplot(data.tmp, aes(methods, MI, fill=group)) + 
    geom_boxplot(outlier.shape = NA, lwd = 0.5)+
    theme_bw()+
    theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1), plot.title = element_text(hjust = 0.5)) +
    ggpubr::stat_compare_means(aes(group=group), label = "p.signif", method = 'wilcox.test',
                               method.args = list(alternative = "less"),
                               label.y = 0.0115)+
    coord_cartesian(ylim=c(0, 0.0125))+
    ylab("Mutual information")+
    ggtitle(paste0("top", per, "% close/distant cell pairs in ", dataset))+
    xlab("")+
    theme(text = element_text(size = 13),
          axis.text.y=element_text(size=10, color = 'black'),
          axis.text.x=element_text(size=10, angle = 45, vjust=1, hjust=1, color = 'black'))
  p1
  
  per <- 20
  data.tmp <- data[which(data$perc==per & data$datasets == dataset),]
  p2 <- ggplot(data.tmp, aes(methods, MI, fill=group)) + 
    geom_boxplot(outlier.shape = NA, lwd = 0.5)+
    theme_bw()+
    theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1), plot.title = element_text(hjust = 0.5)) +
    ggpubr::stat_compare_means(aes(group=group), label = "p.signif", method = 'wilcox.test',
                               method.args = list(alternative = "less"),
                               label.y = 0.0066)+
    coord_cartesian(ylim=c(0,0.007))+
    ylab("Mutual information")+
    ggtitle(paste0("top", per, "% close/distant cell pairs in ", dataset))+
    xlab("")+
    theme(text = element_text(size = 13),
          axis.text.y=element_text(size=10, color = 'black'),
          axis.text.x=element_text(size=10, angle = 45, vjust=1, hjust=1, color = 'black'))
  p2
  
  per <- 30
  data.tmp <- data[which(data$perc==per & data$datasets == dataset),]
  p3 <- ggplot(data.tmp, aes(methods, MI, fill=group)) + 
    geom_boxplot(outlier.shape = NA, lwd = 0.5)+
    theme_bw()+
    theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1), plot.title = element_text(hjust = 0.5)) +
    ggpubr::stat_compare_means(aes(group=group), label = "p.signif", method = 'wilcox.test',
                               method.args = list(alternative = "less"),
                               label.y = 0.0038)+
    coord_cartesian(ylim=c(0,0.004))+
    ylab("Mutual information")+
    ggtitle(paste0("top", per, "% close/distant cell pairs in ", dataset))+
    xlab("")+
    theme(text = element_text(size = 13),
          axis.text.y=element_text(size=10, color = 'black'),
          axis.text.x=element_text(size=10, angle = 45, vjust=1, hjust=1, color = 'black'))
  p3
  
  per <- 40
  data.tmp <- data[which(data$perc==per & data$datasets == dataset),]
  p4 <- ggplot(data.tmp, aes(methods, MI, fill=group)) + 
    geom_boxplot(outlier.shape = NA, lwd = 0.5)+
    theme_bw()+
    theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1), plot.title = element_text(hjust = 0.5)) +
    ggpubr::stat_compare_means(aes(group=group), label = "p.signif", method = 'wilcox.test',
                               method.args = list(alternative = "less"),
                               label.y = 0.0028)+
    coord_cartesian(ylim=c(0,0.003))+
    ylab("Mutual information")+
    ggtitle(paste0("top", per, "% close/distant cell pairs in ", dataset))+
    xlab("")+
    theme(text = element_text(size = 13),
          axis.text.y=element_text(size=10, color = 'black'),
          axis.text.x=element_text(size=10, angle = 45, vjust=1, hjust=1, color = 'black'))
  p4
  
  cowplot::plot_grid(p1, p2, p3, p4, nrow = 2, labels = LETTERS[5:8], label_size = 20)
  
  
  ###########
  ## CK361 ##
  ###########
  
  dataset <- 'CK361'
  per <- 10
  data.tmp <- data[which(data$perc==per & data$datasets == dataset),]
  p1 <- ggplot(data.tmp, aes(methods, MI, fill=group)) + 
    geom_boxplot(outlier.shape = NA, lwd = 0.5)+
    theme_bw()+
    theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1), plot.title = element_text(hjust = 0.5)) +
    ggpubr::stat_compare_means(aes(group=group), label = "p.signif", method = 'wilcox.test',
                               method.args = list(alternative = "less"),
                               label.y = 0.024)+
    coord_cartesian(ylim=c(0, 0.025))+
    ylab("Mutual information")+
    ggtitle(paste0("top", per, "% close/distant cell pairs in ", dataset))+
    xlab("")+
    theme(text = element_text(size = 13),
          axis.text.y=element_text(size=10, color = 'black'),
          axis.text.x=element_text(size=10, angle = 45, vjust=1, hjust=1, color = 'black'))
  p1
  
  per <- 20
  data.tmp <- data[which(data$perc==per & data$datasets == dataset),]
  p2 <- ggplot(data.tmp, aes(methods, MI, fill=group)) + 
    geom_boxplot(outlier.shape = NA, lwd = 0.5)+
    theme_bw()+
    theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1), plot.title = element_text(hjust = 0.5)) +
    ggpubr::stat_compare_means(aes(group=group), label = "p.signif", method = 'wilcox.test',
                               method.args = list(alternative = "less"),
                               label.y = 0.0150)+
    coord_cartesian(ylim=c(0,0.0155))+
    ylab("Mutual information")+
    ggtitle(paste0("top", per, "% close/distant cell pairs in ", dataset))+
    xlab("")+
    theme(text = element_text(size = 13),
          axis.text.y=element_text(size=10, color = 'black'),
          axis.text.x=element_text(size=10, angle = 45, vjust=1, hjust=1, color = 'black'))
  p2
  
  per <- 30
  data.tmp <- data[which(data$perc==per & data$datasets == dataset),]
  p3 <- ggplot(data.tmp, aes(methods, MI, fill=group)) + 
    geom_boxplot(outlier.shape = NA, lwd = 0.5)+
    theme_bw()+
    theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1), plot.title = element_text(hjust = 0.5)) +
    ggpubr::stat_compare_means(aes(group=group), label = "p.signif", method = 'wilcox.test',
                               method.args = list(alternative = "less"),
                               label.y = 0.0095)+
    coord_cartesian(ylim=c(0,0.01))+
    ylab("Mutual information")+
    ggtitle(paste0("top", per, "% close/distant cell pairs in ", dataset))+
    xlab("")+
    theme(text = element_text(size = 13),
          axis.text.y=element_text(size=10, color = 'black'),
          axis.text.x=element_text(size=10, angle = 45, vjust=1, hjust=1, color = 'black'))
  p3
  
  per <- 40
  data.tmp <- data[which(data$perc==per & data$datasets == dataset),]
  p4 <- ggplot(data.tmp, aes(methods, MI, fill=group)) + 
    geom_boxplot(outlier.shape = NA, lwd = 0.5)+
    theme_bw()+
    theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1), plot.title = element_text(hjust = 0.5)) +
    ggpubr::stat_compare_means(aes(group=group), label = "p.signif", method = 'wilcox.test',
                               method.args = list(alternative = "less"),
                               label.y = 0.0065)+
    coord_cartesian(ylim=c(0,0.007))+
    ylab("Mutual information")+
    ggtitle(paste0("top", per, "% close/distant cell pairs in ", dataset))+
    xlab("")+
    theme(text = element_text(size = 13),
          axis.text.y=element_text(size=10, color = 'black'),
          axis.text.x=element_text(size=10, angle = 45, vjust=1, hjust=1, color = 'black'))
  p4
  
  cowplot::plot_grid(p1, p2, p3, p4, nrow = 2, labels = LETTERS[9:12], label_size = 20)
  
  ###########
  ## CK362 ##
  ###########
  
  dataset <- 'CK362'
  per <- 10
  data.tmp <- data[which(data$perc==per & data$datasets == dataset),]
  p1 <- ggplot(data.tmp, aes(methods, MI, fill=group)) + 
    geom_boxplot(outlier.shape = NA, lwd = 0.5)+
    theme_bw()+
    theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1), plot.title = element_text(hjust = 0.5)) +
    ggpubr::stat_compare_means(aes(group=group), label = "p.signif", method = 'wilcox.test',
                               method.args = list(alternative = "less"),
                               label.y = 0.024)+
    coord_cartesian(ylim=c(0, 0.025))+
    ylab("Mutual information")+
    ggtitle(paste0("top", per, "% close/distant cell pairs in ", dataset))+
    xlab("")+
    theme(text = element_text(size = 13),
          axis.text.y=element_text(size=10, color = 'black'),
          axis.text.x=element_text(size=10, angle = 45, vjust=1, hjust=1, color = 'black'))
  p1
  
  per <- 20
  data.tmp <- data[which(data$perc==per & data$datasets == dataset),]
  p2 <- ggplot(data.tmp, aes(methods, MI, fill=group)) + 
    geom_boxplot(outlier.shape = NA, lwd = 0.5)+
    theme_bw()+
    theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1), plot.title = element_text(hjust = 0.5)) +
    ggpubr::stat_compare_means(aes(group=group), label = "p.signif", method = 'wilcox.test',
                               method.args = list(alternative = "less"),
                               label.y = 0.0150)+
    coord_cartesian(ylim=c(0,0.0155))+
    ylab("Mutual information")+
    ggtitle(paste0("top", per, "% close/distant cell pairs in ", dataset))+
    xlab("")+
    theme(text = element_text(size = 13),
          axis.text.y=element_text(size=10, color = 'black'),
          axis.text.x=element_text(size=10, angle = 45, vjust=1, hjust=1, color = 'black'))
  p2
  
  per <- 30
  data.tmp <- data[which(data$perc==per & data$datasets == dataset),]
  p3 <- ggplot(data.tmp, aes(methods, MI, fill=group)) + 
    geom_boxplot(outlier.shape = NA, lwd = 0.5)+
    theme_bw()+
    theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1), plot.title = element_text(hjust = 0.5)) +
    ggpubr::stat_compare_means(aes(group=group), label = "p.signif", method = 'wilcox.test',
                               method.args = list(alternative = "less"),
                               label.y = 0.0099)+
    coord_cartesian(ylim=c(0,0.01))+
    ylab("Mutual information")+
    ggtitle(paste0("top", per, "% close/distant cell pairs in ", dataset))+
    xlab("")+
    theme(text = element_text(size = 13),
          axis.text.y=element_text(size=10, color = 'black'),
          axis.text.x=element_text(size=10, angle = 45, vjust=1, hjust=1, color = 'black'))
  p3
  
  per <- 40
  data.tmp <- data[which(data$perc==per & data$datasets == dataset),]
  p4 <- ggplot(data.tmp, aes(methods, MI, fill=group)) + 
    geom_boxplot(outlier.shape = NA, lwd = 0.5)+
    theme_bw()+
    theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1), plot.title = element_text(hjust = 0.5)) +
    ggpubr::stat_compare_means(aes(group=group), label = "p.signif", method = 'wilcox.test',
                               method.args = list(alternative = "less"),
                               label.y = 0.0068)+
    coord_cartesian(ylim=c(0,0.007))+
    ylab("Mutual information")+
    ggtitle(paste0("top", per, "% close/distant cell pairs in ", dataset))+
    xlab("")+
    theme(text = element_text(size = 13),
          axis.text.y=element_text(size=10, color = 'black'),
          axis.text.x=element_text(size=10, angle = 45, vjust=1, hjust=1, color = 'black'))
  p4
  
  cowplot::plot_grid(p1, p2, p3, p4, nrow = 2, labels = LETTERS[1:4], label_size = 20)
  
  ###########
  ## CK162 ##
  ###########
  
  dataset <- 'CK162'
  per <- 10
  data.tmp <- data[which(data$perc==per & data$datasets == dataset),]
  p1 <- ggplot(data.tmp, aes(methods, MI, fill=group)) + 
    geom_boxplot(outlier.shape = NA, lwd = 0.5)+
    theme_bw()+
    theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1), plot.title = element_text(hjust = 0.5)) +
    ggpubr::stat_compare_means(aes(group=group), label = "p.signif", method = 'wilcox.test',
                               method.args = list(alternative = "less"),
                               label.y = 0.019)+
    coord_cartesian(ylim=c(0, 0.02))+
    ylab("Mutual information")+
    ggtitle(paste0("top", per, "% close/distant cell pairs in ", dataset))+
    xlab("")+
    theme(text = element_text(size = 13),
          axis.text.y=element_text(size=10, color = 'black'),
          axis.text.x=element_text(size=10, angle = 45, vjust=1, hjust=1, color = 'black'))
  p1
  
  per <- 20
  data.tmp <- data[which(data$perc==per & data$datasets == dataset),]
  p2 <- ggplot(data.tmp, aes(methods, MI, fill=group)) + 
    geom_boxplot(outlier.shape = NA, lwd = 0.5)+
    theme_bw()+
    theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1), plot.title = element_text(hjust = 0.5)) +
    ggpubr::stat_compare_means(aes(group=group), label = "p.signif", method = 'wilcox.test',
                               method.args = list(alternative = "less"),
                               label.y = 0.015)+
    coord_cartesian(ylim=c(0,0.0155))+
    ylab("Mutual information")+
    ggtitle(paste0("top", per, "% close/distant cell pairs in ", dataset))+
    xlab("")+
    theme(text = element_text(size = 13),
          axis.text.y=element_text(size=10, color = 'black'),
          axis.text.x=element_text(size=10, angle = 45, vjust=1, hjust=1, color = 'black'))
  p2
  
  per <- 30
  data.tmp <- data[which(data$perc==per & data$datasets == dataset),]
  p3 <- ggplot(data.tmp, aes(methods, MI, fill=group)) + 
    geom_boxplot(outlier.shape = NA, lwd = 0.5)+
    theme_bw()+
    theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1), plot.title = element_text(hjust = 0.5)) +
    ggpubr::stat_compare_means(aes(group=group), label = "p.signif", method = 'wilcox.test',
                               method.args = list(alternative = "less"),
                               label.y = 0.0095)+
    coord_cartesian(ylim=c(0,0.01))+
    ylab("Mutual information")+
    ggtitle(paste0("top", per, "% close/distant cell pairs in ", dataset))+
    xlab("")+
    theme(text = element_text(size = 13),
          axis.text.y=element_text(size=10, color = 'black'),
          axis.text.x=element_text(size=10, angle = 45, vjust=1, hjust=1, color = 'black'))
  p3
  
  per <- 40
  data.tmp <- data[which(data$perc==per & data$datasets == dataset),]
  p4 <- ggplot(data.tmp, aes(methods, MI, fill=group)) + 
    geom_boxplot(outlier.shape = NA, lwd = 0.5)+
    theme_bw()+
    theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1), plot.title = element_text(hjust = 0.5)) +
    ggpubr::stat_compare_means(aes(group=group), label = "p.signif", method = 'wilcox.test',
                               method.args = list(alternative = "less"),
                               label.y = 0.0068)+
    coord_cartesian(ylim=c(0,0.007))+
    ylab("Mutual information")+
    ggtitle(paste0("top", per, "% close/distant cell pairs in ", dataset))+
    xlab("")+
    theme(text = element_text(size = 13),
          axis.text.y=element_text(size=10, color = 'black'),
          axis.text.x=element_text(size=10, angle = 45, vjust=1, hjust=1, color = 'black'))
  p4
  
  cowplot::plot_grid(p1, p2, p3, p4, nrow = 2, labels = LETTERS[5:9], label_size = 20)
  
  ###########
  ## CK368 ##
  ###########
  
  dataset <- 'CK368'
  per <- 10
  data.tmp <- data[which(data$perc==per & data$datasets == dataset),]
  p1 <- ggplot(data.tmp, aes(methods, MI, fill=group)) + 
    geom_boxplot(outlier.shape = NA, lwd = 0.5)+
    theme_bw()+
    theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1), plot.title = element_text(hjust = 0.5)) +
    ggpubr::stat_compare_means(aes(group=group), label = "p.signif", method = 'wilcox.test',
                               method.args = list(alternative = "less"),
                               label.y = 0.014)+
    coord_cartesian(ylim=c(0, 0.015))+
    ylab("Mutual information")+
    ggtitle(paste0("top", per, "% close/distant cell pairs in ", dataset))+
    xlab("")+
    theme(text = element_text(size = 13),
          axis.text.y=element_text(size=10, color = 'black'),
          axis.text.x=element_text(size=10, angle = 45, vjust=1, hjust=1, color = 'black'))
  p1
  
  per <- 20
  data.tmp <- data[which(data$perc==per & data$datasets == dataset),]
  p2 <- ggplot(data.tmp, aes(methods, MI, fill=group)) + 
    geom_boxplot(outlier.shape = NA, lwd = 0.5)+
    theme_bw()+
    theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1), plot.title = element_text(hjust = 0.5)) +
    ggpubr::stat_compare_means(aes(group=group), label = "p.signif", method = 'wilcox.test',
                               method.args = list(alternative = "less"),
                               label.y = 0.0095)+
    coord_cartesian(ylim=c(0,0.01))+
    ylab("Mutual information")+
    ggtitle(paste0("top", per, "% close/distant cell pairs in ", dataset))+
    xlab("")+
    theme(text = element_text(size = 13),
          axis.text.y=element_text(size=10, color = 'black'),
          axis.text.x=element_text(size=10, angle = 45, vjust=1, hjust=1, color = 'black'))
  p2
  
  per <- 30
  data.tmp <- data[which(data$perc==per & data$datasets == dataset),]
  p3 <- ggplot(data.tmp, aes(methods, MI, fill=group)) + 
    geom_boxplot(outlier.shape = NA, lwd = 0.5)+
    theme_bw()+
    theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1), plot.title = element_text(hjust = 0.5)) +
    ggpubr::stat_compare_means(aes(group=group), label = "p.signif", method = 'wilcox.test',
                               method.args = list(alternative = "less"),
                               label.y = 0.0070)+
    coord_cartesian(ylim=c(0,0.0075))+
    ylab("Mutual information")+
    ggtitle(paste0("top", per, "% close/distant cell pairs in ", dataset))+
    xlab("")+
    theme(text = element_text(size = 13),
          axis.text.y=element_text(size=10, color = 'black'),
          axis.text.x=element_text(size=10, angle = 45, vjust=1, hjust=1, color = 'black'))
  p3
  
  per <- 40
  data.tmp <- data[which(data$perc==per & data$datasets == dataset),]
  p4 <- ggplot(data.tmp, aes(methods, MI, fill=group)) + 
    geom_boxplot(outlier.shape = NA, lwd = 0.5)+
    theme_bw()+
    theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1), plot.title = element_text(hjust = 0.5)) +
    ggpubr::stat_compare_means(aes(group=group), label = "p.signif", method = 'wilcox.test',
                               method.args = list(alternative = "less"),
                               label.y = 0.0048)+
    coord_cartesian(ylim=c(0,0.005))+
    ylab("Mutual information")+
    ggtitle(paste0("top", per, "% close/distant cell pairs in ", dataset))+
    xlab("")+
    theme(text = element_text(size = 13),
          axis.text.y=element_text(size=10, color = 'black'),
          axis.text.x=element_text(size=10, angle = 45, vjust=1, hjust=1, color = 'black'))
  p4
  
  cowplot::plot_grid(p1, p2, p3, p4, nrow = 2, labels = LETTERS[10:13], label_size = 20)
  
  #############
  ## Slide14 ##
  #############
  
  dataset <- 'Slide14'
  per <- 10
  data.tmp <- data[which(data$perc==per & data$datasets == dataset),]
  p1 <- ggplot(data.tmp, aes(methods, MI, fill=group)) + 
    geom_boxplot(outlier.shape = NA, lwd = 0.5)+
    theme_bw()+
    theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1), plot.title = element_text(hjust = 0.5)) +
    ggpubr::stat_compare_means(aes(group=group), label = "p.signif", method = 'wilcox.test',
                               method.args = list(alternative = "less"),
                               label.y = 0.0024)+
    coord_cartesian(ylim=c(0, 0.0025))+
    ylab("Mutual information")+
    ggtitle(paste0("top", per, "% close/distant cell pairs in ", dataset))+
    xlab("")+
    theme(text = element_text(size = 13),
          axis.text.y=element_text(size=10, color = 'black'),
          axis.text.x=element_text(size=10, angle = 45, vjust=1, hjust=1, color = 'black'))
  p1
  
  per <- 20
  data.tmp <- data[which(data$perc==per & data$datasets == dataset),]
  p2 <- ggplot(data.tmp, aes(methods, MI, fill=group)) + 
    geom_boxplot(outlier.shape = NA, lwd = 0.5)+
    theme_bw()+
    theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1), plot.title = element_text(hjust = 0.5)) +
    ggpubr::stat_compare_means(aes(group=group), label = "p.signif", method = 'wilcox.test',
                               method.args = list(alternative = "less"),
                               label.y = 0.0014)+
    coord_cartesian(ylim=c(0,0.0015))+
    ylab("Mutual information")+
    ggtitle(paste0("top", per, "% close/distant cell pairs in ", dataset))+
    xlab("")+
    theme(text = element_text(size = 13),
          axis.text.y=element_text(size=10, color = 'black'),
          axis.text.x=element_text(size=10, angle = 45, vjust=1, hjust=1, color = 'black'))
  p2
  
  per <- 30
  data.tmp <- data[which(data$perc==per & data$datasets == dataset),]
  p3 <- ggplot(data.tmp, aes(methods, MI, fill=group)) + 
    geom_boxplot(outlier.shape = NA, lwd = 0.5)+
    theme_bw()+
    theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1), plot.title = element_text(hjust = 0.5)) +
    ggpubr::stat_compare_means(aes(group=group), label = "p.signif", method = 'wilcox.test',
                               method.args = list(alternative = "less"),
                               label.y = 0.00070)+
    coord_cartesian(ylim=c(0,0.00075))+
    ylab("Mutual information")+
    ggtitle(paste0("top", per, "% close/distant cell pairs in ", dataset))+
    xlab("")+
    theme(text = element_text(size = 13),
          axis.text.y=element_text(size=10, color = 'black'),
          axis.text.x=element_text(size=10, angle = 45, vjust=1, hjust=1, color = 'black'))
  p3
  
  per <- 40
  data.tmp <- data[which(data$perc==per & data$datasets == dataset),]
  p4 <- ggplot(data.tmp, aes(methods, MI, fill=group)) + 
    geom_boxplot(outlier.shape = NA, lwd = 0.5)+
    theme_bw()+
    theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1), plot.title = element_text(hjust = 0.5)) +
    ggpubr::stat_compare_means(aes(group=group), label = "p.signif", method = 'wilcox.test',
                               method.args = list(alternative = "less"),
                               label.y = 0.00048)+
    coord_cartesian(ylim=c(0,0.0005))+
    ylab("Mutual information")+
    ggtitle(paste0("top", per, "% close/distant cell pairs in ", dataset))+
    xlab("")+
    theme(text = element_text(size = 13),
          axis.text.y=element_text(size=10, color = 'black'),
          axis.text.x=element_text(size=10, angle = 45, vjust=1, hjust=1, color = 'black'))
  p4
  
  cowplot::plot_grid(p1, p2, p3, p4, nrow = 2, labels = LETTERS[1:4], label_size = 20)
}



############################
## Supplementary Figure 7 ##
############################
if(F){
  setwd('~/2-CCC-benchmark-202204-202208/CCC-benchmark-202212/')
  result <- readRDS('./Data/Step10_Visualization/cage_eval_updata.rds')
  result_4k <- result$pbmc4k
  result_4k$methods <- factor(result_4k$methods, 
                              levels = c("CellPhoneDB2", "CellPhoneDB3", "CellTalker", 
                                         "Connectome", "ICELLNET",  "NATMI", "iTALK", 
                                         "scConnect", "SingleCellSignalR", "CellChat", 
                                         "RNAMagnet", "PyMINEr",
                                         "scSeqComm", "NicheNet","CytoTalk",
                                         "scMLnet", "Domino", "CellCall"))
  p1 <- ggplot(result_4k, aes(methods, AUPRC, fill=methods)) + 
    geom_boxplot(lwd = 0.5, position=position_dodge(width=0.8))+ #outlier.shape = NA, 
    theme_bw()+
    ylab("AUPRC")+
    xlab("")+
    theme(text = element_text(size = 15),
          plot.title = element_text(size = 16, hjust = 0.5, face = 'bold'),
          axis.text.y=element_text(size=13, color = 'black'),
          axis.text.x=element_text(size=13, angle = 45, vjust=1, hjust=1, color = 'black'),
          legend.position="none")+
    ggtitle('PBMC 4K')
  p1
  
  p12 <- ggplot(result_4k, aes(methods, AUROC, fill=methods)) + 
    geom_boxplot(lwd = 0.5, position=position_dodge(width=0.8))+ #outlier.shape = NA, 
    theme_bw()+
    ylab("AUROC")+
    xlab("")+
    theme(text = element_text(size = 15),
          plot.title = element_text(size = 16, hjust = 0.5, face = 'bold'),
          axis.text.y=element_text(size=13, color = 'black'),
          axis.text.x=element_text(size=13, angle = 45, vjust=1, hjust=1, color = 'black'),
          legend.position="none")+
    ggtitle('PBMC 4K')
  p12
  
  
  result_6k <- result$pbmc6k
  result_6k$methods <- factor(result_6k$methods, 
                              levels = c("CellPhoneDB2", "CellPhoneDB3", "CellTalker", 
                                         "Connectome", "ICELLNET",  "NATMI", "iTALK", 
                                         "scConnect", "SingleCellSignalR", "CellChat", 
                                         "RNAMagnet", "PyMINEr",
                                         "scSeqComm", "NicheNet","CytoTalk",
                                         "scMLnet", "Domino", "CellCall"))
  p2 <- ggplot(result_6k, aes(methods, AUPRC, fill=methods)) + 
    geom_boxplot(lwd = 0.5, position=position_dodge(width=0.8))+ #outlier.shape = NA, 
    theme_bw()+
    ylab("AUPRC")+
    xlab("")+
    theme(text = element_text(size = 15),
          plot.title = element_text(size = 16, hjust = 0.5, face = 'bold'),
          axis.text.y=element_text(size=13, color = 'black'),
          axis.text.x=element_text(size=13, angle = 45, vjust=1, hjust=1, color = 'black'),
          legend.position="none")+
    ggtitle('PBMC 6K')
  p22 <- ggplot(result_6k, aes(methods, AUROC, fill=methods)) + 
    geom_boxplot(lwd = 0.5, position=position_dodge(width=0.8))+ #outlier.shape = NA, 
    theme_bw()+
    ylab("AUROC")+
    xlab("")+
    theme(text = element_text(size = 15),
          plot.title = element_text(size = 16, hjust = 0.5, face = 'bold'),
          axis.text.y=element_text(size=13, color = 'black'),
          axis.text.x=element_text(size=13, angle = 45, vjust=1, hjust=1, color = 'black'),
          legend.position="none")+
    ggtitle('PBMC 6K')
  
  
  result_8k <- result$pbmc8k
  result_8k$methods <- factor(result_8k$methods, 
                              levels = c("CellPhoneDB2", "CellPhoneDB3", "CellTalker", 
                                         "Connectome", "ICELLNET",  "NATMI", "iTALK", 
                                         "scConnect", "SingleCellSignalR", "CellChat", 
                                         "RNAMagnet", "PyMINEr",
                                         "scSeqComm", "NicheNet","CytoTalk",
                                         "scMLnet", "Domino", "CellCall"))
  p3 <- ggplot(result_8k, aes(methods, AUPRC, fill=methods)) + 
    geom_boxplot(lwd = 0.5, position=position_dodge(width=0.8))+ #outlier.shape = NA, 
    theme_bw()+
    ylab("AUPRC")+
    xlab("")+
    theme(text = element_text(size = 15),
          plot.title = element_text(size = 16, hjust = 0.5, face = 'bold'),
          axis.text.y=element_text(size=13, color = 'black'),
          axis.text.x=element_text(size=13, angle = 45, vjust=1, hjust=1, color = 'black'),
          legend.position="none")+
    ggtitle('PBMC 8K')
  p32 <- ggplot(result_8k, aes(methods, AUROC, fill=methods)) + 
    geom_boxplot(lwd = 0.5, position=position_dodge(width=0.8))+ #outlier.shape = NA, 
    theme_bw()+
    ylab("AUROC")+
    xlab("")+
    theme(text = element_text(size = 15),
          plot.title = element_text(size = 16, hjust = 0.5, face = 'bold'),
          axis.text.y=element_text(size=13, color = 'black'),
          axis.text.x=element_text(size=13, angle = 45, vjust=1, hjust=1, color = 'black'),
          legend.position="none")+
    ggtitle('PBMC 8K')
  
  
  result <- rbind(result_4k, result_6k, result_8k)
  p4 <- ggplot(result, aes(methods, AUPRC, fill=methods)) + 
    geom_boxplot(lwd = 0.5, position=position_dodge(width=0.8))+ #outlier.shape = NA, 
    theme_bw()+
    ylab("AUPRC")+
    xlab("")+
    theme(text = element_text(size = 15),
          plot.title = element_text(size = 16, hjust = 0.5, face = 'bold'),
          axis.text.y=element_text(size=13, color = 'black'),
          axis.text.x=element_text(size=13, angle = 45, vjust=1, hjust=1, color = 'black'),
          legend.position="none")+
    ggtitle('Total')
  p42 <- ggplot(result, aes(methods, AUROC, fill=methods)) + 
    geom_boxplot(lwd = 0.5, position=position_dodge(width=0.8))+ #outlier.shape = NA, 
    theme_bw()+
    ylab("AUROC")+
    xlab("")+
    theme(text = element_text(size = 15),
          plot.title = element_text(size = 16, hjust = 0.5, face = 'bold'),
          axis.text.y=element_text(size=13, color = 'black'),
          axis.text.x=element_text(size=13, angle = 45, vjust=1, hjust=1, color = 'black'),
          legend.position="none")+
    ggtitle('Total')
  
  cowplot::plot_grid(p1, p12,  labels = LETTERS[1:2], label_size = 20)
  cowplot::plot_grid(p2, p22, labels = LETTERS[3:4], label_size = 20)
  cowplot::plot_grid(p3, p32, labels = LETTERS[5:6], label_size = 20)
  
  cowplot::plot_grid(p4, p42, labels = LETTERS[7:8], label_size = 20)
  
}


############################
## Supplementary Figure 8 ##
############################

if(F){
  setwd('~/2-CCC-benchmark-202204-202208/CCC-benchmark-202212/')
  result <- readRDS('./Data/Step10_Visualization/proteomic_eval_updata.rds')
  result_4k <- result[[1]]
  result_4k$methods <- factor(result_4k$methods, 
                              levels = c("CellPhoneDB2", "CellPhoneDB3", "CellTalker", 
                                         "Connectome", "ICELLNET",  "NATMI", "iTALK", 
                                         "scConnect", "SingleCellSignalR", "CellChat", 
                                         "RNAMagnet", "PyMINEr",
                                         "scSeqComm", "NicheNet","CytoTalk",
                                         "scMLnet", "Domino", "CellCall"))
  p1 <- ggplot(result_4k, aes(methods, AUPRC, fill=methods)) + 
    geom_boxplot(lwd = 0.5, position=position_dodge(width=0.8))+ #outlier.shape = NA, 
    theme_bw()+
    ylab("AUPRC")+
    xlab("")+
    theme(text = element_text(size = 15),
          plot.title = element_text(size = 16, hjust = 0.5, face = 'bold'),
          axis.text.y=element_text(size=13, color = 'black'),
          axis.text.x=element_text(size=13, angle = 45, vjust=1, hjust=1, color = 'black'),
          legend.position="none")+
    ggtitle('PBMC 4K')
  p1
  
  p12 <- ggplot(result_4k, aes(methods, AUROC, fill=methods)) + 
    geom_boxplot(lwd = 0.5, position=position_dodge(width=0.8))+ #outlier.shape = NA, 
    theme_bw()+
    ylab("AUROC")+
    xlab("")+
    theme(text = element_text(size = 15),
          plot.title = element_text(size = 16, hjust = 0.5, face = 'bold'),
          axis.text.y=element_text(size=13, color = 'black'),
          axis.text.x=element_text(size=13, angle = 45, vjust=1, hjust=1, color = 'black'),
          legend.position="none")+
    ggtitle('PBMC 4K')
  p12
  
  
  result_6k <- result[[2]]
  result_6k$methods <- factor(result_6k$methods, 
                              levels = c("CellPhoneDB2", "CellPhoneDB3", "CellTalker", 
                                         "Connectome", "ICELLNET",  "NATMI", "iTALK", 
                                         "scConnect", "SingleCellSignalR", "CellChat", 
                                         "RNAMagnet", "PyMINEr",
                                         "scSeqComm", "NicheNet","CytoTalk",
                                         "scMLnet", "Domino", "CellCall"))
  p2 <- ggplot(result_6k, aes(methods, AUPRC, fill=methods)) + 
    geom_boxplot(lwd = 0.5, position=position_dodge(width=0.8))+ #outlier.shape = NA, 
    theme_bw()+
    ylab("AUPRC")+
    xlab("")+
    theme(text = element_text(size = 15),
          plot.title = element_text(size = 16, hjust = 0.5, face = 'bold'),
          axis.text.y=element_text(size=13, color = 'black'),
          axis.text.x=element_text(size=13, angle = 45, vjust=1, hjust=1, color = 'black'),
          legend.position="none")+
    ggtitle('PBMC 6K')
  p22 <- ggplot(result_6k, aes(methods, AUROC, fill=methods)) + 
    geom_boxplot(lwd = 0.5, position=position_dodge(width=0.8))+ #outlier.shape = NA, 
    theme_bw()+
    ylab("AUROC")+
    xlab("")+
    theme(text = element_text(size = 15),
          plot.title = element_text(size = 16, hjust = 0.5, face = 'bold'),
          axis.text.y=element_text(size=13, color = 'black'),
          axis.text.x=element_text(size=13, angle = 45, vjust=1, hjust=1, color = 'black'),
          legend.position="none")+
    ggtitle('PBMC 6K')
  
  
  result_8k <- result[[3]]
  result_8k$methods <- factor(result_8k$methods, 
                              levels = c("CellPhoneDB2", "CellPhoneDB3", "CellTalker", 
                                         "Connectome", "ICELLNET",  "NATMI", "iTALK", 
                                         "scConnect", "SingleCellSignalR", "CellChat", 
                                         "RNAMagnet", "PyMINEr",
                                         "scSeqComm", "NicheNet","CytoTalk",
                                         "scMLnet", "Domino", "CellCall"))
  p3 <- ggplot(result_8k, aes(methods, AUPRC, fill=methods)) + 
    geom_boxplot(lwd = 0.5, position=position_dodge(width=0.8))+ #outlier.shape = NA, 
    theme_bw()+
    ylab("AUPRC")+
    xlab("")+
    theme(text = element_text(size = 15),
          plot.title = element_text(size = 16, hjust = 0.5, face = 'bold'),
          axis.text.y=element_text(size=13, color = 'black'),
          axis.text.x=element_text(size=13, angle = 45, vjust=1, hjust=1, color = 'black'),
          legend.position="none")+
    ggtitle('PBMC 8K')
  p32 <- ggplot(result_8k, aes(methods, AUROC, fill=methods)) + 
    geom_boxplot(lwd = 0.5, position=position_dodge(width=0.8))+ #outlier.shape = NA, 
    theme_bw()+
    ylab("AUROC")+
    xlab("")+
    theme(text = element_text(size = 15),
          plot.title = element_text(size = 16, hjust = 0.5, face = 'bold'),
          axis.text.y=element_text(size=13, color = 'black'),
          axis.text.x=element_text(size=13, angle = 45, vjust=1, hjust=1, color = 'black'),
          legend.position="none")+
    ggtitle('PBMC 8K')
  
  
  result <- rbind(result_4k, result_6k, result_8k)
  p4 <- ggplot(result, aes(methods, AUPRC, fill=methods)) + 
    geom_boxplot(lwd = 0.5, position=position_dodge(width=0.8))+ #outlier.shape = NA, 
    theme_bw()+
    ylab("AUPRC")+
    xlab("")+
    theme(text = element_text(size = 15),
          plot.title = element_text(size = 16, hjust = 0.5, face = 'bold'),
          axis.text.y=element_text(size=13, color = 'black'),
          axis.text.x=element_text(size=13, angle = 45, vjust=1, hjust=1, color = 'black'),
          legend.position="none")+
    ggtitle('Total')
  p42 <- ggplot(result, aes(methods, AUROC, fill=methods)) + 
    geom_boxplot(lwd = 0.5, position=position_dodge(width=0.8))+ #outlier.shape = NA, 
    theme_bw()+
    ylab("AUROC")+
    xlab("")+
    theme(text = element_text(size = 15),
          plot.title = element_text(size = 16, hjust = 0.5, face = 'bold'),
          axis.text.y=element_text(size=13, color = 'black'),
          axis.text.x=element_text(size=13, angle = 45, vjust=1, hjust=1, color = 'black'),
          legend.position="none")+
    ggtitle('Total')
  
  cowplot::plot_grid(p1, p12,  labels = LETTERS[1:2], label_size = 20)
  cowplot::plot_grid(p2, p22, labels = LETTERS[3:4], label_size = 20)
  cowplot::plot_grid(p3, p32, labels = LETTERS[5:6], label_size = 20)
  
  cowplot::plot_grid(p4, p42, labels = LETTERS[7:8], label_size = 20)
}

############################
## Supplementary Figure 9 ##
############################


scales::show_col(pal_aaas(palette = "default", alpha = 0.6)(10))
mycolors_aaas <- pal_aaas(palette = "default", alpha = 0.6)(10)

mycolor_software <- mycolors_aaas[c(1:5)]
names(mycolor_software) <- c('CytoTalk', "NicheNet", "MISTy", "HoloNet", "stMLnet")
scales::show_col(mycolor_software)


result_index <- readRDS('./Data/Step8_LRTBenchResult/result_index.rds')
result_index <- do.call(rbind, result_index)
result_index <- tibble::rownames_to_column(result_index, 'datasets')
result_index$datasets <- gsub('\\.[0-9]+', '', result_index$datasets)
result_index <- result_index[which(result_index$methods %in% 
                                     c('CytoTalk', "NicheNet", "MISTy", "HoloNet", "stMLnet")), ]
result_index$methods <- factor(result_index$methods,
                               levels = c('CytoTalk', "NicheNet", "MISTy", "HoloNet", "stMLnet"))
result_index <- result_index[result_index$datasets %in% 
                               c("CID4465", "CID44971", "control_P7", "control_P8",
                                 'UKF260_T_ST', 'UKF266_T_ST', 'UKF334_T_ST', 'UKF243_T_ST'), ]
result_index <- result_index[which(!grepl('TAM_RB', result_index$celllines)), ]

rename_dataset <- function(x){
  switch(EXPR = x,
         "GSE120268_AXL_KO" = 'AXL_BC',
         "GSE157680_NRP1_MDA-MB-231_KO" = 'NRP1_BC',        
         "GSE15893_CXCR4_KO" = 'CXCR4_BC',
         "GSE15893_CXCL12_Treatment" = 'CXCL12_BC',
         "GSE160990_TGFB1_Treatment" = 'TGFB1_BC',
         "GSE36051_DLL4_MCF7_Treatment" = 'DLL4_BC(1)',
         "GSE36051_DLL4_MDA-MB-231_Treatment" = 'DLL4_BC(2)',
         "GSE36051_JAG1_MDA-MB-231_Treatment" = 'JAG1_BC',
         "GSE65398_IGF1_Treatment" = 'IGF1_BC(1)',     
         "GSE7561_IGF1_Treatment" = 'IGF1_BC(2)', 
         "GSE69104_CSF1R_KO_TAM_RB" = 'CSF1R_RB_TAM',
         "GSE69104_CSF1R_KO_TAM_EP" = 'CSF1R_EP_TAM',
         "GSE116414_FGFR1_KO" = 'FGFR1_GBM',
         'GSE140145_ALK_KO' = 'ALK_GBM',
         "GSE206947_EFNB2_Treatment" = 'EFNB2_CF',
         "GSE181575_TGFB1_Treatment" = 'TGFB1_CF',
         "GSE123018_TGFB1_Treatment_24hrs" = 'TGFB1_24h_CF',
         "GSE123018_TGFB1_Treatment_2hrs" = 'TGFB1_2h_CF',     
         "GSE123018_TGFB1_Treatment_45mins" = 'TGFB1_45m_CF', 
         "GSE123018_TGFB1_Treatment_6hrs" = 'TGFB1_6h_CF'
         
  )
}

df_plot <- result_index
df_plot$celllines <- lapply(df_plot$celllines,rename_dataset) %>% unlist()
df_plot$AUROC_label <- round(df_plot$AUROC, digits = 2)
df_plot$AUPRC_label <- round(df_plot$AUPRC, digits = 2)

df_plot$celllines <- factor(df_plot$celllines,
                            levels = c('AXL_BC', 'NRP1_BC', 'CXCR4_BC',
                                       'CXCL12_BC', 'TGFB1_BC', 'DLL4_BC(1)',
                                       'DLL4_BC(2)', 'JAG1_BC', 'IGF1_BC(1)',     
                                       'IGF1_BC(2)', 'CSF1R_EP_TAM', 'FGFR1_GBM',  
                                       'ALK_GBM', 'EFNB2_CF', 'TGFB1_CF', 'TGFB1_45m_CF',
                                       'TGFB1_2h_CF', 'TGFB1_6h_CF','TGFB1_24h_CF'))

df_plot <- lapply(unique(df_plot$celllines), fzunction(ds){
  #print(ds)
  df <- df_plot[df_plot$celllines == ds,]
  if(sum(df$AUROC, na.rm = TRUE)==0){
    df = NULL
  }else{
    df = df
  }
  
})
delet_celllines <- which(lapply(df_plot, is.null) %>% unlist())
if(length(delet_celllines)>0) df_plot <- df_plot[-delet_celllines]
df_plot <- df_plot %>% do.call('rbind',.)

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

### barplot —— BreastCancer
if(T){
  
  ##############
  ## CID44971 ## AUPRC
  ############## 
  
  df_plot_sub <- df_plot[which(df_plot$celltype == 'CancerEpithelial' & df_plot$datasets == 'CID44971'), ]
  df_plot_sub$AUPRC_label[is.na(df_plot_sub$AUPRC_label)] <- 0
  p1 <- ggplot(df_plot_sub, aes(x=celllines,y=AUPRC_label,fill=methods)) + 
    geom_bar(stat="identity",position=position_dodge(0.75),width=0.6) +
    scale_y_continuous(expand = c(0, 0.01)) + 
    labs(y = 'AUPRC') +
    coord_cartesian(ylim=c(-0.02,1.1))+
    scale_fill_manual(values = mycolor_software, 
                      labels=c('CytoTalk', "NicheNet", "MISTy", "HoloNet", "stMLnet")) + 
    theme_bar() + #coord_flip() +
    geom_text(aes(label=AUPRC_label),position=position_dodge(0.9),vjust = -0.5,size = 4)+
    ggtitle('CID44971')+
    theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
          axis.text.x = element_text(angle=20, vjust=1, hjust=1))
  
  p1
  
  #############
  ## CID4465 ##
  ############# AUPRC
  df_plot_sub <- df_plot[which(df_plot$celltype == 'CancerEpithelial' & df_plot$datasets == 'CID4465'), ]
  df_plot_sub$AUPRC_label[is.na(df_plot_sub$AUPRC_label)] <- 0
  p2 <- ggplot(df_plot_sub, aes(x=celllines,y=AUPRC_label,fill=methods)) + 
    geom_bar(stat="identity",position=position_dodge(0.75),width=0.6) +
    scale_y_continuous(expand = c(0, 0.01)) + 
    labs(y = 'AUPRC') +
    coord_cartesian(ylim=c(-0.02,1.1))+
    scale_fill_manual(values = mycolor_software, 
                      labels=c('CytoTalk', "NicheNet", "MISTy", "HoloNet", "stMLnet")) + 
    theme_bar() + #coord_flip() +
    geom_text(aes(label=AUPRC_label),position=position_dodge(0.9),vjust = -0.5,size = 4)+
    ggtitle('CID4465')+
    theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
          axis.text.x = element_text(angle=20, vjust=1, hjust=1))
  p2
}
#cowplot::plot_grid(p1, p2, nrow = 2)

### barplot —— cardiac fibroblasts
if(T){
  
  ################
  ## control_P7 ##
  ################ AUPRC
  
  df_plot_sub <- df_plot[which(df_plot$celltype == 'Fibroblast' & df_plot$datasets == 'control_P7'), ]
  df_plot_sub$AUPRC_label[is.na(df_plot_sub$AUPRC_label)] <- 0
  p3 <- ggplot(df_plot_sub, aes(x=celllines,y=AUPRC_label,fill=methods)) + 
    geom_bar(stat="identity",position=position_dodge(0.75),width=0.6) +
    scale_y_continuous(expand = c(0, 0.01)) + 
    labs(y = 'AUPRC') +
    coord_cartesian(ylim=c(-0.02,1.1))+
    scale_fill_manual(values = mycolor_software, 
                      labels=c('CytoTalk', "NicheNet", "MISTy", "HoloNet", "stMLnet")) + 
    theme_bar() + #coord_flip() +
    geom_text(aes(label=AUPRC_label),position=position_dodge(0.9),vjust = -0.5,size = 4)+
    ggtitle('control_P7')+
    theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
          axis.text.x = element_text(angle=20, vjust=1, hjust=1))
  p3
  
  ################
  ## control_P8 ##
  ################ AUPRC
  
  df_plot_sub <- df_plot[which(df_plot$celltype == 'Fibroblast' & df_plot$datasets == 'control_P8'), ]
  df_plot_sub$AUPRC_label[is.na(df_plot_sub$AUPRC_label)] <- 0
  p4 <- ggplot(df_plot_sub, aes(x=celllines,y=AUPRC_label,fill=methods)) + 
    geom_bar(stat="identity",position=position_dodge(0.75),width=0.6) +
    scale_y_continuous(expand = c(0, 0.01)) + 
    labs(y = 'AUPRC') +
    coord_cartesian(ylim=c(-0.02,1.1))+
    scale_fill_manual(values = mycolor_software, 
                      labels=c('CytoTalk', "NicheNet", "MISTy", "HoloNet", "stMLnet")) + 
    theme_bar() + #coord_flip() +
    geom_text(aes(label=AUPRC_label),position=position_dodge(0.9),vjust = -0.5,size = 4) +
    ggtitle('control_P8')+
    theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
          axis.text.x = element_text(angle=20, vjust=1, hjust=1))
  p4
}

### barplot —— GBM
if(T){
  
  #################
  ## UKF243_T_ST ##
  ################# AUPRC
  
  df_plot_sub <- df_plot[which(df_plot$celltype %in% c('macrophages', 'malignant') & df_plot$datasets == 'UKF243_T_ST'), ]
  df_plot_sub$AUPRC_label[is.na(df_plot_sub$AUPRC_label)] <- 0
  p5 <- ggplot(df_plot_sub, aes(x=celllines,y=AUPRC_label,fill=methods)) + 
    geom_bar(stat="identity",position=position_dodge(0.75),width=0.6) +
    scale_y_continuous(expand = c(0, 0.01)) + 
    labs(y = 'AUPRC') +
    coord_cartesian(ylim=c(-0.02,1.1))+
    scale_fill_manual(values = mycolor_software, 
                      labels=c('CytoTalk', "NicheNet", "MISTy", "HoloNet", "stMLnet")) + 
    theme_bar() + #coord_flip() +
    geom_text(aes(label=AUPRC_label),position=position_dodge(0.9),vjust = -0.5,size = 4)+
    ggtitle('UKF243_T_ST')+
    theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
          axis.text.x = element_text(angle=20, vjust=1, hjust=1))
  p5
  
  
  #################
  ## UKF260_T_ST ##
  ################# AUPRC
  
  df_plot_sub <- df_plot[which(df_plot$celltype %in% c('macrophages', 'malignant') & df_plot$datasets == 'UKF260_T_ST'), ]
  df_plot_sub$AUPRC_label[is.na(df_plot_sub$AUPRC_label)] <- 0
  p6 <- ggplot(df_plot_sub, aes(x=celllines,y=AUPRC_label,fill=methods)) + 
    geom_bar(stat="identity",position=position_dodge(0.75),width=0.6) +
    scale_y_continuous(expand = c(0, 0.01)) + 
    labs(y = 'AUPRC') +
    coord_cartesian(ylim=c(-0.02,1.1))+
    scale_fill_manual(values = mycolor_software, 
                      labels=c('CytoTalk', "NicheNet", "MISTy", "HoloNet", "stMLnet")) + 
    theme_bar() + #coord_flip() +
    geom_text(aes(label=AUPRC_label),position=position_dodge(0.9),vjust = -0.5,size = 4)+
    ggtitle('UKF260_T_ST') +
    theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
          axis.text.x = element_text(angle=20, vjust=1, hjust=1))
  p6
  
  #################
  ## UKF266_T_ST ##
  ################# AUPRC
  
  df_plot_sub <- df_plot[which(df_plot$celltype %in% c('macrophages', 'malignant') & df_plot$datasets == 'UKF266_T_ST'), ]
  df_plot_sub$AUPRC_label[is.na(df_plot_sub$AUPRC_label)] <- 0
  p7 <- ggplot(df_plot_sub, aes(x=celllines,y=AUPRC_label,fill=methods)) + 
    geom_bar(stat="identity",position=position_dodge(0.75),width=0.6) +
    scale_y_continuous(expand = c(0, 0.01)) + 
    labs(y = 'AUPRC') +
    coord_cartesian(ylim=c(-0.02,1.1))+
    scale_fill_manual(values = mycolor_software, 
                      labels=c('CytoTalk', "NicheNet", "MISTy", "HoloNet", "stMLnet")) + 
    theme_bar() + #coord_flip() +
    geom_text(aes(label=AUPRC_label),position=position_dodge(0.9),vjust = -0.5,size = 4)+
    ggtitle('UKF266_T_ST') +
    theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
          axis.text.x = element_text(angle=20, vjust=1, hjust=1))
  p7
  
  #################
  ## UKF334_T_ST ##
  ################# AUPRC
  
  df_plot_sub <- df_plot[which(df_plot$celltype %in% c('macrophages', 'malignant') & df_plot$datasets == 'UKF334_T_ST'), ]
  df_plot_sub$AUPRC_label[is.na(df_plot_sub$AUPRC_label)] <- 0
  p8 <- ggplot(df_plot_sub, aes(x=celllines,y=AUPRC_label,fill=methods)) + 
    geom_bar(stat="identity",position=position_dodge(0.75),width=0.6) +
    scale_y_continuous(expand = c(0, 0.01)) + 
    labs(y = 'AUPRC') +
    coord_cartesian(ylim=c(-0.02,1.1))+
    scale_fill_manual(values = mycolor_software, 
                      labels=c('CytoTalk', "NicheNet", "MISTy", "HoloNet", "stMLnet")) + 
    theme_bar() + #coord_flip() +
    geom_text(aes(label=AUPRC_label),position=position_dodge(0.9),vjust = -0.5,size = 4)+
    ggtitle('UKF334_T_ST') +
    theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
          axis.text.x = element_text(angle=20, vjust=1, hjust=1))
  p8
  
  
}

cowplot::plot_grid(p1, p2, p3, p4, ncol = 2, labels = LETTERS[1:4], label_size = 20)
cowplot::plot_grid(p5, p6, p7, p8, ncol = 2, labels = LETTERS[5:8], label_size = 20)


#############################
## Supplementary Figure 10 ##
#############################

if(T){
  setwd('~/2-CCC-benchmark-202204-202208/CCC-benchmark-202212/')
  
  result_index <- readRDS('./Data/Step8_LRTBenchResult/result_index.rds')
  result_index <- do.call(rbind, result_index)
  result_index <- tibble::rownames_to_column(result_index, 'datasets')
  result_index$datasets <- gsub('\\.[0-9]+', '', result_index$datasets)
  result_index <- result_index[which(result_index$methods %in% 
                                       c('CytoTalk', "NicheNet", "MISTy", "HoloNet", "stMLnet")), ]
  
  result_index$methods <- factor(result_index$methods,
                                 levels = c('CytoTalk', "NicheNet", "MISTy", "HoloNet", "stMLnet"))
  
  tmp <- result_index[result_index$datasets %in% 
                        c("CID4465", "CID44971", "control_P7", "control_P8", 
                          'UKF260_T_ST', 'UKF266_T_ST', 'UKF334_T_ST', 'UKF243_T_ST'), ]
  tmp <- tmp[which(!grepl('TAM_RB', tmp$celllines)), ]
  
  
  p1 <- ggpubr::ggviolin(tmp, "methods", "AUROC", color = "methods", add =  c("jitter", "boxplot"))+
    theme(legend.position = "none")+
    theme(axis.text.x = element_text(angle=20, vjust=1, hjust=1), plot.title = element_text(hjust = 0.5))
  
  
  methods <- as.character(unique(result_index$methods))
  stat_result <- lapply(methods, function(method1){
    result1 <- result_index[which(result_index$methods==method1), ]
    tmp <- lapply(methods, function(method2){
      result2 <- result_index[which(result_index$methods==method2), ]
      wilcox.test(result1$AUROC, result2$AUROC, alternative = "greater")$p.value
    })
    tmp <- do.call(cbind, tmp)
    colnames(tmp) <- methods
    tmp
  })
  stat_result <- do.call(rbind, stat_result)
  rownames(stat_result) <- methods
  stat_binary <- ifelse(stat_result<0.05, 1, 0)
  order_names <- c('CytoTalk', "NicheNet", "MISTy", "HoloNet", "stMLnet")
  stat_binary <- stat_binary[order_names, order_names]
  
  p2 <- pheatmap::pheatmap(stat_binary,cluster_rows = FALSE, cluster_cols = FALSE,
                           treeheight_row = 0,
                           treeheight_col = 0,
                           display_numbers = TRUE, number_color = 'white',
                           fontsize=20,
                           #annotation_col = anno,
                           number_format = "%.0f",
                           border_color=NA,
                           angle_col = "45",
                           color = colorRampPalette(c("#B4BFB9", "#6499A7"))(2),
                           #main = paste0('SI in ', data),
                           #na_col = 'black', 
                           silent = FALSE
  )
  cowplot::plot_grid(p1, p2$gtable, labels = LETTERS[1:2], label_size = 20)
}

if(T){
  scales::show_col(pal_aaas(palette = "default", alpha = 0.6)(10))
  mycolors_aaas <- pal_aaas(palette = "default", alpha = 0.6)(10)
  
  mycolor_software <- mycolors_aaas[c(1:5)]
  names(mycolor_software) <- c('CytoTalk', "NicheNet", "MISTy", "HoloNet", "stMLnet")
  scales::show_col(mycolor_software)
  
  
  result_index <- readRDS('./Data/Step8_LRTBenchResult/result_index.rds')
  result_index <- do.call(rbind, result_index)
  result_index <- tibble::rownames_to_column(result_index, 'datasets')
  result_index$datasets <- gsub('\\.[0-9]+', '', result_index$datasets)
  result_index <- result_index[which(result_index$methods %in% 
                                       c('CytoTalk', "NicheNet", "MISTy", "HoloNet", "stMLnet")), ]
  result_index$methods <- factor(result_index$methods,
                                 levels = c('CytoTalk', "NicheNet", "MISTy", "HoloNet", "stMLnet"))
  result_index <- result_index[result_index$datasets %in% 
                                 c("CID4465", "CID44971", "control_P7", "control_P8",
                                   'UKF260_T_ST', 'UKF266_T_ST', 'UKF334_T_ST', 'UKF243_T_ST'), ]
  result_index <- result_index[which(!grepl('TAM_RB', result_index$celllines)), ]
  
  rename_dataset <- function(x){
    switch(EXPR = x,
           "GSE120268_AXL_KO" = 'AXL_BC',
           "GSE157680_NRP1_MDA-MB-231_KO" = 'NRP1_BC',        
           "GSE15893_CXCR4_KO" = 'CXCR4_BC',
           "GSE15893_CXCL12_Treatment" = 'CXCL12_BC',
           "GSE160990_TGFB1_Treatment" = 'TGFB1_BC',
           "GSE36051_DLL4_MCF7_Treatment" = 'DLL4_BC(1)',
           "GSE36051_DLL4_MDA-MB-231_Treatment" = 'DLL4_BC(2)',
           "GSE36051_JAG1_MDA-MB-231_Treatment" = 'JAG1_BC',
           "GSE65398_IGF1_Treatment" = 'IGF1_BC(1)',     
           "GSE7561_IGF1_Treatment" = 'IGF1_BC(2)', 
           "GSE69104_CSF1R_KO_TAM_RB" = 'CSF1R_RB_TAM',
           "GSE69104_CSF1R_KO_TAM_EP" = 'CSF1R_EP_TAM',
           "GSE116414_FGFR1_KO" = 'FGFR1_GBM',
           'GSE140145_ALK_KO' = 'ALK_GBM',
           "GSE206947_EFNB2_Treatment" = 'EFNB2_CF',
           "GSE181575_TGFB1_Treatment" = 'TGFB1_CF',
           "GSE123018_TGFB1_Treatment_24hrs" = 'TGFB1_24h_CF',
           "GSE123018_TGFB1_Treatment_2hrs" = 'TGFB1_2h_CF',     
           "GSE123018_TGFB1_Treatment_45mins" = 'TGFB1_45m_CF', 
           "GSE123018_TGFB1_Treatment_6hrs" = 'TGFB1_6h_CF'
           
    )
  }
  
  df_plot <- result_index
  df_plot$celllines <- lapply(df_plot$celllines,rename_dataset) %>% unlist()
  df_plot$AUROC_label <- round(df_plot$AUROC, digits = 2)
  df_plot$AUROC_label <- round(df_plot$AUROC, digits = 2)
  
  df_plot$celllines <- factor(df_plot$celllines,
                              levels = c('AXL_BC', 'NRP1_BC', 'CXCR4_BC',
                                         'CXCL12_BC', 'TGFB1_BC', 'DLL4_BC(1)',
                                         'DLL4_BC(2)', 'JAG1_BC', 'IGF1_BC(1)',     
                                         'IGF1_BC(2)', 'CSF1R_EP_TAM', 'FGFR1_GBM',  
                                         'ALK_GBM', 'EFNB2_CF', 'TGFB1_CF', 'TGFB1_45m_CF',
                                         'TGFB1_2h_CF', 'TGFB1_6h_CF','TGFB1_24h_CF'))
  
  df_plot <- lapply(unique(df_plot$celllines), function(ds){
    #print(ds)
    df <- df_plot[df_plot$celllines == ds,]
    if(sum(df$AUROC, na.rm = TRUE)==0){
      df = NULL
    }else{
      df = df
    }
    
  })
  delet_celllines <- which(lapply(df_plot, is.null) %>% unlist())
  if(length(delet_celllines)>0) df_plot <- df_plot[-delet_celllines]
  df_plot <- df_plot %>% do.call('rbind',.)
  
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
  
  ### barplot —— BreastCancer
  if(T){
    
    ##############
    ## CID44971 ## AUROC
    ############## 
    
    df_plot_sub <- df_plot[which(df_plot$celltype == 'CancerEpithelial' & df_plot$datasets == 'CID44971'), ]
    df_plot_sub$AUROC_label[is.na(df_plot_sub$AUROC_label)] <- 0
    p1 <- ggplot(df_plot_sub, aes(x=celllines,y=AUROC_label,fill=methods)) + 
      geom_bar(stat="identity",position=position_dodge(0.75),width=0.6) +
      scale_y_continuous(expand = c(0, 0.01)) + 
      labs(y = 'AUROC') +
      coord_cartesian(ylim=c(-0.02,1.1))+
      scale_fill_manual(values = mycolor_software, 
                        labels=c('CytoTalk', "NicheNet", "MISTy", "HoloNet", "stMLnet")) + 
      theme_bar() + #coord_flip() +
      geom_text(aes(label=AUROC_label),position=position_dodge(0.9),vjust = -0.5,size = 4)+
      ggtitle('CID44971')+
      theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
            axis.text.x = element_text(angle=20, vjust=1, hjust=1))
    
    p1
    
    #############
    ## CID4465 ##
    ############# AUROC
    df_plot_sub <- df_plot[which(df_plot$celltype == 'CancerEpithelial' & df_plot$datasets == 'CID4465'), ]
    df_plot_sub$AUROC_label[is.na(df_plot_sub$AUROC_label)] <- 0
    p2 <- ggplot(df_plot_sub, aes(x=celllines,y=AUROC_label,fill=methods)) + 
      geom_bar(stat="identity",position=position_dodge(0.75),width=0.6) +
      scale_y_continuous(expand = c(0, 0.01)) + 
      labs(y = 'AUROC') +
      coord_cartesian(ylim=c(-0.02,1.1))+
      scale_fill_manual(values = mycolor_software, 
                        labels=c('CytoTalk', "NicheNet", "MISTy", "HoloNet", "stMLnet")) + 
      theme_bar() + #coord_flip() +
      geom_text(aes(label=AUROC_label),position=position_dodge(0.9),vjust = -0.5,size = 4)+
      ggtitle('CID4465')+
      theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
            axis.text.x = element_text(angle=20, vjust=1, hjust=1))
    p2
  }
  #cowplot::plot_grid(p1, p2, nrow = 2)
  
  ### barplot —— cardiac fibroblasts
  if(T){
    
    ################
    ## control_P7 ##
    ################ AUROC
    
    df_plot_sub <- df_plot[which(df_plot$celltype == 'Fibroblast' & df_plot$datasets == 'control_P7'), ]
    df_plot_sub$AUROC_label[is.na(df_plot_sub$AUROC_label)] <- 0
    p3 <- ggplot(df_plot_sub, aes(x=celllines,y=AUROC_label,fill=methods)) + 
      geom_bar(stat="identity",position=position_dodge(0.75),width=0.6) +
      scale_y_continuous(expand = c(0, 0.01)) + 
      labs(y = 'AUROC') +
      coord_cartesian(ylim=c(-0.02,1.1))+
      scale_fill_manual(values = mycolor_software, 
                        labels=c('CytoTalk', "NicheNet", "MISTy", "HoloNet", "stMLnet")) + 
      theme_bar() + #coord_flip() +
      geom_text(aes(label=AUROC_label),position=position_dodge(0.9),vjust = -0.5,size = 4)+
      ggtitle('control_P7')+
      theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
            axis.text.x = element_text(angle=20, vjust=1, hjust=1))
    p3
    
    ################
    ## control_P8 ##
    ################ AUROC
    
    df_plot_sub <- df_plot[which(df_plot$celltype == 'Fibroblast' & df_plot$datasets == 'control_P8'), ]
    df_plot_sub$AUROC_label[is.na(df_plot_sub$AUROC_label)] <- 0
    p4 <- ggplot(df_plot_sub, aes(x=celllines,y=AUROC_label,fill=methods)) + 
      geom_bar(stat="identity",position=position_dodge(0.75),width=0.6) +
      scale_y_continuous(expand = c(0, 0.01)) + 
      labs(y = 'AUROC') +
      coord_cartesian(ylim=c(-0.02,1.1))+
      scale_fill_manual(values = mycolor_software, 
                        labels=c('CytoTalk', "NicheNet", "MISTy", "HoloNet", "stMLnet")) + 
      theme_bar() + #coord_flip() +
      geom_text(aes(label=AUROC_label),position=position_dodge(0.9),vjust = -0.5,size = 4) +
      ggtitle('control_P8')+
      theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
            axis.text.x = element_text(angle=20, vjust=1, hjust=1))
    p4
  }
  
  ### barplot —— GBM
  if(T){
    
    #################
    ## UKF243_T_ST ##
    ################# AUROC
    
    df_plot_sub <- df_plot[which(df_plot$celltype %in% c('macrophages', 'malignant') & df_plot$datasets == 'UKF243_T_ST'), ]
    df_plot_sub$AUROC_label[is.na(df_plot_sub$AUROC_label)] <- 0
    p5 <- ggplot(df_plot_sub, aes(x=celllines,y=AUROC_label,fill=methods)) + 
      geom_bar(stat="identity",position=position_dodge(0.75),width=0.6) +
      scale_y_continuous(expand = c(0, 0.01)) + 
      labs(y = 'AUROC') +
      coord_cartesian(ylim=c(-0.02,1.1))+
      scale_fill_manual(values = mycolor_software, 
                        labels=c('CytoTalk', "NicheNet", "MISTy", "HoloNet", "stMLnet")) + 
      theme_bar() + #coord_flip() +
      geom_text(aes(label=AUROC_label),position=position_dodge(0.9),vjust = -0.5,size = 4)+
      ggtitle('UKF243_T_ST')+
      theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
            axis.text.x = element_text(angle=20, vjust=1, hjust=1))
    p5
    
    
    #################
    ## UKF260_T_ST ##
    ################# AUROC
    
    df_plot_sub <- df_plot[which(df_plot$celltype %in% c('macrophages', 'malignant') & df_plot$datasets == 'UKF260_T_ST'), ]
    df_plot_sub$AUROC_label[is.na(df_plot_sub$AUROC_label)] <- 0
    p6 <- ggplot(df_plot_sub, aes(x=celllines,y=AUROC_label,fill=methods)) + 
      geom_bar(stat="identity",position=position_dodge(0.75),width=0.6) +
      scale_y_continuous(expand = c(0, 0.01)) + 
      labs(y = 'AUROC') +
      coord_cartesian(ylim=c(-0.02,1.1))+
      scale_fill_manual(values = mycolor_software, 
                        labels=c('CytoTalk', "NicheNet", "MISTy", "HoloNet", "stMLnet")) + 
      theme_bar() + #coord_flip() +
      geom_text(aes(label=AUROC_label),position=position_dodge(0.9),vjust = -0.5,size = 4)+
      ggtitle('UKF260_T_ST') +
      theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
            axis.text.x = element_text(angle=20, vjust=1, hjust=1))
    p6
    
    #################
    ## UKF266_T_ST ##
    ################# AUROC
    
    df_plot_sub <- df_plot[which(df_plot$celltype %in% c('macrophages', 'malignant') & df_plot$datasets == 'UKF266_T_ST'), ]
    df_plot_sub$AUROC_label[is.na(df_plot_sub$AUROC_label)] <- 0
    p7 <- ggplot(df_plot_sub, aes(x=celllines,y=AUROC_label,fill=methods)) + 
      geom_bar(stat="identity",position=position_dodge(0.75),width=0.6) +
      scale_y_continuous(expand = c(0, 0.01)) + 
      labs(y = 'AUROC') +
      coord_cartesian(ylim=c(-0.02,1.1))+
      scale_fill_manual(values = mycolor_software, 
                        labels=c('CytoTalk', "NicheNet", "MISTy", "HoloNet", "stMLnet")) + 
      theme_bar() + #coord_flip() +
      geom_text(aes(label=AUROC_label),position=position_dodge(0.9),vjust = -0.5,size = 4)+
      ggtitle('UKF266_T_ST') +
      theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
            axis.text.x = element_text(angle=20, vjust=1, hjust=1))
    p7
    
    #################
    ## UKF334_T_ST ##
    ################# AUROC
    
    df_plot_sub <- df_plot[which(df_plot$celltype %in% c('macrophages', 'malignant') & df_plot$datasets == 'UKF334_T_ST'), ]
    df_plot_sub$AUROC_label[is.na(df_plot_sub$AUROC_label)] <- 0
    p8 <- ggplot(df_plot_sub, aes(x=celllines,y=AUROC_label,fill=methods)) + 
      geom_bar(stat="identity",position=position_dodge(0.75),width=0.6) +
      scale_y_continuous(expand = c(0, 0.01)) + 
      labs(y = 'AUROC') +
      coord_cartesian(ylim=c(-0.02,1.1))+
      scale_fill_manual(values = mycolor_software, 
                        labels=c('CytoTalk', "NicheNet", "MISTy", "HoloNet", "stMLnet")) + 
      theme_bar() + #coord_flip() +
      geom_text(aes(label=AUROC_label),position=position_dodge(0.9),vjust = -0.5,size = 4)+
      ggtitle('UKF334_T_ST') +
      theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
            axis.text.x = element_text(angle=20, vjust=1, hjust=1))
    p8
    
    
  }
  
  cowplot::plot_grid(p1, p2, p3, p4, ncol = 2, labels = LETTERS[3:6], label_size = 20)
  cowplot::plot_grid(p5, p6, p7, p8, ncol = 2, labels = LETTERS[7:10], label_size = 20)
}
