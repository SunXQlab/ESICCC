rm(list = ls())
library(Seurat)
library(Connectome)
set.seed(123)

load("./step1_data_process/se_data/step2_scRNA_10X.RData")

connectome.genes <- union(Connectome::ncomms8866_human$Ligand.ApprovedSymbol,
                          Connectome::ncomms8866_human$Receptor.ApprovedSymbol)
genes <- connectome.genes[connectome.genes %in% rownames(se.sc)]
se.sc <- ScaleData(se.sc,features = genes)
sc.con <- CreateConnectome(se.sc,species = 'human', calculate.DOR = TRUE)

result.Connectome <- FilterConnectome(sc.con, max.p = 0.05, min.pct = 0.05, remove.na = T)
save(result.Connectome, file = './step2_sc_tools_benchmark/14-Connectome/result.Connectome.RData')
