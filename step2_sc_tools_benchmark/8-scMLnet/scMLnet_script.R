rm(list = ls())
library(Seurat)
library(scMLnet)
set.seed(123)

# scMLnet v0.2.0 + database[scMLnet] + CID4465
# pacakge: https://github.com/SunXQlab/scMLnet2.0

load("./step1_data_process/se_data/step2_scRNA_10X.RData")

GCMat<- GetAssayData(se.sc, "counts", "RNA")
BarCluTable <- data.frame(Barcode = colnames(se.sc), Cluster = se.sc$celltype_major)
types <- unique(BarCluTable$Cluster)
rm(se.sc)

LigRecLib <- read.table("./step2_sc_tools_benchmark/8-scMLnet/database/LigRec.txt", header = T)
colnames(LigRecLib)[2:3] <- c("source", "target")
TFTarLib <- read.table("./step2_sc_tools_benchmark/8-scMLnet/database/TFTargetGene.txt", header = T)
colnames(TFTarLib)[1:2] <- c("source", "target")
RecTFLib <- read.table("./step2_sc_tools_benchmark/8-scMLnet/database/RecTF.txt", header = T)
colnames(RecTFLib)[1:2] <- c("source", "target")

result.scMLnet <- list()
for (LigClu in types) {
  subtypes <- types[-which(types == LigClu)]
  for (RecClu in subtypes) {
    netList <- tryCatch(RunMLnet(data = GCMat, BarCluTable = BarCluTable, 
                                 RecClu = RecClu, LigClu = LigClu,
                                 LigRec.DB = LigRecLib, 
                                 TFTG.DB = TFTarLib,
                                 RecTF.DB = RecTFLib),
                        error=function(e){NA}
    )
    list.names <- paste(LigClu, RecClu, sep = "_")
    result.scMLnet[[list.names]] <- netList
  }
}

save(result.scMLnet, file = "./step2_sc_tools_benchmark/8-scMLnet/result.scMLnet.RData")