rm(list = ls())

setwd("./step1_data_process/")

library(Seurat)
library(Giotto)
library(tibble)
library(stringr)
set.seed(123)

# Create Seurat Object for ST data
expr.data <- Read10X("./dataset/ST.Data/filtered_count_matrices/CID4465_filtered_count_matrix/", 
                     gene.column = 1)

# Load the image data
img <- Read10X_Image(image.dir = './dataset/ST.Data/spatial/CID4465_spatial/')
img <- img[colnames(expr.data)]
coords <- img@coordinates
coords <- coords[,4:5]
rm(img)

# Load the celltype metadata from stereoscope
meta.tmp <- read.table("./result/stereoscope_result_CID4465/CID4465_spatial_cnt/W.2022-05-15135505.473576.tsv", 
                       header = TRUE, sep = "\t")
meta.tmp <- column_to_rownames(meta.tmp, var = "X")
colnames(meta.tmp) <- str_replace_all(colnames(meta.tmp), "\\.", "_")
identical(rownames(meta.tmp), colnames(expr.data))

meta.tmp$celltype <- apply(meta.tmp, 1, function(x){
  colnames(meta.tmp)[which.max(x)]
})

meta.tmp$celltype <- stringr::str_replace_all(meta.tmp$celltype, "Cancer_Epithelial", "Cancer Epithelial")
meta.tmp$celltype <- stringr::str_replace_all(meta.tmp$celltype, "B_cells", "B-cells")
meta.tmp$celltype <- stringr::str_replace_all(meta.tmp$celltype, "T_cells", "T-cells")

## workdir
temp_dir = './result/Giotto_result/'
myinstructions = createGiottoInstructions(save_dir = temp_dir,
                                          save_plot = FALSE,
                                          show_plot = FALSE)

### Create a Giotto object 
gio <- createGiottoObject(raw_exprs = expr.data,
                             spatial_locs = coords,
                             instructions = myinstructions)

identical(rownames(meta.tmp), colnames(gio@raw_exprs))
gio <- addCellMetadata(gobject = gio,
                       new_metadata = meta.tmp$celltype, 
                       vector_name = "celltype")

### normalize
gio <- normalizeGiotto(gio)
gio <- addStatistics(gobject = gio)

### create a spatial Delaunay network (default)
gio = createSpatialNetwork(gobject = gio, method = "Delaunay")

### select top 25th highest expressing genes
gene_metadata = fDataDT(gio)
lapply(seq(0.1,1,0.05), quantile, x = gene_metadata$mean_expr_det, na.rm = T) %>% unlist()
high_expressed_genes = gene_metadata[mean_expr_det > 0.6962901]$gene_ID

### identify ICGs
CPGscoresHighGenes =  findICG(gobject = gio, 
                              selected_genes = high_expressed_genes,
                              spatial_network_name = 'Delaunay_network',
                              cluster_column = 'celltype',
                              diff_test = 'permutation', 
                              adjust_method = 'fdr',
                              nr_permutations = 500,
                              do_parallel = T, cores = 6)

### filter ICGs
CPGscoresFilt = filterICG(CPGscoresHighGenes, direction = "both")
table(CPGscoresFilt$CPGscores$spec_int[CPGscoresFilt$CPGscores$cell_type=="Cancer Epithelial"])
ICGs_list = lapply(unique(CPGscoresFilt$CPGscores$cell_type), function(x){
  y=CPGscoresFilt$CPGscores[CPGscoresFilt$CPGscores$cell_type==x,]
  z=lapply(unique(y$int_cell_type), function(t){
    y$genes[y$int_cell_type==t]
  })
  names(z)=unique(y$int_cell_type)
  z
})
names(ICGs_list) = unique(CPGscoresFilt$CPGscores$cell_type)
str(ICGs_list)

saveRDS(ICGs_list, file = "./result/Giotto_result/CID4465_bc_ICGs.rds")
saveRDS(gio, file = "./result/Giotto_result/CID4465_giotto_bc.rds")

