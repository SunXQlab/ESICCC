rm(list = ls())

# library package
if(T){
  # MISTy
  library(mistyR)
  library(future)
  
  # Seurat
  library(Seurat)
  library(Giotto)
  
  # data manipulation
  library(Matrix)
  library(tibble)
  library(dplyr)
  library(purrr)
  
  # normalization
  library(sctransform)
  
  # resource
  library(progeny)
  
  # setup parallel execution
  options(future.globals.maxSize = 1024^3)
  plan(multisession)
}
set.seed(123)
source('./step4_tools_intra_benchmark/MistyR/code.R')

# load ST data
load('./step4_tools_intra_benchmark/stMLnet/input/CID4465_seruat.RData')

# run normalization
sct.data <- vst(GetAssayData(
  object = se.st,
  slot = "counts",
  assay = "Spatial"
),
verbosity = 0
)

se.st[["SCT"]] <- CreateAssayObject(data = sct.data$y)

gene.expression <- GetAssayData(se.st, assay = "SCT")
coverage <- rowSums(gene.expression > 0) / ncol(gene.expression)
slide.markers <- names(which(coverage >= 0.05))

Databases <- readRDS('./step4_tools_intra_benchmark/MistyR/prior_knowledge/Databases.rds')
ligands <- Databases$LigRec.DB$source %>% unique() %>% .[. %in% slide.markers]
receptors <- Databases$LigRec.DB$target %>% unique() %>% .[. %in% slide.markers]

icgs <- readRDS('./step1_data_process/result/Giotto_result/CID4465_bc_ICGs.rds')

celltypes <- se.st$celltype %>% as.character()

geometry <- se.st@meta.data %>% select(., row,col)

#########
## run ##
#########

ct <- "Cancer Epithelial"
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

# Define specific properties for each view
view.properties <- list(
  "main" = ifelse(celltypes == ct,1,0),
  "ligand" = ifelse(celltypes != ct,1,0)
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

#########
## Run ##
#########

spot.ids = NULL
out.alias = paste0("step4_tools_intra_benchmark/MistyR/result/results_",ct,"_paraview_10")

# Extracting data
view.data <- map(view.assays,
                 extract_seurat_data,
                 geometry = geometry,
                 visium.slide = se.st
)
str(view.data,max.level = 1)

# Adding all spots ids in case they are not defined
if (is.null(spot.ids)) {
  spot.ids <- rownames(view.data[[1]])
}

# First filter the features from the data
view.data.filt <- map2(view.data, view.features, filter_data_features)
str(view.data.filt,max.level = 1)
view.data.filt[[1]][1:4,1:4]
view.data.filt[[2]][1:4,1:4]

# specific properties: celltype
view.data.spec <- map2(view.data.filt, view.properties, add_specific_properties)
str(view.data.spec,max.level = 1)
view.data.spec[[1]][1:10,1:4]
view.data.spec[[2]][1:10,1:4]

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
  rename(regulator=Predictor,target=Target,score=Importance,type=view)
misty_score$type <- gsub('_10','',misty_score$type)
misty_score$type <- gsub('intra','receptor',misty_score$type)

head(misty_score)
saveRDS(misty_score,paste0("step4_tools_intra_benchmark/MistyR/result",'/misty_score.rds'))
