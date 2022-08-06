rm(list = ls())
library(CytoTalk)

# create conda environment
if(F){
  library(reticulate)  # To install and call Python modules from R.
  conda_create(envname = "r_reticulate_CytoTalk", python_version = "3.7.3")  # Create a new Conda environment to facilitate the Python module installation.
  conda_install(envname = "r_reticulate_CytoTalk", "pybind11")  # Install two necessary Python modules for correctly compiling and using the "pcst_fast" Python module.
  conda_install(envname = "r_reticulate_CytoTalk", "numpy")
  conda_install(envname = "r_reticulate_CytoTalk", "git+https://github.com/fraenkel-lab/pcst_fast.git", pip = TRUE) # To install the "pcst_fast" module.
}

setwd("./step4_tools_intra_benchmark/CytoTalk/")

source("./read_matrix_with_meta.R")

fpath.mat <- paste0("./input/CID4465_norm_count.txt")
fpath.meta <- paste0("./input/CID4465_meta.txt")
lst.sc <- read_matrix_with_meta(fpath.mat, fpath.meta)
celltype <- unique(lst.sc$cell_types)

result_CytoTalk <- list()
for (sender in celltype) {
  reciever <- "Cancer Epithelial"
  for(reciever in celltype){
    print(paste(sender, reciever, sep = "_"))
    results <- tryCatch(CytoTalk::run_cytotalk(lst.sc, sender, reciever, cores = 100, cutoff_a = 0.05, cutoff_b = 0.05),
                        error=function(e){NA}
    )
    cp <- paste(sender, reciever, sep = "_")
    result_CytoTalk[[cp]] <- results
  }
}

save(result_CytoTalk, file = "./result/result_CytoTalk.RData")
