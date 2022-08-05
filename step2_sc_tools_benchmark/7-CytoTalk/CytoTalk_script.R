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

setwd("./step2_sc_tools_benchmark/7-CytoTalk/")

source("./read_matrix_with_meta.R")

id <- "CID4465" 
fpath.mat <- paste0("../3-CellPhoneDB/input/", id,"_cellphonedb_count.txt")
fpath.meta <- paste0("../3-CellPhoneDB/input/", id,"_cellphonedb_meta.txt")
lst.sc <- read_matrix_with_meta(fpath.mat, fpath.meta)
celltype <- unique(lst.sc$cell_types)

result_CytoTalk <- list()
for (sender in celltype) {
  for (reciever in celltype) {
    results <- tryCatch(CytoTalk::run_cytotalk(lst.sc, sender, reciever, cores = 10, cutoff_a = 0.05, cutoff_b = 0.05),
                        error=function(e){NA}
    )
    cp <- paste(sender, reciever, sep = "-")
    result_CytoTalk[[cp]] <- results
  }
}

save(result_CytoTalk, file = './result.CytoTalk.RData')
