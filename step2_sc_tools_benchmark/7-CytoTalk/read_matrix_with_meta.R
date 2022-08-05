errorifnot <- function (expr, msg) 
{
  if (!expr) 
    stop(format_message(msg))
}

new_named_list <- function (mat, cell_types) 
{
  list(mat = mat, cell_types = cell_types)
}

read_matrix_with_meta <- function (fpath_mat, fpath_meta, auto_transform = TRUE) 
{
  #mat <- vroom_sparse_with_rownames(fpath_mat)
  mat <- read.table(fpath_mat)
  mat <- as.matrix(mat)
  #meta <- vroom_with_rownames(fpath_meta)
  meta <- read.table(fpath_meta, header = TRUE, sep = "\t")
  meta <- tibble::column_to_rownames(meta, var = "V1")
  index <- match(rownames(meta), colnames(mat))
  is(TRUE)
  errorifnot(!any(is.na(index)), "meta file does not match matrix colnames")
  cell_types <- meta[index, 1]
  mat <- check_count_data(mat, auto_transform)
  new_named_list(mat, cell_types)
}