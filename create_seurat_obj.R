#' create Seurat object 
#' 
#' @description create a custom seurat object from the files downloaded from the broad insititute single cell database
#' 
create_seurat_obj <- function(path) {
   matrix_path <- list.files(path, ".*raw.mtx", full.names = TRUE)
   barcodes_path <- list.files(path, ".*.barcodes.tsv", full.names = TRUE)
   features_path <-  list.files(path, ".*.features.tsv", full.names = TRUE)
  
   sample_name <- matrix_path |> basename()  |> substr(1,6)

   counts_matrix <- ReadMtx(
     mtx = matrix_path,
     cells = barcodes_path,
     features = features_path,
     cell.column = 1,     # Usually barcodes are in the first column
     feature.column = 2   # Change to 1 if your gene file only has 1 column (IDs instead of Symbols)
   )
  # Create the Seurat object
   obj <- CreateSeuratObject(counts = counts_matrix, project = sample_name)
  
  # Tag with custom metadata
   obj$Tissue <- sample_name
  
   return(obj)
  
}
