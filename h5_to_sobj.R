# Data subset from GSE204765 (h5ad format)
library(anndata)
library(Seurat)
library(Matrix)

wd <- "~/projects/aderm_multiome/samples_to_use"

# Converts h5ad to Seurat RDS
# Note: Assumes data$X contains the count matrix (raw, filtered, or normalized)
convert_h5ad_to_seurat <- function(input_path, output_path, min_features = 100, min_cells = 3) {
  adata <- anndata::read_h5ad(input_path)
  
  # Uses Matrix::t for sparse-friendly transposing to save memory
  counts <- Matrix::t(adata$X)
  
  obj <- CreateSeuratObject( 
    counts       = counts,
    meta.data    = adata$obs, 
    min.features = min_features,
    min.cells    = min_cells
  )
  
  saveRDS(obj, output_path)
  message("Successfully converted and saved: ", basename(output_path))
  return(invisible(obj))
}

# Identify all .h5ad files for conversion
h5ad_files <- list.files(
  path       = wd,
  pattern    = "\\.h5ad$",
  full.names = TRUE
)

# Batch process files; skips existing .rds files to prevent redundant processing of samples from previous runs
lapply(h5ad_files, function(f) {
  out_path <- sub("\\.h5ad$", ".rds", f)
  
  if (!file.exists(out_path)) {
    convert_h5ad_to_seurat(f, out_path)
  } else {
    message("Skipping existing file: ", basename(out_path))
  }
})

# See script [2] for sanity check of converted objects and initial quality control
