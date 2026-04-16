# =============================================================================
# Helper: Convert Seurat object to AnnData (.h5ad) for Python/Scanpy
# Usage: source("helpers/sobjV5_to_anndata.R")
# NOTE: For Seurat v5 objects 
# =============================================================================
seuratV5_to_anndata <- function(sobj,
                                output_path,
                                assay     = "RNA") {
  
  library(zellkonverter)
  library(SingleCellExperiment)
  library(Seurat)
  
  # 1. Setup paths
  output_path <- tools::file_path_sans_ext(output_path)
  h5ad_path   <- paste0(output_path, ".h5ad")
  
  # 2. Pre-processing for Seurat v5
  cat("  Joining layers and converting to standard Assay...\n")
  sobj <- JoinLayers(sobj, assay = assay)
  
  # Force Assay5 to standard Assay for zellkonverter stability
  sobj[[assay]] <- as(sobj[[assay]], "Assay") 
  
  # 3. Metadata Sanitation 
  cat("  Sanitizing 63 metadata columns...\n")
  # Remove lists/complex objects and convert factors to strings
  sobj@meta.data[] <- lapply(sobj@meta.data, function(x) {
    if (is.list(x)) return(NULL) # Remove nested lists
    if (is.factor(x) || is.logical(x)) return(as.character(x))
    return(x)
  })
  
  # 4. Convert to SingleCellExperiment
  cat("  Converting Seurat -> SCE...\n")
  sce <- as.SingleCellExperiment(sobj, assay = assay)
  
  # 5. Handle file overwriting manually
  if (file.exists(h5ad_path)) {
    cat("  Removing existing file...\n")
    file.remove(h5ad_path)
  }
  
  # 6. Write to h5ad
  cat(sprintf("  Writing .h5ad -> %s\n", h5ad_path))
  writeH5AD(sce, file = h5ad_path)
  
  cat("Success!\n")
  invisible(h5ad_path)
}

#TODO: add subset argument to function, null if null
# Usage example:
# sobj_kc <- subset(sobj_kc_clean, Cell.type %in% c("Cornified keratinocytes", "Keratinocytes"))
seuratV5_to_anndata(sobj_kc_clean, output_path = "../04_downstream_analysis/pseudotime/sobj_kc_clean")


