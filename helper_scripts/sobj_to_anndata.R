# =============================================================================
# Helper: Convert Seurat object to AnnData (.h5ad) for Python/Scanpy
# Usage: source("helpers/sobj_to_anndata.R")
# NOTE: Won't work with Seurat5 objects; use sobjV5_to_anndata.R instead !!
# =============================================================================
# TO-DO; test this works lol

seurat_to_anndata <- function(sobj,
                              output_path,
                              assay    = "RNA",
                              overwrite = TRUE) {
  
  library(SeuratDisk)
  
  # Validate output path has no extension
  output_path <- tools::file_path_sans_ext(output_path)
  h5seurat_path <- paste0(output_path, ".h5seurat")
  h5ad_path     <- paste0(output_path, ".h5ad")
  
  # Ensure output directory exists
  out_dir <- dirname(output_path)
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
    cat(sprintf("Created output directory: %s\n", out_dir))
  }
  
  cat(sprintf("Converting Seurat object to AnnData...\n"))
  cat(sprintf("  Cells : %d\n", ncol(sobj)))
  cat(sprintf("  Genes : %d\n", nrow(sobj)))
  cat(sprintf("  Assay : %s\n", assay))
  
  # Step 1: Save as .h5seurat
  cat(sprintf("  Saving .h5seurat -> %s\n", h5seurat_path))
  SaveH5Seurat(sobj,
               filename  = h5seurat_path,
               assay     = assay,
               overwrite = overwrite)
  
  # Step 2: Convert to .h5ad
  cat(sprintf("  Converting to .h5ad -> %s\n", h5ad_path))
  Convert(h5seurat_path,
          dest      = "h5ad",
          overwrite = overwrite)
  
  # Optionally clean up intermediate .h5seurat
  if (file.exists(h5seurat_path)) {
    file.remove(h5seurat_path)
    cat(sprintf("  Cleaned up intermediate .h5seurat\n"))
  }
  
  if (file.exists(h5ad_path)) {
    cat(sprintf("Done. AnnData saved to: %s\n", h5ad_path))
  } else {
    warning("Conversion may have failed — .h5ad file not found at expected path.")
  }
  
  invisible(h5ad_path)  # return path invisibly for use in pipelines
}


# Usage example:
# sobj_kc <- subset(sobj_clean, Cell.type %in% c("Cornified keratinocytes", "Keratinocytes"))
# seurat_to_anndata(sobj_kc, output_path = "../04_downstream_analysis/pseudotime/sobj_kc")
