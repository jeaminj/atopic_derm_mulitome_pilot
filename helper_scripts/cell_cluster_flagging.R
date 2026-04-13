# =============================================================================
# Helper: Flag misassigned cells based on UMAP outlier detection
# Usage: source("helpers/flag_misassigned_by_umap.R")
# =============================================================================

flag_misassigned_by_umap <- function(sobj,
                                     reduction = "umap_harmony",
                                     cell_type_col = "Cell.type",
                                     sd_threshold = 3) {
  
  # Extract UMAP coordinates and dynamically get column names
  umap_coords <- as.data.frame(Embeddings(sobj, reduction = reduction))
  dim1 <- colnames(umap_coords)[1]  # e.g. "umap_harmony_1"
  dim2 <- colnames(umap_coords)[2]  # e.g. "umap_harmony_2"
  cat(sprintf("Using UMAP dimensions: %s, %s\n", dim1, dim2))
  
  umap_coords[[cell_type_col]] <- sobj[[cell_type_col, drop = TRUE]]
  cell_types <- unique(umap_coords[[cell_type_col]])
  
  is_outlier <- rep(FALSE, ncol(sobj))
  names(is_outlier) <- colnames(sobj)
  
  cat("UMAP outlier detection per cell type:\n")
  for (ct in cell_types) {
    ct_cells  <- umap_coords[[cell_type_col]] == ct
    ct_coords <- umap_coords[ct_cells, c(dim1, dim2)]
    
    mean_u1 <- mean(ct_coords[[dim1]])
    mean_u2 <- mean(ct_coords[[dim2]])
    sd_u1   <- sd(ct_coords[[dim1]])
    sd_u2   <- sd(ct_coords[[dim2]])
    
    ct_outliers <- ct_cells & (
      abs(umap_coords[[dim1]] - mean_u1) > sd_threshold * sd_u1 |
        abs(umap_coords[[dim2]] - mean_u2) > sd_threshold * sd_u2
    )
    
    cat(sprintf("  %-35s : %d outliers flagged\n", ct, sum(ct_outliers)))
    is_outlier[ct_outliers] <- TRUE
  }
  
  sobj$is_misassigned <- is_outlier
  sobj$removal_status <- ifelse(is_outlier, "To be removed", "Keep")
  
  cat(sprintf("\nTotal flagged: %d / %d (%.1f%%)\n",
              sum(is_outlier), ncol(sobj),
              100 * sum(is_outlier) / ncol(sobj)))
  
  return(sobj)
}

