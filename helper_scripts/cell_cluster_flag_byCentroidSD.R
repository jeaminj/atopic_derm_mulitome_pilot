#' Flag Cells Based on UMAP Geometric Outliers
#' 
#' @description
#' Identifies cells that are spatially distant from their assigned cell type 
#' centroid in UMAP space using a Standard Deviation (SD) threshold. 
#' Useful for flagging potential misannotations, high-ambient noise, or 
#' doublet-derived "islands."
#'
#' @param sobj A Seurat object.
#' @param reduction The name of the reduction to use (default: "umap_harmony").
#' @param cell_type_col Metadata column containing categorical labels (default: "Cell.type").
#' @param sd_threshold Numeric; number of SDs from the mean to define an outlier (default: 3).
#' 
#' @details 
#' This function calculates the centroid (mean) and SD for each dimension 
#' of the UMAP per cell type. Cells exceeding the threshold in either 
#' dimension are flagged.
#' 
#' @note 
#' CAUTION: In differentiation systems (e.g., Keratinocytes), this may 
#' incorrectly flag "bridge" cells or transitional states that are 
#' biologically real but spatially distant from the main cluster centroid.
#'
#' @return The Seurat object with 'is_misassigned' and 'removal_status' added to metadata.
#' @export

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
    # effectively ignores NAs
    ct_idx <- which(!is.na(umap_coords[[cell_type_col]]) & umap_coords[[cell_type_col]] == ct)
    # Safety check: Ensure we have enough cells to calculate SD
    if (length(ct_idx) < 3) {
      cat(sprintf("  %-35s : Too few cells, skipping\n", ct))
      next
    }
    # Create a logical mask for the full object based on indices
    ct_cells <- rep(FALSE, nrow(umap_coords))
    ct_cells[ct_idx] <- TRUE
    
    ct_coords <- umap_coords[ct_idx, c(dim1, dim2)]
    
    # na.rm = TRUE is a backup safety measure
    mean_u1 <- mean(ct_coords[[dim1]], na.rm = TRUE)
    mean_u2 <- mean(ct_coords[[dim2]], na.rm = TRUE)
    sd_u1   <- sd(ct_coords[[dim1]], na.rm = TRUE)
    sd_u2   <- sd(ct_coords[[dim2]], na.rm = TRUE)
    
    # Identify outliers while ensuring NAs don't return TRUE
    outlier_idx <- ct_idx[
      abs(umap_coords[ct_idx, dim1] - mean_u1) > sd_threshold * sd_u1 |
        abs(umap_coords[ct_idx, dim2] - mean_u2) > sd_threshold * sd_u2
    ]
    
    # Remove any potential NAs that might have crept into the index
    outlier_idx <- outlier_idx[!is.na(outlier_idx)]
    
    cat(sprintf("  %-35s : %d outliers flagged\n", ct, length(outlier_idx)))
    is_outlier[outlier_idx] <- TRUE
  }
  
  sobj$is_misassigned <- is_outlier
  # Flag any NAs as "Unknown" so they don't accidentally get kept or removed without a look
  sobj$removal_status <- ifelse(is.na(sobj[[cell_type_col, drop = TRUE]]), 
                                "Missing Annotation", 
                                ifelse(is_outlier, "Far from centroid", "Within 4 SD"))
  
  cat(sprintf("\nTotal flagged: %d / %d (%.1f%%)\n",
              sum(is_outlier), ncol(sobj),
              100 * sum(is_outlier) / ncol(sobj)))
  
  return(sobj)
}

