#' Advanced QC Outlier Detection
#' @param sobj Seurat object
#' @param nmads Number of Median Absolute Deviations for thresholds (default 3)
#' @param mt_pattern Regex for mitochondrial genes (default "^MT-")
#' @return Seurat object with QC flags in metadata: sobj$qc_outlier_reason
calculate_qc_outliers <- function(sobj, 
                                  nmads = 3, 
                                  mt_pattern = "^MT-") {
  
  suppressPackageStartupMessages({
    library(Seurat)
    library(Matrix)
    library(scater)
  })
  
  cat("Running Adaptive QC Detection...\n")
  
  # 1. Ensure MT percentage is calculated
  if (!"percent.mt" %in% colnames(sobj@meta.data)) {
    sobj[["percent.mt"]] <- PercentageFeatureSet(sobj, pattern = mt_pattern)
  }
  
  # 2. Extract metrics for outlier detection
  # We use log-scale for counts/features as they are usually log-normally distributed
  metrics <- data.frame(
    log_count = log10(sobj$nCount_RNA),
    log_feature = log10(sobj$nFeature_RNA),
    pct_mt = sobj$percent.mt
  )
  
  # 3. Define Outliers using MADs (Median Absolute Deviations)
  # High counts/features = potential doublets
  # Low counts/features = potential empty droplets/debris
  # High MT = dying cells
  is_low_feat  <- isOutlier(metrics$log_feature, nmads = nmads, type = "lower")
  is_high_feat <- isOutlier(metrics$log_feature, nmads = nmads, type = "higher")
  is_high_mt   <- isOutlier(metrics$pct_mt, nmads = nmads, type = "higher")
  
  # 4. Final flagging logic
  sobj$qc_outlier_reason <- "Pass QC metrics"
  sobj$qc_outlier_reason[is_low_feat]  <- "Low Features"
  sobj$qc_outlier_reason[is_high_feat] <- "Potential Doublet"
  sobj$qc_outlier_reason[is_high_mt]   <- "High Mitochondrial"
  
  # Multi-hit logic (if a cell fails two checks)
  multi_hit <- (is_low_feat + is_high_feat + is_high_mt) > 1
  sobj$qc_outlier_reason[multi_hit] <- "Multiple QC Failures"
  
  sobj$is_qc_outlier <- sobj$qc_outlier_reason != "Pass QC metrics"
  
  # Summary report
  stats <- table(sobj$qc_outlier_reason)
  cat("\nQC Summary Table:\n")
  print(stats)
  
  total_flagged <- sum(sobj$is_qc_outlier)
  cat(sprintf("\nTotal cells flagged for removal: %d (%.2f%%)\n", 
              total_flagged, (total_flagged/ncol(sobj))*100))
  
  return(sobj)
}