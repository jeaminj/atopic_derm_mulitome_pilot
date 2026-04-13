# Note: Due to low # of samples, the authors design is not followed as the model will fail to converge
# Instead a simpler design is used
#

library(Seurat)
library(nebula)
library(dplyr)

# --- Subset to Keratinocytes only ---
sobj_kc <- subset(sobj_clean, Cell.type == "Keratinocytes")

# --- Helper function for NB regression DE (diffxpy translation) ---
run_nb_de_nebula <- function(sobj,
                             comparison,
                             group_col   = "disease_status",
                             subject_col = "ind_id") {
  
  # Extract raw counts
  counts <- GetAssayData(sobj, assay = "RNA", layer = "counts")
  
  # Simplified metadata — disease_status only
  meta <- sobj@meta.data %>%
    dplyr::select(all_of(c(subject_col, group_col))) %>%
    dplyr::mutate(
      disease_status = factor(.data[[group_col]],
                              levels = c(comparison$ident.2, comparison$ident.1))
    )
  
  subject_id <- meta[[subject_col]]
  
  # Simplified design matrix ~ 1 + disease_status only
  design <- model.matrix(~ 1 + disease_status, data = meta)
  
  cat(sprintf("Running NEBULA NB regression: %s vs %s (%d cells, %d genes, %d donors)\n",
              comparison$ident.1, comparison$ident.2,
              ncol(counts), nrow(counts), length(unique(subject_id))))
  
  nebula_out <- nebula(
    counts,
    id      = subject_id,
    pred    = design,
    model   = "NBGMM", # simpler model than NBLMM and more stable with few donors
    verbose = TRUE
  )
  
  # Dynamically detect disease_status coefficient columns
  logfc_col <- grep("logFC.*disease_status", colnames(nebula_out$summary), value = TRUE)
  pval_col  <- grep("^p_.*disease_status",   colnames(nebula_out$summary), value = TRUE)
  se_col    <- grep("^se_.*disease_status",  colnames(nebula_out$summary), value = TRUE)
  
  results <- nebula_out$summary %>%
    dplyr::rename(
      log2FC = all_of(logfc_col),
      pval   = all_of(pval_col),
      se     = all_of(se_col)
    ) %>%
    dplyr::mutate(
      padj      = p.adjust(pval, method = "BH"),
      direction = ifelse(log2FC > 0,
                         paste0("Up in ", comparison$ident.1),
                         paste0("Down in ", comparison$ident.1))
    ) %>%
    dplyr::arrange(padj)
  
  return(results)
}

# --- Run: Healthy vs Lesional (KC only) ---
de_kc_lesional_nebula <- run_nb_de_nebula(
  kc_lesional,
  comparison = list(ident.1 = "AD", ident.2 = "Healthy"),
  group_col  = "disease_status"
)

# --- Summary ---
cat("Lesional DE hits (FDR < 0.05):", sum(de_kc_lesional_nebula$padj < 0.05, na.rm = TRUE), "\n")

# --- Compare nebula vs DESeq2 results ---
compare_methods <- de_kc_lesional %>%
  dplyr::select(gene, log2FC_deseq = log2FC, padj_deseq = padj) %>%
  dplyr::inner_join(
    de_kc_lesional_nebula %>% dplyr::select(gene, log2FC_nebula = log2FC, padj_nebula = padj),
    by = "gene"
  )
