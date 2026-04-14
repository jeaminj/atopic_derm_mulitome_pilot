library(DESeq2)
library(dplyr)

out_dir="~/projects/aderm_multiome/04_downstream_analysis/diff_expression/"

run_pseudobulk_de <- function(sobj,
                              comparison,
                              group_col   = "disease_status",
                              subject_col = "ind_id") {
  
  # --- Aggregate counts per donor (pseudobulk) ---
  counts <- GetAssayData(sobj, assay = "RNA", layer = "counts")
  meta   <- sobj@meta.data
  donors <- unique(meta[[subject_col]])
  
  # Sum raw counts across all cells per donor
  pb_counts <- sapply(donors, function(d) {
    cells <- rownames(meta)[meta[[subject_col]] == d]
    Matrix::rowSums(counts[, cells, drop = FALSE])
  })
  colnames(pb_counts) <- donors
  storage.mode(pb_counts) <- "integer"  # coerce to integer without altering values
  
  # --- Build donor-level metadata ---
  pb_meta <- meta %>%
    dplyr::select(all_of(c(subject_col, group_col))) %>%
    dplyr::group_by(across(all_of(subject_col))) %>%
    dplyr::summarise(
      disease_status = unique(.data[[group_col]]),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      disease_status = factor(disease_status,
                              levels = c(comparison$ident.2, comparison$ident.1))
    ) %>%
    as.data.frame()
  rownames(pb_meta) <- pb_meta[[subject_col]]
  pb_meta <- pb_meta[donors, ]
  
  cat(sprintf("Running pseudobulk NB DE: %s vs %s\n", comparison$ident.1, comparison$ident.2))
  cat(sprintf("  %d donors | %d genes | %d cells aggregated\n",
              length(donors), nrow(pb_counts), ncol(counts)))
  
  # --- DESeq2 (NB + Wald test) ---
  dds <- DESeqDataSetFromMatrix(
    countData = as.matrix(pb_counts),
    colData   = pb_meta,
    design    = ~ disease_status
  )
  
  # Filter lowly expressed genes (min 10 counts across donors)
  dds <- dds[rowSums(counts(dds)) >= 10, ]
  
  # Run DESeq2 (NB model + Wald test)
  dds <- DESeq(dds, test = "Wald", fitType = "parametric", quiet = TRUE)
  
  # Extract results for disease_status contrast
  res <- results(dds,
                 contrast = c("disease_status", comparison$ident.1, comparison$ident.2),
                 alpha    = 0.05) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("gene") %>%
    dplyr::rename(log2FC = log2FoldChange) %>%
    dplyr::mutate(
      direction = ifelse(log2FC > 0,
                         paste0("Up in ", comparison$ident.1),
                         paste0("Down in ", comparison$ident.1))
    ) %>%
    dplyr::arrange(padj)
  
  return(res)
}

# --- Run: Healthy vs Lesional (Keratinocyte cells only) ---
sobj_kc <- subset(sobj_clean, Cell.type %in% c("Cornified keratinocytes", "Keratinocytes"))

de_kc_lesional <- run_pseudobulk_de(
  sobj_kc,
  comparison = list(ident.1 = "AD", ident.2 = "Healthy"),
  group_col  = "disease_status"
)


# --- Quick results summary ---
cat("Lesional DE hits (FDR < 0.05):", sum(de_kc_lesional$padj < 0.05, na.rm = TRUE), "\n")

# --- Save results table to file 

write.table(de_kc_lesional, 
            file = paste0(out_dir, "deseq2_pseudobulk_DE_results_ADvsHealthy_KeratinocytesONLY.tsv"),
            col.names = TRUE,
            row.names = FALSE,
            sep = "\t",
            quote = FALSE)
