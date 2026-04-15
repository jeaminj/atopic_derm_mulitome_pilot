# TO DO: Functionize this fully to import as helper script

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

# --- Run: Healthy vs Lesional (All Keratinocyte cells) ---
sobj_kc_all <- subset(sobj_clean, Cell.type %in% c("Cornified keratinocytes", "Keratinocytes"))

de_kc_all_lesional <- run_pseudobulk_de(
  sobj_kc_all,
  comparison = list(ident.1 = "AD", ident.2 = "Healthy"),
  group_col  = "disease_status"
)

# --- Quick results summary ---
cat("Lesional DE hits (FDR < 0.05):", sum(de_kc_lesional$padj < 0.05, na.rm = TRUE), "\n")

# --- Save results table to file 

write.table(de_kc_lesional, 
            file = paste0(out_dir, "deseq2_pseudobulk_DE_results_ADvsHealthy_AllKeratinocytes.tsv"),
            col.names = TRUE,
            row.names = FALSE,
            sep = "\t",
            quote = FALSE)


# Run: Healthy vs Lesional (KC subtypes (KC1, KC2, KC3, KC4, KC5, cKC)) 
# Subset to exclude HF and cycling cells
sobj_kc_clean <- subset(sobj_kc_all, 
                        subset = Cell.type.granular.v2 %in% c("KC 1", "KC 2", "KC 3", "KC 4", "KC 5", "Cornified KC 1", "Cornified KC 2"))

# Merge Cornified KC 1 and 2 into one "cKC" group
sobj_kc_clean$subset_group <- as.character(sobj_kc_clean$Cell.type.granular.v2)
sobj_kc_clean$subset_group[sobj_kc_clean$subset_group %in% c("Cornified KC 1", "Cornified KC 2")] <- "cKC"

# Verify the groups
table(sobj_kc_clean$subset_group)

# Define the final targets
target_subsets <- c("KC 1", "KC 2", "KC 3", "KC 4", "KC 5", "cKC")

# Run DE using a list to store results
de_results_list <- list()

for (sub in target_subsets) {
  cat("Processing pseudobulk DE for:", sub, "\n")
  
  # Narrow down to the specific cell state
  current_subset <- subset(sobj_kc_clean, subset = subset_group == sub)
  
  # Run pseudobulk function
  de_results_list[[sub]] <- run_pseudobulk_de(
    current_subset,
    comparison = list(ident.1 = "AD", ident.2 = "Healthy"),
    group_col  = "disease_status"
  )
}

# Loop through the list and write each to a TSV
for (subset_name in names(de_results_list)) {
  
  # Clean filename (replace spaces with underscores)
  clean_name <- gsub(" ", "_", subset_name)
  file_path <- paste0(out_dir, "deseq2_pseudobulk_DE_ADvsHealthy_", clean_name, ".tsv")
  
  cat("Writing:", file_path, "\n")
  
  write.table(de_results_list[[subset_name]], 
              file = file_path,
              col.names = TRUE,
              row.names = FALSE,
              sep = "\t",
              quote = FALSE)
}
