# kickoff script to the helper function that ranks DE results, runs GSEA, and outputs results/plots
# any specifications (gene set, pval_threshold, species) defined in the function call

source("run_gsea.R") # helper script

de_res_dir <- "~/projects/aderm_multiome/04_downstream_analysis/diff_expression/"

# List all DE result files
list_of_de_res <- list.files(de_res_dir, 
                             pattern = "\\.tsv$",  # adjust if different extension
                             full.names = TRUE)

# Preview what was found
cat("DE result files found:\n")
print(basename(list_of_de_res))

# Run GSEA on each file, using the filename (sans extension) as comparison_name
gsea_results <- lapply(list_of_de_res, function(f) {
  
  comparison_name <- tools::file_path_sans_ext(basename(f))
  cat(sprintf("\nProcessing: %s\n", comparison_name))
  
  de_results <- read.table(f, header = TRUE, sep = "\t")
  
  run_gsea(
    de_results      = de_results,
    comparison_name = comparison_name,
    species         = "Homo sapiens",
    gene_sets       = c("C2"),
    subset          = "CP:REACTOME", # subcollection that authors used
    rank_metric     = "log2FC_pval",
    pval_col        = "pvalue",
    log2fc_col      = "log2FC",
    gene_col        = "gene",
    nperm           = 1000,
    pval_threshold  = 0.05,
    out_dir         = "keratinocyte_results/"
  )
})

# todo: 
# run pseudo de and gsea on other celltypes 

