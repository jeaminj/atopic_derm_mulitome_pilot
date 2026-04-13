# will also act as the output dir
wd="~/projects/aderm_multiome/03_cluster_qc"
setwd(wd)

library(Seurat)
library(clustree)
library(ggplot2)
library(patchwork)
source("../helper_scripts/cell_cluster_flagging.R") # will use to clean the final umap 

# Test multiple resolutions to find the most stable clustering
# *** The data was uploaded with annotations already, so this step is not really necessary ***
# Can still be useful
sobj_integrated <- FindClusters(sobj_integrated, 
                                resolution = seq(0.1, 1.0, by = 0.1),
                                verbose = FALSE)
SaveSeuratRds(sobj_clean, "sobj_cleaned.rds")
# Visualize cluster stability across resolutions
clustree(sobj_integrated, prefix = "RNA_snn_res.")
ggsave("clustree_postIntegration.png", height = 10, width = 10, dpi = 300)


# PLOTS --- using the provided cell type annotations ---

# Assigning colors to match colors from reference paper
cell_type_cols <- c(
  "Keratinocytes"           = "#4472C4",
  "Cornified keratinocytes" = "#E07B39",
  "Sweat gland"             = "#1F7A1F",
  "Fibroblasts"             = "#C00000",
  "Pericyte/SMC"            = "#7030A0",
  "VEC"                     = "#843C0C",
  "LEC"                     = "#FF99CC",
  "Melanocytes"             = "#92A400",
  "Schwann"                 = "#00B0B0",
  "DC"                      = "#9DC3E6",
  "Macrophages"             = "#F4B183",
  "Neutrophils"             = "#92D050",
  "Mast"                    = "#FF7F7F",
  "Plasma"                  = "#C9956C",
  "T/NK"                    = "#FFB6C1"
)
# also matching the legend order as the paper has
cell_type_order <- c(
  "Keratinocytes",
  "Cornified keratinocytes",
  "Sweat gland",
  "Fibroblasts",
  "Pericyte/SMC",
  "VEC",
  "LEC",
  "Melanocytes",
  "Schwann",
  "DC",
  "Macrophages",
  "Neutrophils",
  "Mast",
  "Plasma",
  "T/NK"
)

sobj_integrated$Cell.type <- factor(sobj_integrated$Cell.type, levels = cell_type_order)


p_celltypes <- DimPlot(sobj_integrated, 
                       reduction = "umap_harmony", 
                       group.by = "Cell.type", 
                       cols = cell_type_cols, 
                       pt.size = 0.3) +
  # Use the legend title as your "Broad cell states" header
  labs(x = "UMAP1", y = "UMAP2", color = "Broad cell types") +
  # override.aes makes the legend circles larger
  guides(color = guide_legend(title.position = "top", 
                              title.hjust = 0.5,
                              override.aes = list(size = 4))) + 
  theme(plot.title = element_blank(), # Ensures the main plot title is gone
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(arrow = arrow(length = unit(0.3, "cm"), ends = "last")),
        axis.line.y = element_line(arrow = arrow(length = unit(0.3, "cm"), ends = "last")),
        # Adds a little breathing room above the legend
        legend.margin = margin(t = 10)) 

p_celltypes # note cells (namely keratinocytes among its non-dominant cluster)
ggsave("celltype_umap_uncleaned.png", width = 10, height = 8, dpi = 300)

# Going to 'clean' the clustering based on distance from dominant cluster centroid:

# Using the imported helper script:
sobj_flagged <- flag_misassigned_by_umap(sobj_integrated,
                                         reduction = "umap_harmony",
                                         cell_type_col = "Cell.type",
                                         sd_threshold = 2)
# --- Visualize before removing ---
p_celltypes <- DimPlot(sobj_integrated,
                       reduction = "umap_harmony",
                       group.by = "Cell.type", 
                       cols = cell_type_cols, pt.size = 0.3) + 
  ggtitle("Broad cell states") +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        panel.grid = element_blank(), panel.border = element_blank(),
        axis.line.x = element_line(arrow = arrow(length = unit(0.3, "cm"), ends = "last")),
        axis.line.y = element_line(arrow = arrow(length = unit(0.3, "cm"), ends = "last")))

p_removal <- DimPlot(sobj_flagged,
                     reduction = "umap_harmony",
                     group.by = "removal_status",
                     cols = c("To be removed" = "red", "Keep" = "lightgrey"),
                     order = "To be removed", pt.size = 0.3) +
  ggtitle("Cells flagged for removal (all cell types)") +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        panel.grid = element_blank(), panel.border = element_blank(),
        axis.line.x = element_line(arrow = arrow(length = unit(0.3, "cm"), ends = "last")),
        axis.line.y = element_line(arrow = arrow(length = unit(0.3, "cm"), ends = "last")))

p_celltypes | p_removal
ggsave("celltype_umap_flagged.png", width = 15, height = 8, dpi = 300)

# --- Remove after visual confirmation ---
cells_to_keep <- colnames(sobj_integrated)[!sobj_flagged$is_misassigned]
sobj_clean <- subset(sobj_integrated, cells = cells_to_keep)

# save cleaned to file
SaveSeuratRds(sobj_clean, "sobj_cleaned.rds")


# Now to revisualize the cleaned result grouped across our metadata ---- 
# --- Shared theme for all UMAP plots ---
umap_theme <- list(
  labs(x = "UMAP1", y = "UMAP2"),
  guides(color = guide_legend(title.position = "top",
                              title.hjust = 0.5,
                              override.aes = list(size = 4))),
  theme(plot.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(arrow = arrow(length = unit(0.3, "cm"), ends = "last")),
        axis.line.y = element_line(arrow = arrow(length = unit(0.3, "cm"), ends = "last")),
        legend.margin = margin(t = 10))
)

# Seting variable colors
disease_col <- c("Healthy"  = "#449867",
                 "Lesional" = "#eca667")

# --- Cell type annotation ---
p_celltype <- DimPlot(sobj_clean, reduction = "umap_harmony", group.by = "Cell.type",
                      cols = cell_type_cols, pt.size = 0.3) +
  labs(color = "Broad cell types") + umap_theme
ggsave("umap_cleaned_byCellType.png", width = 10, height = 8, dpi = 300)

# --- Biological variables ---
p_disease <- DimPlot(sobj_clean, reduction = "umap_harmony", group.by = "disease_status",
                     cols = disease_cols, pt.size = 0.3) +
  labs(color = "Disease status") + umap_theme

p_lesional <- DimPlot(sobj_clean, reduction = "umap_harmony", group.by = "lesional",
                      cols = lesional_cols, pt.size = 0.3) +
  labs(color = "Lesional") + umap_theme

p_disease_lesional <- DimPlot(sobj_clean, reduction = "umap_harmony", cols = disease_col,
                              group.by = "disease_lesional", pt.size = 0.2) +
  labs(color = "Disease / Lesional") + umap_theme
ggsave("umap_cleaned_byDisease.png", width = 10, height = 8, dpi = 300)

# --- Technical / batch variables ---
p_batch <- DimPlot(sobj_clean, reduction = "umap_harmony", group.by = "batch",
                   pt.size = 0.3) +
  labs(color = "Batch") + umap_theme

p_chemistry <- DimPlot(sobj_clean, reduction = "umap_harmony", group.by = "chemistry_10x",
                       cols = chem_cols, pt.size = 0.3) +
  labs(color = "10x Chemistry") + umap_theme

p_channel <- DimPlot(sobj_clean, reduction = "umap_harmony", group.by = "channel",
                     pt.size = 0.3) +
  labs(color = "Channel") + umap_theme

p_sample <- DimPlot(sobj_clean, reduction = "umap_harmony", group.by = "sample_name",
                    cols = sample_cols, pt.size = 0.3) +
  labs(color = "Sample") + umap_theme

p_donor <- DimPlot(sobj_clean, reduction = "umap_harmony", group.by = "ind_id",
                   cols = shared_cols, pt.size = 0.3) +
  labs(color = "Donor ID") + umap_theme

# --- Covariates ---
p_sex <- DimPlot(sobj_clean, reduction = "umap_harmony", group.by = "sex",
                 pt.size = 0.3) +
  labs(color = "Sex") + umap_theme

p_phase <- DimPlot(sobj_clean, reduction = "umap_harmony", group.by = "phase",
                   cols = phase_cols, pt.size = 0.3) +
  labs(color = "Cell cycle phase") + umap_theme

# --- Combined layout ---
wrap_plots(
  p_disease, p_lesional, p_disease_lesional,
  p_batch, p_chemistry, p_channel, p_sample,
  p_donor, p_sex, p_phase,
  ncol = 3
) +
  plot_annotation(title = "Post-cleaning UMAP overview",
                  theme = theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold")))

ggsave("umap_cleaned_byVariables_overview.png", width = 10, height = 8, dpi = 300 )

# End of script [4], see script [5] for differential expression, trajectories, and other downstream analysis
