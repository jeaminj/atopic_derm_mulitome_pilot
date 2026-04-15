# Script [2]: QC - atopic dermatitis-mutliome project
library(Seurat)
library(ggplot2)
library(patchwork)

wd="~/projects/aderm_multiome/02_norm_integ"
setwd(wd)

sobj <- LoadSeuratRds("../01_qc/sobj_merged_filtered.rds")

# ---
# In Seurat5, merged objects store each sample as a separate layer within the RNA assay.
# JoinLayers collapses these into a single counts matrix, which is required
# for downstream steps like NormalizeData and FindVariableFeatures to operate
# across all cells simultaneously rather than per-sample.
sobj_layersjoined <- JoinLayers(sobj,
                                assay = "RNA",
                                layers = "counts")

# Log-normalize counts to account for differences in sequencing depth across cells.
# Each cell's counts are divided by total counts, multiplied by 10,000 (default scale.factor),
# then log1p-transformed. Result stored in the "data" layer.
sobj_layersjoined <- NormalizeData(sobj_layersjoined, 
                                   normalization.method = "LogNormalize", 
                                   verbose = FALSE)

# Uncomment line below and rerun normalization if needed
#options(future.globals.maxSize = 8000 * 1024^2) # increased RAM limit to 8gb

# Identify the top 2000 highly variable features using variance-stabilizing transformation (VST).
# These genes drive biological variability and will be used for PCA and downstream clustering.
sobj_layersjoined <- FindVariableFeatures(sobj_layersjoined, 
                                          selection.method = "vst", 
                                          nfeatures = 2000, verbose = FALSE)

# Inspect the top 15 most variable genes for a quick sanity check
top15 <- head(VariableFeatures(sobj_layersjoined), 15)
cat("Top 15 variable genes:", paste(top15, collapse = ", "), "\n")

# Plot variable features, highlighting the top 15 with labels
top_vFeatures_plot <- VariableFeaturePlot(sobj_layersjoined)
top_vFeatures_plot <- LabelPoints(plot = top_vFeatures_plot, 
                                  points = top15, repel = TRUE,
                                  xnudge = 0.1, ynudge = 0,
                                  max.overlaps = 20)
top_vFeatures_plot
ggsave("02a_top_vFeatures.png", width = 10, height = 8, dpi = 300)

# Scale data to zero mean and unit variance across cells.
# This ensures no single gene dominates PCA due to high expression magnitude.
# By default, only the variable features identified above are scaled.
sobj_layersjoined <- ScaleData(sobj_layersjoined, verbose = FALSE)

# Reduce dimensionality using PCA on the scaled variable features.
sobj_layersjoined <- RunPCA(sobj_layersjoined, npcs = 50, verbose = FALSE)

# Visualize variance explained per PC to determine how many PCs to retain.
elbow_plt <- ElbowPlot(sobj_layersjoined, ndims = 50)

ggsave("02b_elbow_plot.png", elbow_plt, width = 10, height = 7, dpi = 300)

# Embed cells in 2D using UMAP for visualization, using the top 30 PCs as input.
# Adjust dims based on ElbowPlot inspection (e.g., 1:20 if elbow is at ~20).
sobj_layersjoined <- RunUMAP(sobj_layersjoined, dims = 1:30, reduction = "pca", verbose = FALSE)

# checkpoint; save normalized sobj to file
SaveSeuratRds(sobj_layersjoined, "sobj_normalized.rds")

# Uncomment if needed 
#sobj_layersjoined <- LoadSeuratRds("02_norm_integ/sobj_normalized.rds")


# This next section is to assess batch effects and determine if integration is necessary
# --- Consistent color palettes for biologically confounded variables ---

# disease_status & lesional share the same color logic:
# AD / Yes (lesional) = same color; Healthy / No (lesional) = same color
disease_cols  <- c("AD" = "#E05C5C",  "Healthy" = "#4DABB5")
lesional_cols <- c("Yes" = "#E05C5C", "No"      = "#4DABB5")

# Shared colors between sample_name and ind_id based on donor identity
shared_cols <- c(
  # - sample_name -
  "040219_1-H1"    = "#2A5783",
  "040219_2-H2"    = "#4DABB5",  # Healthy2
  "101918H51-H51"  = "#A66BB5",  # Healthy5
  "101818L1-L1"    = "#E6A316",
  "052519L1-L1"    = "#E05C5C",  # BCH04
  "100219L-L"      = "#6B9E3F",  # BCH05
  
  # - ind_id -
  "Healthy1" = "#2A5783",
  "Healthy2" = "#4DABB5",
  "Healthy5" = "#A66BB5",
  "BCH01"    = "#E6A316",
  "BCH05"    = "#6B9E3F",
  "BCH04"    = "#E05C5C"
)
# Remaining variables with independent palettes
batch_cols    <- setNames(scales::hue_pal()(length(unique(sobj_layersjoined$batch))),
                          sort(unique(sobj_layersjoined$batch)))
sample_cols   <- setNames(c("#E07B39","#6B9E3F","#4DABB5","#A66BB5"),
                          sort(unique(sobj_layersjoined$sample_name)))
chem_cols     <- c("v2" = "#E07B54", "v3" = "#45AABF")
phase_cols    <- c("G1" = "#E07B54", "G2M" = "#3A9E6F", "S" = "#7B7BBF")

# --- DimPlots grouped by variables where batch effects could arise ---
p1 <- DimPlot(sobj_layersjoined, group.by = "batch",        cols = batch_cols)   + ggtitle("batch")
p2 <- DimPlot(sobj_layersjoined, group.by = "sample_name", cols = shared_cols) + ggtitle("sample_name")
p3 <- DimPlot(sobj_layersjoined, group.by = "chemistry_10x",cols = chem_cols)    + ggtitle("chemistry_10x")
p4 <- DimPlot(sobj_layersjoined, group.by = "channel")                           + ggtitle("channel")
p5 <- DimPlot(sobj_layersjoined, group.by = "ind_id",      cols = shared_cols) + ggtitle("ind_id")                            + ggtitle("ind_id")
p6 <- DimPlot(sobj_layersjoined, group.by = "phase",        cols = phase_cols)   + ggtitle("phase")
p7 <- DimPlot(sobj_layersjoined, group.by = "sex")                               + ggtitle("sex")
p8 <- DimPlot(sobj_layersjoined, group.by = "disease_status", cols = disease_cols) + ggtitle("disease_status")
p9 <- DimPlot(sobj_layersjoined, group.by = "lesional",       cols = lesional_cols) + ggtitle("lesional")

wrap_plots(list(p1, p2, p3, p4, p5, p6, p7, p8, p9), ncol = 3)

ggsave("02c_umaps_before_integ.png", width = 12, height = 9, dpi = 300)
# ---

# Integration ---------------------
# We will use Harmony
library(harmony)

sobj_integrated <- RunHarmony(
  sobj_layersjoined,
  group.by.vars = "channel",  # replicating how authors batch corrrected
  reduction = "pca",
  reduction.save = "harmony",
  verbose = FALSE
)

# Using harmony reduction instead of pca for downstream steps
sobj_integrated <- RunUMAP(sobj_integrated,
                           reduction = "harmony",
                           dims = 1:50,
                           n.neighbors = 15,        # match k=15 from sc.pp.neighbors
                           reduction.name = "umap_harmony",
                           verbose = FALSE)
# lets check the plots now:
p1 <- DimPlot(sobj_integrated, reduction = "umap_harmony", group.by = "batch",        cols = batch_cols)   + ggtitle("batch")
p2 <- DimPlot(sobj_integrated, reduction = "umap_harmony", group.by = "sample_name", cols = shared_cols) + ggtitle("sample_name")
p3 <- DimPlot(sobj_integrated, reduction = "umap_harmony", group.by = "chemistry_10x",cols = chem_cols)    + ggtitle("chemistry_10x")
p4 <- DimPlot(sobj_integrated, reduction = "umap_harmony", group.by = "channel")                           + ggtitle("channel")
p5 <- DimPlot(sobj_integrated, reduction = "umap_harmony", group.by = "ind_id",      cols = shared_cols) + ggtitle("ind_id")                            + ggtitle("ind_id")
p6 <- DimPlot(sobj_integrated, reduction = "umap_harmony", group.by = "phase",        cols = phase_cols)   + ggtitle("phase")
p7 <- DimPlot(sobj_integrated, reduction = "umap_harmony", group.by = "sex")                               + ggtitle("sex")
p8 <- DimPlot(sobj_integrated, reduction = "umap_harmony", group.by = "disease_status", cols = disease_cols) + ggtitle("disease_status")
p9 <- DimPlot(sobj_integrated, reduction = "umap_harmony", group.by = "lesional",       cols = lesional_cols) + ggtitle("lesional")

wrap_plots(list(p1, p2, p3, p4, p5, p6, p7, p8, p9), ncol = 3)

ggsave("02d_umaps_postIntegration.png", width = 12, height = 9, dpi = 300)

# Once integration is confirmed to have worked, continue with standard workflow
# --- k-NN graph with k=15 on 50 PCs ---
# Equivalent to sc.pp.neighbors(adata, n_neighbors=15, n_pcs=50)
# Note: we use the harmony reduction since integration was already performed
sobj_integrated <- FindNeighbors(sobj_integrated,
                                 reduction = "harmony",
                                 dims = 1:50,        # 50 PCs as per methods
                                 k.param = 15,       # k=15 as per methods
                                 verbose = FALSE)

# --- Leiden clustering ---
# Equivalent to sc.tl.leiden
# Seurat's default is Louvain (algorithm=1); Leiden requires algorithm=4
# install.packages("leidenAlg") if not already installed
set.seed(67)
sobj_integrated <- FindClusters(sobj_integrated,
                                resolution = 0.4,
                                algorithm = 4,       # Leiden
                                random.seed = 67,
                                verbose = FALSE)

SaveSeuratRds(sobj_integrated, "sobj_integrated.rds")

# --- End of Integration. See script [4] for Resolution assessment, annotation, and more plots --- 
