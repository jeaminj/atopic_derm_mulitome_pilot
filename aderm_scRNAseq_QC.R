# Script [2]: QC sanity (processed data downloaded) - atopic dermatitis-mutliome project
library(ggplot2)
library(Seurat)
library(SeuratObject)

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

sample_dir="~/projects/aderm_multiome/samples_to_use"
wd="~/projects/aderm_multiome/01_qc"
setwd(wd)

# 1. Identify and name the sample files
file_path_list <- list.files(path = sample_dir, pattern = "\\.rds$", full.names = TRUE)
names(file_path_list) <- gsub("\\.rds$", "", basename(file_path_list))

# check for sanity
print(file_path_list)

# 2. Read into a named list
sobj_list <- lapply(file_path_list, readRDS)

# 3. Merge samples
merged_sobj <- merge(
  x = sobj_list[[1]], # We use the first object as the base
  y = sobj_list[-1],  # and add the rest
  add.cell.ids = names(sobj_list), 
  project = "aderm_sc"
)

# 4. QC Visualizations ~~~~~~~~~~~~
# Sanity check on the state of the sample data and level of processing done
authors_qc_vlns <- VlnPlot(merged_sobj, 
                           features = c("nCount_RNA", "nFeature_RNA", "mt_frac"),
                           ncol = 3,
                           pt.size = 0.2, # Usually cleaner for many samples
                           layer = "counts")

ggsave(file.path(wd, "01a_authors_qc_vlns.png"), authors_qc_vlns, 
       width = 10, height = 7, dpi = 300)

# nCount vs nFeature should show strong positive correlation
p1 <- FeatureScatter(merged_sobj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  theme(legend.position = "none")

# nCount vs mt_frac ensures high-mito cells aren't outliers in count depth
p2 <- FeatureScatter(merged_sobj, feature1 = "nCount_RNA", feature2 = "mt_frac") +
  theme(legend.position = "none")

correlation_plots <- p1 + p2
ggsave(file.path(wd, "01b_qc_correlations.png"), correlation_plots, width = 12, height = 6)


cell_counts <- as.data.frame(table(merged_sobj@meta.data$orig.ident))
cell_ids <- sapply(strsplit(Cells(merged_sobj), "_"), `[`, 1)
cell_counts <- as.data.frame(table(cell_ids))
p3 <- ggplot(cell_counts, aes(x = cell_ids, y = Freq, fill = cell_ids)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Cell Counts per Sample", x = "Sample", y = "Number of Cells") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
ggsave(file.path(wd, "01c_cell_counts_per_sample.png"), p3, width = 8, height = 6)

# Stacked Bar: Cell Type distribution per Sample
metadata <- merged_sobj@meta.data

p_ct_stack <- ggplot(metadata, aes(x = ind_id, fill = Cell.type)) +
  geom_bar(position = "fill") + 
  scale_fill_manual(values = cell_type_cols) +
  theme_minimal() +
  labs(title = "Cell Type Proportion per Sample", y = "Proportion") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("01d_celltype_proportions_custom_colors.png", p_ct_stack, width = 12, height = 7)


# Some more sanity checking is performed after normalization and reducing dimensions, but
# this is a good point to stop and save to file the now merged Seurat object:
SaveSeuratRds(merged_sobj, paste0(wd,"/","sobj_merged_filtered.rds"))


# See script [3] for continuation ~~~
