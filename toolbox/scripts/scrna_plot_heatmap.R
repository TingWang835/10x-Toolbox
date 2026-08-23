# Redirect output and errors to log
log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file, type = "output")
sink(log_file, type = "message")

library(DESeq2)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(readr)
library(tidyr)
library(tibble)

# 1. Load Inputs & Parameters
counts_file    <- snakemake@input[["counts_csv"]]
meta_file      <- snakemake@input[["metadata_csv"]]
degs_file      <- snakemake@input[["deseq2_filtered_csv"]]
output_pdf     <- snakemake@output[["heatmap_pdf"]]
top_n          <- as.numeric(snakemake@params[["top_n_degs"]])
raw_conditions <- snakemake@params[["valid_conditions"]]

# Parse valid_conditions
if (length(raw_conditions) == 1 && grepl("[,\\s]", raw_conditions)) {
  valid_conditions <- unlist(strsplit(raw_conditions, "[\\s,]+"))
} else {
  valid_conditions <- as.character(raw_conditions)
}

counts <- read.csv(counts_file, row.names = 1, check.names = FALSE)
meta   <- read.csv(meta_file, check.names = FALSE)
degs   <- read.csv(degs_file)

# Filter metadata and counts by valid conditions
meta <- meta %>%
  filter(
    pseudobulk_id %in% colnames(counts),
    condition %in% valid_conditions
  )
rownames(meta) <- meta$pseudobulk_id
counts <- counts[, rownames(meta)]

# 2. Extract Top N DEGs per Cell Type
top_genes <- degs %>%
  filter(!is.na(padj)) %>%
  group_by(cell_type) %>%
  slice_min(order_by = padj, n = top_n, with_ties = FALSE) %>%
  pull(gene) %>%
  unique()

if (length(top_genes) == 0) {
  warning("No significant DEGs found across cell types to plot. Creating empty placeholder PDF.")
  pdf(output_pdf)
  plot.new()
  text(0.5, 0.5, "No Significant DEGs Found")
  dev.off()
  sink(type = "message")
  sink(type = "output")
  q(save = "no")
}

# 3. Compute VST Transformation Globally
dds <- DESeqDataSetFromMatrix(
  countData = round(counts),
  colData   = meta,
  design    = ~ condition
)
dds <- dds[rowSums(counts(dds)) >= 10, ]
vst_mat <- assay(vst(dds, blind = FALSE))

# 4. Collapse Technical Replicates by (cell_type, condition)
vst_sub <- vst_mat[rownames(vst_mat) %in% top_genes, , drop = FALSE]
vst_df  <- as.data.frame(vst_sub)
vst_df$gene <- rownames(vst_df)

vst_long <- vst_df %>%
  pivot_longer(-gene, names_to = "pseudobulk_id", values_to = "vst_val") %>%
  left_join(meta %>% select(pseudobulk_id, cell_type, condition), by = "pseudobulk_id") %>%
  group_by(gene, cell_type, condition) %>%
  summarize(mean_vst = mean(vst_val), .groups = "drop") %>%
  mutate(group_id = paste(cell_type, condition, sep = "___"))

# Reconstruct aggregated matrix
mean_vst_mat <- vst_long %>%
  select(gene, group_id, mean_vst) %>%
  pivot_wider(names_from = group_id, values_from = mean_vst) %>%
  column_to_rownames("gene") %>%
  as.matrix()

# 5. Build Column Annotations & Color Scales
col_meta <- data.frame(group_id = colnames(mean_vst_mat)) %>%
  separate(group_id, into = c("cell_type", "condition"), sep = "___", remove = FALSE)

unique_conds   <- unique(col_meta$condition)
palette_colors <- c("#2b5c8f", "#d95f02", "#7570b3", "#e7298a", "#17b7dfd8","#7de28e02", "#3605bb9f", "#d1b410")
# palette for marking conditions
cond_colors    <- setNames(palette_colors[seq_along(unique_conds)], unique_conds)

col_annot <- HeatmapAnnotation(
  `Cell Type` = col_meta$cell_type,
  Condition   = col_meta$condition,
  col = list(Condition = cond_colors),
  annotation_name_side = "left"
)

# Row Z-score normalization
z_mat <- t(scale(t(mean_vst_mat)))
z_mat[is.na(z_mat)] <- 0

col_fun <- colorRamp2(c(-2, 0, 2), c("#fde725ff", "#2a788eff", "#440154ff"))

# 6. Render Grouped Heatmap
pdf(output_pdf, width = 12, height = 10)

ht <- Heatmap(
  z_mat,
  name                  = "Z-Score",
  col                   = col_fun,
  top_annotation        = col_annot,
  column_split          = col_meta$cell_type,
  cluster_columns       = FALSE,
  cluster_rows          = TRUE,
  show_column_names     = FALSE,
  show_row_names        = length(top_genes) <= 100,
  row_names_gp          = gpar(fontsize = 8),
  column_title_gp       = gpar(fontsize = 8, fontface = "bold"),
  column_title_rot      = 45
)

draw(ht, merge_legend = TRUE)
dev.off()

sink(type = "message")
sink(type = "output")