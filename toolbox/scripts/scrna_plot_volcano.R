# Redirect output and errors to log
log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file, type = "output")
sink(log_file, type = "message")

library(ggplot2)
library(dplyr)
library(readr)

# 1. Load Inputs & Parameters
full_degs_file <- snakemake@input[["deseq2_full_csv"]]
output_pdf     <- snakemake@output[["volcano_pdf"]]
padj_cutoff    <- as.numeric(snakemake@params[["padj_cutoff"]])
log2fc_cutoff  <- as.numeric(snakemake@params[["log2fc_cutoff"]])
raw_conditions <- snakemake@params[["valid_conditions"]]

# Parse valid_conditions
if (length(raw_conditions) == 1 && grepl("[,\\s]", raw_conditions)) {
  valid_conditions <- unlist(strsplit(raw_conditions, "[\\s,]+"))
} else {
  valid_conditions <- as.character(raw_conditions)
}

full_degs <- read.csv(full_degs_file)

# Match contrast column against valid conditions (e.g., "400min_vs_300min")
if ("contrast" %in% colnames(full_degs) && nrow(full_degs) > 0) {
  expected_contrasts <- c(
    paste(valid_conditions, collapse = "_vs_"),
    paste(rev(valid_conditions), collapse = "_vs_")
  )
  full_degs <- full_degs %>%
    filter(contrast %in% expected_contrasts | is.na(contrast))
}


if (nrow(full_degs) == 0) {
  warning("Full DEG table is empty after filtering. Generating empty placeholder PDF.")
  pdf(output_pdf)
  plot.new()
  text(0.5, 0.5, "No DEG Data Available for Active Conditions")
  dev.off()
  sink(type = "message")
  sink(type = "output")
  q(save = "no")
}

# 2. Annotate Significance Status
df_plot <- full_degs %>%
  filter(!is.na(pvalue), !is.na(log2FoldChange)) %>%
  mutate(
    padj_clean = ifelse(padj < 1e-50, 1e-50, padj),
    neg_log10_padj = -log10(padj_clean),
    significance = case_when(
      padj < padj_cutoff & log2FoldChange >= log2fc_cutoff  ~ "Up-regulated",
      padj < padj_cutoff & log2FoldChange <= -log2fc_cutoff ~ "Down-regulated",
      TRUE                                                  ~ "Not Significant"
    )
  )

sig_colors <- c(
  "Up-regulated"    = "#d95f02",
  "Down-regulated"  = "#2b5c8f",
  "Not Significant" = "grey70"
)

# 3. Generate Faceted Volcano Plot across Cell Types
p <- ggplot(df_plot, aes(x = log2FoldChange, y = neg_log10_padj, color = significance)) +
  geom_point(alpha = 0.6, size = 1.2) +
  scale_color_manual(values = sig_colors) +
  geom_vline(xintercept = c(-log2fc_cutoff, log2fc_cutoff), linetype = "dashed", color = "black", alpha = 0.5) +
  geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed", color = "black", alpha = 0.5) +
  facet_wrap(~ cell_type, scales = "free_y") +
  labs(
    title = "Pseudobulk Differential Expression across Cell Types",
    subtitle = sprintf("Conditions: %s | Cutoffs: padj < %g, |log2FC| >= %g", 
                       paste(valid_conditions, collapse = " vs "), padj_cutoff, log2fc_cutoff),
    x = "log2 Fold Change",
    y = "-log10 (Adjusted p-value)",
    color = "Status"
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

# 4. Export PDF Plot
pdf(output_pdf, width = 14, height = 11)
print(p)
dev.off()

sink(type = "message")
sink(type = "output")