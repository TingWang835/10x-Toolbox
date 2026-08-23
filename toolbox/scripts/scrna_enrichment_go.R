# Redirect output and errors to log
log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file, type = "output")
sink(log_file, type = "message")

library(clusterProfiler)
library(org.Ce.eg.db)
library(dplyr)
library(readr)
library(ggplot2)

# 1. Load Inputs & Parameters
degs_file    <- snakemake@input[["deseq2_filtered_csv"]]
counts_file  <- snakemake@input[["counts_csv"]]
output_csv   <- snakemake@output[["go_csv"]]
output_pdf   <- snakemake@output[["go_dotplot_pdf"]]
p_adj_cutoff <- as.numeric(snakemake@params[["go_padj_cutoff"]])
org_db       <- snakemake@params[["org_db"]]
geneid_type  <- snakemake@params[["geneid_type"]]

library(org_db, character.only = TRUE)

degs   <- read.csv(degs_file)
counts <- read.csv(counts_file, row.names = 1, check.names = FALSE)

# Define Universe: Genes passing low-count filter in expression matrix
gene_universe <- rownames(counts)

if (nrow(degs) == 0) {
  warning("No DEGs provided. Creating empty outputs.")
  write_csv(data.frame(), output_csv)
  pdf(output_pdf); plot.new(); text(0.5, 0.5, "No DEGs Available"); dev.off()
  sink(type = "message"); sink(type = "output"); q(save = "no")
}

# 2. Perform GO Enrichment per Cell Type
go_results_list <- list()

for (ct in unique(degs$cell_type)) {
  ct_degs <- degs %>% filter(cell_type == ct) %>% pull(gene) %>% unique()
  
  if (length(ct_degs) < 5) next 
  # Skip cell types with too few DEGs
  
  ego <- enrichGO(
    gene          = ct_degs,
    universe      = gene_universe,
    OrgDb         = get(org_db),
    keyType       = geneid_type,
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = p_adj_cutoff,
    qvalueCutoff  = 0.05
  )
  
  if (!is.null(ego) && nrow(ego@result) > 0) {
    res_df <- as.data.frame(ego) %>%
      filter(p.adjust < p_adj_cutoff) %>%
      mutate(cell_type = ct)
    go_results_list[[ct]] <- res_df
  }
}

# 3. Export Combined Table & Generate Dotplot
if (length(go_results_list) > 0) {
  all_go_df <- bind_rows(go_results_list)
  write_csv(all_go_df, output_csv)
  
  # Parse GeneRatio fraction string into a numeric decimal
  all_go_df <- all_go_df %>%
    mutate(GeneRatio_num = sapply(GeneRatio, function(x) eval(parse(text = x))))
  
  # Select Top 5 GO Terms per Cell Type for plotting
  top_go <- all_go_df %>%
    group_by(cell_type) %>%
    slice_min(order_by = p.adjust, n = 5, with_ties = FALSE)
  
  p <- ggplot(top_go, aes(x = GeneRatio_num, y = reorder(Description, GeneRatio_num), color = p.adjust, size = Count)) +
    geom_point() +
    scale_color_gradient(low = "red", high = "blue") +
    scale_x_continuous(n.breaks = 3) + # Now works correctly on continuous numeric data
    facet_wrap(~ cell_type, scales = "free_y", ncol = 2) +
    labs(
      title = "GO Biological Process Enrichment per Cell Type",
      x = "Gene Ratio",
      y = NULL,
      color = "Adjusted p-value",
      size = "Count"
    ) +
    theme_bw(base_size = 10) +
    theme(
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  pdf(output_pdf, width = 16, height = 12)
  print(p)
  dev.off()
} else {
  write_csv(data.frame(), output_csv)
  pdf(output_pdf); plot.new(); text(0.5, 0.5, "No Enriched GO Terms Found"); dev.off()
}

sink(type = "message")
sink(type = "output")