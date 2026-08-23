# Redirect output and errors to Snakemake log
log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file, type = "output")
sink(log_file, type = "message")

library(DESeq2)
library(dplyr)
library(readr)

# Load input files and parameters
counts_file      <- snakemake@input[["counts_csv"]]
meta_file        <- snakemake@input[["metadata_csv"]]
analysis         <- as.character(snakemake@params[["analysis"]]) 
raw_conditions   <- snakemake@params[["valid_conditions"]]
padj_cutoff      <- as.numeric(snakemake@params[["padj_cutoff"]])
log2fc_cutoff    <- as.numeric(snakemake@params[["log2fc_cutoff"]])

# Separate raw_conditions into correct format
if (length(raw_conditions) == 1 && grepl("[,\\s]", raw_conditions)) {
  valid_conditions <- unlist(strsplit(raw_conditions, "[\\s,]+"))
} else {
  valid_conditions <- as.character(raw_conditions)
}

# 1. Resolve Active Conditions based on Analysis Mode
if (analysis == "wald") {
  if (length(valid_conditions) != 2) {
    stop(paste0(
      "For 'wald' (pair-wise) analysis, valid_conditions must contain exactly 2 elements: [ref_cond, target_cond].\n",
      "Received: ", paste(valid_conditions, collapse = ", ")
    ))
  }
  ref_cond    <- valid_conditions[1]
  target_cond <- valid_conditions[2]
  active_conditions <- c(ref_cond, target_cond)
  
} else if (analysis == "lrt") {
  if (length(valid_conditions) < 3) {
    stop("For 'lrt' analysis, valid_conditions must contain 3 or more time points.")
  }
  active_conditions <- valid_conditions
  
} else {
  stop("Invalid analysis type specified. Must be 'wald' or 'lrt'.")
}

# 2. Read and Align Data
counts <- read.csv(counts_file, row.names = 1, check.names = FALSE)
meta   <- read.csv(meta_file, check.names = FALSE)

meta <- meta[meta$pseudobulk_id %in% colnames(counts), ]
rownames(meta) <- meta$pseudobulk_id

all_results_list <- list()
unique_cell_types <- unique(meta$cell_type)

# 3. Process Per-Cell-Type Differential Expression
for (ct in unique_cell_types) {
  
  # Subset metadata to current cell type and allowed conditions
  sub_meta <- meta %>%
    filter(
      cell_type == ct,
      condition %in% active_conditions
    )
  
  cond_count <- table(sub_meta$condition)
  
  # Logging progress for active cell type
  message(sprintf("\n--- Processing cell type: %s ---", ct))
  message("Sample counts per condition:")
  print(cond_count)
  
  # Ensure sufficient biological replicates (>= 2 samples per condition)
  if (length(cond_count) < length(active_conditions) || any(cond_count < 2)) {
    message(sprintf("Skipping %s: Insufficient replicates per condition (minimum 2 required).", ct))
    next
  }
  
  # Align counts matrix columns with filtered metadata rows
  sub_counts <- counts[, rownames(sub_meta)]
  
  # Set condition factor levels
  sub_meta$condition <- factor(sub_meta$condition, levels = active_conditions)
  
  # Build DESeq2 Object
  dds <- DESeqDataSetFromMatrix(
    countData = round(sub_counts),
    colData   = sub_meta,
    design    = ~ condition
  )
  
  # Filter low-count transcripts across replicates
  dds <- dds[rowSums(counts(dds)) >= 10, ]
  
  # Run Model and Extract Results
  if (analysis == "lrt") {
    dds <- DESeq(dds, test = "LRT", reduced = ~ 1, quiet = TRUE)
    res <- results(dds)
    res_df <- as.data.frame(res)
    res_df$contrast <- paste(active_conditions, collapse = "_vs_")
    
  } else { # "wald"
    dds <- DESeq(dds, test = "Wald", quiet = TRUE)
    res <- results(dds, contrast = c("condition", target_cond, ref_cond))
    res_df <- as.data.frame(res)
    res_df$contrast <- paste0(target_cond, "_vs_", ref_cond)
  }
  
  res_df$gene      <- rownames(res_df)
  res_df$cell_type <- ct
  
  all_results_list[[ct]] <- res_df
}


# 4. Export Combined and Filtered DEG Tables
dir.create(dirname(snakemake@output[["deseq2_full_csv"]]), recursive = TRUE, showWarnings = FALSE)

if (length(all_results_list) > 0) {
  full_results <- bind_rows(all_results_list)
  
  # Filter hits by significance thresholds
  if (analysis == "wald") {
    filtered_results <- full_results %>%
      filter(
        !is.na(padj),
        !is.na(log2FoldChange),
        padj < padj_cutoff,
        abs(log2FoldChange) >= log2fc_cutoff
      )
  } else {
    filtered_results <- full_results %>%
      filter(
        !is.na(padj),
        padj < padj_cutoff
      )
  }
} else {
  warning("No cell types met the criteria (>=2 replicates per condition). Generating empty output files.")
  # Create empty dataframe structure with expected columns
  full_results <- data.frame(
    baseMean=numeric(), log2FoldChange=numeric(), lfcSE=numeric(),
    stat=numeric(), pvalue=numeric(), padj=numeric(),
    contrast=character(), gene=character(), cell_type=character()
  )
  filtered_results <- full_results
}


# Always write both files to satisfy Snakemake output targets
write_csv(full_results, snakemake@output[["deseq2_full_csv"]])
write_csv(filtered_results, snakemake@output[["deseq2_filtered_csv"]])

# Close log sinks
sink(type = "message")
sink(type = "output")