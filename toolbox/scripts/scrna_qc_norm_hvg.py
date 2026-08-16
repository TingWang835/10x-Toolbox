import sys
import scanpy as sc

# Redirect standard error to Snakemake log file
sys.stderr = open(snakemake.log[0], "w")

# 1. Read AnnData
adata = sc.read_h5ad(snakemake.input["h5ad"])

# 2. Extract QC attributes
if "chromosome" not in adata.var.columns:
    adata.var["chromosome"] = "unknown"

# Flag mitochondrial features (checking chromosome tags and common C. elegans prefixes)
mito_chroms = {"chrM", "MT", "M", "MtDNA"}
adata.var["mt"] = adata.var["chromosome"].isin(mito_chroms) | adata.var["gene_symbol"].str.startswith(("ctc-", "nd-", "cox-", "atp-", "cyb-"))

# Calculate cell and gene QC metrics
sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
)

# 3. Perform Quality Filtering
sc.pp.filter_cells(adata, min_genes=snakemake.params["min_genes"])
sc.pp.filter_genes(adata, min_cells=snakemake.params["min_cells"])

# Slice and create an explicit copy to prevent AnnData view warnings
valid_cells = (
    (adata.obs["n_genes_by_counts"] < snakemake.params["max_genes"])
    & (adata.obs["pct_counts_mt"] < snakemake.params["max_pct_mito"])
)
adata = adata[valid_cells, :].copy()

# 4. Preserve Raw Unnormalized Counts
# Store sparse expression array in layers before log-normalization
adata.layers["counts"] = adata.X.copy()

# 5. Normalization & Log Transformation
sc.pp.normalize_total(adata, target_sum=float(snakemake.params["target_sum"]))
sc.pp.log1p(adata)

# 6. Highly Variable Gene (HVG) Selection
sc.pp.highly_variable_genes(
    adata,
    n_top_genes=int(snakemake.params["n_top_genes"]),
    flavor=str(snakemake.params["flavor"]).lower(),
    subset=False  # Keep all genes in matrix; flags HVGs in adata.var['highly_variable']
)

# 7. Write final processed AnnData target
adata.write_h5ad(snakemake.output["normalized_h5ad"])