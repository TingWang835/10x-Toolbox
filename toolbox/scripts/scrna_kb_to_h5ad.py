import scanpy as sc


# 1. Load the matrix (kb_python outputs cells x genes directly)
adata = sc.read_mtx(snakemake.input["matrix"])

# 2. Read barcode and gene metadata
with open(snakemake.input["barcodes"]) as f:
    barcodes = [line.strip() for line in f]

with open(snakemake.input["genes"]) as f:
    genes = [line.strip() for line in f]

# 3. Assign names directly without transposing
adata.obs_names = barcodes
adata.var_names = genes
adata.var_names_make_unique()

# 4. Add metadata and save
adata.obs["sample_id"] = snakemake.wildcards["run"]
adata.obs["aligner"] = "kb_python"

adata.write_h5ad(snakemake.output["h5ad"])