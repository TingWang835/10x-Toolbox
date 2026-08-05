import sys
import scanpy as sc


# Redirect stderr to log file
sys.stderr = open(snakemake.log[0], "w")

adata = sc.read_h5ad(snakemake.input["h5ad"])

# 1. Preserve raw counts layer
adata.layers["counts"] = adata.X.copy()

# 2. Total count normalization & Log1p transformation
sc.pp.normalize_total(adata, target_sum=snakemake.params["target_sum"])
sc.pp.log1p(adata)

# 3. Highly Variable Gene (HVG) selection
sc.pp.highly_variable_genes(
    adata,
    n_top_genes=snakemake.params["n_top_genes"],
    flavor=snakemake.params["flavor"],
    subset=False  # Flags genes in adata.var['highly_variable'] without dropping non-HVGs
)

# 4. Save processed AnnData
adata.write_h5ad(snakemake.output["normalized_h5ad"])