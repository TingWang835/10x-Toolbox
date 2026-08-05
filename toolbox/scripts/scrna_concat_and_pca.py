import sys
import numpy as np
import scanpy as sc
import anndata as ad

# Redirect stderr/stdout to log file
sys.stderr = open(snakemake.log[0], "w")
sys.stdout = sys.stderr

input_files = snakemake.input["h5ads"]
print(f"Loading {len(input_files)} samples sequentially for concatenation...")

# Extract sample keys
keys = [path.split("/")[-1].replace("_normalized.h5ad", "") for path in input_files]

# Option A: Generator passed directly to concat (memory efficient)
combined_adata = ad.concat(
    (sc.read_h5ad(path) for path in input_files),
    keys=keys,
    label="batch",
    index_unique="-",
    merge="unique"
)

print(f"Concatenated Dataset: {combined_adata.n_obs} cells across {len(keys)} samples.")

# -------------------------------------------------------------
# Rest of the PCA pipeline (HVGs + Sparse Scaling)
# -------------------------------------------------------------

n_top = int(snakemake.params["n_top_genes"])
sc.pp.highly_variable_genes(
    combined_adata, 
    n_top_genes=n_top, 
    batch_key="batch", 
    flavor=snakemake.params["flavor"]
)

# Save unscaled normalized data in raw slot for downstream differential expression
combined_adata.raw = combined_adata

# Filter to highly variable genes for PCA computation
combined_adata_hvg = combined_adata[:, combined_adata.var["highly_variable"]].copy()

# Zero-center=False maintains sparse matrix representation
z_cap = float(snakemake.params["z_cap"]) if snakemake.params["z_cap"] else None
sc.pp.scale(combined_adata_hvg, max_value=z_cap, zero_center=False)

# Run PCA on scaled HVG matrix
n_pcs = int(snakemake.params["n_pc"])
sc.tl.pca(
    combined_adata_hvg, 
    svd_solver=snakemake.params["svd_solver"], 
    n_comps=n_pcs
)

# Transfer cell embeddings (obsm) and metadata (uns) back to main object
combined_adata.obsm["X_pca"] = combined_adata_hvg.obsm["X_pca"]
combined_adata.uns["pca"] = combined_adata_hvg.uns["pca"]

# Transfer gene loadings (varm) across all genes
n_total_genes = combined_adata.n_vars
pcs_matrix = np.zeros((n_total_genes, n_pcs))
hvg_mask = combined_adata.var["highly_variable"].values
pcs_matrix[hvg_mask, :] = combined_adata_hvg.varm["PCs"]
combined_adata.varm["PCs"] = pcs_matrix

# Write concatenated object
combined_adata.write_h5ad(snakemake.output["concat_h5ad"])
print("Merged PCA dataset saved successfully!")