import sys
import anndata as ad
import numpy as np
import scanpy as sc

sys.stderr = open(snakemake.log[0], "w")
sys.stdout = sys.stderr

input_files = snakemake.input["h5ads"]
keys = [path.split("/")[-1].replace("_normalized.h5ad", "") for path in input_files]
# strip the sample id as "keys" from path of normalized h5ad
# this keys will be used in adata.obs["batch"]

# Read first object to retain global .var metadata structure (gene_symbol, chromosome, etc.)
first_adata = sc.read_h5ad(input_files[0])
var_meta = first_adata.var.copy()

# Concatenate datasets across cells (obs)
# merge="only_same" keeps var columns that are identical across all inputs
combined_adata = ad.concat(
    (sc.read_h5ad(path) for path in input_files),
    keys=keys,
    label="batch",
    index_unique="-",
    merge="same",
)

# Restore rich gene metadata if dropped during concat
for col in var_meta.columns:
    if col not in combined_adata.var.columns:
        combined_adata.var[col] = var_meta[col]

# Identify HVGs across batches
n_top = int(snakemake.params["n_top_genes"])
sc.pp.highly_variable_genes(
    combined_adata,
    n_top_genes=n_top,
    batch_key="batch",
    flavor=snakemake.params["flavor"],
)

# Preserve log-normalized expression matrix for downstream differential expression
combined_adata.raw = combined_adata

# Subset to HVGs for PCA computation
combined_adata_hvg = combined_adata[
    :, combined_adata.var["highly_variable"]
].copy()

# Scale without centering to maintain sparse representation
z_cap = float(snakemake.params["z_cap"]) if snakemake.params["z_cap"] else None
sc.pp.scale(combined_adata_hvg, max_value=z_cap, zero_center=False)

# Run PCA
n_pcs = int(snakemake.params["n_pc"])
sc.tl.pca(
    combined_adata_hvg, svd_solver=snakemake.params["svd_solver"], n_comps=n_pcs
)

# Transfer PCA back
combined_adata.obsm["X_pca"] = combined_adata_hvg.obsm["X_pca"]
combined_adata.uns["pca"] = combined_adata_hvg.uns["pca"]

# Reconstruct gene loadings (varm) across total features
n_total_genes = combined_adata.n_vars
pcs_matrix = np.zeros((n_total_genes, n_pcs))
hvg_mask = combined_adata.var["highly_variable"].values
pcs_matrix[hvg_mask, :] = combined_adata_hvg.varm["PCs"]
combined_adata.varm["PCs"] = pcs_matrix

combined_adata.write_h5ad(snakemake.output["concat_h5ad"])