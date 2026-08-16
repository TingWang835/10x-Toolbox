import sys
import scanpy as sc

# Redirect standard error and standard output to Snakemake log file
sys.stderr = open(snakemake.log[0], "w")
sys.stdout = sys.stderr

# Loading concatenated AnnData matrix
adata = sc.read_h5ad(snakemake.input["h5ad"])

batch_key = "batch"
method = str(snakemake.params["method"]).lower()
use_rep = "X_pca"

# Ensure batch_key exists in metadata
if batch_key not in adata.obs:
    raise KeyError(
        f"Specified batch_key 'adata.obs[{batch_key}]' not found in adata.obs! "
        f"Available keys: {list(adata.obs.columns)}"
    )


# =============================================================================
# 1. Batch Correction Step
# =============================================================================
if method == "harmony":
    import harmonypy as hm

    # harmonypy accepts PCs x Cells matrix shape (n_obs, n_pcs)
    pca_mat = adata.obsm["X_pca"]
    ho = hm.run_harmony(
        data_mat=pca_mat,
        meta_data=adata.obs,
        vars_use=[batch_key],
        max_iter_harmony=10,
        random_state=0,
    )
    
    # Store corrected PCs
    adjusted_pcs = ho.Z_corr
    if adjusted_pcs.shape[0] != adata.n_obs:
        adjusted_pcs = adjusted_pcs.T

    adata.obsm["X_pca_harmony"] = adjusted_pcs
    use_rep = "X_pca_harmony"

elif method == "bbknn":
    import bbknn

    # BBKNN constructs graph directly on PCA representation
    bbknn.bbknn(adata, batch_key=batch_key, use_rep="X_pca")
    use_rep = "X_pca"

elif method == "scvi":
    import scvi
    import torch

    # Prepare scVI object using unnormalized counts layer
    scvi_adata = adata.copy()
    if "counts" in adata.layers:
        scvi_adata.X = adata.layers["counts"].copy()

    # Subset features to HVGs for faster and robust VAE training
    if "highly_variable" in adata.var:
        scvi_adata = scvi_adata[:, adata.var["highly_variable"]].copy()

    scvi.model.SCVI.setup_anndata(scvi_adata, batch_key=batch_key)
    model = scvi.model.SCVI(
        scvi_adata, n_latent=30, n_layers=2, gene_likelihood="nb"
    )

    use_gpu = torch.cuda.is_available()
    print(f"PyTorch CUDA Available: {use_gpu}")

    model.train(
        max_epochs=200,
        early_stopping=True,
        accelerator="gpu" if use_gpu else "cpu",
        devices=1 if use_gpu else None,
    )

    # Store latent representation
    adata.obsm["X_scVI"] = model.get_latent_representation(scvi_adata)
    use_rep = "X_scVI"

else:
    raise ValueError(f"Unsupported batch correction method: {method}")

# =============================================================================
# 2. Neighbors & Spatial Embeddings Step
# =============================================================================
# BBKNN constructs its own neighbor graph internally during integration
if method != "bbknn":
    print(f"Computing neighbor graph using representation '{use_rep}'...")
    sc.pp.neighbors(
        adata,
        use_rep=use_rep,
        n_neighbors=int(snakemake.params["n_neighbors"]),
    )

embed = str(snakemake.params["embed"]).lower()
min_dist = float(snakemake.params["min_dist"])

if embed in ["umap", "umap&tsne"]:
    print("Calculating UMAP coordinates...")
    sc.tl.umap(adata, min_dist=min_dist)

if embed in ["tsne", "umap&tsne"]:
    print("Calculating t-SNE coordinates...")
    sc.tl.tsne(adata, use_rep=use_rep, n_jobs=snakemake.threads)

# =============================================================================
# 3. Save Final Integrated AnnData Target
# =============================================================================
print(f"Writing integrated AnnData object to {snakemake.output['embedding_h5ad']}...")
adata.write_h5ad(snakemake.output["embedding_h5ad"])
print("Batch correction and spatial embeddings calculation complete!")