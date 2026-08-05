import sys
import scanpy as sc
import numpy as np

# Log stdout/stderr to Snakemake log file
sys.stderr = open(snakemake.log[0], "w")
sys.stdout = sys.stderr

print("Loading merged AnnData matrix...")
adata = sc.read_h5ad(snakemake.input["h5ad"])

batch_key = snakemake.params["batch_key"]
method = snakemake.params["method"]

# Check for batch_key existence
if batch_key not in adata.obs:
    raise KeyError(f"Specified batch_key '{batch_key}' not found in adata.obs! Available keys: {list(adata.obs.columns)}")

print(f"Executing batch correction method '{method}' using batch key '{batch_key}'...")

if method == "harmony":
    import harmonypy as hm

    print("Running Harmony batch correction (direct API)...")
    
    # 1. Run Harmony directly on PCA embeddings
    # harmonypy expects PCs x Cells matrix shape (n_pcs, n_obs)
    pca_mat = adata.obsm["X_pca"].T
    
    ho = hm.run_harmony(
        data_mat=pca_mat,
        meta_data=adata.obs,
        vars_use=[batch_key],
        max_iter_harmony=10,
        random_state=0
    )
    
    # 2. Extract integrated PCs and ensure shape is (n_obs, n_pcs)
    # ho.Z_corr shape in modern harmonypy is (n_obs, n_pcs)
    adjusted_pcs = ho.Z_corr
    if adjusted_pcs.shape[0] != adata.n_obs:
        adjusted_pcs = adjusted_pcs.T

    adata.obsm["X_pca_harmony"] = adjusted_pcs

    # 3. Recompute neighbors and UMAP using Harmony components
    sc.pp.neighbors(adata, use_rep="X_pca_harmony", key_added="harmony_neighbors")
    sc.tl.umap(adata, neighbors_key="harmony_neighbors")

elif method == "bbknn":
    import bbknn
    bbknn.bbknn(adata, batch_key=batch_key, use_rep="X_pca")
    sc.tl.umap(adata)

elif method == "scvi":
    import scvi
    import torch

    use_gpu = torch.cuda.is_available()
    print(f"PyTorch CUDA Available: {use_gpu}")
    if use_gpu:
        print(f"Using GPU device: {torch.cuda.get_device_name(0)}")


    if adata.raw is not None:
        print("Extracting raw count matrix for scVI training...")
        scvi_adata = adata.raw.to_adata()
    else:
        scvi_adata = adata.copy()

    # Retain HVGs for faster and cleaner scVI training
    if "highly_variable" in adata.var:
        scvi_adata = scvi_adata[:, adata.var["highly_variable"]].copy()

    # Setup AnnData for scVI
    scvi.model.SCVI.setup_anndata(scvi_adata, batch_key=batch_key)
    
    model = scvi.model.SCVI(
        scvi_adata, 
        n_latent=30, 
        n_layers=2,
        gene_likelihood="nb"
    )
    
    print("Training scVI Variational Autoencoder...")
    model.train(
        max_epochs=200, 
        early_stopping=True,
        accelerator="gpu" if use_gpu else "cpu",
        devices=1 if use_gpu else None
    )
    
    # Store latent representation back into main adata object
    adata.obsm["X_scVI"] = model.get_latent_representation(scvi_adata)
    sc.pp.neighbors(adata, use_rep="X_scVI", key_added="scvi_neighbors")
    sc.tl.umap(adata, neighbors_key="scvi_neighbors")

else:
    raise ValueError(f"Unsupported batch correction method: {method}")

print(f"Writing batch-corrected dataset to {snakemake.output['batch_h5ad']}...")
adata.write_h5ad(snakemake.output["batch_h5ad"])
print("Batch correction step complete!")