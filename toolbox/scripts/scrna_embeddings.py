import sys
import scanpy as sc

# Redirect log output
sys.stderr = open(snakemake.log[0], "w")
sys.stdout = sys.stderr

adata = sc.read_h5ad(snakemake.input["h5ad"])

use_rep = snakemake.params["use_rep"]
embed = str(snakemake.params["embed"]).lower()
n_neighbors = int(snakemake.params["n_neighbors"])
min_dist = float(snakemake.params["min_dist"])


# 1. Compute Nearest Neighbor Graph
sc.pp.neighbors(
    adata, 
    use_rep=use_rep, 
    n_neighbors=n_neighbors, 
    key_added="neighbors"
)

# 2. Compute Requested Embedding(s)
if embed == "umap":
    sc.tl.umap(
        adata, 
        min_dist=min_dist,  
        neighbors_key="neighbors"
    )

elif embed == "tsne":
    sc.tl.tsne(
        adata, 
        use_rep=use_rep, 
        n_jobs=snakemake.threads
    )

else:
    sc.tl.umap(
        adata, 
        min_dist=min_dist,  
        neighbors_key="neighbors"
    )
    sc.tl.tsne(
        adata, 
        use_rep=use_rep, 
        n_jobs=snakemake.threads
    )

# 3. Save updated AnnData with stored spatial embeddings
print(f"Writing updated AnnData to {snakemake.output['embedding_h5ad']}...")
adata.write_h5ad(snakemake.output["embedding_h5ad"])
print("Dimension reduction coordinates calculated and stored successfully!")