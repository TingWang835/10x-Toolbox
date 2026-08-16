import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc

# Redirect stdout/stderr to log file
sys.stderr = open(snakemake.log[0], "w")
sys.stdout = sys.stderr

adata = sc.read_h5ad(snakemake.input["h5ad"])
markers_df = pd.read_csv(snakemake.input["markers_csv"])

resolution = float(snakemake.params["resolution"])
cluster_key = snakemake.params["cluster_key"]
top_n = int(snakemake.params["top_n_markers"])

# 1. Community Detection (Leiden Clustering)
sc.tl.leiden(
    adata, 
    resolution=resolution, 
    key_added=cluster_key
)

# 2. Parse Marker Dictionary & Match Case
# build markers dictionary in cell-type: genes (lower case)
marker_dict = {}
for _, row in markers_df.iterrows():
    cell_type = row["cell_type"].strip()
    genes = [g.strip().lower() for g in str(row["markers"]).split(";") if g.strip()]
    marker_dict[cell_type] = genes

# build gene list from adata.var_names (lower case)
var_lower = {gene.lower(): gene for gene in adata.var_names}

# 3. Calculate ScType Scores across Cells
cell_scores = pd.DataFrame(index=adata.obs_names)

# evaluate genes in markers dictionary exist in gene list, if not: pass
for cell_type, genes in marker_dict.items():
    valid_genes = [var_lower[g] for g in genes if g in var_lower]
    if not valid_genes:
        continue
    
    # Extract matrix slice, extract cells with the matching genes to marker.csv, form np array by filling 0s in empty anndata space
    sub_matrix = adata[:, valid_genes].X
    if not isinstance(sub_matrix, np.ndarray):
        sub_matrix = sub_matrix.toarray()
    
    # Standardize gene expression across cells
    scaled_expr = (sub_matrix - np.mean(sub_matrix, axis=0)) / (np.std(sub_matrix, axis=0) + 1e-6)
    cell_scores[cell_type] = np.mean(scaled_expr, axis=1)

# 4. Aggregate Scores per Leiden Cluster & Map Cell Types
cluster_scores = cell_scores.groupby(adata.obs[cluster_key]).mean()
cluster_scores.to_csv(snakemake.output["cluster_scores_csv"])

# Assign top scoring lineage label to each cluster
cluster_to_celltype = cluster_scores.idxmax(axis=1).to_dict()
adata.obs["cell_type"] = adata.obs[cluster_key].map(cluster_to_celltype).astype("category")

# 5. Compute Differential Markers for Report & DotPlot Validation
print("Computing cluster marker genes (Wilcoxon Rank-Sum)...")
sc.tl.rank_genes_groups(
    adata, 
    groupby="cell_type", 
    method="wilcoxon", 
    key_added="rank_cell_type"
)

# Export full marker table to CSV
marker_df = sc.get.rank_genes_groups_df(adata, group=None, key="rank_cell_type")
marker_df.to_csv(snakemake.output["markers_csv"], index=False)

# Collect top N markers for plot display
groups = adata.obs["cell_type"].cat.categories
result = adata.uns["rank_cell_type"]
plot_markers = []

for group in groups:
    genes = [gene for gene in result["names"][group][:top_n]]
    plot_markers.extend(genes)

unique_plot_markers = list(dict.fromkeys(plot_markers))

# 6. Generate Outputs & Plots
print("Generating diagnostic plots...")

# A. Dual UMAP (Leiden Clusters vs. Assigned Lineages)
fig, axes = plt.subplots(1, 2, figsize=(12, 5))
sc.pl.umap(adata, color=cluster_key, title=f"Leiden Clusters ({cluster_key})", show=False, ax=axes[0])
sc.pl.umap(adata, color="cell_type", title="Annotated Lineages", show=False, ax=axes[1])
plt.tight_layout()
fig.savefig(snakemake.output["umap_plot"], bbox_inches="tight", dpi=300)
plt.close()

# B. Marker Expression Verification DotPlot
dp = sc.pl.dotplot(
    adata,
    var_names=unique_plot_markers,
    groupby="cell_type",
    standard_scale="var",
    figsize=(10, 6),
    show=False
)
dp.savefig(snakemake.output["dotplot"], bbox_inches="tight", dpi=300)
plt.close()

# 7. Write Final Annotated Object
print(f"Saving annotated AnnData to {snakemake.output['h5ad_annotated']}...")
adata.write_h5ad(snakemake.output["h5ad_annotated"])
print("Annotation module completed successfully!")