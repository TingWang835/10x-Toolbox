import sys
import pandas as pd
import scanpy as sc

# Redirect standard output and error to log
sys.stderr = open(snakemake.log[0], "w")
sys.stdout = sys.stderr

adata = sc.read_h5ad(snakemake.input["h5ad"])

cluster_key = snakemake.params["cluster_key"]
cluster_map = snakemake.params["cluster_mapping"]
canonical_markers = snakemake.params["canonical_markers"]
top_n = snakemake.params["top_n_markers"]

# 1. Differential Expression Analysis (Wilcoxon Rank-Sum)
print(f"Finding marker genes per cluster using group key: '{cluster_key}'...")
sc.tl.rank_genes_groups(adata, groupby=cluster_key, method="wilcoxon")

# Export top differential markers to CSV for reference
markers_df = pd.DataFrame(adata.uns["rank_genes_groups"]["names"]).head(top_n)
markers_df.to_csv(snakemake.output["marker_csv"], index=False)
print(f"Saved top {top_n} cluster markers to {snakemake.output['marker_csv']}")

# 2. Map Clusters to Manual Cell Type Labels
print("Mapping cluster IDs to biological cell types...")
# Ensure cluster categories match dictionary keys as strings
adata.obs[cluster_key] = adata.obs[cluster_key].astype(str)
adata.obs["cell_type"] = adata.obs[cluster_key].map(cluster_map).fillna("Unknown")

# 3. Generate Annotated UMAP
print("Generating annotated UMAP...")
ax_umap = sc.pl.umap(
    adata, 
    color=["cell_type", cluster_key], 
    title=["Cell Types", f"Clusters ({cluster_key})"],
    show=False
)
ax_umap[0].figure.savefig(snakemake.output["umap_plot"], bbox_inches="tight")

# 4. Generate Marker Verification Dotplot
print("Generating canonical marker dotplot...")
# Filter canonical markers to those actually present in the dataset
valid_markers = [g for g in canonical_markers if g in adata.var_names]

if valid_markers:
    dotplot_dict = sc.pl.dotplot(
        adata, 
        var_names=valid_markers, 
        groupby="cell_type", 
        show=False
    )
    dotplot_dict["mainplot_ax"].figure.savefig(
        snakemake.output["dotplot"], bbox_inches="tight"
    )
else:
    print("Warning: None of the provided canonical markers were found in adata.var_names. Creating empty placeholder.")
    with open(snakemake.output["dotplot"], "w") as f:
        f.write("")

# 5. Save Updated AnnData
print("Saving annotated AnnData matrix...")
adata.write_h5ad(snakemake.output["h5ad_annotated"])
print(f"Successfully finished manual annotation outputting to {snakemake.output['h5ad_annotated']}")