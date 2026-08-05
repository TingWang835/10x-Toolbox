import sys
import scanpy as sc
import celltypist
from celltypist import models

# Redirect output and error streams to log
sys.stderr = open(snakemake.log[0], "w")
sys.stdout = sys.stderr

print("Loading dataset...")
adata = sc.read_h5ad(snakemake.input["h5ad"])

model_name = snakemake.params["model_name"]
majority_voting = snakemake.params["use_majority_voting"]
cluster_key = snakemake.params["cluster_key"]

# 1. Download/Load CellTypist Model
models.download_models(force_update=False)
model = models.Model.load(model=model_name)

# 2. Run Automated Annotation
# CellTypist expects log1p normalized data (target_sum=10,000)
predictions = celltypist.annotate(
    adata, 
    model=model, 
    majority_voting=majority_voting
)

# Convert predictions back into the main AnnData object
adata = predictions.to_adata(insert_labels=True, insert_conf=True)

# Define primary annotation key based on majority voting setting
annot_key = "majority_voting" if majority_voting and "majority_voting" in adata.obs else "predicted_labels"
adata.obs["auto_cell_type"] = adata.obs[annot_key]


# 3. Output Visualizations

# A. Cell Type UMAP
ax_umap = sc.pl.umap(
    adata, 
    color=["auto_cell_type", cluster_key], 
    title=["Auto Cell Types", f"Clusters ({cluster_key})"],
    show=False
)
ax_umap[0].figure.savefig(snakemake.output["umap_plot"], bbox_inches="tight")

# B. Prediction Confidence Map
ax_conf = sc.pl.umap(
    adata, 
    color="conf_score", 
    title="CellTypist Prediction Confidence",
    cmap="viridis",
    show=False
)
ax_conf.figure.savefig(snakemake.output["probability_plot"], bbox_inches="tight")

# 4. Save Output Matrix
adata.write_h5ad(snakemake.output["h5ad_annotated"])