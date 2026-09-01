import os
import sys
import numpy as np
import pandas as pd
import scanpy as sc



# Redirect stdout/stderr to Snakemake log file
sys.stderr = open(snakemake.log[0], "w")
sys.stdout = sys.stderr

# Load inputs and configuration parameters
adata = sc.read_h5ad(snakemake.input["h5ad"])
markers_df = pd.read_csv(snakemake.input["markers_csv"])

resolution = float(snakemake.params["resolution"])
cluster_key = str(snakemake.params["cluster_key"])
threshold = float(snakemake.params["score_threshold"])

os.makedirs(
    os.path.dirname(snakemake.output["cluster_scores_csv"]), exist_ok=True
)

# =============================================================================
# 3. Create obs["condition"] column based on sample_name
# =============================================================================
# Load runinfo metadata table
runinfo = pd.read_csv(snakemake.input["runinfo"])

# Build fallback lookup mapping for both Run and SampleName identifiers
run_to_cond = dict(zip(runinfo["SampleName"], runinfo["condition"]))

# Assign condition metadata column
adata.obs["condition"] = (
    adata.obs["sample_name"].map(run_to_cond).astype("category")
)

# Verify no NaN values were introduced during mapping
if adata.obs["condition"].isna().any():
    print(
        "Warning: Some sample_ids could not be mapped to a condition in runinfo!"
    )


# =============================================================================
# Assign Condition Metadata
# =============================================================================
# Map SampleName (or Run) to condition
run_to_cond = dict(zip(runinfo["SampleName"], runinfo["condition"]))

# Fallback check to match either sample_name or batch keys in obs
sample_key = "sample_name" if "sample_name" in adata.obs else "batch"

if sample_key in adata.obs:
    adata.obs["condition"] = (
        adata.obs[sample_key].map(run_to_cond).astype("category")
    )
    if adata.obs["condition"].isna().any():
        print("Warning: Unmapped samples detected when assigning conditions!")


# =============================================================================
# 1. Community Detection (Leiden Clustering)
# =============================================================================
print(f"Running Leiden clustering at resolution {resolution}...")
sc.tl.leiden(
    adata,
    resolution=resolution,
    key_added=cluster_key,
    flavor="igraph",
)

# =============================================================================
# 2. Parse Marker Sets & Build Mapping Tables
# =============================================================================
symbol_to_var = {}
if "gene_symbol" in adata.var.columns:
    for var_id, sym in zip(adata.var_names, adata.var["gene_symbol"]):
        if pd.notna(sym):
            symbol_to_var[str(sym).strip().lower()] = var_id
else:
    for var_id in adata.var_names:
        symbol_to_var[str(var_id).strip().lower()] = var_id

marker_dict = {}
for _, row in markers_df.iterrows():
    cell_type = str(row["cell_type"]).strip()
    raw_genes = [
        g.strip().lower() for g in str(row["markers"]).split(";") if g.strip()
    ]
    valid_var_ids = [symbol_to_var[g] for g in raw_genes if g in symbol_to_var]
    if valid_var_ids:
        marker_dict[cell_type] = valid_var_ids

# =============================================================================
# 3. Vectorized ScType Score Calculation
# =============================================================================
cell_scores = pd.DataFrame(index=adata.obs_names)

for cell_type, var_ids in marker_dict.items():
    sub_matrix = adata[:, var_ids].X
    if hasattr(sub_matrix, "toarray"):
        sub_matrix = sub_matrix.toarray()

    means = np.mean(sub_matrix, axis=0)
    stds = np.std(sub_matrix, axis=0) + 1e-6
    scaled_expr = (sub_matrix - means) / stds

    cell_scores[cell_type] = np.mean(scaled_expr, axis=1)

# =============================================================================
# 4. Cluster Aggregation & Annotation Mapping
# =============================================================================
cluster_scores = cell_scores.groupby(adata.obs[cluster_key]).mean()
cluster_scores.to_csv(snakemake.output["cluster_scores_csv"])

max_scores = cluster_scores.max(axis=1)
top_lineages = cluster_scores.idxmax(axis=1)

cluster_to_celltype = {}
for cluster_id, lineage in top_lineages.items():
    if max_scores[cluster_id] >= threshold:
        cluster_to_celltype[cluster_id] = lineage
    else:
        cluster_to_celltype[cluster_id] = "Unknown"

adata.obs["cell_type"] = (
    adata.obs[cluster_key].map(cluster_to_celltype).astype("category")
)


# =============================================================================
# 5. Differential Marker Gene Identification
# =============================================================================
print("Computing lineage marker genes (Wilcoxon Rank-Sum)...")
sc.tl.rank_genes_groups(
    adata, groupby="cell_type", method="wilcoxon", key_added="rank_cell_type"
)

# Export full marker table
marker_df = sc.get.rank_genes_groups_df(adata, group=None, key="rank_cell_type")
marker_df.to_csv(snakemake.output["markers_csv"], index=False)


# =============================================================================
# 6. Save Annotated AnnData Object
# =============================================================================
# Clean string types in main AnnData object and adata.raw directly before saving
for attr in [adata.obs, adata.var]:
    attr.index = attr.index.astype(str)
    for col in attr.columns:
        if isinstance(attr[col].dtype, (pd.StringDtype, pd.ArrowDtype)):
            attr[col] = attr[col].astype(object)

if adata.raw is not None:
    adata.raw.var.index = adata.raw.var.index.astype(str)
    for col in adata.raw.var.columns:
        if isinstance(
            adata.raw.var[col].dtype, (pd.StringDtype, pd.ArrowDtype)
        ):
            adata.raw.var[col] = adata.raw.var[col].astype(object)

print(f"Saving annotated dataset to {snakemake.output['h5ad_annotated']}...")
adata.write_h5ad(snakemake.output["h5ad_annotated"])
print("Cell-type annotation step complete!")