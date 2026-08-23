import os
import sys
import numpy as np
import pandas as pd
import scanpy as sc

sys.stderr = open(snakemake.log[0], "w")
sys.stdout = sys.stderr

adata = sc.read_h5ad(snakemake.input["h5ad"])

cell_type_key = str(snakemake.params["cell_type_key"])
sample_key = str(snakemake.params["sample_key"])
run_key = str(snakemake.params.get("run_key"))
cond_key = str(snakemake.params["cond_key"])
min_cells = int(snakemake.params["min_cells_per_group"])

os.makedirs(os.path.dirname(snakemake.output["counts_csv"]), exist_ok=True)

# Select unnormalized raw counts
if adata.raw is not None:
    counts_matrix = adata.raw.X
    var_names = adata.raw.var_names
else:
    counts_matrix = adata.X
    var_names = adata.var_names

# Combine sample, run, and cell type into unique pseudobulk group key
# e.g., "GSM3618670___run_0___Seam Cells"
adata.obs["pseudobulk_group"] = (
    adata.obs[sample_key].astype(str)
    + "___"
    + adata.obs[run_key].astype(str)
    + "___"
    + adata.obs[cell_type_key].astype(str)
)

pb_counts = {}
pb_meta = []

# Aggregate raw counts for ALL groups present
for group_id, group_obs in adata.obs.groupby("pseudobulk_group"):
    cell_type = group_obs[cell_type_key].iloc[0]
    sample_id = group_obs[sample_key].iloc[0]
    run_id = group_obs[run_key].iloc[0]
    condition = group_obs[cond_key].iloc[0]

    if pd.isna(cell_type) or cell_type == "Unknown" or len(group_obs) < min_cells:
        continue

    sub_indices = group_obs.index
    sub_matrix = counts_matrix[adata.obs_names.get_indexer(sub_indices)]

    # Matrix sum and flatten
    summed_counts = np.array(sub_matrix.sum(axis=0)).flatten()

    pb_counts[group_id] = summed_counts
    pb_meta.append(
        {
            "pseudobulk_id": group_id,
            "sample_id": sample_id,
            "run_id": run_id,
            "cell_type": cell_type,
            "condition": condition,
            "n_cells": len(group_obs),
        }
    )

counts_df = pd.DataFrame(pb_counts, index=var_names)
metadata_df = pd.DataFrame(pb_meta)

counts_df.to_csv(snakemake.output["counts_csv"])
metadata_df.to_csv(snakemake.output["metadata_csv"], index=False)