import sys
import anndata as ad
import scanpy as sc

# Redirect standard output and error to log
sys.stderr = open(snakemake.log[0], "w")

input_files = snakemake.input["h5ads"]

# Read all run-level AnnData objects
adatas = [sc.read_h5ad(f) for f in input_files]

# Concatenate along the cell axis (obs), joining variables (var) that match
merged_adata = ad.concat(
    adatas,
    axis=0,
    join="outer",
    label="run_id",
    keys=[snakemake.wildcards.get("run", f"run_{i}") for i in range(len(input_files))],
    index_unique="-"
)

# Preserve sample-level metadata
merged_adata.obs["sample_name"] = snakemake.wildcards["sample"]
merged_adata.obs["aligner"] = snakemake.wildcards["aligner"]

# Make variable names unique if duplicate gene symbols exist
merged_adata.var_names_make_unique()

# Write merged AnnData
merged_adata.write_h5ad(snakemake.output["merged_h5ad"])