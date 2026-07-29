import scanpy as sc

# Load kb-python matrix and metadata files
adata = sc.read_mtx(snakemake.input["matrix"])

with open(snakemake.input["barcodes"]) as f:
    barcodes = [line.strip() for line in f]

with open(snakemake.input["genes"]) as f:
    genes = [line.strip() for line in f]

# Orient matrix to standard (cells x genes) format
adata = adata.T
adata.obs_names = barcodes
adata.var_names = genes
adata.var_names_make_unique()

adata.obs["sample_id"] = snakemake.wildcards["run"]
adata.obs["aligner"] = "kb_python"

adata.write_h5ad(snakemake.output["h5ad"])