import os
import scanpy as sc

filtered_dir = os.path.dirname(snakemake.input["matrix"])

adata = sc.read_10x_mtx(
    filtered_dir,
    var_names="gene_symbols",
    make_unique=True
)

adata.obs["sample_id"] = snakemake.wildcards["run"]
adata.obs["aligner"] = "starsolo"

adata.write_h5ad(snakemake.output["h5ad"])