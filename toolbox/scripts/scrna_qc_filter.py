import sys
import gffutils
import scanpy as sc

# Redirect standard error to Snakemake log
sys.stderr = open(snakemake.log[0], "w")

# 1. Load GTF into a gffutils database
db = gffutils.create_db(
    snakemake.input["gtf"],
    dbfn=":memory:",  # Keeps index in RAM for fast execution
    disable_infer_transcripts=True,
    disable_infer_genes=True,
    force=True
)

# 2. Extract WBGene -> Chromosome mapping
chrom_map = {gene.id: gene.seqid for gene in db.features_of_type("gene")}

# 3. Read AnnData and map chromosomes
adata = sc.read_h5ad(snakemake.input["h5ad"])

# Map to adata.var['gene_ids'] (which stores the WBGene IDs)
adata.var["chrom"] = adata.var_names.map(chrom_map).fillna("unknown")

# Flag features originating from specified chromosomes
adata.var["mt"] = adata.var["chrom"].isin(["chrM", "MT", "M"])
# adata.var["chrx"] = adata.var["chrom"] == "chrx" # in case need to filter cells with > y% of chrx genes

# Calculate QC metrics based on flagged features
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

# 4. Perform Quality Control Filtering
sc.pp.filter_cells(adata, min_genes=snakemake.params["min_genes"])
sc.pp.filter_genes(adata, min_cells=snakemake.params["min_cells"])
adata = adata[adata.obs["n_genes_by_counts"] < snakemake.params["max_genes"], :]
adata = adata[adata.obs["pct_counts_mt"] < snakemake.params["max_pct_mito"], :]
# adata = adata[adata.obs["pct_counts_chrx"] < snakemake.params["max_pct_chrx"], :] # in case need to filter cells with  of chrx genes > y%

# Save filtered dataset
adata.write_h5ad(snakemake.output["filtered_h5ad"])