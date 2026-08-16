import os
import sys
import gffutils
import pandas as pd
import scanpy as sc

aligner = snakemake.params.get("aligner", "").lower()

# 1. Load matrix according to aligner type
if aligner in ["starsolo", "star"]:
    matrix_dir = os.path.dirname(snakemake.input.matrix)

    adata = sc.read_10x_mtx(
        matrix_dir,
        var_names="gene_ids",
        cache=False,
    )
    adata.obs_names_make_unique()
    adata.var_names_make_unique()

elif aligner in ["kb", "kb_python", "kallisto"]:
    adata = sc.read_mtx(snakemake.input.matrix)

    barcodes = pd.read_csv(snakemake.input.barcodes, header=None)[0].values
    genes = pd.read_csv(snakemake.input.genes, header=None)[0].values

    adata.obs_names = barcodes
    adata.var_names = genes

    adata.obs_names_make_unique()
    adata.var_names_make_unique()

else:
    sys.exit(
        f"Error: Invalid or unsupported aligner '{aligner}'. Expected 'starsolo' or 'kb_python'."
    )

# 2. Connect to or build gffutils database safely
if hasattr(snakemake.input, "db") and os.path.exists(snakemake.input.db):
    db = gffutils.FeatureDB(snakemake.input.db)
else:
    gtf_path = getattr(snakemake.input, "gtf", None)
    if not gtf_path:
        sys.exit(
            "Error: Neither a valid 'db' file nor a 'gtf' input file was provided to the rule."
        )

    db_path = f"{gtf_path}.db"

    if os.path.exists(db_path):
        try:
            db = gffutils.FeatureDB(db_path)
        except Exception:
            db = gffutils.create_db(
                gtf_path,
                dbfn=db_path,
                force=True,
                disable_infer_genes=True,
                disable_infer_transcripts=True,
                merge_strategy="merge",
            )
    else:
        db = gffutils.create_db(
            gtf_path,
            dbfn=db_path,
            force=True,
            disable_infer_genes=True,
            disable_infer_transcripts=True,
            merge_strategy="merge",
        )

gene_records = {}
target_gene_types = {"gene", "pseudogene", "ncRNA_gene"}

# 3. Extract gene metadata from DB
for feature in db.all_features():
    g_id = feature.attributes.get("gene_id", [None])[0] or feature.id

    if not isinstance(g_id, str) or not g_id.startswith("WBGene"):
        continue

    is_gene_level = (
        feature.featuretype in target_gene_types
        or feature.featuretype.endswith("_gene")
    )

    if g_id not in gene_records or is_gene_level:
        symbol = (
            feature.attributes.get("gene_name", [None])[0]
            or (
                feature.attributes.get("Name", [None])[0]
                if is_gene_level
                else None
            )
            or feature.attributes.get("locus", [None])[0]
            or g_id
        )

        raw_biotype = (
            feature.attributes.get("gene_biotype", [None])[0]
            or feature.attributes.get("gene_type", [None])[0]
            or feature.featuretype
        )
        biotype = "protein_coding" if raw_biotype == "gene" else raw_biotype

        gene_records[g_id] = {
            "chromosome": feature.seqid,
            "start": feature.start,
            "end": feature.end,
            "strand": feature.strand,
            "gene_symbol": symbol,
            "biotype": biotype,
        }

gene_meta = pd.DataFrame.from_dict(gene_records, orient="index")

# 4. Map metadata into adata.var
for col in ["chromosome", "start", "end", "strand", "gene_symbol", "biotype"]:
    if col in gene_meta.columns:
        adata.var[col] = adata.var_names.map(gene_meta[col])

# Fallbacks for unannotated loci
adata.var["gene_symbol"] = adata.var["gene_symbol"].fillna(
    adata.var_names.to_series()
)
adata.var["biotype"] = adata.var["biotype"].fillna("unannotated")

# 5. Save output
adata.write_h5ad(snakemake.output.h5ad)