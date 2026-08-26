import pandas as pd
import scanpy as sc
from pathlib import Path
import shutil

# Ingestion mode loaded from config
INGEST_MODE = config.get("POST_FASTQ_INGEST", "concat_h5ad").lower()

# -----------------------------------------------------------------------------
# Path 1: Pre-concatenated / Batched H5AD
# -----------------------------------------------------------------------------
rule ingest_concat_h5ad:
    input:
        raw_concat = f"{READS_DIR}/rawdata/pca_concat_batch_embed.h5ad"
    output:
        h5ad = f"{READS_DIR}/h5ad/{ALIGNER}/pca_concat_batch_embed.h5ad",
        csv = f"{READS_DIR}/h5ad/{ALIGNER}/local_concat_runinfo.csv"
    run:
        Path(output.h5ad).parent.mkdir(parents=True, exist_ok=True)
        # 1. Copy file to standard workflow location
        shutil.copyfile(input.raw_concat, output.h5ad)
        
        # 2. Extract quick cell count runinfo from metadata
        adata = sc.read_h5ad(output.h5ad, backed="r")
        n_cells = adata.n_obs
        n_genes = adata.n_vars
        batches = adata.obs["batch"].unique().tolist() if "batch" in adata.obs else ["Combined"]
        
        run_df = pd.DataFrame({
            "Run": batches,
            "SampleName": batches,
            "Condition": "Batched_Input",
            "CellCount": n_cells // len(batches),
            "GeneCount": n_genes,
            "Platform": "H5AD_CONCAT"
        })
        run_df.to_csv(output.csv, index=False)
        adata.file.close()

# -----------------------------------------------------------------------------
# Path 2: Individual H5AD Files
# -----------------------------------------------------------------------------
rule ingest_indiv_h5ad:
    input:
        h5ad_files = [str(p) for p in Path(f"{READS_DIR}/rawdata").glob("*.h5ad")]
    output:
        csv = f"{READS_DIR}/h5ad/{ALIGNER}/local_indiv_runinfo.csv"
    run:
        dest_dir = Path(f"{READS_DIR}/h5ad/{ALIGNER}/individual")
        dest_dir.mkdir(parents=True, exist_ok=True)
        runinfo_data = []
        
        for file_path in input.h5ad_files:
            p = Path(file_path)
            sample_id = p.stem
            
            # Copy to individual target folder
            dest = dest_dir / p.name
            shutil.copyfile(p, dest)
            
            # Peek at metadata efficiently
            adata = sc.read_h5ad(p, backed="r")
            condition = adata.obs["condition"].iloc[0] if "condition" in adata.obs else "Unknown"
            
            runinfo_data.append({
                "Run": sample_id,
                "SampleName": sample_id,
                "Condition": condition,
                "CellCount": adata.n_obs,
                "GeneCount": adata.n_vars,
                "Platform": "H5AD_INDIVIDUAL"
            })
            adata.file.close()
            
        pd.DataFrame(runinfo_data).to_csv(output.csv, index=False)

# -----------------------------------------------------------------------------
# Path 3: Raw H5 / MTX Conversion to H5AD
# -----------------------------------------------------------------------------
rule ingest_raw_h5_mtx:
    input:
        raw_dir = f"{READS_DIR}/rawdata"
    output:
        csv = f"{READS_DIR}/h5ad/{ALIGNER}/local_h5mtx_runinfo.csv"
    run:
        out_indiv = Path(f"{READS_DIR}/h5ad/{ALIGNER}/individual")
        out_indiv.mkdir(parents=True, exist_ok=True)
        
        runinfo_data = []
        
        def convert_to_h5ad(src_path, sample_id):
            if str(src_path).endswith(".h5"):
                adata = sc.read_10x_h5(src_path)
            else:
                adata = sc.read_10x_mtx(src_path, make_unique=True)
                
            out_file = out_indiv / f"{sample_id}.h5ad"
            adata.write_h5ad(out_file, compression="gzip")
            
            return {
                "Run": sample_id,
                "SampleName": sample_id,
                "Condition": "Unknown",
                "CellCount": adata.n_obs,
                "GeneCount": adata.n_vars,
                "Platform": "RAW_H5_MTX"
            }

        # Convert MTX directories
        for mtx in Path(input.raw_dir).rglob("matrix.mtx*"):
            sample_dir = mtx.parent
            sample_id = sample_dir.name if sample_dir.name != "filtered_feature_bc_matrix" else sample_dir.parent.name
            runinfo_data.append(convert_to_h5ad(sample_dir, sample_id))
            
        # Convert standalone .h5 files
        for h5 in Path(input.raw_dir).glob("*.h5"):
            runinfo_data.append(convert_to_h5ad(h5, h5.stem))
            
        pd.DataFrame(runinfo_data).to_csv(output.csv, index=False)

# -----------------------------------------------------------------------------
# Dynamic Checkpoint & Router to Unify Output Files
# -----------------------------------------------------------------------------
def get_post_fastq_runinfo(wildcards):
    if INGEST_MODE == "concat_h5ad":
        return f"{READS_DIR}/h5ad/{ALIGNER}/local_concat_runinfo.csv"
    elif INGEST_MODE == "indiv_h5ad":
        return f"{READS_DIR}/h5ad/{ALIGNER}/local_indiv_runinfo.csv"
    elif INGEST_MODE == "raw_h5_mtx":
        return f"{READS_DIR}/h5ad/{ALIGNER}/local_h5mtx_runinfo.csv"
    else:
        raise ValueError(f"Invalid POST_FASTQ_INGEST option: '{INGEST_MODE}'")

rule standardize_post_fastq_runinfo:
    input:
        get_post_fastq_runinfo
    output:
        f"{READS_DIR}/runinfo_scrna.csv"
    shell:
        """
        mv {input} {output}
        """