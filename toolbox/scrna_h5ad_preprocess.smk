# =============================================================================
# Helper function
# =============================================================================
def get_concat_and_pca_input(wildcards):
    """prepare input list for Anndata concat and pca."""
    runinfo_df = pd.read_csv(f"{READS_DIR}/runinfo_scrna.csv")
    samples = runinfo_df["SampleName"].astype(str).str.strip().unique().tolist()
    scrna_aligner = config.get("SCRNA_ALIGNER", "starsolo").lower()

    h5ads = expand("{rddir}/h5ad/{aln}/grouped/{s}_normalized.h5ad",
        rddir=READS_DIR, aln=scrna_aligner, s=samples)
    return h5ads


# =============================================================================
# Rules
# =============================================================================
rule h5ad_qc_filter:
    """
    Calculates QC metrics, filters out low-quality cells based on gene counts,
    mitochondrial content, and ribo percentage, then saves the filtered AnnData.
    """
    input:
        unpack(get_refs),
        h5ad = f"{READS_DIR}/h5ad/{{aligner}}/grouped/{{sample}}_grouped.h5ad"
    output:
        filtered_h5ad = temp(f"{READS_DIR}/h5ad/{{aligner}}/grouped/{{sample}}_qc_filtered.h5ad")
    log:
        f"{LOG_DIR}/preprocess/qc_filter/{{aligner}}/{{sample}}.log"
    params:
        min_genes = config.get("QC_MIN_GENES", 200),
        max_genes = config.get("QC_MAX_GENES", 6000),
        min_cells = config.get("QC_MIN_CELLS", 3),
        max_pct_mito = config.get("QC_MAX_PCT_MITO", 15.0),
        # max_pct_chrx = config.get("QC_MAX_PCT_CHRX", 15.0)
    conda:
        "../env/scrna_h5ad_preprocess.yaml"
    script:
        "scripts/scrna_qc_filter.py"


rule normalize_and_hvg:
    """
    Normalizes expression counts, log-transforms the matrix, 
    and identifies highly variable genes (HVGs).
    """
    input:
        h5ad = f"{READS_DIR}/h5ad/{{aligner}}/grouped/{{sample}}_qc_filtered.h5ad"
    output:
        normalized_h5ad = temp(f"{READS_DIR}/h5ad/{{aligner}}/grouped/{{sample}}_normalized.h5ad")
    log:
        f"{LOG_DIR}/preprocess/normalize_hvg/{{aligner}}/{{sample}}.log"
    params:
        target_sum = config.get("NORM_TARGET_SUM", 1e4),
        n_top_genes = config.get("HVG_N_TOP_GENES", 2000),
        flavor = config.get("HVG_FLAVOR", "seurat").lower()
    conda:
        "../env/scrna_h5ad_preprocess.yaml"
    script:
        "scripts/scrna_normalize_and_hvg.py"


rule scrna_concat_and_pca:
    input:
        h5ads = get_concat_and_pca_input
    output:
        concat_h5ad = f"{READS_DIR}/h5ad/{{aligner}}/pca_concat.h5ad"
    params:
        n_pc = config.get("N_PC", 50),
        z_cap  =  config.get("Z_CAP", 10),
        svd_solver = config.get("SVD_SOLVER", "auto").lower(),
        n_top_genes = config.get("HVG_N_TOP_GENES", 2000),
        flavor = config.get("HVG_FLAVOR", "seurat").lower()
    log:
        f"{LOG_DIR}/preprocess/concat_pca_{{aligner}}.log"
    conda:
        "../env/scrna_h5ad_preprocess.yaml"
    resources:
        mem_mb = 15000
    script:
        "scripts/scrna_concat_and_pca.py"


rule scrna_batch_correction:
    """
    Applies selected batch correction algorithm (Harmony, BBKNN, or scVI) 
    to the combined single-cell AnnData matrix.
    """
    input:
        h5ad = f"{READS_DIR}/h5ad/{{aligner}}/pca_concat.h5ad"
    output:
        batch_h5ad = f"{READS_DIR}/h5ad/{{aligner}}/pca_concat_batch.h5ad"
    log:
        f"{LOG_DIR}/preprocess/batch_correction_{{aligner}}.log"
    params:
        method = config.get("BATCH_METHOD", "harmony").lower(),
        batch_key = config.get("BATCH_KEY", "batch")
    conda:
        "../env/scrna_h5ad_preprocess.yaml"
    threads: 8
    resources:
        mem_mb = 15000,
        gpu = 1
    script:
        "scripts/scrna_batch.py"



batch_method = config.get("BATCH_METHOD", "harmony").lower()
rep_map = {
    "harmony": "X_pca_harmony",
    "bbknn": "X_pca",
    "scvi": "X_scVI"
}
use_rep = rep_map.get(batch_method, "X_pca")

rule scrna_embeddings:
    """
    Computes UMAP and t-SNE coordinates on batch-corrected embeddings
    and stores spatial layouts directly inside the AnnData object.
    """
    input:
        h5ad = f"{READS_DIR}/h5ad/{{aligner}}/pca_concat_batch.h5ad"
    output:
        embedding_h5ad = f"{READS_DIR}/h5ad/{{aligner}}/pca_concat_embed.h5ad"
    log:
        f"{LOG_DIR}/preprocess/embeddings_{{aligner}}.log"
    params:
        use_rep = use_rep,
        embed = config.get("EMBED", "umap"),
        n_neighbors = config.get("UMAP_N_NEIGHBORS", 15),
        min_dist = config.get("UMAP_MIN_DIST", 0.5)
    conda:
        "../env/scrna_h5ad_preprocess.yaml"
    threads: 8
    script:
        "scripts/scrna_embeddings.py"