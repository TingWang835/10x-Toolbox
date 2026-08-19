# =============================================================================
# Helper function
# =============================================================================
def get_concat_input(wildcards):
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
rule qc_normalize_hvg:
    """
    Filters cells/genes, normalizes expression, log-transforms, 
    and identifies HVGs in a single in-memory pass.
    """
    input:
        h5ad=f"{READS_DIR}/h5ad/{{aligner}}/grouped/{{sample}}_grouped.h5ad",
    output:
        normalized_h5ad=temp(
            f"{READS_DIR}/h5ad/{{aligner}}/grouped/{{sample}}_normalized.h5ad"
        ),
    log:
        f"{LOG_DIR}/preprocess/qc_norm_hvg/{{aligner}}/{{sample}}.log",
    params:
        min_genes=config.get("QC_MIN_GENES", 200),
        max_genes=config.get("QC_MAX_GENES", 6000),
        min_cells=config.get("QC_MIN_CELLS", 3),
        max_pct_mito=config.get("QC_MAX_PCT_MITO", 15.0),
        target_sum=config.get("NORM_TARGET_SUM", 1e4),
        n_top_genes=config.get("HVG_N_TOP_GENES", 2000),
        flavor=config.get("HVG_FLAVOR", "seurat").lower()
    conda:
        "../env/scrna_h5ad_preprocess.yaml"
    script:
        "scripts/scrna_qc_norm_hvg.py"


rule scrna_concat_and_pca:
    input:
        h5ads=get_concat_input,
    output:
        concat_h5ad=temp(f"{READS_DIR}/h5ad/{{aligner}}/pca_concat.h5ad"),
    params:
        n_pc=config.get("N_PC", 50),
        z_cap=config.get("Z_CAP", 10),
        svd_solver=config.get("SVD_SOLVER", "auto").lower(),
        n_top_genes=config.get("HVG_N_TOP_GENES", 2000),
        flavor=config.get("HVG_FLAVOR", "seurat").lower()
    log:
        f"{LOG_DIR}/preprocess/concat_pca_{{aligner}}.log",
    conda:
        "../env/scrna_h5ad_preprocess.yaml"
    resources:
        mem_mb=15000,
    script:
        "scripts/scrna_concat_and_pca.py"


rule scrna_batch_and_embed:
    """
    Applies batch correction (Harmony/BBKNN/scVI) and computes 
    UMAP/t-SNE embeddings before saving final .h5ad.
    """
    input:
        h5ad=f"{READS_DIR}/h5ad/{{aligner}}/pca_concat.h5ad",
        runinfo = f"{READS_DIR}/runinfo_scrna.csv"
    output:
        embedding_h5ad=f"{READS_DIR}/h5ad/{{aligner}}/pca_concat_batch_embed.h5ad",
    log:
        f"{LOG_DIR}/preprocess/batch_embed_{{aligner}}.log",
    params:
        method=config.get("BATCH_METHOD", "harmony").lower(),
        embed=config.get("EMBED", "umap&tsne"),
        n_neighbors=config.get("UMAP_N_NEIGHBORS", 15),
        min_dist=config.get("UMAP_MIN_DIST", 0.5),
    conda:
        "../env/scrna_h5ad_preprocess.yaml"
    threads: 8
    resources:
        mem_mb=15000,
        gpu=1,
    script:
        "scripts/scrna_batch_embed.py"