rule h5ad_merge:
    """
    Merges individual run-level AnnData (.h5ad) files into a single sample or dataset AnnData object.
    Supports outputs from any aligner (starsolo, kb_python, etc.).
    """
    input:
        h5ads = lambda wc: [
            f"{READS_DIR}/h5ad/{wc.aligner}/{run}.h5ad" 
            for run in get_runs_for_sample(wc.sample)
        ]
    output:
        merged_h5ad = f"{READS_DIR}/h5ad_merged/{{aligner}}/{{sample}}.h5ad"
    log:
        f"{LOG_DIR}/h5ad_merge/{{aligner}}/{{sample}}.log"
    conda:
        "../env/anndata_merge.yaml"
    script:
        "/scripts/anndata_merge.py"