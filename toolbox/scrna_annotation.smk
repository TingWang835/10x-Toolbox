rule scrna_annotation_cals:
    """
    Performs Leiden clustering, scores cluster markers against embryonic 
    lineages using ScType-style gene sets, maps cell types, and outputs annotated h5ad.
    """
    input:
        h5ad=f"{READS_DIR}/h5ad/{{aligner}}/pca_concat_batch_embed.h5ad",
        markers_csv=config.get("ANNOTATION_MARKER", "databases/annotation_marker/celegans_embryo_markers.csv")
    output:
        h5ad_annotated=f"{READS_DIR}/annotation/{{aligner}}/annotated.h5ad",
        cluster_scores_csv=f"{READS_DIR}/annotation/{{aligner}}/reports/cluster_scores.csv",
        markers_csv=f"{READS_DIR}/annotation/{{aligner}}/reports/marker_genes.csv"
    log:
        f"{LOG_DIR}/annotation/scrna_annotation_{{aligner}}.log"
    params:
        resolution=config.get("RESOLUTION", 0.8),
        cluster_key=config.get("CLUSTER_KEY", "leiden"),
        score_threshold=config.get("SCORE_THRESHOLD", 0.05)
    conda:
        "../env/scrna_annotation.yaml"
    threads: 8
    script:
        "scripts/scrna_annotation_sctype.py"


rule scrna_annotation_plots:
    """
    Generates standalone UMAP, t-SNE, and ordered marker dotplots 
    from the annotated .h5ad object and precomputed marker genes CSV.
    """
    input:
        h5ad=f"{READS_DIR}/annotation/{{aligner}}/annotated.h5ad",
        markers_csv=f"{READS_DIR}/annotation/{{aligner}}/reports/marker_genes.csv"
    output:
        plots_dir=directory(f"{READS_DIR}/annotation/{{aligner}}/plots")
    log:
        f"{LOG_DIR}/annotation/scrna_annotation_plots_{{aligner}}.log"
    params:
        cluster_key=config.get("CLUSTER_KEY", "leiden"),
        top_n_markers=config.get("TOP_N_MARKERS", 3),
        embed=config.get("EMBED", "umap&tsne"),
        umap_tsne_width = config.get("UMAP_TSNE_WIDTH", 12),
        umap_tsne_tall = config.get("UMAP_TSNE_TALL", 8),
    conda:
        "../env/scrna_annotation.yaml"
    threads: 2
    script:
        "scripts/scrna_annotation_plots.py"