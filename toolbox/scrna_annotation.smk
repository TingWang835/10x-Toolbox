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


rule scrna_pseudobulk_prep:
    """
    Aggregates single-cell raw counts by cell_type and sample ID
    to generate pseudobulk count matrices and metadata for DESeq2/edgeR.
    """
    input:
        h5ad=f"{READS_DIR}/annotation/{{aligner}}/annotated.h5ad"
    output:
        counts_csv=f"{READS_DIR}/annotation/{{aligner}}/pseudobulk/pseudobulk_counts.csv",
        metadata_csv=f"{READS_DIR}/annotation/{{aligner}}/pseudobulk/pseudobulk_metadata.csv"
    log:
        f"{LOG_DIR}/annotation/scrna_pseudobulk_prep_{{aligner}}.log"
    params:
        cell_type_key=config.get("CELL_TYPE_KEY", "cell_type"),
        sample_key=config.get("SAMPLE_KEY", "sample_name"),
        run_key=config.get("RUN_KEY", "run_id"),
        cond_key=config.get("CONDITION_KEY", "condition"),
        min_cells_per_group=config.get("MIN_CELLS_PER_GROUP", 10)
    conda:
        "../env/scrna_annotation.yaml"
    threads: 2
    script:
        "scripts/scrna_pseudobulk_prep.py"


rule scrna_pseudobulk_deseq2:
    """
    Run DESeq2 differential expression.
    """
    input:
        counts_csv=f"{READS_DIR}/annotation/{{aligner}}/pseudobulk/pseudobulk_counts.csv",
        metadata_csv=f"{READS_DIR}/annotation/{{aligner}}/pseudobulk/pseudobulk_metadata.csv"
    output:
        deseq2_full_csv=f"{READS_DIR}/annotation/{{aligner}}/pseudobulk/deseq2_full_results.csv",
        deseq2_filtered_csv=f"{READS_DIR}/annotation/{{aligner}}/pseudobulk/deseq2_filtered_results.csv"
    log:
        f"{LOG_DIR}/annotation/scrna_deseq2_{{aligner}}.log"
    params:
        padj_cutoff=config.get("PADJ_CUTOFF", 0.05),
        log2fc_cutoff=config.get("LOG2FC_CUTOFF", 0.25),
        analysis = config.get("DESEQ_ANALYSIS", "wald").lower(),
        valid_conditions = config.get("COMPARE_CONDITIONS", [])
    conda:
        "../env/scrna_r_env.yaml" 
    threads: 4
    script:
        "scripts/scrna_pseudobulk_deseq2.R"


rule scrna_plot_heatmap:
    input:
        counts_csv=f"{READS_DIR}/annotation/{{aligner}}/pseudobulk/pseudobulk_counts.csv",
        metadata_csv=f"{READS_DIR}/annotation/{{aligner}}/pseudobulk/pseudobulk_metadata.csv",
        deseq2_filtered_csv=f"{READS_DIR}/annotation/{{aligner}}/pseudobulk/deseq2_filtered_results.csv"
    output:
        heatmap_pdf=f"{READS_DIR}/annotation/{{aligner}}/pseudobulk/plots/top_degs_heatmap.pdf"
    params:
        top_n_degs=config.get("TOP_N_DEGS", 10),
        valid_conditions = config.get("COMPARE_CONDITIONS", [])
    conda: 
        "../env/scrna_r_env.yaml"
    log:
        f"{READS_DIR}/logs/annotation/scrna_heatmap_{{aligner}}.log"
    script:
        "scripts/scrna_plot_heatmap.R"


rule scrna_plot_volcano:
    input:
        deseq2_full_csv=f"{READS_DIR}/annotation/{{aligner}}/pseudobulk/deseq2_full_results.csv"
    output:
        volcano_pdf=f"{READS_DIR}/annotation/{{aligner}}/pseudobulk/plots/all_celltypes_volcano.pdf"
    params:
        padj_cutoff=config.get("PADJ_CUTOFF", 0.05),
        log2fc_cutoff=config.get("LOG2FC_CUTOFF", 0.25),
        valid_conditions = config.get("COMPARE_CONDITIONS", []) 
    conda: 
        "../env/scrna_r_env.yaml"
    log:
        f"{READS_DIR}/logs/annotation/scrna_volcano_{{aligner}}.log"
    script:
        "scripts/scrna_plot_volcano.R"


rule scrna_enrichment_go:
    input:
        deseq2_filtered_csv=f"{READS_DIR}/annotation/{{aligner}}/pseudobulk/deseq2_filtered_results.csv",
        counts_csv=f"{READS_DIR}/annotation/{{aligner}}/pseudobulk/pseudobulk_counts.csv"
    output:
        go_csv=f"{READS_DIR}/annotation/{{aligner}}/pseudobulk/go_enrichment_results.csv",
        go_dotplot_pdf=f"{READS_DIR}/annotation/{{aligner}}/pseudobulk/plots/go_enrichment_dotplot.pdf"
    params:
        go_padj_cutoff=config.get("GO_PADJ_CUTOFF", 0.05),
        org_db=config.get("ORG_DB", "org.Ce.eg.db"),
        geneid_type=config.get("GENEID_TYPE", "SYMBOL").upper()
    conda: 
        "../env/scrna_r_env.yaml"
    log:
        f"{READS_DIR}/logs/annotation/scrna_go_{{aligner}}.log"
    script:
        "scripts/scrna_enrichment_go.R"