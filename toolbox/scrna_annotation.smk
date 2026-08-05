rule scrna_manual_annotation:
    """
    Find marker genes per cluster, map clusters to cell types,
    and output annotated AnnData and diagnostic plots.
    """
    input:
        h5ad = f"{READS_DIR}/h5ad/{{aligner}}/pca_concat_batch.h5ad"
    output:
        h5ad_annotated = f"{READS_DIR}/annotation/h5ad/{{aligner}}_annotated.h5ad",
        umap_plot = f"{READS_DIR}/annotation/plots/{{aligner}}_celltypes_umap.png",
        dotplot = f"{READS_DIR}/annotation/plots/{{aligner}}_markers_dotplot.png",
        marker_csv = f"{READS_DIR}/annotation/tables/{{aligner}}_cluster_markers.csv"
    log:
        f"{LOG_DIR}/annotation/manual_annotation_{{aligner}}.log"
    params:
        cluster_key = config.get("COMM_DETECT", "leiden"),
        # Dictionary mapping cluster IDs (as strings) to cell type labels
        cluster_mapping = config.get("CLUSTER_MAP", {
            "0": "T Cells",
            "1": "B Cells",
            "2": "Monocytes",
            "3": "NK Cells"
        }),
        # Canonical marker genes to display in the verification dotplot
        canonical_markers = config.get("CANONICAL_MARKERS", [
            "CD3D", "CD3E",  # T Cells
            "CD19", "MS4A1", # B Cells
            "CD14", "LYZ",   # Monocytes
            "NCAM1", "KLRB1" # NK Cells
        ]),
        top_n_markers = config.get("TOP_N_MARKERS", 10)
    conda:
        "../env/scrna_preprocess.yaml"
    threads: 4
    script:
        "scripts/scrna_manual_annotate.py"




rule scrna_auto_annotation:
    """
    Automatically annotate cell types using pre-trained CellTypist models.
    """
    input:
        h5ad = f"{READS_DIR}/h5ad/{{aligner}}/pca_concat_batch.h5ad"
    output:
        h5ad_annotated = f"{READS_DIR}/annotation/h5ad/{{aligner}}_auto_annotated.h5ad",
        umap_plot = f"{READS_DIR}/annotation/plots/{{aligner}}_auto_celltypes_umap.png",
        probability_plot = f"{READS_DIR}/annotation/plots/{{aligner}}_celltypist_conf_umap.png"
    log:
        f"{LOG_DIR}/annotation/auto_annotation_{{aligner}}.log"
    params:
        model_name = config.get("CELLTYPIST_MODEL", "Immune_All_Low.pkl"),
        use_majority_voting = config.get("USE_MAJORITY_VOTING", True),
        cluster_key = config.get("COMM_DETECT", "leiden")
    conda:
        "../env/scrna_preprocess.yaml" # Ensure celltypist is installed in this env
    threads: 4
    script:
        "scripts/scrna_auto_annotate.py"