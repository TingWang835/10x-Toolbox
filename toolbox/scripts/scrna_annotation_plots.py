import os
import sys
import matplotlib.pyplot as plt
import gc
import pandas as pd
import scanpy as sc
import seaborn as sns

# Redirect stdout/stderr to log file
sys.stderr = open(snakemake.log[0], "w")
sys.stdout = sys.stderr

# Load annotated AnnData object and calculated differential expression markers
adata = sc.read_h5ad(snakemake.input["h5ad"])
marker_df = pd.read_csv(snakemake.input["markers_csv"])

cluster_key = str(snakemake.params["cluster_key"])
top_n = int(snakemake.params["top_n_markers"])
embed = str(snakemake.params["embed"]).lower()

plots_dir = str(snakemake.output["plots_dir"])
os.makedirs(plots_dir, exist_ok=True)

# =============================================================================
# Collect top N markers from precalculated CSV
# =============================================================================
categories = [
    c for c in adata.obs["cell_type"].cat.categories if c != "Unknown"
]
plot_markers = []

for cat in categories:
    # Filter CSV for group markers sorted by logfoldchanges / pvals
    sub_df = marker_df[marker_df["group"] == cat]
    cat_genes = list(sub_df["names"].head(top_n))
    plot_markers.extend(cat_genes)

unique_plot_markers = [
    g for g in dict.fromkeys(plot_markers) if g in adata.var_names
]

# Map WBGene IDs -> gene_symbol strings
has_symbols = "gene_symbol" in adata.var.columns
if has_symbols:
    var_to_symbol = adata.var["gene_symbol"].to_dict()
    unique_plot_symbols = [
        str(var_to_symbol.get(g, g))
        for g in unique_plot_markers
        if pd.notna(var_to_symbol.get(g, g))
    ]
    unique_plot_symbols = list(dict.fromkeys(unique_plot_symbols))
else:
    unique_plot_symbols = unique_plot_markers

# =============================================================================
# Diagnostic Plot Generation
# =============================================================================

# -----------------------------------------------------------------------------
# A. Standalone UMAP Plot
# -----------------------------------------------------------------------------
umap_tsne_figsz = (snakemake.params["umap_tsne_width"], snakemake.params["umap_tsne_tall"])

# assigning palettes by cluster numbers
n_celltypes = len(adata.obs["cell_type"].cat.categories)
celltype_palette = sns.color_palette("husl", n_colors=n_celltypes).as_hex()

n_leiden = len(adata.obs[cluster_key].cat.categories)
leiden_palette = sns.color_palette("husl", n_colors=n_leiden).as_hex()

if embed in ["umap", "umap&tsne"] and "X_umap" in adata.obsm:
    # 1. Overall Cell Type UMAP
    fig, ax = plt.subplots(figsize=umap_tsne_figsz)
    sc.pl.umap(
        adata,
        color="cell_type",
        title="Annotated Lineages (All Conditions)",
        palette=celltype_palette,
        legend_loc="right margin",
        show=False,
        ax=ax,
    )
    ax.set_aspect("equal", adjustable="box")
    plt.savefig(
        os.path.join(plots_dir, "umap_celltypes_all.png"),
        dpi=300,
        bbox_inches=None,
    )
    plt.close(fig)

    # 2. Overall Leiden Cluster UMAP
    fig, ax = plt.subplots(figsize=umap_tsne_figsz)
    sc.pl.umap(
        adata,
        color=cluster_key,
        title=f"Leiden Clusters ({cluster_key})",
        legend_loc="right margin",
        palette=leiden_palette,
        show=False,
        ax=ax,
    )
    ax.set_aspect("equal", adjustable="box")
    plt.savefig(
        os.path.join(plots_dir, "umap_leiden_all.png"),
        dpi=300,
        bbox_inches=None,
    )
    plt.close(fig)

    # 3. Individual Single-Panel UMAPs per Condition
    if "condition" in adata.obs.columns:
        conditions = adata.obs["condition"].cat.categories
        for cond in conditions:
            sub_adata = adata[adata.obs["condition"] == cond]

            fig, ax = plt.subplots(figsize=umap_tsne_figsz)
            sc.pl.umap(
                sub_adata,
                color="cell_type",
                title=f"Annotated Lineages ({cond})",
                palette=celltype_palette,
                legend_loc="right margin",
                show=False,
                ax=ax,
            )
            ax.set_aspect("equal", adjustable="box")
            out_path = os.path.join(plots_dir, f"umap_celltypes_cond_{cond}.png")
            fig.savefig(out_path, bbox_inches=None, dpi=300)
        
        # Free memory allocations explicitly
    plt.close(fig)
    del sub_adata
    gc.collect()

# -----------------------------------------------------------------------------
# B. Standalone t-SNE Plot
# -----------------------------------------------------------------------------
if embed in ["tsne", "umap&tsne"] and "X_tsne" in adata.obsm:
    # 1. Overall Cell Type t-SNE
    fig, ax = plt.subplots(figsize=umap_tsne_figsz)
    sc.pl.tsne(
        adata,
        color="cell_type",
        title="Annotated Lineages (All Conditions)",
        palette=celltype_palette,
        legend_loc="right margin",
        show=False,
        ax=ax,
    )
    ax.set_aspect("equal", adjustable="box")
    plt.savefig(
        os.path.join(plots_dir, "tsne_celltypes_all.png"),
        dpi=300,
        bbox_inches=None, 
    )
    plt.close(fig)

    # 2. Overall Leiden Cluster t-SNE
    fig, ax = plt.subplots(figsize=umap_tsne_figsz)
    sc.pl.tsne(
        adata,
        color=cluster_key,
        title=f"Leiden Clusters ({cluster_key})",
        palette=leiden_palette,
        legend_loc="right margin",
        show=False,
        ax=ax,
    )
    ax.set_aspect("equal", adjustable="box")
    plt.savefig(
        os.path.join(plots_dir, "tsne_leiden_all.png"),
        dpi=300,
        bbox_inches=None,
    )
    plt.close(fig)

    # 3. Individual Single-Panel t-SNEs per Condition
    if "condition" in adata.obs.columns:
        conditions = adata.obs["condition"].cat.categories
        for cond in conditions:
            sub_adata = adata[adata.obs["condition"] == cond]

            fig, ax = plt.subplots(figsize=umap_tsne_figsz)
            sc.pl.tsne(
                sub_adata,
                color="cell_type",
                title=f"Annotated Lineages ({cond})",
                palette=celltype_palette,
                legend_loc="right margin",
                show=False,
                ax=ax,
            )
            ax.set_aspect("equal", adjustable="box")
            out_path = os.path.join(plots_dir, f"tsne_celltypes_cond_{cond}.png")
            fig.savefig(out_path, bbox_inches=None, dpi=300)
        
        # Free memory allocations explicitly
    plt.close(fig)
    del sub_adata
    gc.collect()

# -----------------------------------------------------------------------------
# C. Ordered / Diagonalized Marker Dot Plot
# -----------------------------------------------------------------------------
if unique_plot_symbols:
    print(f"Generating ordered dotplot for {len(unique_plot_symbols)} marker symbols...")
    
    if has_symbols:
        valid_symbols = [
            s for s in unique_plot_symbols if s in adata.var["gene_symbol"].values
        ]
    else:
        valid_symbols = [s for s in unique_plot_symbols if s in adata.var_names]

    if valid_symbols:
        dp = sc.pl.dotplot(
            adata,
            var_names=valid_symbols,
            groupby="cell_type",
            gene_symbols="gene_symbol" if has_symbols else None,
            use_raw=False,
            standard_scale="var",
            dendrogram=True,
            figsize=(12, 7),
            show=False,
        )

        dotplot_path = os.path.join(plots_dir, "markers_dotplot.png")
        if hasattr(dp, "savefig"):
            dp.savefig(dotplot_path, bbox_inches="tight", dpi=300)
        else:
            plt.savefig(dotplot_path, bbox_inches="tight", dpi=300)
        plt.close()

print("Plot generation complete!")