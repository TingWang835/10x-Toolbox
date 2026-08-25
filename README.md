# 10x-Toolbox
10x Toolbox is a semi-automatic analyzing workflow for single cell sequencing. It was designed to simplify complecated Bioinfo workflows into short commands, while utilizing a project specific config.yaml to record and guarantee repeatability. 

* Uses Starsolo and kb_python as aligner, supports virtually all UMI-based 3' and 5' scRNA-seq protocols.

* Runtime environment: Python, R, shell

* Execution manager: Snakemake, Conda

* Accept input format: fastq, h5, MTX, h5ad (1:1 sample or concatenated)

* Main outputs: 
    - h5ad (PCA, tSNE, UMAP embedded); 
    - tSNE, UMAP.png;
    - Pseudobulk DESeq2 DEG_results.csv
    - Heatmap, Volcano_plot.pdf
    - Enrichment result and dotplot




##### Sample tSNE plots using Parker et al (2019) C.elegans dataset 

<img width="1092" height="906" alt="Screenshot from 2026-08-23 17-37-52" src="https://github.com/user-attachments/assets/c2bc8fb8-e8a7-4e70-8afc-7de9486939e5" />

##### Sample Heapmap using Parker et al (2019) dataset (300min vs 400min after bleach C.elegans embryo cells)

## Workflow Architecture
<img width="763" height="962" alt="10x Toolbox scRNA flow chart drawio" src="https://github.com/user-attachments/assets/ed05578d-a5ae-4269-925e-516ee8946868" />
