# 10x-Toolbox
10x Toolbox is a semi-automatic analyzing workflow for single cell sequencing. It was designed to simplify complicated Bioinfo workflows into short commands, while utilizing a project specific config.yaml to record and guarantee repeatability. 

* Uses Starsolo and kb_python as aligner, supports virtually all UMI-based 3' and 5' scRNA-seq protocols.

* Runtime environment: Python, R, shell

* Execution manager: Snakemake, Conda


## Workflow Architecture
10x Toolbox accepts input format: `fastq`, `h5`, `MTX`, `h5ad` (1:1 sample or concatenated) \
And outputs analyzed plots and files including: \
    - PCA, tSNE, UMAP embedded h5ad  \ 
    - tSNE, UMAP.png \
    - Pseudobulk DESeq2 DEG_results.csv \
    - Heatmap, Volcano_plot.pdf \
    - Enrichment result and dotplot 

Hereunder is a simplified flow chart for the workflow's Architecture
<img width="763" height="962" alt="10x Toolbox scRNA flow chart drawio" src="https://github.com/user-attachments/assets/aa2c82d6-d4e7-4276-a6d6-94154aa80b89" />

## Reproducibility
<img width="692" height="285" alt="Untitled Diagram drawio" src="https://github.com/user-attachments/assets/7d03508c-3a8e-40db-b9a2-5395c6336545" />


10x Toolbox applies reproducibility controls by using 3 separate files:
1. `reads/project_name/runinfo.csv`
   - Records experiment Run_number, Sample_name and other metadata.
   - Controls which sample/run should be processed in the workflow.
2. `reads/project_name/config.yaml`
   - Records project specific config variables for each module.
   - Controls results by applying the same analysis variables when repeated.
   - Allows instance switching between projects.
3. `env/module.yaml`
   - Contains env.yaml for each individual module.
   - controls package versions.


## Output Showcase

<img width="5465" height="2399" alt="tsne_celltypes_combined_grid" src="https://github.com/user-attachments/assets/d561da8f-c388-49a9-9a0a-4bf22c2c4edb" />

##### Sample tSNE plots using Parker et al (2019) C.elegans dataset 

<img width="1092" height="906" alt="Screenshot from 2026-08-23 17-37-52" src="https://github.com/user-attachments/assets/c2bc8fb8-e8a7-4e70-8afc-7de9486939e5" />

##### Sample Heapmap using Parker et al (2019) dataset (300min vs 400min after bleach C.elegans embryo cells)




