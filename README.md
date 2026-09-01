# 10x-Toolbox
10x Toolbox is a semi-automatic analyzing workflow for single cell sequencing. It is my personal project to simplify complicated Bioinformatics workflows into short commands, while learning in practice of the repective analytic skills. 

* Uses Starsolo and kb_python as aligner, supports virtually all UMI-based 3' and 5' scRNA-seq protocols.

* Runtime environment: Python, R, shell.

* Execution manager: Snakemake, Conda.

* Tested on Linux Mint system.


### Bookmarks 
1. Introduction \
[1.1 Workflow Architecture](#1.1-Workflow-Architecture) \
[1.2 Rsult Showcase](#1.2-Result-Showcase) 

2. Setup \
[2.1 Miniconda 3](#2.1-Miniconda-3) \
[2.2 Install Snakemake](#2.2-Install-Snakemake) \
[2.3 Setup Jupyter Kernel](#2.3-Setup-Jupyter-Kernel)

3. scRNA-Seq Analysis \
[3.1 Commands](#3.1-Commands) \
[3.2 Reference & Databases](#3.2-References-&-Databases) \
[3.3 Download Reference](#3.3-Download-Reference) \
[3.4 Online or Local Data](#3.4-Online-or-Local-Data)\
[3.5 QC](#3.5-QC) \
[3.6 Alignment](#3.6-Alignment) \
[3.7 Preprocess](#3.7-Preprocess) \
[3.8 Annotation & Enrichment Analysis](#3.8-Annotation-&-Enrichment-Analysis) 






# 1. Introduction
<a id="1.1-Workflow-Architecture"></a>

## 1.1 Workflow Architecture

### Core Structures
Core structure of the 10x Toolbox workflow was constructed by 4 major file types: 
1. `Snakefile` : Terminal control, assigns functions for each short command, coordinates downstream `smks` and `scripts`. 
2. Modular `smk` files : carries modulated snakemake rules, micromanage input, output, `env` and log path for `shells` code blocks or Python and R `scripts`. 
3. Python, R, `scripts` : Carries actually codes for a particular function. 
4. `Env.yaml` : manages packages, versions and dependencies via conda.

### Reproducibility Controls
<img width="692" height="285" alt="Untitled Diagram drawio" src="https://github.com/user-attachments/assets/7d03508c-3a8e-40db-b9a2-5395c6336545" />


The workflow guarantee a reproducible outcome via a 3-way control: 
#### 1. `reads/project_name/config.yaml` 
Config.yaml serves as a project specific switch board for all variables used in the workflow. By keeping the variables unchanged, it guarantees the project being analyzed under the same conditions. Moreover, it also serves as an easy access for changing variables for the analysis process. Instead of going through the code "jungles" to set the right numbers or terms, a simple modification in config.yaml will get the job done. 

#### 2. `reads/project_name/runinfo.csv` 
runinfo.csv is a common SRA download file recording the run/sample numbers and their metadata. This workflow was designed to rely on the file to control samples for analysis, in simple term, what you see on runinfo is what will get analyzed. By filtering the list of runs and samples, the workflow parse only the wanted one for downstream process. 

Since the workflow heavy relies on runinfo for sample control, a runinfo.csv will be automatically generated when local files are ingested. For more detail please see chapter [3.4 Online or Local Data](#3.4-Online-or-Local-Data)

#### 3. `envs/Env.yaml`
Env.yaml under the envs folder controls information about the packages, version and their dependencies for each individual module, this info will be passed into Conda for managing suitable environment. Under certain circumstance e.g. re-analyzing a older dataset using a particular package version, the version number can be specified to guarantee a repeatable outcome.


### Inputs & outputs
10x Toolbox accepts input format: `fastq`, `h5`, `MTX`, `h5ad` (1:1 sample or concatenated) 

And outputs analyzed plots and files including: \
    - PCA, tSNE, UMAP embedded h5ad \
    - tSNE, UMAP.png \
    - Pseudobulk DESeq2 DEG_results.csv \
    - Heatmap, Volcano_plot.pdf \
    - Enrichment result and dotplot 

Simplified flow chart for the workflow's Architecture
<img width="763" height="962" alt="10x Toolbox scRNA flow chart drawio" src="https://github.com/user-attachments/assets/aa2c82d6-d4e7-4276-a6d6-94154aa80b89" />



<a id="1.2-Result-Showcase"></a>

## 1.2 Result Showcase

<img width="2732" height="1199" alt="tsne_celltypes_combined_grid" src="https://github.com/user-attachments/assets/d561da8f-c388-49a9-9a0a-4bf22c2c4edb" />

##### Sample tSNE plots using Packer et al (2019, PMID: 31488706) C.elegans dataset 

<img width="1092" height="906" alt="Screenshot from 2026-08-23 17-37-52" src="https://github.com/user-attachments/assets/c2bc8fb8-e8a7-4e70-8afc-7de9486939e5" />

##### Sample Heatmap using Packer et al (2019, PMID: 31488706) dataset (300min vs 400min after bleach C.elegans embryo cells)



# 2. Setup
<a id="2.1-Miniconda-3"></a>

## 2.1 Miniconda 3
   1. Download latest miniconda 3
   ```bash
      wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
   ```
   2. Install 
   ```bash
      bash Miniconda3-latest-Linux-x86_64.sh
   ```
   3. Reload shell
   ```bash
      source ~/.bashrc
   ```



<a id="2.2-Install-Snakemake"></a>

## 2.2 Install Snakemake

   1. Download and unpack github package under your working directory.
   2. Make sure snakemake_install.yml is in the working dir.
   3. CD terminal to working dir and run:
   ```bash
      conda env create -f snakemake_install.yml -n new_env_name
   ```
   4. Run following command in directory terminal (first time only), to permit execution of run.sh and runinfo_cond.sh on your linux system.
   ```bash
      chmod +x run.sh
      chmod +x runinfo_cond.sh
   ``` 
   run.sh is a little bash shortcut to run terminal command:
   ```bash 
   # it cuts down boiler plate heavy code:
   snakemake --config your_PRJNAME --use-conda --cores 12 --printshellcmds download_refs
   
   # into:
   ./run.sh your_PRJNAME download_refs
   ```
   run.sh also register `$EXTRA_FLAGS` so other snakemake flags can be appended.\
   e.g. if need to refresh datetime of current command's input files 
   ```bash
   ./run.sh your_PRJNAME download_refs --touch
   ```
   runinfo_cond.sh is another bash tool to automatically add condition based on SampleName in runinfo.csv, for more details please refer to [What is runinfo_cond.sh?](#What-is-runinfo_cond.sh?)



<a id="2.3-Setup-Jupyter-Kernel"></a>

## 2.3 Setup Jupyter Kernel
ipynb is a useful tool in parsing h5ad files, it needs a kernel containing scanpy. \
To avoid installing scanpy in base kernel, a link between snakemake env and new kernel can be connected. 

   1. After installing snakemake env, activate it by command
   ```bash
   conda activate your_snakemake_env
   ```
   2. Create a Jupyter Kernel connection your current activated env
   ```bash
   python -m ipykernel install --user --name your_env_name --display-name "Jupyter_display_name"
   # your_env_name will be the name stored in your computer
   # "Jupyter_display_name" will be the name displayed in you ipynb kernel selection
   ````
   3. Find query_h5ad.ipynb in directory, change its running kernel to Jupyter_display_name
   

# 3. scRNA-Seq Analysis
<a id="3.1-Commands"></a>

## 3.1 Commands
Here is a list of short commands and their functions. this list can be reviewed any time with terminal command `./ run.sh your_PRJNAME note`. 

#### Universal Functions 
   `download_refs`      - Download fa, gff3 and gtf files 
  
#### scRNAseq Analysis 
   `scrna_runinfo`      - Prepare runinfo.csv, move local ingests to respective paths\
   `scrna_fastqs`       - Download SRA fastq dataset \
   `scrna_qc`           - run fastQC and MultiQC\
   `scrna_align`        - Alignment with Starsolo or kb_python, convert and group individual MTX to h5ad\
   `scrna_preprocess`   - Preprocess include QC, filter, normalize, HVG, PCA, concat, batch correction, embed\
   `scran_annotation`   - Characterize clusters based on marker + plot characterized UMAP and/or tSNE\
      `scrna_anno_cal`  - Subdivision of scrna_annotation, run only characterization\
      `scrna_anno_plot` - Subdivision of scrna_annotation, run only plotting\
   `scrna_pseudobulk`   - Prep pseudobulk, Compute DEG by DESeq2\
   `scrna_degplot`      - Generate Heatmap and Volcano plots\
   `scrna_enrich`       - Perform functional enrichment study\
   `scrna_all`          - Run all comands in consequence

<a id="3.2-References-&-Databases"></a>

## 3.2 References & Databases
This Workflow stores `Reference` and `Database` in well sorted files separated from `Reads` and `Results`, making them sharable bewteen different projects.

1. `References` 
* NCBI: `fasta`, `gff3` download from NCBI website based on Ref variables in config file.  Chromosome names translated to USCS standard, `gtf` file generated form `gff3`.Stored under `./refs/ncbi/species/assembly/release/acc`.

* Ensembl: `fasta`, `gff3` download from Ensmebl/Wormbase website based on Ref variables in config file.  Chromosome names translated to USCS standard, `gtf` file generated form `gff3`. Stored under `./refs/ensembl/species/assembly/release`.


2. `Databases`\
scRNAseq project stores 2 types of files in ./databases folder
* 10x whitelist: barcode wihtelist.txt files for cell ranger and starsolo aligners. file name can be set in config file.

* annotation marker list: this is a csv list recording gene marker of each cell type. used in annotation and cluster identification process. Allocate the file under `./databases/annotation_marker`

3. `Reads` & `Results`\
10x Toolbox tights both `reads` and `results` to their project, in order for snakemake to find the right path, every project should be allocated in their specific main folder, which shares the same name as config variable "PRJNAME" under`./reads/your_PRJNAME`. For more details please see chapter [3.4 Online or Local Data](#3.4-Online-or-Local-Data).


<a id="3.3-Download-Reference"></a>

### 3.3 Download Reference
  Enter the following config variables:
   ```yaml
   REF: 
      SOURCE: "ensembl"
      SPECIES: "caenorhabditis_elegans"
      ASSEMBLY: "WBcel235" # e.g. ncbi: GRCh38, ensembl: R64-2-1
      RELEASE: "38" # ncbi: p14 (patch #), ensemble: 83
      ACC: "" # leave empty for ensembl
   ```
   
   In your snakemake env terminal enter:
   ```bash
   ./run.sh your_PRJNAME download_ref
  ```
  
   fasta, gff3 files will be downloaded, chromosome name translated to USCS standard, gtf file will be generated from gff3. 

<a id="3.4-Online-or-Local-Data"></a>

### 3.4 Online or Local Data
Because the workflow accepts various local or online data types, it is necessary to determine data source type in config variable `DATASOURCE`. 
* If `sra` was chosen, `PRJNUMBER` is needed to get runinfo.csv. 
* For `sra` and `local_fastq`, the `SUFFIX` variables are needed to correctly identify their contents
```yaml
# choose Data source 
# sra: download runinfo.csv based on PRJNUMBER from SRA website, use it as guide;
# local_fastq: running fastq files from reads/{PRJNAME}/rawdata, generate runinfo.csv using file name (illumina naming standard);
# local_concat_h5ad: for local concat and batch corrected h5ad, 
# local_indiv_h5ad: for local individual h5ad
# local_raw_h5_mtx: for local MTX or h5 file ingestion.
DATASOURCE: "sra"

#PRJN number (optional) only used when DATASOURCE = sra
PRJNUMBER: "PRJNA523834"

# Assign Fastq suffix for correctly identificating the content (e.g. _1.fastq.gz for CB + UMI _2.fastq.gz for cDNA)
CB_UMI_SUFFIX: "1"  # Illumina default 1
CDNA_SUFFIX: "2"    # Illumina default 2
```
   After filling in variables,  run following command In your snakemake env terminal:
   ```bash
   ./run.sh your_PRJNAME scrna_runinfo
  ```

#### Online dataset
If DATASOURCE == `sra`, a runinfo.csv will be downloaded from SRA website. columns will be simplified to include `Run, SampleName, Experiment, SampleName, LibraryLayout, Platform`. Filter the list for target run/samples for download.

After runinfo.csv is ready, use the following terminal command to start download via prefetch and fasterq-dump:
   ```bash
   ./run.sh your_PRJNAME scrna_fastqs
  ```

#### Local dataset
If DATASOURCE ==`local_x` , allocate local files into `reads/your_PRJNAME/rawdata` before running `scrna_runinfo` command.
The workflow will copy the allocated files and automatically direct them into the respective path :
1. `local_fastq`: copy files end with *fastq.gz to `reads/your_PRJNAME/fastqs`, build runinfo.csv by parsing file names (suitable for Illumina naming format).

2. `local_concat_h5ad`: copy file ends with *.h5ad to `reads/your_PRJNAME/h5ad`, change name to pca_concat_batch_embed.h5ad, parse metadata in obs and var matrix to generate runinfo.csv.

3. `local_indiv_h5ad`: copy files end with *.h5ad to `reads/your_PRJNAME/h5ad/grouped`, parse metadata in obs and var matrix to generate runinfo.csv.

4. `local_raw_h5_mtx`: 
   * For raw MTX files, the workflow accepts both Cell Ranger and GEO download arrangement as shown below. the script will parse sample name from the mother folder (Cell Ranger) or from the filename (GEO) automatically. MTX files will be converted to h5ad and output to `reads/your_PRJNAME/h5ad/grouped`.
   ```
   Cell Ranger style:
   reads/your_PRJNAME/rawdata/
   ├── sample_1/
   │   ├── matrix.mtx.gz
   │   ├── barcodes.tsv.gz
   │   └── features.tsv.gz
   └── sample_2/
       ├── matrix.mtx.gz
       ├── barcodes.tsv.gz
       └── features.tsv.gz
    ```

   ```
   GEO standard:
   reads/your_PRJNAME/rawdata/
   ├── sample_1_matrix.mtx.gz
   ├── sample_1_barcodes.tsv.gz
   ├── sample_1_genes.tsv.gz
   ├── sample_2_matrix.mtx.gz
   ├── sample_2_barcodes.tsv.gz
   └── sample_2_genes.tsv.gz
   ```
   * Raw h5 files will be converted to h5ad and output to `reads/your_PRJNAME/h5ad/grouped`, a runinfo.csv will be generated using filenames.


<a id="3.5-QC"></a>

### 3.5 QC (fastqs exclusive)
   If the dataset is in fastq format (Online or local), fastQC and multiQC can be run to check data quality.

   No specific config variable needed.

   In your snakemake env terminal enter:
   ```bash
   ./run.sh your_PRJNAME scrna_qc
  ```

  Outputs Location:  `reads/your_PRJNAME/qc/` 


<a id="3.6-Alignment"></a>

### 3.6 Alignment (fastqs exclusive)
   After qualified by QC, the fastqs are ready for alignment.

   The following config variables can be set for either kb_python or Starsolo aligner:
   ```yaml
   # choose rna alinger f(starsolo, kb_python )
   SCRNA_ALIGNER: "kb_python"

   # ------------------------------------------------
   #  (optional) Star aligner specific variable
   # ------------------------------------------------
   10X_VERSION: "v2"
   V2_WHITELIST_NAME: "737K-august-2016.txt"
   V3_WHITELIST_NAME: "3M-february-2018.txt"

   # suffix array index bases
   # Organism Class           |   Genome Size   |   Recommended sa Value
   # Small Eukaryotes / Yeast |     ~12 Mb      |      10 (Default choice)
   # Insects / Worms          | 100 Mb – 300 Mb |      12
   # Large Mammals / Plants   |     ~3 Gb+      |      14 (STAR maximum limit)
   # Bacteria                 |  ~4 Mb – 5 Mb   |       8
   STAR_SA: "12"

   # Cell Barcode and UMI decoder type 
   SOLO_TYPE: "CB_UMI_Simple" #choose from "CB_UMI_Simple", "CB_UMI_Complex", "SmartSeq" (case sensitive)

   # Starsolo cell filtering variables
   SOLO_CELL_FILTER: "CellRanger2.2"
   SOLO_EXPECTED_CELLS: 3000  # typical range from 500~10000
   SOLO_MAX_PERCENTILE: 0.99
   SOLO_MAX_MIN_RATIO: 10

   # ------------------------------------------------
   #(optional) Kallisto variables
   # ------------------------------------------------
   # determine CB and UMI length, decide CD + UMI correction
   KB_TECHNOLOGY: "10XV2" # e.g., 10XV1, 10XV2, 10XV3
   ```

   Run the following terminal command:
   ```bash
   ./run.sh your_PRJNAME scrna_align
  ```
   
   `scrna_align` is a multitasking command, its function includes:
   * Generating index files for selected aligner (stored respective `refs`folder).
   * (Starsolo only) Download 10x whitelists (stored in `databases/10x_whitelists`).
   * Run alignment.
   * Parse the aligned mtx file and convert to h5ads.
   * Grouped and merge individual run/lane h5ads by SampleName column in runinfo.csv, outputs to: `reads/your_PRJNAME/h5ad/grouped`.



<a id="3.7-Preprocess"></a>

### 3.7 Preprocess (fastqs, local_indiv_h5ad, local_raw_h5_mtx)
Preprocess step run further QC, HVG and normalized on the grouped files before concatenated them into one unified h5ad file, followed by PCA, batch correction and UMAP and/or tSNE embedding.


Fill the following config variables to adjust preprocess criterias.
```yaml
# QC filter variables
QC_MIN_GENES: 200       # cells unique gene expression < 200
QC_MAX_GENES: 6000      # cells unique gene expression > 6000
QC_MIN_CELLS: 3         # cell number < 3
QC_MAX_PCT_MITO: 15.0   # cells with unique Mitochrondria gene >15.0% (potentially damanaged)

# normalize and HVG selection variables
NORM_TARGET_SUM: 10000      # normalized to counts per: 10000 (default) or None (Scanpy take medium across samples)
HVG_N_TOP_GENES: 2000        # select the top 2000 highly variable genes for noise reduction
HVG_FLAVOR: "seurat"         # option: "seurat", "cell_ranger" | "seurat_v3" and "pearson_residuals" needs un-logged input

# Concat and PCA variables
Z_CAP: 10             # Z-Score cap to prevent extreme expression bias
N_PC: 50              # Compute the top N number of Priciple Component
SVD_SOLVER: "arpack"  # Choose from "arpack", "random" (for large dataset), "auto"

# Batch correction variables
BATCH_METHOD: "harmony"  # Options: "harmony", "bbknn", "scvi"

# Cluster plots variables 
EMBED: "umap&tsne"  # Option "tsne", "umap", "umap&tsne" choose which 2D projection data to enbed in h5ad
N_NEIGHBOR: 15
```


Run terminal command:
```bash
./run.sh your_PRJNAME scrna_preprocess
```
For storage saving purpose, the workflow will remove the grouped h5ads after preproces was successful. Backup the grouped h5ad files if needed, otherwise the finish output can be found at `reads/your_PRJNAME/h5ad/pca_concat_batch_embed.h5ad`.


<a id="3.8-Annotation-&-Enrichment-Analysis"></a>

### 3.8 Annotation & Enrichment Analysis (All Data Source)

In some dataset, sample names were occupied by serial number and therefore does not carry readable information about the treatment conditions. For more reader friendly results and plots without disturbing the original metadata, this workflow adds a `condition` column in runinfo.csv and match the cells in preprocessed h5ad to create `adata.obs["condition"]` at the very beginning of annotation. 

<a id="What-is-runinfo_cond.sh?"></a>

#### What is runinfo_cond.sh?
To aid mass addition of `condition` based on `SampleName` when dealing with enormous dataset, I have coded a little bash tool runinfo_cond.sh:
1. Open runinfo_cond.sh in a code reader.
2. Modify the code starting at line 15 to create a dictionary for conditions and sample names, hereunder is an exmample for Packer et al 2019 dataset.
```bash
BEGIN {
    OFS=","
    # --- DEFINE YOUR CATEGORIES HERE ---
    # cond[sample_name]="treatment_condition"
    cond["GSM3618670"]="300min"
    cond["GSM3618671"]="400min"
    cond["GSM3618672"]="500min_1"
    cond["GSM3618673"]="500min_2"
    cond["GSM3618674"]="mix1"
    cond["GSM3618675"]="mix2"
    cond["GSM3618676"]="mix3"
    # -------------------------------------
}  
```
3. Run terminal command, and the conditions will be added to your runinfo.csv:
```bash
./runinfo_cond.sh your_PRJNAME scrna
```


### Annotation and Enrichment Analysis
The following config variables are shared in a number for snakemake short commands due to overlapping requirements between package configurations. These short commands includes:
   * `sncrna_annotation`
   * Sub-commands `scrna_anno_cal` and `scrna_anno_plot` 
   * `scrna_pseudobulk` 
   * `scrna_degplot`
   * `scrna_enrich`
```yaml
ANNOTATION_MARKER: databases/annotation_marker/celegans_embryo_markers.csv  # locate marker list file under "databases/annotation_marker/your_marker_name.csv"

CLUSTER_KEY: "leiden"      # "leiden" (default) better at producing well connected clusters, using "louvain" only in reproducing historic data analysis.
RESOLUTION: 0.2            # low res: 0.1-0.3 , high res:1.0-2.0+ determine cluster separatation
TOP_N_MARKERS: 3
SCORE_THRESHOLD: 0.05
UMAP_MIN_DIST: 0.5

# UMAP tSNE plot size
UMAP_TSNE_WIDTH: 4
UMAP_TSNE_TALL: 4

# pseudobulk preperation variable 
CELL_TYPE_KEY: "cell_type" # for finding the correct cell type metadata in e.g. adata.obs["cell_type"]
SAMPLE_KEY: "sample_name"  # same as above 
RUN_KEY: "run_id"          # same as above 
CONDITION_KEY: "condition" # same as above 
MIN_CELLS_PER_GROUP: 10    # filter subgroups with < X number of cells


# DESeq2 variables
DESEQ_ANALYSIS: "wald" # for pair-wise: "wald", for ANOVA: "lrt"
COMPARE_CONDITIONS: ["300min","400min"] # option1: for wald analysis, list EXACTLY 2 conditions
                                                         # option2: for lrt analysis list all selected conditions [control, treatment1, treatment2]
PADJ_CUTOFF: 0.05
LOG2FC_CUTOFF: 0.25

# Heatmap variables
TOP_N_DEGS: 10 # show top X number of DEGs from each cell type in Heatmap.(ranked by P-adj)

# Enrichment variables
GO_PADJ_CUTOFF: 0.05
GENEID_TYPE: "WORMBASE"
# choose from:
  # ENSEMBL
  # WORMBASE 
  # ENTREZID (ncbi gene_ids)
  # SYMBOL (standard gene short names)
```


### 3.8.1 Annotation
The workflow only support sctype to identify cell types from gene marker list at the moment (more features will be added later). A gene marker list in csv format can be created or downloaded and stored under the path determined by config variable `ANNOTATION_MARKER`. Hereunder is an example I have used for embronic C.elegans cells:

```csv
cell_type,markers
Early Progenitors (Cleavage),pie-1;skn-1;pal-1;mex-3
Pharynx Muscle,myo-2;pha-4;ceh-22
Pharynx Epithelium / Marginal,pha-1;ceh-34
Body Wall Muscle,hlh-1;myo-3;unc-120;pat-3
Differentiating Neurons (Pan),rab-3;unc-104;snt-1
Cholinergic Neurons,unc-17;cho-1
GABAergic Neurons,unc-25;unc-47
Glial Support Cells,vap-1;itx-1
Hypodermis / Cuticle Synthesis,dpy-7;col-12;lin-26;noah-1
Seam Cells,scm-1;nhr-25
Intestine / Endoderm,elt-2;ges-1;ifb-2
Excretory Canal / Cell,sulp-4;lim-4
Coelomocytes,cup-4;unc-122
Germline (P4/Z2/Z3),nos-1;nos-2;pgl-1;fbf-1
```

Run terminal command
```bash
./run.sh your_PRJNAME scrna_annotation
```
Because the annotation calculation is a CPU heavy process, it is not ideal to re-run if only for adjusting plot color or size. Therefore separate sub commands can be used:

```bash
./run.sh your_PRJNAME scrna_anno_cal # re-run sctype calculation
./run.sh your_PRJNAME scrna_anno_plot # re-run UMAP and/or tSNE plots
```

Outputs of annotation, DEG, heatmap volcano plots and enrichment analysis is summarized in chatper [3.8.final Output Summary](#3.8.final-Output-Summary).


### 3.8.2 Pseudobulk
For reducing noise, the workflow pseudobulks cells in the same sctype clusters to amplify the differently expressed genes (DEG), and use DESeq2 to quantify the differences.

Because various dataset may store cell metadata under differently name columns, the config provide flexible variables to help parsing the right cell type, sample name and run id. After setting all the right variables, run termian command:

```bash
./run.sh your_PRJNAME scrna_pseudobulk
```

Outputs of annotation, DEG, heatmap volcano plots and enrichment analysis is summarized in chatper [3.8.final Output Summary](#3.8.final-Output-Summary).

### 3.8.3 Heatmap, Volcano plot, Enrichment
Heatmap, Volcano plot, Enrichment result and dotplot can be generated from DEG results from DESeq2. Using the following 2 terminal commands to produce or re-run the step after changes:
```bash
./run.sh your_PRJNAME scrna_degplots
./run.sh your_PRJNAME scrna_enrichment
```

<a id="3.8.final-Output-Summary"></a>

### 3.8.final Output Summary
Here is the summary of all the results generate in chapter 3.8. For convenience, the outputs aree all located under `reads/your_PRJNAME/annotation/aligner` (separated by aligner sub-folders for cross comparison).

These outputs are separated into 3 major folders:


1. `Reports`: contains annotation results of statistic between leiden clusters vs cell type clusters and individual genes vs cell type clusters.
2. `umap_tsne`: contains tSNE, UMAP and marker gene dotplot.png.
3. `pseudobulk`: contains counts, metadatas, DESeq2 results and enrichment results.csv for pseudobulk method.
   * sub-folder `plots`: heatmap, volcano plot and enrichment dotplot is located here.


That is all the function included in 10x Toolbox, Thank you for your patient and happy Bioinformating!

