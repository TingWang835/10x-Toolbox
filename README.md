# 10x-Toolbox
10x Toolbox is a semi-automatic analyzing workflow for single cell sequencing. It is my personal project to simplify complicated Bioinfo workflows into short commands, while utilizing a project specific config.yaml to record and guarantee repeatability. 

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
[3.2 Reference & Databases](#3.2-Reference-&-Databases) \
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
2. Modular `smk` files : carries modulated snakemake rules, micromanage input, output, `env` and log path for `shell` or `scripts`. 
3. Python, R, `scripts` : Carries actually codes for a particular function. 
4. `Env.yaml` : manages packages, versions and dependencies via conda.

### Reproducibility Controls
<img width="692" height="285" alt="Untitled Diagram drawio" src="https://github.com/user-attachments/assets/7d03508c-3a8e-40db-b9a2-5395c6336545" />


The workflow guarantee a reproducible outcome via a 3-way control: 
#### 1. `reads/project_name/config.yaml` 
Config.yaml serves as a project specific switch board for all variables used in the workflow. By keeping the variables unchanged, it guarantees the project being analyzed under the same conditions. On the other hand, it also serves as an easy access for changing variables for the analysis process. Instead of going through the code "jungles" to set the right numbers or terms, a simple modification in config.yaml will get the job done. 

#### 2. `reads/project_name/runinfo.csv` 
Runinfo.csv is a common SRA download file recording the run/sample numbers and their metadata. This workflow was designed to rely on the file to control samples for analysis, in simple term, what you see on runinfo is what will get analyzed. By filtering the list of runs and samples, the workflow parse only the wanted one for downstream process. \

Since the workflow heavy relies on runinfo for sample control, a runinfo.csv will be automatically generated when local files are ingested. For more detail please see chapter [3.4 Online or Local Data](#3.4-Online-or-Local-Data)

#### 3. `envs/Env.yaml`
Env.yaml under the envs folder controls information about the packages, version and their dependencies for each individual module, where this info will be passed into Conda for managing the most suitable environment. Under certain circumstance e.g. re-analyzing a older dataset using a particular package version, the version number can be specified to guarantee a repeatable outcome.


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
   4. Run following command in directory terminal (first time only), to permit execution of run.sh and runinfo_cond.sh.
   ```bash
      chmod +x run.sh
      chmod +x runinfo_cond.sh
   ``` 




<a id="2.3-Setup-Jupyter-Kernel"></a>

## 2.3 Setup Jupyter Kernel
ipynb is a useful tool for parsing h5ad files, it needs a kernel containing scanpy. \
To avoid installing scanpy in base kernel, a link between snakemake env and new kernel can be connected. 

   1. After installing snakemake env, activate it by command
   ```bash
   conda activate your_snakemake_env
   ```
   2. Create a Jupyter Kernel connection your current activated env
   ```bash
   python -m ipykernel install --user --name your_env_name --display-name "Jupyter_display_name"
   # your_env_name will be the name in your computer
   # "Jupyter_display_name" will be the name displayed in you ipynb kernel selection
   ````
   3. Find query_h5ad.ipynb, change its running kernel to Jupyter_display_name
   

# 3. scRNA-Seq Analysis
<a id="3.1-Commands"></a>

## 3.1 Commands
Here is a list of short commands and their functions. this list can be reviewed any time with terminal command `./ run.sh your_PRJNAME note`. 

#### Universal Functions 
   `snakemake download_refs`      - Download fa, gff3 and gtf files 
  
#### scRNAseq Analysis 
   `snakemake scrna_runinfo`      - Prepare runinfo.csv, move local ingests to respective paths\
   `snakemake scrna_fastqs`       - Download SRA fastq dataset \
   `snakemake scrna_qc`           - run fastQC and MultiQC\
   `snakemake scrna_align`        - Alignment with Starsolo or kb_python, convert and group individual MTX to h5ad\
   `snakemake scrna_preprocess`   - Preprocess include QC, filter, normalize, HVG, PCA, concat, batch correction, embed\
   `snakemake scran_annotation`   - Characterize clusters based on marker + plot characterized UMAP and/or tSNE\
      `snakemake scrna_anno_cal`  - Subdivision of scrna_annotation, run only characterization\
      `snakemake_scrna_anno_plot` - Subdivision of scrna_annotation, run only plotting\
   `snakemake scrna_pseudobulk`   - Prep pseudobulk, Compute DEG by DESeq2\
   `snakemake scrna_degplot`      - Generate Heatmap and Volcano plots\
   `snakemake scrna_enrich`       - Perform functional enrichment study\
   `snakemake scrna_all`          - Run all comands in consequence

<a id="3.2-References-Databases-Reads-&-Results"></a>

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
10x Toolbox tights both `reads` and `results` to their project, in order for snakemake to find the right path, every project should be allocated in their specific main folder, which shares the same name as config variable "PRJNAME" under`./reads/your_projname`. For more deatils please see chapter [3.4 Online or Local Data](#3.4-Online-or-Local-Data).


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
Data source type is determined by  config variable `DATASOURCE`. \
If `sra` was chosen, `PRJNUMBER` is needed to get runinfo.csv. \
For `sra` and `local_fastq`, the `SUFFIX` variables is needed to correctly identify their contents
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
   After filling in variables,  runn following command In your snakemake env terminal:
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
If DATASOURCE ==`local_x` , allocate local files into `reads/your_PRJNAME/rawdata` before running the command above.
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
   If the dataset are in fastq format (Online or local), fastQC and multiQC can be run to check data quality.

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

At the last step of preprocess, the script will add `condition` column as a metadata in adata.obs matrix. This will be the sample condition (e.g. control_1, treatment_1) printed on future plots and results.

A quick way to add the sample conditions based on sample name:
1. Open runinfo_cond.sh in code code reader.
2. Modify the code at line 15 to create a dictionary for adding the conditions, hereunder is an exmample for Packer et al 2019 dataset.
```bash
BEGIN {
    OFS=","
    # --- DEFINE YOUR CATEGORIES HERE ---
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

Now the runinfo.csv is ready to for preprocess steps, the following config variables can be determined.
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


Finally, run terminal command:
```bash
./run.sh your_PRJNAME scrna_preprocess
```
For storage saving purpose, the workflow will remove the grouped h5ads after preproces was successful. Copy the grouped h5ad files to a separate location if needed, otherwise the finish output can be found at `reads/your_PRJNAME/h5ad/pca_concat_batch_embed.h5ad`.


<a id="3.8-Annotation-&-Enrichment-Analysis"></a>

### 3.8 Annotation & Enrichment Analysis (All Data Source)








