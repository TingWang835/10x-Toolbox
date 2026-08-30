# 10x-Toolbox
10x Toolbox is a semi-automatic analyzing workflow for single cell sequencing. It is my personal project to simplify complicated Bioinfo workflows into short commands, while utilizing a project specific config.yaml to record and guarantee repeatability. 

* Uses Starsolo and kb_python as aligner, supports virtually all UMI-based 3' and 5' scRNA-seq protocols.

* Runtime environment: Python, R, shell.

* Execution manager: Snakemake, Conda.

* Tested on Linux Mint system.


### Bookmarks 
1. Introduction \
[1.1 Workflow Architecture](#1.1-Workflow-Architecture) \
[1.2 Reproducibility](#1.2-Reproducibility) \
[1.3 Rsult Showcase](#1.3-Result-Showcase) 

2. Setup \
[2.1 Miniconda 3](#2.1-Miniconda-3) \
[2.1 Install Snakemake](#2.1-Install-Snakemake) \
[2.2 Setup Jupyter Kernel](#2.2-Setup-Jupyter-Kernel)

3. 




[2.2 Download Reference](#2.2-Download-Reference) \
[2.3 QC](#2.3-QC) \
DNA and VCF related functions \
[2.4 DNA Aligner](#2.4-DNA-Aligner) \
[2.5 VCF](#2.5-VCF) \
[2.6 DNA rigidity score](#2.6-DNA-rigidity-score) \
RNAseq related functions \
[2.7 RNA Aligner](#2.7-RNA-aligner) \
[2.8 RNA expression and tools compare](#2.8-RNA-expression) \
[2.9 Expression report analysis Heatmap and function enrichment](#2.9-expression-report-analysis-heatmap-and-function-enrichment) \
Chip_Seq related functions \
[2.10 Chip Aligner](#2.10-Chip-aligner) \
[2.11 Chip Peak Caller](#2.11-Chip-Peak-Caller) \
[2.12 Chip Annotation and Enrichment Analysis](#2.12-Chip-Annotation-and-Enrichment-Analysis) \



# Introduction
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

Since the workflow heavy relies on runinfo for sample control, a runinfo.csv will be automatically generated when local files are ingested. For more detail please see chapter [3.3 Getting Started](#3.3-Getting-Started)

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



<a id="1.3-Result-Showcase"></a>
## 1.3 Result Showcase

<img width="2732" height="1199" alt="tsne_celltypes_combined_grid" src="https://github.com/user-attachments/assets/d561da8f-c388-49a9-9a0a-4bf22c2c4edb" />

##### Sample tSNE plots using Packer et al (2019, PMID: 31488706) C.elegans dataset 

<img width="1092" height="906" alt="Screenshot from 2026-08-23 17-37-52" src="https://github.com/user-attachments/assets/c2bc8fb8-e8a7-4e70-8afc-7de9486939e5" />

##### Sample Heatmap using Packer et al (2019, PMID: 31488706) dataset (300min vs 400min after bleach C.elegans embryo cells)

# 2 Setup
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
<a id="2.1-Install-Snakemake"></a>
## 2.1 Install Snakemake

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

<a id="2.2-Setup-Jupyter-Kernel"></a>
## 2.2 Setup Jupyter Kernel
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
   

# 3 10x Toolbox User Manual
<a id="3.1-Commands"></a>
## 3.1 Commands
Here is a list of short commands and their functions. this list can be reviewed any time with terminal command `./ run.sh your_PRJNAME note`. 

#### Universal Functions 
   `snakemake download_refs`      - Download fa, gff3 and gtf files 
  
#### scRNAseq Analysis 
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
10x Toolbox tights both `reads` and `results` to their project, in order for snakemake to find the right path, every project should be allocated in their specific main folder, which shares the same name as config variable "PRJNAME" under`./reads/your_projname`. For more deatils please see chapter [3.3 Getting Started](#3.3-Getting-Started).


## 3.3 Getting Started









