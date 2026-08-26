import os
import pandas as pd

# =============================================================================
# Setup and Configuration
# =============================================================================
PRJNAME = config.get("PRJNAME") 
if not PRJNAME:
    raise ValueError("ERROR: Specify PRJNAME via --config")

# Load project-specific configuration from the project directory
configfile: f"reads/{PRJNAME}/config.yaml"

# --- Mandatory Configuration Checks ---
if "REF" not in config or not isinstance(config["REF"], dict):
    raise ValueError("\n\n[CONFIGURATION ERROR]\nMissing 'REF' block structure in config.yaml!\n")

mandatory_ref_vars = {
    "SOURCE": "choose from ncbi or ensembl",
    "SPECIES": "Latin binomial name (e.g., saccharomyces_cerevisiae)",
    "ASSEMBLY": "Major assembly build designation (e.g., R64-2-1 or GRCh38)"
}

for var, description in mandatory_ref_vars.items():
    if not config["REF"].get(var):
        raise ValueError(
            f"\n\n[CONFIGURATION ERROR]\n"
            f"Variable 'REF -> {var}' is missing!\n"
            f"Description: {description}\n"
        )

# =============================================================================
# Unified Global Path & Context Variables
# =============================================================================
# REF block variables
ref_cfg  = config["REF"]
source   = ref_cfg["SOURCE"].lower()
species  = ref_cfg["SPECIES"].lower()
assembly = ref_cfg["ASSEMBLY"]
release  = ref_cfg.get("RELEASE", "")
acc = ref_cfg.get("ACC") or "GCF_000848505.1"

# Name variations for structural renaming rules
species_low = species.replace(" ", "_")
species_cap = species_low.capitalize()

# Paths
REFS_DIR = f"refs/{source}/{species_low}/{assembly}/{release}"
READS_DIR = f"reads/{PRJNAME}"
LOG_DIR = f"reads/{PRJNAME}/logs"


# Sample file patterns
cb_umi_suffix = config.get("CB_UMI_SUFFIX", "1")
cdna_suffix = config.get("CDNA_SUFFIX", "2")




# =============================================================================
# Helper Functions
# =============================================================================
def get_refs(wildcards):
    """
    Dynamically routes pipeline requirements to the isolated database path layout.
    """
    if source == "ensembl":
        return {
            "fasta": f"{REFS_DIR}/{species_cap}.{assembly}.{release}.fa",
            "gff": f"{REFS_DIR}/{species_cap}.{assembly}.{release}.gff3",
            "gtf": f"{REFS_DIR}/{species_cap}.{assembly}.{release}.gtf",
            "fai": f"{REFS_DIR}/{species_cap}.{assembly}.{release}.fa.fai"
        }
    elif source == "ncbi":
        return {
            "fasta": f"{REFS_DIR}/{acc}.fa",
            "gff": f"{REFS_DIR}/{acc}.gff3",
            "gtf": f"{REFS_DIR}/{acc}.gtf",
            "fasta": f"{REFS_DIR}/{acc}.fa.fai"
        }
    else:
        raise ValueError(f"Unsupported REF_SOURCE target: {source}")



# =============================================================================
# scRNA-seq Helper Functions
# =============================================================================
def get_scrna_runinfo(wildcards):
    """
    Triggers the appropriate metadata checkpoint and returns the filtered 
    scRNA run accessions based on the configured datasource.
    """
    datasource = config["DATASOURCE"].lower()

    if datasource == "sra":
        # Evaluates SRA metadata checkpoint
        checkpoints.fetch_sra_metadata.get(**wildcards)
        csv_path = f"{READS_DIR}/sra_runinfo.csv"

    elif datasource == "local_fastq":
        # Evaluates local Illumina processing checkpoint
        checkpoints.process_illumina_local.get(**wildcards)
        csv_path = f"{READS_DIR}/local_runinfo.csv"

    else:
        raise ValueError(f"Invalid DATASOURCE '{datasource}' specified in config. Expected 'sra' or 'local_fastq'.")

    # Read dynamically created metadata CSV
    return csv_path

def get_scrna_fastq(wildcards):
    """Generates QC report targets, merge by multiqc and forces fastq trimming."""
    runinfo_df = pd.read_csv(f"{READS_DIR}/runinfo_scrna.csv")
    runs = runinfo_df["Run"].astype(str).str.strip().unique().tolist()
    return [
        ancient(f"{READS_DIR}/fastqs/{r}_{pair}.fastq.gz")
        for r in runs
        for pair in [1, 2]
    ]

def get_scrna_qc(wildcards):
    """Run fastqc and multiqc"""
    runinfo_df = pd.read_csv(f"{READS_DIR}/runinfo_scrna.csv")
    runs = runinfo_df["Run"].astype(str).str.strip().unique().tolist()

    fastqc_r1 = expand("{rddir}/qc/cb_umi/{r}_{sfb}_fastqc.zip",
                    rddir=READS_DIR, r=runs, sfb=cb_umi_suffix)
    fastqc_r2 = expand("{rddir}/qc/cdna/{r}_{sfc}_fastqc.zip",
                    rddir=READS_DIR, r=runs, sfc=cdna_suffix)
    multiqc_r1 = f"{READS_DIR}/qc/multiqc_cb_umi.html"
    multiqc_r2 = f"{READS_DIR}/qc/multiqc_cdna.html"
  
    return fastqc_r1 + fastqc_r2 + [multiqc_r1, multiqc_r2]

def get_scrna_align(wildcards):
    """Align fastq.gz with index to generate bam and summary report."""
    runinfo_df = pd.read_csv(f"{READS_DIR}/runinfo_scrna.csv")
    runs = runinfo_df["Run"].astype(str).str.strip().unique().tolist()

    scrna_aligner = config.get("SCRNA_ALIGNER", "starsolo").lower()
    
    if scrna_aligner in ["starsolo", "star"]:
        align = expand("{rddir}/aligner/scrna/starsolo/{r}/Aligned.sortedByCoord.out.bam",
                       rddir=READS_DIR, r=runs)
        return align

    elif scrna_aligner == "kb_python":
        align = expand("{rddir}/aligner/scrna/kb/{r}/counts_filtered/cells_x_genes.mtx",
                       rddir=READS_DIR, r=runs)
        return align
    else:
        raise ValueError(f"Invalid aligner '{scrna_aligner}' specified in config.")


def get_scrna_clean_convert_group(wildcards):
    """Returns all expected sample-level merged .h5ad targets."""
    scrna_aligner = config.get("SCRNA_ALIGNER", "starsolo").lower()
    
    # Handle STAR alias
    if scrna_aligner == "star":
        scrna_aligner = "starsolo"
        
    if scrna_aligner not in ["starsolo", "kb_python"]:
        raise ValueError(f"Invalid aligner '{scrna_aligner}' specified in config.")
        
    runinfo_df = pd.read_csv(f"{READS_DIR}/runinfo_scrna.csv")
    samples = runinfo_df["SampleName"].astype(str).str.strip().unique().tolist()
    
    # DAG steps merged h5ad > trigger individual h5ad conversion > trigger clean geneid
    return expand(
        "{rddir}/h5ad/{aln}/grouped/{s}_grouped.h5ad",
        rddir=READS_DIR, aln=scrna_aligner, s=samples)


def get_scrna_h5ad_preprocess(wildcards):
    """Run QC, filter, normalize, hvg selection, Concat, PCA, batch correction and embed UMAP and/or tSNE for .h5ad files"""
    scrna_aligner = config.get("SCRNA_ALIGNER", "starsolo").lower()
    embed = f"{READS_DIR}/h5ad/{scrna_aligner}/pca_concat_batch_embed.h5ad"
    return embed

def get_scrna_anno_cal(wildcards):
    scrna_aligner = config.get("SCRNA_ALIGNER", "starsolo").lower()
    return f"{READS_DIR}/annotation/{scrna_aligner}/annotated.h5ad"

def get_scrna_anno_plot(wildcards):
    scrna_aligner = config.get("SCRNA_ALIGNER", "starsolo").lower()
    return f"{READS_DIR}/annotation/{scrna_aligner}/plots"

def get_scrna_pseudobulk(wildcards):
    scrna_aligner = config.get("SCRNA_ALIGNER", "starsolo").lower()
    return f"{READS_DIR}/annotation/{scrna_aligner}/pseudobulk/deseq2_full_results.csv"

def get_scrna_heatmap_volcano(wildcards):
    scrna_aligner = config.get("SCRNA_ALIGNER", "starsolo").lower()
    heatmap = f"{READS_DIR}/annotation/{scrna_aligner}/pseudobulk/plots/top_degs_heatmap.pdf"
    volcano = f"{READS_DIR}/annotation/{scrna_aligner}/pseudobulk/plots/all_celltypes_volcano.pdf"
    return heatmap, volcano

def get_scrna_enrich(wildcards):
    scrna_aligner = config.get("SCRNA_ALIGNER", "starsolo").lower()
    return f"{READS_DIR}/annotation/{scrna_aligner}/pseudobulk/go_enrichment_results.csv"

# =============================================================================
# Include Modular Rule Files
# =============================================================================
include: "toolbox/get_refs.smk"
include: "toolbox/scrna_getdata_sra.smk"
include: "toolbox/scrna_getdata_local.smk"
include: "toolbox/scrna_aligner.smk" 
include: "toolbox/scrna_h5ad_preprocess.smk" 
include: "toolbox/scrna_annotation.smk" 



# =============================================================================
# Terminal Rules
# =============================================================================
rule note:
    input:
    run:
        print("\n" + "="*50)
        print("THIS SNAKEFILE SERVES AS TERMINAL FOR BIOINFORMATIC TOOLBOX")
        print("Please specify a target rule:")
        print("  snakemake download_refs     - Download fa, gff3 and gtf files")
        print("  snakemake scrna_getdata     - Acquire scRNA-Seq dataset")
        print("  snakemake qc                - run QC")
        print("RNAseq Analysis")
        print("  snakemake rna_align         - Run RNA alignment")
        print("  snakemake rna_exp           - Run RNA expression analysis")
        print("  snakemake exp_tools_compare - (optional) Compare expression analysis tools deseq2 vs edger")
        print("  snakemake rna_report        - Evaluate expression data with pca, MA, p-value histagram and volcano plots,")
        print("                                generate csv for list of significantly different expressed genes")
        print("  snakemake rna_enrich        - Generate heatmap and perform functional enrichment study")
        print("  snakemake rna_all           - Gets samples and refs, run QC, rna align, exp, report and enrich")
        print("="*50 + "\n")


rule download_refs:
    """Utility run-target that safely maps upstream conversion assets."""
    input: 
        lambda wildcards: list(get_refs(wildcards).values())



# =============================================================================
# scRNA-Seq Analysis
# =============================================================================
rule scrna_runinfo:
    """Download runinfo.csv or create one from fastq files under reads/project_name/rawdata"""
    input:
        get_scrna_runinfo
    output:
        f"{READS_DIR}/runinfo_scrna.csv"
    shell:
        """
        mv {input} {output}
        """

rule scrna_fastqs:
    """Download run files from SRA (or use local fastq files), run fastqc and multiqc"""
    input: get_scrna_fastq

rule scrna_qc:
    """Run Fastqc and Multiqc on fastq files"""
    input: get_scrna_qc

rule scrna_align:
    """Align RNAseq reads using starsolo or kb_python."""
    input: 
        get_scrna_align,
        get_scrna_clean_convert_group


rule scrna_preprocess:
    """Run through all preprocessing procedures including:
    QC, Filter, Normalize, HVG, PCA, Concat, Batch Correction, Embed."""
    input: get_scrna_h5ad_preprocess

rule scrna_annotation:
    """Run Annotation calculations and plots. """
    input: 
        get_scrna_anno_cal,
        get_scrna_anno_plot

rule scrna_anno_cal:
    """ Run only the calculations for annitation. """
    input: get_scrna_anno_cal

rule scrna_anno_plot:
    """ Run only the plots for annitation. """
    input: get_scrna_anno_plot

rule scrna_pseudobulk:
    """Prepare dataset for a pseudobulk clustering by cell_types, run_id and sample_name / condition,
       Run Deseq2 for DEG analysis"""
    input: get_scrna_pseudobulk

rule scrna_degplots:
    """Generate heatmap and volcano plots"""
    input: get_scrna_heatmap_volcano

rule scrna_enrich:
    """Run Enrichment analysis"""
    input: get_scrna_enrich

rule scrna_all:
    """A shortcut to run a full course of RNAseq with qc, aligner, exp, report and enrich."""
    input: 
        get_scrna_qc,
        get_scrna_align,
        get_scrna_clean_convert_group,
        get_scrna_h5ad_preprocess,
        get_scrna_anno_cal,
        get_scrna_anno_plot,
        get_scrna_pseudobulk,
        get_scrna_heatmap_volcano,
        get_scrna_enrich





