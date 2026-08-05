localrules: 
# =============================================================================
# Check Points
# =============================================================================





# =============================================================================
# Download Sample data
# =============================================================================
rule parse_scrna_metadata:
    """
    Parses global runinfo.csv to filter single-cell RNA-seq runs,
    maps GSM accessions to SRR runs, and allows sample filtering.
    """
    input:
        csv = f"{READS_DIR}/runinfo.csv"
    output:
        csv = f"{READS_DIR}/runinfo_scrna.csv"
    log:
        f"{LOG_DIR}/getdata/parse_scrna_metadata.log"
    run:
        import pandas as pd

        df = pd.read_csv(input.csv)

        if "LibraryStrategy" in df.columns:
            df = df[df["LibraryStrategy"].str.upper() == "RNA-SEQ"]

        cols_to_keep = [c for c in ["Run", "SampleName", "Experiment", "LibraryLayout", "Platform"] if c in df.columns]
        df_filtered = df[cols_to_keep]

        df_filtered.to_csv(output.csv, index=False)

rule prefetch_sra:
    """
    Prefetches raw .sra container from NCBI into the sra directory.
    """
    output:
        sra_dir = temp(directory(f"{READS_DIR}/sra/{{run}}"))
    log:
        f"{LOG_DIR}/prefetch/{{run}}.log"
    conda:
        "../env/scrna_getdata.yaml"
    threads: 2
    resources:
        mem_mb = 4000
    shell:
        """
        exec 2> {log}
        mkdir -p {READS_DIR}/sra
        
        # Point output directory to the parent 'sra' folder. 
        # prefetch automatically creates the '{wildcards.run}' subfolder inside it.
        prefetch {wildcards.run} \
            --max-size u \
            -O {READS_DIR}/sra
        """

rule parse_scrna_fastqs:
    """
    Parse .sra using fasterq-dump (multi-threaded),
    uses a dedicated temp directory for fast disk I/O, and compresses outputs via pigz.
    """
    input:
        sra_dir = f"{READS_DIR}/sra/{{run}}"
    output:
        r1 = f"{READS_DIR}/fastqs/{{run}}_1.fastq.gz",
        r2 = f"{READS_DIR}/fastqs/{{run}}_2.fastq.gz"
    log:
        f"{LOG_DIR}/download_fastq/{{run}}.log"
    conda:
        "../env/scrna_getdata.yaml"
    threads: 4
    resources:
        mem_mb = 8000
    shell:
        """
        exec 2> {log}
        mkdir -p {READS_DIR}/fastqs
        
        # Create a temporary working directory for fasterq-dump cache files
        TMP_DIR=$(mktemp -d -p {READS_DIR}/fastqs tmp_{wildcards.run}_XXXXXX)

        # Pass the downloaded local SRA directory directly
        fasterq-dump {input.sra_dir} \
            --split-files \
            --include-technical \
            --threads {threads} \
            --outdir {READS_DIR}/fastqs \
            --temp $TMP_DIR

        # Remove temp scratch dir
        rm -rf $TMP_DIR
        
        # Compress with pigz
        pigz -f -p {threads} {READS_DIR}/fastqs/{wildcards.run}_*.fastq
        """


# =============================================================================
# Quality Control
# =============================================================================
# -------------------------------------------------------------------
# FastQC Rules with Config-Driven Read Suffixes
# -------------------------------------------------------------------
rule fastqc_cb_umi:
    input:
        r1 = f"{READS_DIR}/fastqs/{{run}}_{cb_umi_suffix}.fastq.gz"
    output:
        html = f"{READS_DIR}/qc/cb_umi/{{run}}_{cb_umi_suffix}_fastqc.html",
        zip  = f"{READS_DIR}/qc/cb_umi/{{run}}_{cb_umi_suffix}_fastqc.zip"
    params:
        outdir = f"{READS_DIR}/qc/cb_umi"
    conda: "../env/qc.yaml"
    threads: 4
    shell:
        """
        fastqc --quiet --threads {threads} --outdir {params.outdir} {input.r1}
        """


rule fastqc_cdna:
    input:
        r2 = f"{READS_DIR}/fastqs/{{run}}_{cdna_suffix}.fastq.gz"
    output:
        html = f"{READS_DIR}/qc/cdna/{{run}}_{cdna_suffix}_fastqc.html",
        zip  = f"{READS_DIR}/qc/cdna/{{run}}_{cdna_suffix}_fastqc.zip"
    params:
        outdir = f"{READS_DIR}/qc/cdna"
    conda: "../env/qc.yaml"
    threads: 4
    shell:
        """
        fastqc --quiet --threads {threads} --outdir {params.outdir} {input.r2}
        """

# -------------------------------------------------------------------
# MultiQC Rules
# -------------------------------------------------------------------
rule multiqc_cb_umi:
    input:
        lambda wildcards: expand(
            "{rddir}/qc/cb_umi/{r}_" + cb_umi_suffix + "_fastqc.zip",
            rddir=READS_DIR, r=get_scrna_runinfo(wildcards)
        )
    output:
        report = f"{READS_DIR}/qc/multiqc_cb_umi.html"
    params:
        indir  = f"{READS_DIR}/qc/cb_umi",
        title  = "Cell Barcode & UMI QC Report"
    conda: "../env/qc.yaml"
    shell:
        """
        multiqc {params.indir} \
            --filename {output.report} \
            --title "{params.title}" \
            --force
        """


rule multiqc_cdna:
    input:
        lambda wildcards: expand(
            "{rddir}/qc/cdna/{r}_" + cdna_suffix + "_fastqc.zip",
            rddir=READS_DIR, r=get_scrna_runinfo(wildcards)
        )
    output:
        report = f"{READS_DIR}/qc/multiqc_cdna.html"
    params:
        indir  = f"{READS_DIR}/qc/cdna",
        title  = "cDNA Transcripts QC Report"
    conda: "../env/qc.yaml"
    shell:
        """
        multiqc {params.indir} \
            --filename {output.report} \
            --title "{params.title}" \
            --force
        """