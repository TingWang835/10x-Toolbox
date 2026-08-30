localrules: 


    
# =============================================================================
# Quality Control
# =============================================================================
# -------------------------------------------------------------------
# FastQC Rules with Config-Driven Read Suffixes
# -------------------------------------------------------------------
rule fastqc_cb_umi:
    input:
        r1 = ancient(f"{READS_DIR}/fastqs/{{run}}_{cb_umi_suffix}.fastq.gz")
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
        r2 = ancient(f"{READS_DIR}/fastqs/{{run}}_{cdna_suffix}.fastq.gz")
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