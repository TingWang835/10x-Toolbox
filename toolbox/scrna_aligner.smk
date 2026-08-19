localrules:  # no local rules here due to heavy ram requirement for RNA alignment

# =============================================================================
# Variables
# =============================================================================
v2_whitelist = config.get("V2_WHITELIST_NAME", "737K-august-2016.txt")
v3_whitelist = config.get("V3_WHITELIST_NAME", "3M-february-2018.txt")

# =============================================================================
# Helper Functions
# =============================================================================
def get_scrna_index_prefix(wildcards):
    """Retrieves the reference index path based on the aligner type."""
    scrna_aligner =  config.get("SCRNA_ALIGNER", "kb_python")
    if wildcards.scrna_aligner == "starsolo":
        return f"{REFS_DIR}/star_index"
    elif wildcards.scrna_aligner == "kb_python":
        return f"{REFS_DIR}/kb_python_index"
    return f"{REFS_DIR}/{wildcards.aligner}"

def get_10x_whitelist(wildcards):
    """Retrieve the correct 10x whitelist file path."""
    version = config.get("10X_VERSION", "v2")
    if version == "v2":
        return f"databases/10x_whitelists/v2/{v2_whitelist}"
    elif version == "v3":
        return f"databases/10x_whitelists/v3/{v3_whitelist}"
    else:
        raise ValueError(f"Unsupported 10x version {version}")

def get_10x_lengths(wildcards):
    """Retrieve CB and UMI lengths based on 10x version."""
    version = config.get("10X_VERSION", "v2")
    if version == "v2":
        return {"umi_len": 10, "cb_len": 16}
    elif version == "v3":
        return {"umi_len": 12, "cb_len": 16}
    else:
        raise ValueError(f"Unsupported 10x version {version}")

def get_h5ads_for_group(wildcards):
    """
    Finds all runs belonging to wildcards.sample and returns
    paths to their individual .h5ad files.
    """
    runinfo_df = pd.read_csv(f"{READS_DIR}/runinfo_scrna.csv")
    
    # Filter runs for current sample
    runs = runinfo_df[runinfo_df["SampleName"] == wildcards.sample]["Run"].tolist()
    
    # Use wildcards.aligner dynamically
    return [
        f"{READS_DIR}/h5ad/{wildcards.aligner}/individuals/{run}.h5ad" 
        for run in runs
    ]

def get_aligner_output(wildcards):
    """Retrieves Matrix, Gene file and Barcodes based on aligner type."""
    scrna_aligner = config.get("SCRNA_ALIGNER", "kb_python").lower()

    if scrna_aligner in {"starsolo", "star"}:
        return {
            "matrix": f"{READS_DIR}/aligner/scrna/starsolo/{wildcards.run}/Solo.out/Gene/filtered/matrix.mtx.gz",
            "barcodes": f"{READS_DIR}/aligner/scrna/starsolo/{wildcards.run}/Solo.out/Gene/filtered/barcodes.tsv.gz",
            "genes": f"{READS_DIR}/aligner/scrna/starsolo/{wildcards.run}/Solo.out/Gene/filtered/features.tsv.gz",
        }
    elif scrna_aligner in {"kb_python", "kb", "kallisto"}:
        return {
            "matrix": f"{READS_DIR}/aligner/scrna/kb/{wildcards.run}/counts_filtered/cells_x_genes.mtx",
            "barcodes": f"{READS_DIR}/aligner/scrna/kb/{wildcards.run}/counts_filtered/cells_x_genes.barcodes.txt",
            "genes": f"{READS_DIR}/aligner/scrna/kb/{wildcards.run}/counts_filtered/cells_x_genes.genes.txt",
        }
    else:
        raise ValueError(f"Unsupported SCRNA_ALIGNER target: {scrna_aligner}")
# =============================================================================
# Index Rules
# =============================================================================
rule star_index:
    input: 
        unpack(get_refs)
    output: 
        directory(f"{REFS_DIR}/star_index")
    log: 
        f"{LOG_DIR}/aligner/scrna/star_index.log"
    conda: 
        "../env/scrna_aligner.yaml"
    params: 
        sa = config.get("STAR_SA", "14") 
    threads: 4
    shell:
        """
        mkdir -p {output}
        STAR --runMode genomeGenerate \
             --genomeDir {output} \
             --genomeFastaFiles {input.fasta} \
             --sjdbGTFfile {input.gtf} \
             --runThreadN {threads} \
             --genomeSAindexNbases {params.sa} > {log} 2>&1
        """

rule kb_index:
    """
    Generates single-cell kallisto index and transcript-to-gene mapping using kb-python.
    """
    input:
        unpack(get_refs)
    output:
        idx = f"{REFS_DIR}/kb_index/transcriptome.idx",
        t2g = f"{REFS_DIR}/kb_index/transcripts_to_genes.txt"
    log:
        f"{LOG_DIR}/aligner/scrna/kb_index.log"
    conda:
        "../env/scrna_aligner.yaml"
    threads: 4
    shell:
        """
        kb ref -i {output.idx} \
               -g {output.t2g} \
               -f1 {REFS_DIR}/kb_index/transcriptome.fa \
               {input.fasta} {input.gtf} > {log} 2>&1
        """

rule download_10x_whitelist:
    """
    Downloads 10x Genomics cell barcode whitelists for v2 or v3 chemistry.
    """
    output:
        v2 = f"databases/10x_whitelists/v2/{v2_whitelist}",
        v3 = f"databases/10x_whitelists/v3/{v3_whitelist}"
    log:
        f"{LOG_DIR}/aligner/scrna/download_whitelists.log"
    shell:
        """
        exec 2> {log}
        mkdir -p $(dirname {output.v2}) $(dirname {output.v3})

        # Download v2 whitelist if not present
        if [ ! -f {output.v2} ]; then
            wget -qO- https://raw.githubusercontent.com/Lab-of-Adaptive-Immunity/cr_whitelists/main/{v2_whitelist} > {output.v2}
        fi

        # Download and decompress v3 whitelist if not present
        if [ ! -f {output.v3} ]; then
            wget -qO- https://raw.githubusercontent.com/Lab-of-Adaptive-Immunity/cr_whitelists/main/{v3_whitelist}.gz | zcat > {output.v3}
        fi
        """

# =============================================================================
# Alignment Rules
# =============================================================================

rule starsolo_align:
    """
    Aligns single-cell paired FASTQs to genome index using STARsolo, 
    performing cell barcode error correction, UMI deduplication, and feature quantification.
    """
    input:
        r1 = ancient(f"{READS_DIR}/fastqs/{{run}}_{cb_umi_suffix}.fastq.gz"),
        r2 = ancient(f"{READS_DIR}/fastqs/{{run}}_{cdna_suffix}.fastq.gz"),
        index = f"{REFS_DIR}/star_index",
        whitelist = get_10x_whitelist
    output:
        bam = f"{READS_DIR}/aligner/scrna/starsolo/{{run}}/Aligned.sortedByCoord.out.bam",
        summary = f"{READS_DIR}/aligner/scrna/starsolo/{{run}}/Solo.out/Gene/Summary.csv",
        matrix = f"{READS_DIR}/aligner/scrna/starsolo/{{run}}/Solo.out/Gene/filtered/matrix.mtx.gz",
        barcodes = f"{READS_DIR}/aligner/scrna/starsolo/{{run}}/Solo.out/Gene/filtered/barcodes.tsv.gz",
        features = f"{READS_DIR}/aligner/scrna/starsolo/{{run}}/Solo.out/Gene/filtered/features.tsv.gz"
    log:
        f"{LOG_DIR}/aligner/scrna/star/{{run}}.log"
    conda:
        "../env/scrna_aligner.yaml"
    threads: 4
    resources:
        mem_mb = 15000
    params:
        out_prefix = f"{READS_DIR}/aligner/scrna/starsolo/{{run}}/", 
        cb_len = lambda wc: get_10x_lengths(wc)["cb_len"],
        umi_len = lambda wc: get_10x_lengths(wc)["umi_len"],
        solo_type = config.get("SOLO_TYPE", "CB_UMI_Simple"),
        barcode_read_len = 0, # to prevent error due to uneven CB + UMI length
        cell_filter_type = config.get("SOLO_CELL_FILTER", "CellRanger2.2"),
        expected_cells = config.get("SOLO_EXPECTED_CELLS", 3000),
        max_percentile = config.get("SOLO_MAX_PERCENTILE", 0.99),
        max_min_ratio = config.get("SOLO_MAX_MIN_RATIO", 10)
    shell:
        """
        exec 2> {log}

        STAR --runThreadN {threads} \
             --genomeDir {input.index} \
             --readFilesIn {input.r2} {input.r1} \
             --readFilesCommand zcat \
             --outFileNamePrefix {params.out_prefix} \
             --outSAMtype BAM SortedByCoordinate \
             --soloType {params.solo_type} \
             --soloCBwhitelist {input.whitelist} \
             --soloCBlen {params.cb_len} \
             --soloUMIlen {params.umi_len} \
             --soloBarcodeReadLength {params.barcode_read_len}\
             --soloCellFilter {params.cell_filter_type} {params.expected_cells} {params.max_percentile} {params.max_min_ratio} \
             --soloFeatures Gene \
             --soloOutFileNames Solo.out/ gene.tsv barcodes.tsv matrix.mtx

        # change gene.tsv to features.tsv for better downstream compatibility
        if [ -f {params.out_prefix}Solo.out/Gene/filtered/gene.tsv ]; then
            mv {params.out_prefix}Solo.out/Gene/filtered/gene.tsv {params.out_prefix}Solo.out/Gene/filtered/features.tsv
        fi

        gzip -f {params.out_prefix}Solo.out/Gene/filtered/*.tsv || true
        gzip -f {params.out_prefix}Solo.out/Gene/filtered/*.mtx || true
        """

rule kb_python_align:
    """
    Aligns single-cell paired FASTQs and quantifies gene expression using kb-python (kallisto + bustools).
    """
    input:
        r1 = ancient(f"{READS_DIR}/fastqs/{{run}}_{cb_umi_suffix}.fastq.gz"),
        r2 = ancient(f"{READS_DIR}/fastqs/{{run}}_{cdna_suffix}.fastq.gz"),
        idx = f"{REFS_DIR}/kb_index/transcriptome.idx",
        t2g = f"{REFS_DIR}/kb_index/transcripts_to_genes.txt"
    output:
        matrix = f"{READS_DIR}/aligner/scrna/kb/{{run}}/counts_filtered/cells_x_genes.mtx",
        barcodes = f"{READS_DIR}/aligner/scrna/kb/{{run}}/counts_filtered/cells_x_genes.barcodes.txt",
        genes = f"{READS_DIR}/aligner/scrna/kb/{{run}}/counts_filtered/cells_x_genes.genes.txt"
    log:
        f"{LOG_DIR}/aligner/scrna/kb/{{run}}.log"
    conda:
        "../env/scrna_aligner.yaml"
    threads: 4
    resources:
        mem_mb = 15000
    params:
        out_dir = f"{READS_DIR}/aligner/scrna/kb/{{run}}",
        technology = config.get("KB_TECHNOLOGY", "10XV2")
    shell:
        """
        exec 2> {log}

        kb count \
            -i {input.idx} \
            -g {input.t2g} \
            -x {params.technology} \
            -o {params.out_dir} \
            -t {threads} \
            --filter bustools \
            {input.r1} {input.r2}
        """



# =============================================================================
# Anndata (.h5ad) conversion
# =============================================================================


ALIGNER_TYPE = config.get("SCRNA_ALIGNER", "kb_python").lower()

rule mtx_to_h5ad:
    input:
        unpack(get_refs),
        unpack(get_aligner_output)
    output:
        h5ad = temp(f"{READS_DIR}/h5ad/{ALIGNER_TYPE}/individuals/{{run}}.h5ad")
    params: 
        aligner = ALIGNER_TYPE
    conda:
        "../env/scrna_aligner.yaml"
    script:
        "scripts/scrna_mtx_to_h5ad.py"


rule h5ad_group:
    """
    Merges individual run-level AnnData (.h5ad) files into a single sample AnnData object.
    """
    input:
        h5ads = get_h5ads_for_group
    output:
        grouped_h5ad = temp(f"{READS_DIR}/h5ad/{{aligner}}/grouped/{{sample}}_grouped.h5ad")
    log:
        f"{LOG_DIR}/h5ad_group/{{aligner}}/{{sample}}.log"
    conda:
        "../env/scrna_aligner.yaml"
    script:
        "scripts/anndata_grouped.py"