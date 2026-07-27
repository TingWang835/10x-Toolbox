import sys
sys.path.append("toolbox/scripts")
from uscs_chromosomes import CHROM_DICTS, get_ucsc_mapping

localrules: download_ensembl_refs, download_ncbi_refs, uscs_translator, gff_to_gtf, samtools_faidx

# =============================================================================
# Helper Function
# =============================================================================
def get_raw_refs(wildcards):
    """
    Dynamically routes pipeline requirements to raw downloaded reference files.
    """
    if source == "ensembl":
        return {
            "fasta": f"{REFS_DIR}/{species_cap}.{assembly}.{release}.raw.fa",
            "gff": f"{REFS_DIR}/{species_cap}.{assembly}.{release}.raw.gff3"
        }
    elif source == "ncbi":
        return {
            "fasta": f"{REFS_DIR}/{acc}.raw.fa",
            "gff": f"{REFS_DIR}/{acc}.raw.gff3"
        }
    else:
        raise ValueError(f"Unsupported REF_SOURCE target: {source}")


# =============================================================================
# Rules
# =============================================================================
rule download_ncbi_refs:
    output:
        fasta = f"{REFS_DIR}/{acc}.raw.fa", 
        gff   = f"{REFS_DIR}/{acc}.raw.gff3" 
    conda: "../env/get_refs.yaml" 
    log: f"{LOG_DIR}/get_refs/download_ncbi_ref.log" 
    shell:
        """
        mkdir -p {REFS_DIR} 
        exec 2> {log}

        TMP_ZIP="{REFS_DIR}/{acc}_dataset.zip"

        # Download zip (include fa and gff) via NCBI datasets API CLI
        datasets download genome accession {acc} \
            --include genome,gff3 \
            --filename $TMP_ZIP

        # Unzip package
        unzip -q -o $TMP_ZIP -d {REFS_DIR}/tmp_out

        # Move fa and gff to path
        mv {REFS_DIR}/tmp_out/ncbi_dataset/data/{acc}/*_genomic.fna {output.fasta}
        mv {REFS_DIR}/tmp_out/ncbi_dataset/data/{acc}/genomic.gff {output.gff}

        # Purge temporary packages
        rm -rf $TMP_ZIP {REFS_DIR}/tmp_out
        """


rule download_ensembl_refs:
    output:
        fasta = f"{REFS_DIR}/{species_cap}.{assembly}.{release}.raw.fa",
        gff   = f"{REFS_DIR}/{species_cap}.{assembly}.{release}.raw.gff3"
    log: f"{LOG_DIR}/get_refs/download_ensembl_ref.log"
    conda: "../env/get_refs.yaml"
    params:
        url_base = f"ftp://ftp.ensemblgenomes.org/pub/metazoa/release-{release}" if species_low == "caenorhabditis_elegans" else f"ftp://ftp.ensembl.org/pub/release-{release}"
    shell:
        """
        exec 2> {log}

        # 1. Stream, decompress, and rename the Genomic FASTA file
        wget -qO- {params.url_base}/fasta/{species_low}/dna/{species_cap}.{assembly}.dna.toplevel.fa.gz | \
            gunzip -c > {output.fasta}

        # 2. Stream, decompress, and rename the GFF3 Annotation file
        wget -qO- {params.url_base}/gff3/{species_low}/{species_cap}.{assembly}.{release}.gff3.gz | \
            gunzip -c > {output.gff}
        """


rule uscs_translator:
    """Translates FASTA and GFF3 files to UCSC chromosome naming convention."""
    input:
        unpack(get_raw_refs)
    output:
        fasta = f"{REFS_DIR}/{species_cap}.{assembly}.{release}.fa" if source == "ensembl" else f"{REFS_DIR}/{acc}.fa",
        gff   = f"{REFS_DIR}/{species_cap}.{assembly}.{release}.gff3" if source == "ensembl" else f"{REFS_DIR}/{acc}.gff3"
    params:
        mapping = get_ucsc_mapping(species_low)
    log:
        f"{LOG_DIR}/get_refs/uscs_translator.log"
    run:
        import os
        
        # 1. Process FASTA File
        with open(input.fasta, "r") as infile, open(output.fasta, "w") as outfile:
            for line in infile:
                if line.startswith(">"):
                    parts = line.strip().split()
                    orig_chr = parts[0].replace(">", "")
                    new_chr = params.mapping.get(orig_chr, orig_chr)
                    
                    rest = " " + " ".join(parts[1:]) if len(parts) > 1 else ""
                    outfile.write(f">{new_chr}{rest}\n")
                else:
                    outfile.write(line)

        # 2. Process GFF3 File
        with open(input.gff, "r") as infile, open(output.gff, "w") as outfile:
            for line in infile:
                if line.startswith("#"):
                    if line.startswith("##sequence-region"):
                        parts = line.strip().split()
                        if len(parts) >= 2:
                            parts[1] = params.mapping.get(parts[1], parts[1])
                            line = " ".join(parts) + "\n"
                    outfile.write(line)
                else:
                    parts = line.split("\t")
                    if len(parts) > 0:
                        parts[0] = params.mapping.get(parts[0], parts[0])
                    outfile.write("\t".join(parts))

        # 3. Clean up raw input files
        if os.path.exists(input.fasta):
            os.remove(input.fasta)
        if os.path.exists(input.gff):
            os.remove(input.gff)


rule gff_to_gtf:
    """Locally converts GFF3 to GTF format via gffread."""
    input:
        gff = lambda wildcards: get_refs(wildcards)["gff"]
    output:
        gtf = f"{REFS_DIR}/{species_cap}.{assembly}.{release}.gtf" if source == "ensembl" else f"{REFS_DIR}/{acc}.gtf"
    log: f"{LOG_DIR}/get_refs/gff_to_gtf.log"
    conda: "../env/get_refs.yaml"
    shell:
        """
        gffread {input.gff} -T -o {output.gtf} > {log} 2>&1
        """


rule samtools_faidx:
    """Generates .fai index tracks on reference fasta structures."""
    input:
        fasta = lambda wildcards: get_refs(wildcards)["fasta"]
    output:
        fai = f"{REFS_DIR}/{species_cap}.{assembly}.{release}.fa.fai" if source == "ensembl" else f"{REFS_DIR}/{acc}.fa.fai"
    log: f"{LOG_DIR}/get_refs/samtools_faidx.log"
    conda: "../env/get_refs.yaml"
    shell:
        """
        samtools faidx {input.fasta} > {log} 2>&1
        """

checkpoint fetch_sra_metadata:
    """
    Fetches SRA runinfo CSV using NCBI Entrez Direct (esearch + efetch).
    """
    output:
        csv = f"{READS_DIR}/runinfo.csv"
    log:
        f"{LOG_DIR}/getdata/fetch_sra_metadata.log"
    conda:
        "../env/get_refs.yaml"
    params:
        prj = config.get("PRJNUMBER")
    shell:
        """
        exec 2> {log}
        mkdir -p {READS_DIR}
        mkdir -p $(dirname {log})

        # Fetch SRA runinfo directly in CSV format
        esearch -db sra -query "{params.prj}" | efetch -format runinfo > {READS_DIR}/raw_runinfo.csv

        # Remove empty lines or header collisions if any
        head -n 1 {READS_DIR}/raw_runinfo.csv > {output.csv}
        grep -v "^Run," {READS_DIR}/raw_runinfo.csv >> {output.csv} || true
        rm {READS_DIR}/raw_runinfo.csv
        """