import re
import csv
import pandas as pd
import scanpy as sc
from pathlib import Path
import shutil

ALIGNER = config.get("SCRNA_ALIGNER", "kb_python")


# =============================================================================
# SRA Download Reads
# =============================================================================
checkpoint fetch_sra_metadata:
    """
    Fetches SRA runinfo CSV using NCBI Entrez Direct (esearch + efetch).
    """
    output:
        csv = f"{READS_DIR}/sra_original_runinfo.csv"
    log:
        f"{LOG_DIR}/getdata/fetch_sra_metadata.log"
    conda:
        "../env/scrna_getdata.yaml"
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


rule parse_scrna_metadata:
    """
    Parses global runinfo.csv to filter single-cell RNA-seq runs,
    maps GSM accessions to SRR runs, and allows sample filtering.
    """
    input:
        csv = f"{READS_DIR}/sra_original_runinfo.csv"
    output:
        csv = f"{READS_DIR}/sra_runinfo.csv"
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
        sra_dir = directory(f"{READS_DIR}/sra/{{run}}")
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
    Parse .sra using fasterq-dump (multi-threaded) and compress with pigz.
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

        TMP_DIR=$(mktemp -d -p {READS_DIR}/fastqs tmp_{wildcards.run}_XXXXXX)

        fasterq-dump {input.sra_dir} \
            --split-files \
            --include-technical \
            --threads {threads} \
            --outdir {READS_DIR}/fastqs \
            --temp $TMP_DIR

        rm -rf $TMP_DIR

        pigz -f -p {threads} {READS_DIR}/fastqs/{wildcards.run}_*.fastq
        """


# =============================================================================
# Local fastq Ingest
# =============================================================================

checkpoint process_illumina_local:
    input:
        raw_files = [str(p) for p in Path(READS_DIR, "rawdata").glob("*.fastq.gz")]
    output:
        csv = f"{READS_DIR}/local_fastq_runinfo.csv"
    conda:
        "../env/scrna_getdata.yaml"
    run:
        RAW_DIR = Path(READS_DIR) / "rawdata"
        FASTQ_DIR = Path(READS_DIR) / "fastq"
        FASTQ_DIR.mkdir(parents=True, exist_ok=True)

        pattern = re.compile(
            r"^(?P<sample>.+)_S[0-9oO]+_L(?P<lane>[0-9oO]+)_(?P<read>[RI][1234])_[0-9oO]+\.fastq\.gz$"
        )
        
        read_map = {"R1": "1", "R2": "2", "I1": "3", "I2": "4"}
        file_groups = {}

        # 1. Group chunk files by Sample + Lane + Read
        for file_path in map(Path, input.raw_files):
            match = pattern.match(file_path.name)
            if not match:
                continue
            
            d = match.groupdict()
            sample_name = d["sample"]
            clean_lane = f"{int(re.sub(r'[oO]', '0', d['lane'])):03d}"
            read_code = read_map[d["read"]]
            
            group_key = (sample_name, clean_lane, read_code)
            file_groups.setdefault(group_key, []).append(file_path)

        # 2. Concatenate or Symlink Chunks & Track Full Run Names
        runs_seen = {} # Maps full_run_id -> sample_name
        
        for (sample_name, lane, read_code), sources in file_groups.items():
            sources.sort()
            
            # Combine Sample and Lane into the full Run ID (e.g., NA12878_L001)
            full_run_id = f"{sample_name}_L{lane}"
            target_file = FASTQ_DIR / f"{full_run_id}_{read_code}.fastq.gz"
            
            if not target_file.exists():
                if len(sources) == 1:
                    target_file.symlink_to(sources[0].resolve())
                else:
                    with open(target_file, "wb") as outfile:
                        for src in sources:
                            with open(src, "rb") as infile:
                                shutil.copyfileobj(infile, outfile)

            runs_seen[full_run_id] = sample_name

        # 3. Export complete Run IDs to local_runinfo.csv
        headers = ["Run", "SampleName", "LibraryLayout", "Platform"]
        runinfo_data = [
            {
                "Run": run_id, 
                "SampleName": sample_name, 
                "LibraryLayout": "PAIRED", 
                "Platform": "ILLUMINA"
            } 
            for run_id, sample_name in sorted(runs_seen.items())
        ]
        
        with open(output.csv, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=headers)
            writer.writeheader()
            writer.writerows(runinfo_data)


# =============================================================================
# Local Batch Corrected h5ad Ingest
# =============================================================================
rule ingest_concat_h5ad:
    input:
        raw_concat = f"{READS_DIR}/rawdata/pca_concat_batch_embed.h5ad"
    output:
        h5ad = f"{READS_DIR}/h5ad/{ALIGNER}/pca_concat_batch_embed.h5ad",
        csv = f"{READS_DIR}/local_concat_runinfo.csv"
    run:
        Path(output.h5ad).parent.mkdir(parents=True, exist_ok=True)
        # 1. Copy file to standard workflow location
        shutil.copyfile(input.raw_concat, output.h5ad)
        
        # 2. Extract quick cell count runinfo from metadata
        adata = sc.read_h5ad(output.h5ad, backed="r")
        n_cells = adata.n_obs
        n_genes = adata.n_vars
        batches = adata.obs["batch"].unique().tolist() if "batch" in adata.obs else ["Combined"]
        
        run_df = pd.DataFrame({
            "Run": batches,
            "SampleName": batches,
            "CellCount": n_cells // len(batches),
            "GeneCount": n_genes,
            "Datasource": "H5AD_CONCAT"
        })
        run_df.to_csv(output.csv, index=False)
        adata.file.close()

# =============================================================================
# Local Individual h5ad Ingest
# =============================================================================
rule ingest_indiv_h5ad:
    input:
        raw_dir = f"{READS_DIR}/rawdata"
    output:
        csv = f"{READS_DIR}/local_indiv_runinfo.csv"
    run:
        dest_dir = Path(f"{READS_DIR}/h5ad/{ALIGNER}/grouped")
        dest_dir.mkdir(parents=True, exist_ok=True)
        runinfo_data = []
        
        # Glob dynamically during rule execution
        h5ad_files = list(Path(input.raw_dir).glob("*.h5ad"))
        
        for p in h5ad_files:
            sample_id = p.stem
            dest = dest_dir / p.name
            shutil.copyfile(p, dest)
            
            adata = sc.read_h5ad(p, backed="r")
            runinfo_data.append({
                "Run": sample_id,
                "SampleName": sample_id,
                "CellCount": adata.n_obs,
                "GeneCount": adata.n_vars,
                "Datasource": "H5AD_INDIVIDUAL"
            })
            adata.file.close()
            
        pd.DataFrame(runinfo_data).to_csv(output.csv, index=False)

# =============================================================================
# Local Raw h5 & MTX Ingest
# =============================================================================
rule ingest_raw_h5_mtx:
    input:
        raw_dir = f"{READS_DIR}/rawdata"
    output:
        csv = f"{READS_DIR}/local_h5mtx_runinfo.csv"
    run:
        out_indiv = Path(f"{READS_DIR}/h5ad/{ALIGNER}/grouped")
        out_indiv.mkdir(parents=True, exist_ok=True)
        
        runinfo_data = []
        
        def convert_to_h5ad(src_path, sample_id):
            src_p = Path(src_path)
        
            if str(src_path).endswith(".h5"):
                adata = sc.read_10x_h5(src_path)
            else:
                # Check for compressed or uncompressed filenames
                matrix_file = src_p / "matrix.mtx.gz" if (src_p / "matrix.mtx.gz").exists() else src_p / "matrix.mtx"
                barcodes_file = src_p / "barcodes.tsv.gz" if (src_p / "barcodes.tsv.gz").exists() else src_p / "barcodes.tsv"
                features_file = src_p / "features.tsv.gz" if (src_p / "features.tsv.gz").exists() else (
                            src_p / "genes.tsv.gz" if (src_p / "genes.tsv.gz").exists() else src_p / "genes.tsv")
            
                # 1. Read sparse matrix (transpose so cells are rows, genes are columns)
                adata = sc.read_mtx(matrix_file)
            
                # 2. Assign cell barcodes
                barcodes_df = pd.read_csv(barcodes_file, header=None, sep="\t")
                adata.obs_names = barcodes_df[0].astype(str).values
            
                # 3. Assign gene names (supports 1, 2, or 3 column features.tsv)
                features_df = pd.read_csv(features_file, header=None, sep="\t")
                if features_df.shape[1] >= 2:
                    adata.var_names = features_df[1].astype(str).values
                    adata.var["gene_ids"] = features_df[0].astype(str).values
                else:
                    adata.var_names = features_df[0].astype(str).values
                
                adata.var_names_make_unique()
            
            out_file = out_indiv / f"{sample_id}.h5ad"
            adata.write_h5ad(out_file, compression="gzip")
        
            return {
                "Run": sample_id,
                "SampleName": sample_id,
                "CellCount": adata.n_obs,
                "GeneCount": adata.n_vars,
                "Datasource": "RAW_H5_MTX"
            }

        # 1. Convert MTX Subdirectories (Layout A)
        for mtx in Path(input.raw_dir).rglob("matrix.mtx*"):
            sample_dir = mtx.parent
            sample_id = sample_dir.name if sample_dir.name not in ["filtered_feature_bc_matrix", "counts"] else sample_dir.parent.name
            runinfo_data.append(convert_to_h5ad(sample_dir, sample_id))

        # 2. Convert Prefixed MTX Files in Root rawdata (Layout B)
        for mtx in Path(input.raw_dir).glob("*_matrix.mtx*"):
            sample_id = mtx.name.split("_matrix.mtx")[0]
            adata = sc.read_10x_mtx(input.raw_dir, prefix=f"{sample_id}_", make_unique=True)
            out_file = out_indiv / f"{sample_id}.h5ad"
            adata.write_h5ad(out_file, compression="gzip")
            runinfo_data.append({
                "Run": sample_id,
                "SampleName": sample_id,
                "CellCount": adata.n_obs,
                "GeneCount": adata.n_vars,
                "Datasource": "RAW_H5_MTX"
            })
            
        # 3. Convert standalone .h5 files
        for h5 in Path(input.raw_dir).glob("*.h5"):
            runinfo_data.append(convert_to_h5ad(h5, h5.stem))
            
        pd.DataFrame(runinfo_data).to_csv(output.csv, index=False)



