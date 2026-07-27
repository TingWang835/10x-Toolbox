# Snakefile
import os
import re
import csv

RAW_DIR = f"{READS_DIR}/rawdata"
LINK_DIR = f"{READS_DIR}/fastq"

rule symlink_illumina:
    input:
        # Detects all .fastq.gz files in raw directory
        raw_files = glob_wildcards(f"{RAW_DIR}/{{filename}}.fastq.gz").filename
    output:
        csv = "local_runinfo.csv"
    run:
        os.makedirs(LINK_DIR, exist_ok=True)
        
        # Illumina regex: SampleName_S<sample>_L<lane>_<read>_<chunk>.fastq.gz
        pattern = re.compile(r"^(?P<sample>.+)_S(?P<s_num>\d+)_L(?P<lane>\d+)_(?P<read>[RI]\d+)_(?P<chunk>\d+)$")
        
        # Read code mapping: R1/R2 -> _1/_2, I1/I2 -> _3/_4 (adjust if needed)
        read_map = {
            "R1": "1",
            "R2": "2",
            "R3": "3",
            "I1": "4",
            "I2": "5"
        }

        samples_seen = set()
        runinfo_data = []

        for fn in input.raw_files:
            match = pattern.match(fn)
            if not match:
                continue  # Skip files not matching standard Illumina layout
            
            d = match.groupdict()
            sample_name = d["sample"]
            read_code = read_map.get(d["read"], d["read"])
            
            # Target symlink name: e.g., S1L001001_1.fastq.gz
            link_basename = f"S{d['s_num']}L{d['lane']}{d['chunk']}_{read_code}.fastq.gz"
            
            target_link = os.path.join(LINK_DIR, link_basename)
            source_file = os.path.abspath(os.path.join(RAW_DIR, f"{fn}.fastq.gz"))

            # Create symlink if it doesn't already exist
            if not os.path.exists(target_link):
                os.symlink(source_file, target_link)

            # Record run info metadata (one entry per sample)
            run_id = f"S{d['s_num']}L{d['lane']}{d['chunk']}"
            if run_id not in samples_seen:
                samples_seen.add(run_id)
                runinfo_data.append({
                    "Run": run_id,
                    "SampleName": sample_name,
                    "LibraryLayout": "PAIRED" if "R2" in read_map else "SINGLE",
                    "Platform": "ILLUMINA"
                })

        # Write output CSV
        headers = ["Run", "SampleName", "LibraryLayout", "Platform"]
        with open(output.csv, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=headers)
            writer.writeheader()
            writer.writerows(runinfo_data)