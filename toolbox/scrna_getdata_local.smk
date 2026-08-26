import re
import csv
import shutil
from pathlib import Path

RAW_DIR = Path(READS_DIR) / "rawdata"
FASTQ_DIR = Path(READS_DIR) / "fastq"

import re
import csv
import shutil
from pathlib import Path

# ... standard imports ...

checkpoint process_illumina_local:
    input:
        raw_files = [str(p) for p in Path(READS_DIR, "rawdata").glob("*.fastq.gz")]
    output:
        csv = f"{READS_DIR}/local_runinfo.csv"
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