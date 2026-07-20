#!/bin/bash
# Bioinformatic Toolbox Master Runner

# Capture Project and Target
PRJ=$1 #capture 1st word after cmd
TARGET=$2 #capture 2nd word 

# Validation: Ensure project and target are provided
if [ -z "$PRJ" ] || [ -z "$TARGET" ]; then
    echo "Usage: bash run.sh [project_name] [target]"
    echo "e.g. bash run.sh ebola_2014 vcf_all"
    echo "Available targets: runinfo, qc, bam, vcf, vcf_all, rigid"
    exit 1
fi

# Execution with your standardized flags
snakemake \
    --config PRJNAME="$PRJ" \
    --use-conda \
    --cores 4 \
    --printshellcmds \
    $EXTRA_FLAGS \
    "$TARGET"