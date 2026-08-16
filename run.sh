#!/bin/bash
# Bioinformatic Toolbox Master Runner

# Capture Project, Target, and any optional trailing flags
PRJ=$1
TARGET=$2
EXTRA_FLAGS=${@:3}

# Validation: Ensure project and target are provided
if [ -z "$PRJ" ] || [ -z "$TARGET" ]; then
    echo "Usage: bash run.sh [project_name] [target] [extra_flags...]"
    echo "e.g. bash run.sh ebola_2014 vcf_all --dry-run"
    echo "Available targets: runinfo, qc, bam, vcf, vcf_all, rigid"
    exit 1
fi

# Execution with your standardized flags
snakemake \
    --config PRJNAME="$PRJ" \
    --use-conda \
    --cores 12 \
    --printshellcmds \
    $EXTRA_FLAGS \
    "$TARGET"