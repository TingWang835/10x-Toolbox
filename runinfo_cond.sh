#!/bin/bash

# Ensure a file was provided
if [ -z "$1" ]; then
    echo "Usage: $0 input.csv"
    exit 1
fi
PRJNAME=$1
TYPE=$2
INPUT_FILE="reads/${PRJNAME}/runinfo_${TYPE}.csv"

OUTPUT_FILE="reads/${PRJNAME}/runinfo_${TYPE}_cat.csv"

awk -F, '
BEGIN {
    OFS=","
    # --- DEFINE YOUR CATEGORIES HERE ---
    # cond[sample_name]="treatment_condition"
    cond["GSM3618670"]="300min"
    cond["GSM3618671"]="400min"
    cond["GSM3618672"]="500min_1"
    cond["GSM3618673"]="500min_2"
    cond["GSM3618674"]="mix1"
    cond["GSM3618675"]="mix2"
    cond["GSM3618676"]="mix3"
    # -------------------------------------
} 
NR==1 {
    print $0, "condition"; next
} 
{
    print $0, (cond[$2] ? cond[$2] : "Unassigned")
}' "$INPUT_FILE" > "$OUTPUT_FILE"

echo "Success! Saved to: $OUTPUT_FILE"