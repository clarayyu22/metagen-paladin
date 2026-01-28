#!/bin/bash
set -e

# Arguments
BAM_FILE=$1
OUTPUT_DIR=$2
folder=$(dirname "$BAM_FILE")
file_name_no_ext=$(basename "$BAM_FILE" .bam)
file_name="$folder/$file_name_no_ext"

echo "Running Rscript to process BAM file:"
echo "  BAM: $BAM_FILE"
echo "  Output: $OUTPUT_DIR"

Rscript scripts/bam_to_tsv.R "$BAM_FILE" "$OUTPUT_DIR"

rm "${file_name}_sorted.bam.bai"
rm "${file_name}_sorted.bam"
