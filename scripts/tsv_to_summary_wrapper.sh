#!/bin/bash
set -e

# Arguments
TSV_FILE=$1
OUTPUT_DIR=$2
folder=$(dirname "$TSV_FILE")
file_name_no_ext=$(basename "$TSV_FILE" .tsv)
file_name="$folder/$file_name_no_ext"

echo "Running Rscript to process TSV file:"
echo "  TSV: $TSV_FILE"
echo "  Output: $OUTPUT_DIR"

Rscript scripts/tsv_to_summary.R "$TSV_FILE" "$OUTPUT_DIR"


