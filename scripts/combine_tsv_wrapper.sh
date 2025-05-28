#!/bin/bash
set -e

INPUT_DIR=$1
OUTPUT_DIR=$2

echo "Running R summary script..."
Rscript scripts/combine_tsv.R "$INPUT_DIR" "$OUTPUT_DIR"
echo "Done."