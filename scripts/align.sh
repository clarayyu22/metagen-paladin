#!/bin/bash
set -e

# Check if we have the correct number of arguments
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <fasta_name> <fastq_path> <output_prefix> <threads>"
    exit 1
fi

FASTA_REF_DIR=$1
FASTQ_PATH=$2
OUT_FILE=$3
THREADS=$4

echo "starting alignment..."

# Check if the output file already exists
if [ ! -f "${OUT_FILE}.sam" ]; then
    paladin align -a -t "$THREADS" "$FASTA_REF_DIR" "$FASTQ_PATH" > "${OUT_FILE}.sam"
    # paladin align -a -t "$THREADS" -o "$OUT_FILE" "$FASTA_REF_DIR" "$FASTQ_PATH"
    
    if [ $? -eq 0 ]; then
        echo "PALADIN alignment completed successfully"
    else
        echo "Alignment failed, output file not created!" >&2
        exit 1
    fi
else
    echo "Output file ${OUT_FILE}.sam already exists, skipping alignment"
fi
