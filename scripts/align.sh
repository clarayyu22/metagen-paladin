#!/bin/bash
set -e

# Check if we have the correct number of arguments
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <fasta_name> <fastq_path> <output_prefix> <threads>"
    exit 1
fi

FASTA_PATH=$1
FASTQ_PATH=$2
OUT_FILE=$3
THREADS=$4

echo "starting alignment..."

# Check if the output file already exists
if [ ! -f "${OUT_FILE}.bam" ]; then
    paladin align -a -t "$THREADS" "$FASTA_PATH" "$FASTQ_PATH" > "${OUT_FILE}.sam"
    samtools view --threads "$THREADS" -S -b "${OUT_FILE}.sam" > "${OUT_FILE}.bam"
    #rm "${OUT_FILE}.sam"
    
    if [ $? -eq 0 ]; then
        echo "PALADIN alignment completed successfully"
    else
        echo "Alignment failed, output file not created!" >&2
        exit 1
    fi
else
    echo "Output file ${OUT_FILE}.bam already exists, skipping alignment"
fi
