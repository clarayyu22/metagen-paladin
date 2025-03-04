#!/bin/bash
# script to run Paladin alignment
set -e

if [ "$#" -ne 6 ]; then
    echo "Usage: $0 <fasta_name> <fastq_path>"
    exit 1
fi

FASTA_DIR=$1
FASTA_NAME=$2
FASTQ_PATH=$3
THREADS=$4
OUTPUT_ALIGN_DIR=$5
DONE=$6

FASTA_PATH="${FASTA_DIR}${FASTA_NAME}/${FASTA_NAME}.fasta.txt"
SEQUENCE_ID=$(basename "$FASTQ_PATH" | sed 's/\..*$//')
OUT_FILE="${OUTPUT_ALIGN_DIR}${FASTA_NAME}_${SEQUENCE_ID}_out"

echo "starting alignment..."
./paladin/paladin align "$FASTA_PATH" "$FASTQ_PATH" -t "$THREADS" -o "$OUT_FILE"

if [ $? -eq 0 ]; then
    touch "$DONE"
    echo "Alignment completed: $OUT_FILE"
else
    echo "Alignment failed, output file not created!" >&2
    exit 1
fi



