#!/bin/bash
set -euo pipefail

FASTA=$1
FASTQ=$2
OUTDIR=$3
NCORES=$4
DONE_MARKER=$5

MAX_RETRIES=3
TIMEOUT_MINUTES=360   # Max time allowed per attempt
SUCCESS=0
ATTEMPT=1

while [[ $ATTEMPT -le $MAX_RETRIES ]]; do
    echo "Attempt $ATTEMPT: running paladin align at $(date)"
    timeout ${TIMEOUT_MINUTES}m bash scripts/align.sh "$FASTA" "$FASTQ" "$OUTDIR" "$NCORES" && SUCCESS=1 && break

    echo "paladin_align failed or timed out on attempt $ATTEMPT"
    sleep $(( RANDOM % 60 + 30 ))  # backoff: wait 30–90 sec
    ATTEMPT=$((ATTEMPT + 1))
done

if [[ $SUCCESS -eq 1 ]]; then
    echo "paladin_align succeeded."
    touch "$DONE_MARKER"
    rm -f "$FASTQ"
else
    echo "paladin_align failed after $MAX_RETRIES attempts."
    exit 1
fi
