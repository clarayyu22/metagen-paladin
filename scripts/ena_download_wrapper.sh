#!/bin/bash
set -euo pipefail

SAMPLE=$1
OUTDIR=$2
NCORES=$3
FORWARD_FASTQ=$4

# Random sleep between 50-70 seconds
SLEEP_TIME=$(( (RANDOM % 20) + 50 ))
echo "Sleeping for $SLEEP_TIME seconds before starting download..."
sleep $SLEEP_TIME

# Check ENA server connection
check_server() {
    echo "Checking ENA FTP server..."
    if ping -c 1 ftp.sra.ebi.ac.uk &>/dev/null; then
        echo "Server is reachable."
        return 0
    else
        echo "Server unreachable. Aborting."
        return 1
    fi
}

# Retry logic
MAX_RETRIES=10
TIMEOUT_SECONDS=600
RETRY_DELAY_MIN=30
RETRY_DELAY_MAX=180
MAX_TOTAL_SECONDS=3600  # One hour max retry window

ATTEMPT=1
SUCCESS=0
START_TIME=$(date +%s)

while [[ $ATTEMPT -le $MAX_RETRIES ]]; do
    CURRENT_TIME=$(date +%s)
    ELAPSED=$((CURRENT_TIME - START_TIME))
    if [[ $ELAPSED -ge $MAX_TOTAL_SECONDS ]]; then
        echo "Exceeded maximum retry duration of 1 hour. Cleaning up and exiting."
        rm -f "$OUTDIR/${SAMPLE}.sra" "$OUTDIR/${SAMPLE}_1.fastq.gz" "$OUTDIR/${SAMPLE}_2.fastq.gz" "$OUTDIR/${SAMPLE}.fastq.gz"
        exit 1
    fi

    echo "Attempt $ATTEMPT (elapsed ${ELAPSED}s): Starting fastq-dl at $(date)"
    check_server || exit 1

    if timeout $TIMEOUT_SECONDS fastq-dl --cpus "$NCORES" -a "$SAMPLE" --outdir "$OUTDIR"; then
        echo "fastq-dl succeeded on attempt $ATTEMPT."
        SUCCESS=1
        break
    else
        echo "fastq-dl failed on attempt $ATTEMPT."
        SLEEP_BACKOFF=$(( (RANDOM % (RETRY_DELAY_MAX - RETRY_DELAY_MIN + 1)) + RETRY_DELAY_MIN ))
        echo "Sleeping for $SLEEP_BACKOFF seconds before retrying..."
        sleep $SLEEP_BACKOFF
        ATTEMPT=$((ATTEMPT + 1))
    fi
done

if [ -f "$OUTDIR/${SAMPLE}.sra" ]; then
    cd "$OUTDIR"
    fastq-dump --split-3 "${SAMPLE}.sra" --gzip
    rm -f "${SAMPLE}.sra"
fi

if [[ $SUCCESS -ne 1 ]]; then
    echo "fastq-dl failed after $MAX_RETRIES attempts. Cleaning up and exiting."
    rm -f "$OUTDIR/${SAMPLE}.sra" "$OUTDIR/${SAMPLE}_1.fastq.gz" "$OUTDIR/${SAMPLE}_2.fastq.gz" "$OUTDIR/${SAMPLE}.fastq.gz"
    exit 1
fi

# Post-download processing
R1="$OUTDIR/${SAMPLE}_1.fastq.gz"
R2="$OUTDIR/${SAMPLE}_2.fastq.gz"
SINGLE="$OUTDIR/${SAMPLE}.fastq.gz"

# Handle layout
if [ -f "$R2" ]; then
    rm -f "$R2"
fi

if [ -f "$R1" ] && [ -f "$SINGLE" ]; then
    rm -f "$SINGLE"
elif [ -f "$SINGLE" ]; then
    cp "$SINGLE" "$FORWARD_FASTQ"
elif [ -f "$R1" ]; then
    # No need to do anything, R1 is already named properly
    :
else
    echo "No suitable forward read file found after download. Cleaning up and failing."
    rm -f "$OUTDIR/${SAMPLE}_1.fastq.gz" "$OUTDIR/${SAMPLE}_2.fastq.gz" "$OUTDIR/${SAMPLE}.fastq.gz"
    exit 1
fi

# Cleanup
if [ -f "$SINGLE" ]; then
    rm -f "$SINGLE"
fi

# Sanity check: Validate final output
if [[ ! -s "$FORWARD_FASTQ" || $(stat -c%s "$FORWARD_FASTQ") -lt 10000 ]]; then
    echo "Forward FASTQ file is too small or empty: $FORWARD_FASTQ"
    rm -f "$FORWARD_FASTQ"
    exit 1
fi

# Success flag
touch "$OUTDIR/ena_done.txt"
echo "Finished successfully. Flag ena_done.txt created."