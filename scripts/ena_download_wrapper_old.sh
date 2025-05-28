#!/bin/bash
set -euo pipefail

SAMPLE=$1
OUTDIR=$2
NCORES=$3
FORWARD_FASTQ=$4

# Random sleep between 40-60 seconds
SLEEP_TIME=$(( (RANDOM % 20) + 40 ))
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

ATTEMPT=1
SUCCESS=0

#mkdir -p "$OUTDIR"

while [[ $ATTEMPT -le $MAX_RETRIES ]]; do
    echo "Attempt $ATTEMPT: Starting fastq-dl at $(date)"
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

if [[ $SUCCESS -ne 1 ]]; then
    echo "fastq-dl failed after $MAX_RETRIES attempts. Touching dummy output."
    touch "$FORWARD_FASTQ"
    exit 0
fi

echo "Finished fastq-dl successfully at $(date)"

# Post-download processing
R1="$OUTDIR/${SAMPLE}_1.fastq.gz"
R2="$OUTDIR/${SAMPLE}_2.fastq.gz"
SINGLE="$OUTDIR/${SAMPLE}.fastq.gz"

if [ -f "$OUTDIR/${SAMPLE}.sra" ]; then
    cd "$OUTDIR"
    fastq-dump --split-3 *.sra --gzip
    rm "${SAMPLE}.sra"
fi

# Handle layout
if [ -f "$R2" ]; then
    rm -f "$R2"
fi

if [ -f "$R1" ] && [ -f "$SINGLE" ]; then
    rm -f "$SINGLE"
elif [ -f "$SINGLE" ]; then
    cp "$SINGLE" "$FORWARD_FASTQ"
elif [ -f "$R1" ]; then
    # Already good
    :
else
    echo "No suitable forward read file found after download" >&2
    touch "$FORWARD_FASTQ"  # fail gracefully
fi

# Cleanup
if [ -f "$SINGLE" ]; then
    rm -f "$SINGLE"
fi
