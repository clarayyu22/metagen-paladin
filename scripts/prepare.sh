#!/bin/bash
set -e

#requires: folder called {FASTA} with {FASTA}.fasta.txt inside
FASTA_PATH=$1
PALADIN_REF_DIR=$2
FASTA_NAME=$(basename $FASTA_PATH)

#run paladin prepare
echo "Copying Fasta into new reference directory $PALADIN_REF_DIR"

if [ ! -f "$PALADIN_REF_DIR/$FASTA_NAME.bwt" ]; then
    echo "Preparing Paladin database at $PALADIN_REF_DIR"
    if [ ! -d "$PALADIN_REF_DIR" ]; then
        echo "Error: Reference directory does not exist!"
        exit 1
    fi
    cp $FASTA_PATH $PALADIN_REF_DIR
    paladin index -r3 "$PALADIN_REF_DIR/$FASTA_NAME"
    #paladin index -r3 $PALADIN_REF_DIR$FASTA_NAME
    # paladin prepare -r1 -f $PALADIN_REF_DIR/$FASTA_NAME
fi

echo "Paladin prepare for ${FASTA_PATH} completed."
