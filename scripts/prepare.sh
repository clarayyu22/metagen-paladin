#!/bin/bash
set -e

#requires: folder called {FASTA} with {FASTA}.fasta.txt inside
FASTA_PATH=$1

#run paladin prepare
echo "Preparing Paladin database at $FASTA_PATH"
paladin prepare -r1 -f $FASTA_PATH
echo "Paladin prepare for ${FASTA_PATH} completed."


