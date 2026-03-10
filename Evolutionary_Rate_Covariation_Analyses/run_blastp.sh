#!/bin/bash
# run_blastp.sh
# This script runs a BLASTp search to identify plastid-targeted proteins

# Define input and output files
QUERY_FILE="query_sequences.fasta"  # <<< Update with your input FASTA file
DB="nr"  # Non-redundant protein database (default for BLASTp)
OUTPUT_FILE="blastp_output.txt"

# Number of threads for parallel processing (adjust as needed)
THREADS=4

# Run BLASTp search
echo "Running BLASTp search against non-redundant protein database..."
blastp -query $QUERY_FILE -db $DB -out $OUTPUT_FILE -outfmt "6 std stitle" -evalue 1e-5 -num_threads $THREADS

# Output message upon completion
echo "BLASTp search complete. Results saved in $OUTPUT_FILE."

# NOTE: Make sure BLAST+ is installed and your input query file is correctly specified.
