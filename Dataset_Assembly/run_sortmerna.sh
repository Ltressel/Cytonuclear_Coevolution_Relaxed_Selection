#!/bin/bash
# run_sortmerna.sh
# Remove rRNA sequences from trimmed RNA-seq reads

# Set the path to SortMeRNA's rRNA database folder
DB_DIR="/path/to/sortmerna/rRNA_databases"  # <<< UPDATE THIS PATH >>>

# Create output directory if it doesn't exist
mkdir -p non_rRNA_reads

# Run SortMeRNA on the trimmed reads
for file in trimmed_reads/*_R1_*.fastq.gz; do
    base=$(basename $file "_R1_*.fastq.gz")
    sortmerna --ref $DB_DIR/*.fasta --reads trimmed_reads/${base}_R1_*.fastq.gz \
              --reads trimmed_reads/${base}_R2_*.fastq.gz --other non_rRNA_reads/ \
              --aligned rRNA_output/ --fastx
done

# Output message when completed
echo "rRNA removal complete. Non-rRNA reads saved in non_rRNA_reads/"

# NOTE: Update the DB_DIR variable to point to your local rRNA databases.
# Download SortMeRNA from: https://github.com/biocore/sortmerna
