#!/bin/bash
# run_fastqc.sh
# Run FastQC on all raw RNA-seq files

# Create output directory if it doesn't exist
mkdir -p fastqc_output

# Run FastQC on all .fastq.gz files
fastqc *.fastq.gz -o fastqc_output/

# Output message when completed
echo "FastQC analysis complete. Results saved in fastqc_output/"

# NOTE: Ensure FastQC is installed and available in your system's PATH.
# You can download FastQC from: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
