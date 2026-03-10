#!/bin/bash
# run_trinity.sh
# Assemble transcriptomes from non-rRNA reads using Trinity

# Create output directory if it doesn't exist
mkdir -p trinity_output

# Run Trinity with paired-end non-rRNA reads
Trinity --seqType fq --left non_rRNA_reads/*_R1_*.fastq.gz \
        --right non_rRNA_reads/*_R2_*.fastq.gz \
        --max_memory 50G --CPU 16 --output trinity_output/

# Output message when completed
echo "Transcriptome assembly complete. Results saved in trinity_output/"

# NOTE: Ensure Trinity is installed and available in your system's PATH.
# You can download Trinity from: https://github.com/trinityrnaseq/trinityrnaseq
# Adjust --max_memory and --CPU based on your system’s resources.
