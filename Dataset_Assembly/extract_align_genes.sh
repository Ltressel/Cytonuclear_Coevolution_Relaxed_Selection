#!/bin/bash
# extract_align_genes.sh
# Extract mt and nuclear genes using BLAST and align sequences with MAFFT

# Define paths for input files and output directories
TRANSCRIPTOME="trinity_output/Trinity.fasta"  # <<< UPDATE WITH YOUR TRANSCRIPTOME FILE >>>
DB="target_gene_db.fasta"                     # <<< UPDATE WITH YOUR BLAST GENE DATABASE >>>
OUT_DIR="blast_output"
ALIGN_DIR="aligned_sequences"

# Create output directories if they don't exist
mkdir -p $OUT_DIR $ALIGN_DIR

# Run BLASTx to identify ribosomal protein genes
blastx -query $TRANSCRIPTOME -db $DB -out $OUT_DIR/blast_output.txt -outfmt 6

# Extract FASTA sequences of BLAST hits
awk '{print $1}' $OUT_DIR/blast_output.txt | xargs samtools faidx $TRANSCRIPTOME > $OUT_DIR/blast_hits.fasta

# Run MAFFT to align the extracted sequences
mafft --auto $OUT_DIR/blast_hits.fasta > $ALIGN_DIR/aligned_sequences.fasta

# Output message when completed
echo "Gene extraction and alignment complete. Results saved in $ALIGN_DIR/"

# NOTE: Update the TRANSCRIPTOME and DB paths with your actual data.
# Ensure BLAST and MAFFT are installed and available in your system's PATH.
# You can download BLAST from: https://blast.ncbi.nlm.nih.gov/ and MAFFT from: https://mafft.cbrc.jp/alignment/software/
