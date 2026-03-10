#!/bin/bash
# run_hyphy.sh
# Run HyPhy to analyze dN/dS using the MG94xHKY85 codon model

# Define input files and model
ALIGNMENT_FILE="aligned_sequences.phy"  # <<< UPDATE THIS WITH YOUR ALIGNED SEQUENCES >>>
TREE_FILE="constraint_tree.tre"         # <<< UPDATE THIS WITH YOUR TREE FILE >>>
MODEL="MG94xHKY85"
OUTPUT_DIR="hyphy_output"
mkdir -p $OUTPUT_DIR

# Run HyPhy with the specified model, alignment, and tree
hyphy $MODEL --alignment $ALIGNMENT_FILE --tree $TREE_FILE --output $OUTPUT_DIR

# Output message when completed
echo "HyPhy analysis complete. Results saved in $OUTPUT_DIR"

# NOTE: Ensure HyPhy is installed and available in your system's PATH.
# You can download HyPhy from: https://www.hyphy.org/
