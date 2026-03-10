#!/bin/bash
# run_relax_hyphy.sh
# This script runs HyPhy (RELAX) to detect relaxed selection in the 50-kb clade

# Define input files
ALIGNMENT_FILE="aligned_sequences.fas"  # <<< Update this with your aligned sequences file
TREE_FILE="constraint_tree.tre"         # <<< Update this with your constraint tree file
OUTPUT_DIR="hyphy_relax_output"
mkdir -p $OUTPUT_DIR

# Run RELAX in HyPhy
echo "Running HyPhy RELAX..."
hyphy relax --alignment $ALIGNMENT_FILE --tree $TREE_FILE --output $OUTPUT_DIR/relax_results.json

# Output message
echo "Relaxed selection detection complete. Check the $OUTPUT_DIR for results."

# NOTE: Make sure that HyPhy is installed and that your alignment and tree files are correctly specified.
