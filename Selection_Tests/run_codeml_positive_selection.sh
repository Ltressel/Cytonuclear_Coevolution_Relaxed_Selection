#!/bin/bash
# run_codeml_positive_selection.sh
# This script runs PAML (CODEML) to detect positive selection using the branch-site model

# Define input files
ALT_CTL_FILE="codeml_alt.ctl"  # Alternative model control file
NULL_CTL_FILE="codeml_null.ctl"  # Null model control file
OUTPUT_DIR="paml_output"
mkdir -p $OUTPUT_DIR

# Run the alternative model (allows dN/dS > 1 in 50-kb clade)
echo "Running alternative model..."
codeml $ALT_CTL_FILE
mv rst paml_output/alt_model_rst
mv 2ML.t paml_output/alt_model_2ML.t
mv 2ML.b paml_output/alt_model_2ML.b

# Run the null model (constrains dN/dS = 1 in the 50-kb clade)
echo "Running null model..."
codeml $NULL_CTL_FILE
mv rst paml_output/null_model_rst
mv 2ML.t paml_output/null_model_2ML.t
mv 2ML.b paml_output/null_model_2ML.b

# Output message
echo "Positive selection detection complete. Check the $OUTPUT_DIR for results."

# NOTE: Make sure that PAML is installed and that your control files (codeml_alt.ctl and codeml_null.ctl)
# are correctly set up with paths to your alignment and tree files.
