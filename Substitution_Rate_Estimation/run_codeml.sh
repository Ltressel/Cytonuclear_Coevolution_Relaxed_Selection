#!/bin/bash
# run_codeml.sh
# Run PAML (CODEML) to estimate dN/dS and perform LRTs

# Define input files and control file
CTL_FILE="codeml.ctl"  # <<< UPDATE THIS WITH YOUR CODEML CONTROL FILE >>>
OUTPUT_DIR="paml_output"
mkdir -p $OUTPUT_DIR

# Run CODEML with the control file
codeml $CTL_FILE

# Output message when completed
echo "PAML analysis complete. Results saved in $OUTPUT_DIR"

# NOTE: Ensure PAML is installed and available in your system's PATH.
# You can download PAML from: http://abacus.gene.ucl.ac.uk/software/paml.html
# Customize the control file (codeml.ctl) as needed for your analysis.
