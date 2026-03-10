# Evolutionary_Rate_Covariation_Analyses

This directory contains scripts and methods for performing evolutionary rate covariation (ERC) analyses. 

### Contents
- ```run_blastp.sh```: Script for running a BLASTp search to identify plastid-targeted proteins.
- ```run_erc_analysis.R```: R script for performing the evolutionary rate covariation (ERC) analysis using Spearman rank correlations and bootstrap resampling.
- ```visualize_erc_results.R```: R script for visualizing the bootstrapped correlations between CpRP and nuclear-encoded gene groups.

### Step-by-Step Breakdown
**1. Run BLASTp to Identify Plastid-Targeted Proteins**

**Script:** ```run_blastp.sh```

**Tool:** [BLAST+](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html#downloadblastdata)

**Description:** This script runs BLASTp against a non-redundant protein database to categorize proteins as plastid-targeted or non-targeted, based on the annotation of the obtained hits.

**Usage:** ```bash run_blastp.sh```

**2. Perform ERC Analysis Using Spearman Rank Correlation and Bootstrap Resampling**

**Script:** ```run_erc_analysis.R```

**Tool:** R

**Description:** This script performs evolutionary rate covariation analysis by calculating Spearman rank correlations between normalized branch lengths of mt and each nuclear-encoded gene set (Glycolysis, Cytonuclear Ribosomal, and Cell Cycle). Correlation estimates are bootstrapped (10,000 iterations) to test for statistical differences.

**Usage:** ```Rscript run_erc_analysis.R```

**3. Visualize ERC Analysis Results**

**Script:** ```visualize_erc_results.R```

**Tool:** R

**Description:** This script visualizes the bootstrapped Spearman correlations between mt and nuclear-encoded gene sets (Glycolysis, Cytonuclear Ribosomal, and Cell Cycle). It generates density plots with the mean correlation value and the 95% confidence interval for each gene pair.

**Usage:** ```Rscript visualize_erc_results.R```
