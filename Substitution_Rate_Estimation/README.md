# Substitution_Rate_Estimation

This directory contains scripts and methods for estimating the nonsynonymous (*dN*) and synonymous (*dS*) substitution rates for genes in the study. The analyses were conducted using **PAML v4.7** and **HyPhy**, with likelihood ratio tests (LRTs) to compare substitution rates between branches.

### Contents
- ```run_codeml.sh```: Script to run PAML (CODEML) for estimating *dN*, *dS*, *dN/dS*, and performing LRTs.
- ```run_hyphy.sh```: Script to run HyPhy with the MG94xHKY85 codon model.
- ```run_wilcoxon_test.R```: R script for performing Wilcoxon rank sum tests on the *dN* and *dS* values.
- ```generate_boxplots.R```: R script for generating boxplots of *dN* and *dS*, and combining them into a single figure.

### Step-by-Step Breakdown
**1. Run PAML (CODEML) for *dN* and *dS* Estimation**

**Script:** ```run_codeml.sh```

**Tool:** [PAML v4.7](https://github.com/abacus-gene/paml)

**Description:** This script runs PAML (CODEML) to estimate *dN* and *dS* for each gene and conduct likelihood ratio tests (LRTs) between models.

**Usage:** ```bash run_codeml.sh```

**2. Run HyPhy for MG94xHKY85 Codon Model Analysis**

**Script:** ```run_hyphy.sh```

**Tool:** [HyPhy](https://www.hyphy.org/)

**Description:** This script runs HyPhy using the MG94xHKY85 codon model to estimate dN and dS and conduct LRTs between clades.

**Usage:** ```bash run_hyphy.sh```

**3. Perform Wilcoxon Rank Sum Test**

**Script:** ```run_wilcoxon_test.R```

**Tool:** R (with necessary statistical packages)

**Description:** This R script performs Wilcoxon rank sum tests on the resulting dN and dS values and adjusts P values using false discovery rate correction (FDR).

**Usage:** ```Rscript run_wilcoxon_test.R```

**4. Generate Boxplots of *dN* and *dS***

**Script:** ```generate_boxplots.R```

**Tool:** R (with necessary packages: ```ggplot2```, ```dplyr```, ```tidyr```, ```patchwork```)

**Description:** This R script generates boxplots for *dN* and *dS* values for the 50-kb legumes and other legumes. It combines the boxplots into a single figure and saves the output as a PNG file.

**Usage:** ```Rscript generate_boxplots.R```
