# Selection_Tests
This directory contains scripts and methods for detecting positive and relaxed selection within the 50-kb clade. The analyses were conducted using **PAML (CODEML)** for positive selection tests and **HyPhy (RELAX)** for relaxed selection tests. Additionally, Python scripts are included for visualization of LRT p-values and *dN/dS* values.

### Contents
- ```run_codeml_positive_selection.sh```: Script for running PAML (CODEML) to detect positive selection using branch-site models.
- ```run_relax_hyphy.sh```: Script for running HyPhy (RELAX) to detect relaxed selection.
- ```p_values_for_LRTs.py```: Python script to visualize the p-values from the LRT results.


### Step-by-Step Breakdown
**1. Run PAML for Positive Selection Detection**

**Script:** ```run_codeml_positive_selection.sh```

**Tool:** [PAML v4.7](https://github.com/abacus-gene/paml)

**Description:** This script runs PAML (CODEML) to test for positive selection using the branch-site model. Two models are compared: the alternative model (which allows for positive selection in the 50-kb clade) and the null model (which constrains dN/dS to 1 within the clade).

**Usage:** ```bash run_codeml_positive_selection.sh```

**2. Run HyPhy for Relaxed Selection Detection**

**Script:** ```run_relax_hyphy.sh```

**Tool:** [HyPhy](https://www.hyphy.org/)

**Description:** This script runs HyPhy (RELAX) to test for relaxed selection in the 50-kb clade. The null model assumes no change in selection intensity (k = 1), while the alternative model allows for relaxed selection (k < 1).

**Usage:** ```bash run_relax_hyphy.sh```

**3. Visualize LRT p-values**

**Script:** ```p_values_for_LRTs.py```

**Tool:** Python (with required libraries such as ```matplotlib```, ```seaborn```)

**Description:** This Python script visualizes the p-values from the LRT results for the positive and relaxed selection tests.

**Usage:** ```python p_values_for_LRTs.py```
