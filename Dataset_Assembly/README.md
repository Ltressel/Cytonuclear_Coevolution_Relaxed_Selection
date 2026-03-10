# Dataset_Assembly
This directory contains scripts and methods for assembling the datasets used in the manuscript. Below is a breakdown of each step in the assembly pipeline including raw read quality control, adapter trimming, rRNA removal, transcriptome assembly, and sequence extraction.

### Contents
- ```run_fastqc.sh```: Quality control of raw RNA-seq reads using FastQC.
- ```run_trimgalore.sh```: Adapter trimming and low-quality read filtering using TrimGalore.
- ```run_sortmerna.sh```: Ribosomal RNA removal using SortMeRNA.
- ```run_trinity.sh```: Transcriptome assembly using Trinity.
- ```extract_align_genes.sh```: Gene extraction and alignment using BLAST and MAFFT.

### Step-by-Step Breakdown
**1. Quality Control**

**Script:** ```run_fastqc.sh```

**Tool:** [FastQC v.0.11.5](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

**Description:** This step checks the quality of raw 150 bp paired-end reads.

**Usage:** ```bash run_fastqc.sh```

**2. Adapter Trimming and Filtering**

**Script:** ```run_trimgalore.sh```

**Tool:** [TrimGalore v.0.5.0](https://github.com/FelixKrueger/TrimGalore)

**Description:** Trims adapters and filters low-quality reads from RNA-seq data.

**Usage:** ```bash run_trimgalore.sh```

**3. Ribosomal RNA Removal**

**Script:** ```run_sortmerna.sh```

**Tool:** [SortMeRNA v.2.1](https://github.com/sortmerna/sortmerna)

**Description:** Removes ribosomal RNA (rRNA) sequences using SortMeRNA with 8 rRNA databases.

**Usage:** ```bash run_sortmerna.sh```

**4. Transcriptome Assembly**

**Script:** ```run_trinity.sh```

**Tool:** [Trinity v.2.8.4](https://github.com/trinityrnaseq/trinityrnaseq)

**Description:** Assembles RNA-seq reads de novo into transcriptomes using Trinity.

**Usage:** ```bash run_trinity.sh```

**5. Gene Extraction and Alignment**

**Script:** ```extract_align_genes.sh```

**Tools:** [BLAST v.2.12.0+](https://blast.ncbi.nlm.nih.gov/Blast.cgi), [MAFFT v.7.490](https://mafft.cbrc.jp/alignment/software/source.html)

**Description:** Extracts genes from transcriptomes and mitogenomes, and aligns sequences.

**Usage:** ```bash extract_align_genes.sh```

**Important Notes for Users**
These scripts are designed to run the dataset assembly pipeline described in the manuscript, but **you will need to customize them with your own specific information** before use. Specifically:
- **File paths:** Update the paths in each script (e.g., input data, output directories, reference databases) to match the locations on your local system.
- **Database locations:** Ensure that reference databases (e.g., SortMeRNA rRNA databases, BLAST target gene databases) are correctly set in the scripts.
- **Software installation:** Make sure that the necessary tools (e.g., FastQC, TrimGalore, SortMeRNA, Trinity, BLAST, MAFFT) are installed and properly added to your system's ```PATH```.
- **Comments:** Sections that require updates are marked with comments like  ```# <<< UPDATE THIS PATH >>>```
