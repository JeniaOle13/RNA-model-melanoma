# Transcriptomic Response of a Model Melanoma to Probiotic Bacteria

## Overview

This repository contains the complete analytical workflow and results for our investigation *Bifidobacterium adolescentis* 150 and *Lacticaseibacillus rhamnosus* K32 influence gene expression patterns in a murine melanoma  B16-F1 using bulk RNA-seq.

🔗 **Access the full study report:**  [https://jeniaole13.github.io/RNA-model-melanoma/](https://jeniaole13.github.io/RNA-model-melanoma/)

## Key Findings

-   Identification of differentially expressed genes (DEGs) following probiotic treatment.
-   Characterization of immune-related pathways modulated by probiotics.
-   Comparative analysis of transcriptional profiles between treatment groups.
-   Potential mechanisms of probiotic-mediated anti- or pro-tumor effects.

## Repository structure
```
RNA-seq-report_files/       # Files for Quarto report
data/
├── htseq/*.counts          # HTSeq counts data
├── metadata.tsv            # Metadata
├── ensbl2geneid.tsv        # Mapping ENSEMBL to SYMBOL file
├── cell_marker_mouse.csv   # Cell Markers 2.0 for *Mus musculus*
└── counts.tsv              # Read counts data
figures/                    # Figures for Quarto report
RNA-seq-report.html         # Quarto HTML file
RNA-seq-report.qmd          # Quarto qmd file
```

## Data Availability
Raw sequencing data is available uned BioProject [PRJNA1214537](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1214537/).

## Contacts
For any questions regarding this work, please contact us at jeniaole01@gmail.com

## Funding
Financial support for this study was provided by the Russian Science Foundation under the grant #22-75-10029 ([https://rscf.ru/project/22-75-10029/](https://rscf.ru/project/22-75-10029/)).
