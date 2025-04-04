# Transcriptomic Response of a Murine Model Melanoma to Probiotic Bacteria

## Overview

This repository contains the complete analytical workflow and results for our investigation *Bifidobacterium adolescentis* 150 and *Lacticaseibacillus rhamnosus* K32 influence gene expression patterns in a murine melanoma  B16-F1 using bulk RNA-seq.

🔗 **Access the full study report:**  [https://jeniaole13.github.io/RNA-model-melanoma/](https://jeniaole13.github.io/RNA-model-melanoma/)

## Key Findings
Combined Effect of *B. adolescentis* 150 and *L. rhamnosus* K32: 
- Upregulated: epithelial-mesenchymal transition (EMT), NF-κB/TNF-α signaling, hypoxia, IL-6/STAT5 signaling.
- Downregulated: E2F targets, G2/M checkpoints.

Effects of *B. adolescentis*  150
- Upregulated: TGF-β signaling, WNT/β-catenin signaling.
- Downregulated: immune cell markers.

Effects of *L. rhamnosus* K32
Upregulated: IFN-γ and IFN-α signaling, IL-6/JAK/STAT3 signaling; immune cell markers.

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
README.md                   # Repository description file
RNA-seq-report.html         # Quarto HTML file
RNA-seq-report.qmd          # Quarto qmd file
```

## Data Availability
Raw sequencing data is available uned BioProject [PRJNA1214537](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1214537/).

## Contacts
For any questions regarding this work, please contact us at jeniaole01@gmail.com

## Funding
Financial support for this study was provided by the Russian Science Foundation under the grant #22-75-10029 ([https://rscf.ru/project/22-75-10029/](https://rscf.ru/project/22-75-10029/)).
