# Bulk RNA-seq – Input Data

This directory contains the expected input files for the bulk RNA-seq pipeline.

⚠️ Raw data are not included in this repository due to size and licensing constraints.

## Expected structure

data/
├── countmatrix/
│ └── expression_Matrix_BRCA.txt
├── genes_info.txt
└── meta_data/
└── TCGA_combined_meta_data_BRCA.txt


## Required files

### 1. expression_Matrix_BRCA.txt
- Gene expression count matrix
- Format: tab-separated
- Dimensions: genes × samples
- Values: raw counts

### 2. genes_info.txt
- Gene annotation table
- Must contain a gene length column:
  - `Length`, `gene_length`, or `length`
- Used for RPK/TPM normalization

### 3. TCGA_combined_meta_data_BRCA.txt
- Sample metadata
- Tab-separated
- Should contain:
  - Sample identifiers
  - Subtype information (e.g., PAM50)

## Data source

The dataset is derived from publicly available cancer transcriptomics resources (e.g., TCGA BRCA).

Data can be downloaded from:
- GDC Data Portal (https://portal.gdc.cancer.gov/)
- cBioPortal
- Other publicly accessible repositories

## Notes

- The folder structure must be preserved for the pipeline to run.
- The script assumes genes and gene metadata are aligned by order.
- No raw data are distributed in this repository.

