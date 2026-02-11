# Single-cell RNA-seq – Input Data

This directory contains the expected input files for the Seurat-based single-cell RNA-seq workflow.

⚠️ Raw data are not included in this repository due to file size limitations.

## Expected structure

data/
├── sc_counts_BRCA.txt (MatrixMarket format)
├── sc_barcodes_BCRA.tsv
├── sc_genes_BCRA.tsv
└── meta_data.txt


## Required files

### 1. sc_counts_BRCA.txt
- Sparse count matrix
- Format: MatrixMarket (.mtx or similar)
- Dimensions: cells × genes

### 2. sc_barcodes_BCRA.tsv
- Cell barcode identifiers
- One barcode per row

### 3. sc_genes_BCRA.tsv
- Gene identifiers
- One gene per row
- Must match matrix column order

### 4. meta_data.txt
- Cell-level metadata
- Must contain either:
  - `cell_id` column (recommended), or
  - Same row order as barcodes

Optional metadata columns:
- cell type annotations
- sample origin
- batch information

## Data source

Single-cell dataset derived from publicly available cancer transcriptomic studies.

Potential sources:
- GEO (Gene Expression Omnibus)
- ArrayExpress
- Single Cell Portal

## Notes

- The Seurat pipeline expects the matrix to be transposed to genes × cells internally.
- Ensure gene and barcode ordering matches the count matrix.
- Large processed objects (e.g., Seurat `.rds` files) are not stored in this repository.
