# Spatial Transcriptomics – Input Data

This directory contains the expected input files for the Seurat-based spatial transcriptomics workflow.

⚠️ Raw spatial datasets and large Seurat objects are not included in this repository due to file size limitations.

---

## Expected structure

data/
├── obj_final_subset.rds
└── interactions.csv


---

## Required files

### 1. obj_final_subset.rds

Processed Seurat spatial object.

- Format: `.rds` (Seurat object)
- Contains:
  - Gene expression data
  - Spatial coordinates
  - Cell metadata (e.g., `celltype_profile`, `orig.ident`)
  - Xenium assay layers (if applicable)

This object serves as the main input for the spatial analysis pipeline.

---

### 2. interactions.csv

Ligand–receptor interaction database.

- Format: `.csv`
- Must contain at least:
  - `ligand` column
  - `receptor` column

Used for ligand–receptor interaction scoring and spatial communication analysis.

---

## Data source

Spatial dataset derived from publicly available cancer spatial transcriptomics studies (e.g., Xenium or similar imaging-based platforms).

Potential sources:

- 10x Genomics Xenium datasets
- GEO (Gene Expression Omnibus)
- ArrayExpress

---

## Notes

- The pipeline assumes the spatial object is already created and optionally subset beforehand.
- The object must contain spatial coordinates accessible via `GetTissueCoordinates()`.
- The metadata column `celltype_profile` is required for cluster vs annotation comparison and ligand–receptor analysis.
- Large intermediate files are not stored in this repository.
