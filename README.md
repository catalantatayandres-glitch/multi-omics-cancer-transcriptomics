# Multi-Omics Cancer Transcriptomics

Reproducible transcriptomics pipelines applied to cancer datasets, integrating bulk RNA-seq, single-cell RNA-seq, and spatial transcriptomics analyses.

This repository demonstrates structured, end-to-end computational workflows for transcriptomic data processing, exploratory analysis, clustering, and intercellular communication modeling.

---

## Overview

Cancer is a multi-scale biological system. This project explores transcriptomic data at three complementary levels:

- **Bulk RNA-seq** — tissue-level gene expression profiling  
- **Single-cell RNA-seq** — cellular heterogeneity and clustering  
- **Spatial transcriptomics** — tissue architecture and cell–cell communication  

Each module is self-contained and fully script-based.

---

## Repository Structure

multi-omics-cancer-transcriptomics/
├── bulk-rnaseq/
├── scrnaseq-seurat/
└── spatial-transcriptomics/


Each module contains:

- `scripts/` – Reproducible R pipeline  
- `data/` – Expected input structure (raw data not included)  
- `results/` – Output description and representative figures  

---

## Methodological Scope

### Bulk RNA-seq
- Gene length normalization (RPK)
- Library size normalization (TPM)
- Technical bias assessment
- PCA-based dimensionality reduction
- Unsupervised clustering
- Subtype comparison

### Single-cell RNA-seq
- Quality control filtering
- Feature selection and scaling
- PCA and UMAP embedding
- Graph-based clustering (Louvain)
- Marker gene identification

### Spatial transcriptomics
- Transcriptomic clustering in spatial context
- Tissue-coordinate visualization
- Cluster vs annotation comparison (Rand Index / ARI)
- Ligand–receptor interaction scoring
- Spatially weighted communication modeling

---

## Reproducibility

All analyses are fully script-based and can be executed using:

```bash
Rscript <module>/scripts/<pipeline>.R
Large raw datasets are not included due to size constraints.

---

## Tools & Technologies

- **R**
- **Seurat**
- **ggplot2**
- **dplyr**
- **matrixStats**
- Spatial transcriptomics (Xenium-compatible Seurat objects)

Author
Andrés Catalán Tatay
MSc Molecular Life Sciences
Focus: Cancer transcriptomics, bioinformatics, multi-scale data integration
