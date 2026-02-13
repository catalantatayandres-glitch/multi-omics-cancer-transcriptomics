# Spatial Transcriptomics Analysis

This module integrates transcriptomic clustering with spatial tissue organization to investigate cellular structure and cell–cell communication in cancer tissue.

The workflow extends single-cell analysis by incorporating spatial coordinates and ligand–receptor interaction modeling.

---

## Objective

Spatial transcriptomics preserves the physical location of cells within the tissue, enabling the study of how transcriptional identity relates to tissue architecture.

This workflow aims to:

- Visualize transcriptomic clusters in spatial context
- Compare clustering to biological annotations
- Identify cell-type marker genes
- Infer cell–cell communication via ligand–receptor interactions
- Incorporate spatial proximity into interaction scoring

---

## Analytical Workflow

### Step 1 — Transcriptomic Preprocessing

- Normalization
- Scaling
- Principal Component Analysis (PCA)
- UMAP embedding

### Step 2 — Clustering

- k-nearest neighbor graph construction
- Louvain community detection
- Cluster visualization in transcriptomic space

### Step 3 — Spatial Visualization

- Projection of clusters onto tissue coordinates
- Molecule-level visualization
- Comparison with annotated cell types

### Step 4 — Cluster Validation

- Rand Index
- Adjusted Rand Index (ARI)
- Agreement between clustering and annotation

### Step 5 — Marker Gene Detection

- Marker genes per cell type
- Biological interpretation of spatial domains

### Step 6 — Cell–Cell Communication Modeling

Ligand–receptor interaction scoring based on:

1. Expression strength (log fold-change product)
2. Spatial proximity weighting:
   - Centroid distance similarity
   - Neighborhood composition correction

---

## Methodological Considerations

- Spatial proximity influences inferred biological interactions.
- Clustering agreement does not guarantee biological equivalence.
- Ligand–receptor scoring represents potential communication, not direct signaling evidence.
- Spatial weighting reduces false-positive long-range interactions.

---

## Folder Structure

spatial-transcriptomics/
├── scripts/ # Spatial analysis pipeline
├── data/ # Expected input structure
└── results/ # Output description and representative figures


See:

- `data/README.md` for required input files
- `results/README.md` for generated outputs

---

## Reproducibility

From repository root:

```bash
Rscript spatial-transcriptomics/scripts/spatial-transcriptomics_pipeline.R
All outputs are generated programmatically.

Large spatial objects are not included in this repository due to file size constraints.

---

## Key Concepts Demonstrated

-Integration of transcriptomics and tissue architecture

-Spatially-aware clustering validation

-Marker gene localization

-Ligand–receptor interaction inference

-Spatial weighting of communication networks

---

## Context

This module represents the spatial layer of a multi-scale cancer transcriptomics framework, complementing:

Bulk RNA-seq (tissue-level expression)

Single-cell RNA-seq (cellular heterogeneity)
