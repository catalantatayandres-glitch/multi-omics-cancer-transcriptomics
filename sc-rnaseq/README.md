# Single-cell RNA-seq Analysis (Seurat Workflow)

This module implements a reproducible single-cell RNA-seq workflow using Seurat, applied to cancer transcriptomic data.

The analysis focuses on cellular heterogeneity, dimensionality reduction, graph-based clustering, and marker gene identification.

---

## Objective

Single-cell RNA-seq enables transcriptomic profiling at cellular resolution, allowing the identification of distinct cell populations within heterogeneous tumor environments.

This workflow aims to:

- Perform rigorous quality control filtering
- Identify transcriptionally variable genes
- Reveal cellular structure via dimensionality reduction
- Detect cluster-specific marker genes
- Support biological annotation of cell populations

---

## Analytical Workflow

### Step 1 — Quality Control

- Number of detected genes (nFeature_RNA)
- Total UMI counts (nCount_RNA)
- Mitochondrial transcript percentage (percent.mt)
- Filtering of low-quality or stressed cells

### Step 2 — Normalization & Feature Selection

- Log-normalization
- Identification of highly variable genes
- Feature scaling

### Step 3 — Dimensionality Reduction

- Principal Component Analysis (PCA)
- Elbow plot inspection
- UMAP embedding

### Step 4 — Graph-based Clustering

- k-nearest neighbor graph construction
- Louvain community detection
- Cluster visualization in UMAP space

### Step 5 — Marker Gene Identification

- Wilcoxon rank-sum test
- Minimum expression threshold (min.pct)
- Log2 fold-change filtering

---

## Methodological Considerations

- Clustering resolution controls population granularity.
- UMAP is based on selected principal components.
- Marker genes are computed after unsupervised clustering.
- Annotation requires biological interpretation beyond statistical thresholds.

---

## Folder Structure

scrnaseq-seurat/
├── scripts/ # Seurat analysis pipeline
├── data/ # Expected input structure
└── results/ # Output description and representative figures


See:

- `data/README.md` for required input files
- `results/README.md` for generated outputs

---

## Reproducibility

From repository root:

```bash
Rscript scrnaseq-seurat/scripts/seurat_pipeline.R
All results are generated programmatically.

Large processed Seurat objects are not included in the repository due to file size constraints.

---

## Key Concepts Demonstrated

-Cellular quality control metrics

-Variance-driven feature selection

-Non-linear dimensionality reduction (UMAP)

-Graph-based community detection

-Marker gene discovery for biological annotation

---

## Context

This module represents the cellular-resolution layer of a multi-scale cancer transcriptomics framework, complemented by:

Bulk RNA-seq analysis

Spatial transcriptomics integration
