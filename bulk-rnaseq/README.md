# Bulk RNA-seq Analysis

This module implements a reproducible exploratory workflow for bulk RNA-seq data derived from cancer transcriptomic datasets.

The analysis focuses on normalization strategies, technical bias assessment, dimensionality reduction, and unsupervised sample clustering.

---

## Objective

Bulk RNA-seq provides tissue-level gene expression measurements across heterogeneous cell populations.  
This workflow aims to:

- Demonstrate proper normalization of count data
- Evaluate technical biases inherent to RNA-seq
- Identify latent sample structure
- Compare unsupervised clusters to known molecular subtypes

---

## Analytical Workflow

### Step 1 — Gene Length & Library Size Normalization

- RPK (Reads Per Kilobase)
- TPM (Transcripts Per Million)
- Validation of normalization effects

### Step 2 — Technical Bias Assessment

- Gene length bias
- Sequencing depth bias
- Mean–variance relationship

These diagnostics illustrate why normalization is necessary for downstream analysis.

### Step 3 — Dimensionality Reduction & Clustering

- PCA on top variable genes
- k-means clustering
- Hierarchical clustering (Ward's method)
- Contingency analysis vs known subtypes

---

## Methodological Considerations

- Normalization is performed explicitly (RPK → TPM) for educational clarity.
- PCA is computed on variance-ranked genes to reduce noise.
- Clustering is exploratory and not intended as a formal subtype classifier.
- Differential expression should be performed on raw counts using appropriate statistical frameworks (e.g., DESeq2, edgeR, limma-voom).

---

## Folder Structure

bulk-rnaseq/
├── scripts/ # Analysis pipeline
├── data/ # Expected input structure
└── results/ # Output description and representative figures


See:

- `data/README.md` for required input structure
- `results/README.md` for generated outputs

---

## Reproducibility

From repository root:

```bash
Rscript bulk-rnaseq/scripts/bulk_pipeline.R
All outputs are generated programmatically.

Large datasets are not included in this repository due to size constraints.

---

## Key Concepts Demonstrated

- Count normalization strategies

- Variance structure of RNA-seq data

- Dimensionality reduction

- Unsupervised clustering

- Integration with sample metadata

---

## Context

This module represents the bulk transcriptomic layer of a multi-scale cancer analysis framework, complemented by:

Single-cell RNA-seq analysis

Spatial transcriptomics integration

