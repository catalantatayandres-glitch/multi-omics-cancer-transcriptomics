#!/usr/bin/env Rscript
# -----------------------------------------------------------------------------
# scRNA-seq (Seurat): QC -> normalization -> PCA/UMAP -> clustering -> markers
# Author: Andrés Catalán Tatay
#
# Pipeline overview:
#   Step 1: Load sparse count matrix + metadata and build Seurat object
#   Step 2: QC metrics + basic filtering
#   Step 3: Standard Seurat workflow (Normalize -> HVGs -> Scale -> PCA -> UMAP)
#   Step 4: Clustering + marker genes
#
# Expected inputs (project root expected):
#   data/sc/sc_counts_BRCA.txt              (MatrixMarket; sparse counts)
#   data/sc/sc_barcodes_BCRA.tsv            (cell barcodes; 1 column)
#   data/sc/sc_genes_BCRA.tsv               (gene IDs; 1 column)
#   data/sc/meta_data.txt                   (cell metadata; optional cell_id column)
#
# Outputs:
#   results/sc_week3/
#     figures/                              (QC + UMAP + PCA plots)
#     cluster_markers.csv
#     seurat_processed.rds
#     sessionInfo.txt
#     qc_thresholds.txt
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
  library(ggplot2)
  library(dplyr)
})

# ---- Reproducibility --------------------------------------------------------
set.seed(1)

# ---- Paths ------------------------------------------------------------------
project_root <- getwd()
data_dir <- file.path(project_root, "data", "sc")

counts_mtx   <- file.path(data_dir, "sc_counts_BRCA.txt")
barcodes_tsv <- file.path(data_dir, "sc_barcodes_BCRA.tsv")
genes_tsv    <- file.path(data_dir, "sc_genes_BCRA.tsv")
meta_path    <- file.path(data_dir, "meta_data.txt")

out_dir <- file.path(project_root, "results", "sc_week3")
fig_dir <- file.path(out_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Early checks (fail fast) -----------------------------------------------
stopifnot(
  "counts_mtx not found"   = file.exists(counts_mtx),
  "barcodes_tsv not found" = file.exists(barcodes_tsv),
  "genes_tsv not found"    = file.exists(genes_tsv),
  "meta_path not found"    = file.exists(meta_path)
)

# ---- Helpers ----------------------------------------------------------------
save_plot <- function(p, filename, width = 6.5, height = 4.5) {
  ggsave(file.path(fig_dir, filename), p, width = width, height = height)
}

# ---- Load matrix ------------------------------------------------------------
message("[Load] Reading MatrixMarket counts...")
counts <- readMM(counts_mtx)

message("[Load] Reading barcodes + genes...")
cell_ids <- read.table(barcodes_tsv, header = FALSE, stringsAsFactors = FALSE)$V1
gene_ids <- read.table(genes_tsv, header = FALSE, stringsAsFactors = FALSE)$V1

# Your current layout: counts is (cells x genes)
stopifnot(
  "Mismatch: nrow(counts) vs barcodes" = (length(cell_ids) == nrow(counts)),
  "Mismatch: ncol(counts) vs genes"    = (length(gene_ids) == ncol(counts))
)
dimnames(counts) <- list(cell_ids, gene_ids)

message("[Load] Reading metadata...")
meta <- read.table(meta_path, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

# Metadata alignment (preferred: meta$cell_id)
if ("cell_id" %in% colnames(meta)) {
  idx <- match(cell_ids, meta$cell_id)
  if (anyNA(idx)) {
    warning("Some cell_ids are missing in meta$cell_id (NAs after match). Keeping those rows as NA.")
  }
  meta <- meta[idx, , drop = FALSE]
} else {
  # Otherwise assume same order if number of rows matches
  stopifnot("meta rows do not match number of cells" = (nrow(meta) == length(cell_ids)))
  rownames(meta) <- cell_ids
}

# ---- Build Seurat object ----------------------------------------------------
message("[Step 1] Creating Seurat object...")
# Seurat expects genes x cells
obj <- CreateSeuratObject(
  counts = t(counts),
  meta.data = meta,
  min.cells = 3,
  min.features = 200
)

# ---- QC metrics -------------------------------------------------------------
message("[Step 2] QC metrics + plots...")
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")

VlnPlot(
  obj,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3,
  pt.size = 0
)
save_plot(p_vln, "qc_violin.png", width = 10, height = 3.8)

# Two QC scatters (saved together as a patchwork-like combined plot)
p_scatter1 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
p_scatter2 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
p_scatter <- p_scatter1 + p_scatter2
save_plot(p_scatter, "qc_scatter.png", width = 10, height = 4.2)

# Filter thresholds (record them for reproducibility)
qc_min_features <- 300
qc_max_features <- 6000
qc_max_mt <- 20

writeLines(
  c(
    paste0("qc_min_features=", qc_min_features),
    paste0("qc_max_features=", qc_max_features),
    paste0("qc_max_mt=", qc_max_mt)
  ),
  con = file.path(out_dir, "qc_thresholds.txt")
)

obj <- subset(
  obj,
  subset = nFeature_RNA > qc_min_features &
    nFeature_RNA < qc_max_features &
    percent.mt < qc_max_mt
)

# ---- Standard Seurat workflow ----------------------------------------------
message("[Step 3] Normalize -> HVGs -> Scale -> PCA...")
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
obj <- ScaleData(obj, features = VariableFeatures(obj))
obj <- RunPCA(obj, features = VariableFeatures(obj))

p_elbow <- ElbowPlot(obj, ndims = 50)
save_plot(p_elbow, "pca_elbow.png", width = 6.5, height = 4)

dims_use <- 1:20
message("[Step 3] UMAP + neighbors + clustering (dims = ", paste(dims_use, collapse = ","), ")...")
obj <- RunUMAP(obj, dims = dims_use)
obj <- FindNeighbors(obj, dims = dims_use)
obj <- FindClusters(obj, resolution = 0.3)

p_umap_cluster <- DimPlot(obj, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
  ggtitle("UMAP: Seurat clusters")
save_plot(p_umap_cluster, "umap_clusters.png", width = 6.5, height = 5)

if ("cell_type" %in% colnames(obj@meta.data)) {
  p_umap_celltype <- DimPlot(obj, reduction = "umap", group.by = "cell_type", label = FALSE) +
    ggtitle("UMAP: provided cell-type labels")
  save_plot(p_umap_celltype, "umap_cell_type.png", width = 7.5, height = 5)
}

# ---- Marker genes -----------------------------------------------------------
message("[Step 4] Finding marker genes (FindAllMarkers)...")
markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers, file.path(out_dir, "cluster_markers.csv"), row.names = FALSE)

# ---- Save outputs -----------------------------------------------------------
message("[Save] Writing Seurat object + session info...")
saveRDS(obj, file.path(out_dir, "seurat_processed.rds"))

sink(file.path(out_dir, "sessionInfo.txt"))
print(sessionInfo())
sink()

message("Done. Results written to: ", out_dir)
