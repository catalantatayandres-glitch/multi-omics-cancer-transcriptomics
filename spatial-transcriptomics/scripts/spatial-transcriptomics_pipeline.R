#!/usr/bin/env Rscript
# -----------------------------------------------------------------------------
# Spatial transcriptomics (Seurat): exploration -> clustering -> markers -> LR interactions
# Author: Andrés Catalán Tatay
#
# Pipeline overview:
#   Step 1: Preprocess (Normalize/Scale/PCA/UMAP) + clustering
#   Step 2: Visualize in transcriptomic space and tissue space
#   Step 3: QC distributions (reads/features) + spatial QC plots
#   Step 4: Compare clustering vs annotation (Rand Index / ARI)
#   Step 5: Ligand–receptor interactions + spatial weighting
#
# Inputs (project root expected):
#   data/spatial/obj_final_subset.rds
#   data/spatial/interactions.csv
#
# Outputs:
#   results/spatial/figures/         (plots)
#   results/spatial/                (tables + RDS objects + sessionInfo)
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(tibble)
})

options(future.globals.maxSize = 100000 * 1024^2)

# ---- Parameters -------------------------------------------------------------
dims_use        <- 1:10
resolution_use  <- 0.5
neighbor_radius <- 25

# ---- Paths ------------------------------------------------------------------
project_root <- getwd()

data_dir <- file.path(project_root, "data", "spatial")
obj_path <- file.path(data_dir, "obj_final_subset.rds")
lr_path  <- file.path(data_dir, "interactions.csv")

out_dir <- file.path(project_root, "results", "spatial")
fig_dir <- file.path(out_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

stopifnot(
  "Spatial object not found" = file.exists(obj_path),
  "LR database not found"    = file.exists(lr_path)
)

save_plot <- function(p, filename, width = 7.5, height = 5) {
  ggsave(filename = file.path(fig_dir, filename), plot = p, width = width, height = height)
}

# Small helper: safe "polychrome" palette (won't crash if unavailable)
safe_cols <- function(n) {
  out <- tryCatch(
    Seurat::DiscretePalette(n, palette = "polychrome"),
    error = function(e) NULL
  )
  if (is.null(out)) out <- Seurat::DiscretePalette(n)  # fallback
  out
}

# ---- Load spatial object ----------------------------------------------------
message("[Load] Reading spatial Seurat object...")
spatial_obj <- readRDS(obj_path)

message("[Load] Object dimensions (features x cells): ", paste(dim(spatial_obj), collapse = " x "))

# =============================================================================
# Step 1: Preprocessing (as scRNA-seq) + clustering
# =============================================================================
message("[Step 1] Normalize -> Scale -> PCA -> UMAP...")

spatial_obj <- NormalizeData(spatial_obj)

# (kept as in your script) Scale all features; for large datasets this can be heavy.
spatial_obj <- ScaleData(spatial_obj)

spatial_obj <- RunPCA(spatial_obj, features = rownames(spatial_obj))
spatial_obj <- RunUMAP(spatial_obj, dims = dims_use)

# UMAP plot (transcriptomic space)
p_umap <- DimPlot(spatial_obj, cols = safe_cols(length(levels(Idents(spatial_obj)))))
save_plot(p_umap, "umap_transcriptomic.png", width = 7.5, height = 5)

# (kept from your script; memory cleanup)
if ("Xenium" %in% names(spatial_obj@assays)) {
  spatial_obj@assays$Xenium@layers$scale.data <- NULL
}
gc()

message("[Step 1] Graph-based clustering (Neighbors -> Clusters)...")
spatial_obj <- FindNeighbors(spatial_obj, dims = dims_use)
spatial_obj <- FindClusters(spatial_obj, resolution = resolution_use)

p_umap_clust <- DimPlot(spatial_obj, label = TRUE, cols = safe_cols(length(levels(Idents(spatial_obj))))) +
  ggtitle("UMAP: transcriptomic clustering")
save_plot(p_umap_clust, "umap_clusters.png", width = 7.5, height = 5)

# Compare clustering vs provided annotation (if exists)
if ("celltype_profile" %in% colnames(spatial_obj@meta.data)) {
  p_umap_annot <- DimPlot(spatial_obj, group.by = "celltype_profile",
                          cols = safe_cols(length(unique(spatial_obj$celltype_profile)))) +
    ggtitle("UMAP: provided cell type annotation")
  save_plot(p_umap_annot, "umap_celltype_profile.png", width = 8.5, height = 5)
} else {
  message("[Step 1] Note: 'celltype_profile' not found in metadata.")
}

# =============================================================================
# Step 2: Spatial visualization (tissue space) + molecules
# =============================================================================
message("[Step 2] Tissue-space plots...")

p_img_clust <- ImageDimPlot(spatial_obj, cols = safe_cols(length(levels(Idents(spatial_obj))))) +
  ggtitle("Tissue: transcriptomic clusters")
save_plot(p_img_clust, "tissue_clusters.png", width = 8, height = 6)

if ("celltype_profile" %in% colnames(spatial_obj@meta.data)) {
  p_img_annot <- ImageDimPlot(spatial_obj, group.by = "celltype_profile",
                              cols = safe_cols(length(unique(spatial_obj$celltype_profile)))) +
    ggtitle("Tissue: provided cell type annotation")
  save_plot(p_img_annot, "tissue_celltype_profile.png", width = 8, height = 6)
}

# Molecules plot (kept from your script)
message("[Step 2] Molecules plot (C1QA, CD4, CD8A) ...")
marker_plot <- ImageDimPlot(
  spatial_obj,
  group.by = "celltype_profile",
  cols = safe_cols(length(unique(spatial_obj$celltype_profile))),
  molecules = c("C1QA", "CD4", "CD8A"),
  fov = "zoom"
)

# Save as PDF for zooming
pdf(file.path(fig_dir, "marker_molecules.pdf"), width = 10, height = 7)
print(marker_plot)
dev.off()

# =============================================================================
# Step 3: QC distributions (reads/features) + spatial QC
# =============================================================================
message("[Step 3] QC distributions...")

# These features are technology-specific (Xenium in your object)
if (all(c("nCount_Xenium", "nFeature_Xenium") %in% colnames(spatial_obj@meta.data))) {
  
  p_vln_count_orig <- VlnPlot(spatial_obj, features = "nCount_Xenium",
                              group.by = "orig.ident", pt.size = 0) + NoLegend()
  save_plot(p_vln_count_orig, "vln_nCount_by_origident.png", width = 10, height = 4)
  
  p_vln_feat_orig <- VlnPlot(spatial_obj, features = "nFeature_Xenium",
                             group.by = "orig.ident", pt.size = 0) + NoLegend()
  save_plot(p_vln_feat_orig, "vln_nFeature_by_origident.png", width = 10, height = 4)
  
  if ("celltype_profile" %in% colnames(spatial_obj@meta.data)) {
    p_vln_count_ct <- VlnPlot(spatial_obj, features = "nCount_Xenium",
                              group.by = "celltype_profile", pt.size = 0) + NoLegend()
    save_plot(p_vln_count_ct, "vln_nCount_by_celltype.png", width = 12, height = 5)
    
    p_vln_feat_ct <- VlnPlot(spatial_obj, features = "nFeature_Xenium",
                             group.by = "celltype_profile", pt.size = 0) + NoLegend()
    save_plot(p_vln_feat_ct, "vln_nFeature_by_celltype.png", width = 12, height = 5)
  }
  
  p_vln_count_cl <- VlnPlot(spatial_obj, features = "nCount_Xenium",
                            group.by = "seurat_clusters", pt.size = 0) + NoLegend()
  save_plot(p_vln_count_cl, "vln_nCount_by_cluster.png", width = 10, height = 4)
  
  p_vln_feat_cl <- VlnPlot(spatial_obj, features = "nFeature_Xenium",
                           group.by = "seurat_clusters", pt.size = 0) + NoLegend()
  save_plot(p_vln_feat_cl, "vln_nFeature_by_cluster.png", width = 10, height = 4)
  
  # In tissue space
  p_img_count <- ImageFeaturePlot(spatial_obj, features = "nCount_Xenium")
  save_plot(p_img_count, "tissue_nCount_Xenium.png", width = 8, height = 6)
  
  p_img_feat <- ImageFeaturePlot(spatial_obj, features = "nFeature_Xenium")
  save_plot(p_img_feat, "tissue_nFeature_Xenium.png", width = 8, height = 6)
  
} else {
  message("[Step 3] Note: nCount_Xenium / nFeature_Xenium not found in meta.data.")
}

# =============================================================================
# Step 4: Rand Index / Adjusted Rand Index
# =============================================================================
message("[Step 4] Rand Index / ARI...")

rand_index <- function(x, y) {
  if (length(x) != length(y)) stop("x and y must have the same length")
  
  n <- length(x)
  total_pairs <- choose(n, 2)
  
  tab <- table(x, y)
  a <- sum(choose(tab, 2))
  c <- sum(choose(rowSums(tab), 2))
  d <- sum(choose(colSums(tab), 2))
  
  b <- total_pairs - (c + d - a)
  rand <- (a + b) / total_pairs
  
  expected_index <- c * d / total_pairs
  max_index <- (c + d) / 2
  ari <- (a - expected_index) / (max_index - expected_index)
  
  list(RandIndex = rand, AdjustedRandIndex = ari)
}

if ("celltype_profile" %in% colnames(spatial_obj@meta.data)) {
  ri <- rand_index(x = spatial_obj$seurat_clusters, y = spatial_obj$celltype_profile)
  write.csv(
    data.frame(RandIndex = ri$RandIndex, AdjustedRandIndex = ri$AdjustedRandIndex),
    file.path(out_dir, "rand_index_clust_vs_celltype.csv"),
    row.names = FALSE
  )
} else {
  message("[Step 4] Skipped: 'celltype_profile' missing.")
}

# =============================================================================
# Step 5: Ligand–receptor interactions + spatial weighting
# =============================================================================
message("[Step 5] Ligand–receptor interactions...")

interactions <- read.csv(lr_path, row.names = 1)

# Marker genes per cell type (requires celltype_profile)
if ("celltype_profile" %in% colnames(spatial_obj@meta.data)) {
  
  Idents(spatial_obj) <- "celltype_profile"
  celltype_markers <- FindAllMarkers(spatial_obj, only.pos = TRUE)
  
  write.csv(celltype_markers, file.path(out_dir, "celltype_markers.csv"), row.names = FALSE)
  
  top5 <- celltype_markers %>%
    group_by(cluster) %>%
    filter(avg_log2FC > 1) %>%
    slice_head(n = 5) %>%
    ungroup()
  
  p_dot <- DotPlot(spatial_obj, features = unique(top5$gene)) + RotatedAxis()
  save_plot(p_dot, "dotplot_top5_markers_per_celltype.png", width = 10, height = 6)
  
  # Keep only ligands/receptors
  overexpressed_lr <- subset(
    celltype_markers,
    gene %in% interactions$ligand | gene %in% interactions$receptor
  )
  write.csv(overexpressed_lr, file.path(out_dir, "overexpressed_ligands_receptors.csv"), row.names = FALSE)
  
  # Interaction scores (product of logFC)
  resulting_interactions <- tibble(sender = character(),
                                   receiver = character(),
                                   ligand = character(),
                                   receptor = character(),
                                   score = numeric())
  
  for (ct in unique(overexpressed_lr$cluster)) {
    sender_tmp  <- subset(overexpressed_lr, cluster == ct)
    ligands_tmp <- intersect(sender_tmp$gene, interactions$ligand)
    
    for (l in ligands_tmp) {
      receptors_tmp <- subset(interactions, ligand == l)$receptor
      receivers_tmp <- subset(overexpressed_lr, gene %in% receptors_tmp)
      
      if (nrow(receivers_tmp) == 0) next
      
      for (i in seq_len(nrow(receivers_tmp))) {
        score_tmp <- subset(sender_tmp, gene == l)$avg_log2FC * receivers_tmp[i, "avg_log2FC"]
        resulting_interactions <- bind_rows(
          resulting_interactions,
          tibble(sender = ct,
                 receiver = receivers_tmp[i, "cluster"],
                 ligand = l,
                 receptor = receivers_tmp[i, "gene"],
                 score = as.numeric(score_tmp))
        )
      }
    }
  }
  
  resulting_interactions$cell_type_pair <- paste(resulting_interactions$sender,
                                                 resulting_interactions$receiver, sep = "-->")
  resulting_interactions$ligand_receptor_pair <- paste(resulting_interactions$ligand,
                                                       resulting_interactions$receptor, sep = "-->")
  
  write.csv(resulting_interactions, file.path(out_dir, "lr_interactions_scores.csv"), row.names = FALSE)
  
  # Dotplot (top interactions per sender)
  top5_int <- resulting_interactions %>%
    group_by(sender) %>%
    arrange(desc(score)) %>%
    slice_head(n = 5) %>%
    ungroup()
  
  p_lr <- ggplot(top5_int) +
    geom_dotplot(aes(x = cell_type_pair, y = ligand_receptor_pair, fill = score),
                 binwidth = 2, stackdir = "center", dotsize = 0.5, method = "histodot") +
    RotatedAxis() + NoLegend() +
    ggtitle("Top LR interactions (score = ligand logFC × receptor logFC)")
  save_plot(p_lr, "lr_dotplot_top5_per_sender.png", width = 12, height = 6)
  
  # --- Spatial weighting (centroids) ----------------------------------------
  message("[Step 5] Spatial weighting using centroid distances...")
  
  coords <- GetTissueCoordinates(spatial_obj)
  
  avg_coords_per_ct <- data.frame(x = numeric(), y = numeric(), celltype = character())
  for (ct in levels(spatial_obj$celltype_profile)) {
    cells_tmp <- rownames(subset(spatial_obj@meta.data, celltype_profile == ct))
    coords_tmp <- subset(coords, cell %in% cells_tmp)
    if (nrow(coords_tmp) == 0) next
    avg_tmp <- colMeans(coords_tmp[, 1:2, drop = FALSE])
    avg_coords_per_ct <- rbind(avg_coords_per_ct,
                               data.frame(x = avg_tmp[1], y = avg_tmp[2], celltype = ct))
  }
  
  p_centroids <- ggplot(avg_coords_per_ct) +
    geom_point(aes(x = x, y = y, color = celltype)) +
    coord_fixed() +
    ggtitle("Cell-type centroids (tissue space)")
  save_plot(p_centroids, "celltype_centroids.png", width = 8, height = 6)
  
  avg_dists <- dist(avg_coords_per_ct[, c("x", "y")])
  
  avg_sim <- 1 / (as.matrix(avg_dists) / max(avg_dists))
  rownames(avg_sim) <- avg_coords_per_ct$celltype
  colnames(avg_sim) <- avg_coords_per_ct$celltype
  diag(avg_sim) <- max(avg_sim[avg_sim != Inf], na.rm = TRUE) + 5
  
  pdf(file.path(fig_dir, "centroid_similarity_heatmap.pdf"), width = 9, height = 8)
  heatmap(avg_sim, Rowv = NA, Colv = NA, scale = "none", revC = TRUE)
  dev.off()
  
  spatial_score <- numeric(nrow(resulting_interactions))
  for (i in seq_len(nrow(resulting_interactions))) {
    sd_tmp <- resulting_interactions$sender[i]
    rc_tmp <- resulting_interactions$receiver[i]
    if (!(sd_tmp %in% rownames(avg_sim) && rc_tmp %in% colnames(avg_sim))) next
    spatial_score[i] <- avg_sim[sd_tmp, rc_tmp] * resulting_interactions$score[i]
  }
  resulting_interactions$spatial_score <- spatial_score
  write.csv(resulting_interactions, file.path(out_dir, "lr_interactions_with_centroid_spatial_score.csv"),
            row.names = FALSE)
  
  # --- Neighborhood-based weighting -----------------------------------------
  message("[Step 5] Neighborhood-based spatial weighting...")
  
  coords_df <- coords
  # Neighborhood matrix: counts of neighboring cell types
  neighborhood_matrix <- matrix(nrow = nrow(coords_df), ncol = length(levels(spatial_obj$celltype_profile)))
  
  for (i in seq_len(nrow(coords_df))) {
    cell_i <- coords_df$cell[i]
    x_min <- coords_df$x[i] - neighbor_radius
    x_max <- coords_df$x[i] + neighbor_radius
    y_min <- coords_df$y[i] - neighbor_radius
    y_max <- coords_df$y[i] + neighbor_radius
    
    cells_in_int <- subset(coords_df, x > x_min & x < x_max & y > y_min & y < y_max)$cell
    neighborhood_matrix[i, ] <- as.vector(table(spatial_obj@meta.data[cells_in_int, "celltype_profile"]))
  }
  
  rownames(neighborhood_matrix) <- coords_df$cell
  colnames(neighborhood_matrix) <- levels(spatial_obj$celltype_profile)
  saveRDS(neighborhood_matrix, file.path(out_dir, "spatial_neighbors.rds"))
  
  # Average neighborhood composition per cell type
  neighborhood_scores_normalized <- matrix(nrow = length(levels(spatial_obj$celltype_profile)),
                                           ncol = length(levels(spatial_obj$celltype_profile)))
  rownames(neighborhood_scores_normalized) <- levels(spatial_obj$celltype_profile)
  colnames(neighborhood_scores_normalized) <- levels(spatial_obj$celltype_profile)
  
  for (ct in levels(spatial_obj$celltype_profile)) {
    cells_ct <- rownames(subset(spatial_obj@meta.data, celltype_profile == ct))
    mat_ct <- neighborhood_matrix[cells_ct, , drop = FALSE]
    neighborhood_scores_normalized[ct, ] <- colSums(mat_ct) / table(spatial_obj$celltype_profile)
  }
  
  pdf(file.path(fig_dir, "neighborhood_scores_normalized_heatmap.pdf"), width = 9, height = 8)
  heatmap(neighborhood_scores_normalized, Rowv = NA, Colv = NA, scale = "none", revC = TRUE)
  dev.off()
  
  spatial_score_improved <- numeric(nrow(resulting_interactions))
  for (i in seq_len(nrow(resulting_interactions))) {
    sd_tmp <- resulting_interactions$sender[i]
    rc_tmp <- resulting_interactions$receiver[i]
    if (!(sd_tmp %in% rownames(neighborhood_scores_normalized) && rc_tmp %in% colnames(neighborhood_scores_normalized))) next
    spatial_score_improved[i] <- neighborhood_scores_normalized[sd_tmp, rc_tmp] * resulting_interactions$score[i]
  }
  
  resulting_interactions$spatial_score_improved <- spatial_score_improved
  write.csv(resulting_interactions, file.path(out_dir, "lr_interactions_with_neighborhood_spatial_score.csv"),
            row.names = FALSE)
  
  p_lr_spatial <- ggplot(resulting_interactions) +
    geom_dotplot(aes(x = cell_type_pair, y = ligand_receptor_pair, fill = spatial_score_improved),
                 binwidth = 2, stackdir = "center", dotsize = 0.5, method = "histodot") +
    RotatedAxis() + NoLegend() +
    ggtitle("LR interactions weighted by neighborhood-based spatial proximity")
  save_plot(p_lr_spatial, "lr_dotplot_spatial_score_improved.png", width = 12, height = 7)
  
} else {
  message("[Step 5] Skipped: 'celltype_profile' missing (required for LR workflow).")
}

# ---- Save final object + session info --------------------------------------
saveRDS(spatial_obj, file.path(out_dir, "spatial_processed.rds"))
sink(file.path(out_dir, "sessionInfo.txt")); print(sessionInfo()); sink()

message("Done. Results written to: ", out_dir)