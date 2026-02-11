#!/usr/bin/env Rscript
# -----------------------------------------------------------------------------
# Bulk RNA-seq QC + normalization (RPK/TPM) + PCA & clustering
# Author: Andrés Catalán Tatay
#
# Pipeline overview:
#   Step 1: Normalize raw counts by gene length and library size (RPK / TPM)
#   Step 2: Assess technical biases (gene length, sequencing depth, mean–variance)
#   Step 3: PCA + unsupervised clustering of samples
#
# Inputs (project root expected):
#   data/bulk/countmatrix/expression_Matrix_BRCA.txt
#   data/bulk/genes_info.txt
#   data/bulk/meta_data/TCGA_combined_meta_data_BRCA.txt
# Optional:
#   data/bulk/countmatrix/expression_Matrix_BRCA_VST-transformation.txt
#
# Outputs:
#   results/bulk/normalized/   (RPK + TPM matrices)
#   results/bulk/figures/      (QC + PCA/cluster plots)
#   results/bulk/             (cluster tables + sample cluster labels + sessionInfo)
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(ggplot2)
  library(matrixStats)
})

# ---- Paths -----------------------------------------------------------------
project_root <- getwd()
data_dir     <- file.path(project_root, "data", "bulk")

counts_path  <- file.path(data_dir, "countmatrix", "expression_Matrix_BRCA.txt")
gene_path    <- file.path(data_dir, "genes_info.txt")
meta_path    <- file.path(data_dir, "meta_data", "TCGA_combined_meta_data_BRCA.txt")
vst_path     <- file.path(data_dir, "countmatrix", "expression_Matrix_BRCA_VST-transformation.txt") # optional

out_dir      <- file.path(project_root, "results", "bulk")
norm_dir     <- file.path(out_dir, "normalized")
fig_dir      <- file.path(out_dir, "figures")
dir.create(norm_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir,  recursive = TRUE, showWarnings = FALSE)

# Fallbacks (useful when running in other environments)
if (!file.exists(gene_path) && file.exists("/mnt/data/genes_info.txt")) {
  gene_path <- "/mnt/data/genes_info.txt"
}
if (!file.exists(meta_path) && file.exists("/mnt/data/TCGA_combined_meta_data_BRCA.txt")) {
  meta_path <- "/mnt/data/TCGA_combined_meta_data_BRCA.txt"
}

# Early sanity checks (fail fast with clear error)
stopifnot(
  "counts_path not found" = file.exists(counts_path),
  "gene_path not found"   = file.exists(gene_path),
  "meta_path not found"   = file.exists(meta_path)
)
# ---- Helpers ---------------------------------------------------------------
read_counts <- function(path) {
  x <- read.table(path, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
  
  # If first column is an ID (non-numeric), move to rownames
  if (!is.numeric(x[[1]])) {
    rownames(x) <- x[[1]]
    x[[1]] <- NULL
  }
  
  m <- as.matrix(x)
  storage.mode(m) <- "numeric"
  m
}

read_gene_info <- function(path) {
  
  x <- read.table(path, header = TRUE, sep = "", stringsAsFactors = FALSE, check.names = FALSE, quote = "", comment.char = "")
  
  # Standardize "Length" column name (Length / gene_length / length)
  if ("Length" %in% colnames(x)) {
    x$Length <- as.numeric(x$Length)
  } else if ("gene_length" %in% colnames(x)) {
    x$Length <- as.numeric(x$gene_length)
  } else if ("length" %in% colnames(x)) {
    x$Length <- as.numeric(x$length)
  } else {
    stop("No length column found (Length/gene_length/length) in gene_info")
  }
  
  # Standardize gene ID (optional; not used for re-mapping in this script)
  if (!("gene_id" %in% colnames(x))) {
    if ("ensembl_gene_id_1" %in% colnames(x)) x$gene_id <- x$ensembl_gene_id_1
    else if ("ensembl_gene_id" %in% colnames(x)) x$gene_id <- x$ensembl_gene_id
  }
  
  return(x)
}

save_plot <- function(p, filename, width = 6.5, height = 4.5) {
  ggsave(file.path(fig_dir, filename), p, width = width, height = height)
}

# ---- Load data --------------------------------------------------------------
message("[Load] Reading count matrix...")
counts_raw <- read_counts(counts_path)   # expected: genes x samples

# If matrix is samples x genes, transpose to genes x samples

if (nrow(counts_raw) < ncol(counts_raw)) {
  message("Detected samples x genes matrix (", nrow(counts_raw), " x ", ncol(counts_raw), "). Transposing to genes x samples.")
  counts_raw <- t(counts_raw)
}

message("[Load] Reading gene metadata...")  
gene_info  <- read_gene_info(gene_path)

message("[Load] Reading sample metadata...")
meta <- read.table(meta_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, quote = "", comment.char = "")

# NOTE:
# gene_info and counts_raw are assumed to come from the same reference / ordering.
# No explicit re-mapping by gene_id is performed (avoids mismatches across references).
stopifnot(
  "Mismatch: counts rows vs gene_info rows" = (nrow(counts_raw) == nrow(gene_info))
)

gene_lengths <- gene_info$Length
names(gene_lengths) <- rownames(counts_raw)

# Guard against impossible / missing lengths (prevents Inf/NA explosions)
if (anyNA(gene_lengths) || any(gene_lengths <= 0)) {
  stop("Invalid gene lengths detected (NA or <= 0). Check gene_info Length column.")
}

# ---- Step 1: Library size & gene length normalization (RPK / TPM) ------------------------------------

message("[Step 1] Computing RPK and TPM...")

# RPK = counts / (length_kb)
length_kb <- gene_lengths / 1000
rpk <- counts_raw / length_kb

# TPM: scale each sample so sum(RPK) = 1e6
scale_factors <- colSums(rpk, na.rm = TRUE) / 1e6
tpm <- sweep(rpk, 2, scale_factors, "/")

# Persist matrices (useful for downstream steps)
write.table(rpk, file.path(norm_dir, "rpk_matrix.tsv"), sep = "\t", quote = FALSE, col.names = NA)
write.table(tpm, file.path(norm_dir, "tpm_matrix.tsv"), sep = "\t", quote = FALSE, col.names = NA)

# ---- Step 2: Quality control and bias assessment ---------------------------------------------------------------

#Analyse 3 different biases: 
# - Gene length
# - Sample depth (library size)
# - Mean-Variance relationship

message("[Step 2] Running QC plots (length, depth, mean–variance)...")

# Avoid log10(0) in plots (do NOT modify matrices)
eps <- 1e-6

# Example gene for bias demo (use DPM1 if present else first gene)
gene_example <- if ("DPM1" %in% rownames(counts_raw)) "DPM1" else rownames(counts_raw)[1]

# Length bias (sample 1): raw vs RPK
len_df <- data.frame(
  gene   = rownames(counts_raw),
  length = pmax(gene_lengths, eps),
  counts = pmax(as.numeric(counts_raw[, 1]), eps),
  rpk    = pmax(as.numeric(rpk[, 1]), eps)
)

p_len_raw <- ggplot(len_df, aes(x = length, y = counts)) +
  geom_point(alpha = 0.35, size = 0.8) +
  scale_x_log10() + scale_y_log10() +
  labs(title = "Length bias (raw counts, sample 1)", x = "Gene length (log10)", y = "Counts (log10)")
save_plot(p_len_raw, "length_bias_raw_sample1.png")

p_len_rpk <- ggplot(len_df, aes(x = length, y = rpk)) +
  geom_point(alpha = 0.35, size = 0.8) +
  scale_x_log10() + scale_y_log10() +
  labs(title = "After RPK (sample 1)", x = "Gene length (log10)", y = "RPK (log10)")
save_plot(p_len_rpk, "length_bias_rpk_sample1.png")

# Depth bias: example gene vs library size, raw vs TPM
depth_df <- data.frame(
  sample    = colnames(counts_raw),
  depth_raw = pmax(colSums(counts_raw), eps),
  gene_raw  = pmax(as.numeric(counts_raw[gene_example, ]), eps),
  gene_tpm  = pmax(as.numeric(tpm[gene_example, ]), eps)
)

p_depth_raw <- ggplot(depth_df, aes(x = depth_raw, y = gene_raw)) +
  geom_point(alpha = 0.6, size = 1) +
  scale_x_log10() + scale_y_log10() +
  labs(title = paste0("Depth bias (raw): ", gene_example, " vs library size"),
       x = "Library size (log10)", y = paste0(gene_example, " counts (log10)"))
save_plot(p_depth_raw, "depth_bias_raw_example_gene.png")

p_depth_tpm <- ggplot(depth_df, aes(x = depth_raw, y = gene_tpm)) +
  geom_point(alpha = 0.6, size = 1) +
  scale_x_log10() + scale_y_log10() +
  labs(title = paste0("After TPM: ", gene_example, " vs library size"),
       x = "Library size (log10)", y = paste0(gene_example, " TPM (log10)"))
save_plot(p_depth_tpm, "depth_bias_tpm_example_gene.png")

# Mean–variance relationship (raw counts)
gene_means <- rowMeans(counts_raw)
gene_vars  <- rowVars(counts_raw)
mv_df <- data.frame(
  mean     = pmax(gene_means, eps),
  variance = pmax(gene_vars,  eps)
)
p_mv <- ggplot(mv_df, aes(x = mean, y = variance)) +
  geom_point(alpha = 0.35, size = 0.8) +
  scale_x_log10() + scale_y_log10() +
  labs(title = "Mean–variance relationship (raw counts)", x = "Mean (log10)", y = "Variance (log10)")
save_plot(p_mv, "mean_variance_raw.png")

# QC sumary
qc_summary <- data.frame(
  n_genes = nrow(counts_raw),
  n_samples = ncol(counts_raw),
  library_size_min = min(colSums(counts_raw)),
  library_size_median = median(colSums(counts_raw)),
  library_size_max = max(colSums(counts_raw)),
  stringsAsFactors = FALSE
)
write.csv(qc_summary, file.path(out_dir, "qc_summary.csv"), row.names = FALSE)

# ---- Step 3: Dimensionality reduction and sample clustering ----------------------------------------------

message("[Step 3] PCA + clustering...")

# Choose expression matrix for PCA:
#   - Prefer VST if available
#   - Otherwise use log2(TPM + 1)

expr <- NULL
if (file.exists(vst_path)) {
  message("[Step 3] Using VST matrix: ", vst_path)
  vst <- as.matrix(read.table(vst_path, header = TRUE, check.names = FALSE))
  # Robustly orient to samples x genes
  if (nrow(vst) == nrow(meta)) {
    expr <- vst
  } else if (ncol(vst) == nrow(meta)) {
    expr <- t(vst)
  } else {
    message("[Step 3] VST dims do not match metadata; falling back to log2(TPM+1).")
    expr <- t(log2(tpm + 1))
  }
} else {
  message("[Step 3] VST not found; using log2(TPM+1).")
  expr <- t(log2(tpm + 1))  # samples x genes
}

# Align metadata if possible (expects a sample_id column)
if ("sample_id" %in% colnames(meta)) {
  target_ids <- rownames(expr)
  if (is.null(target_ids) || any(is.na(target_ids))) {
    target_ids <- colnames(counts_raw)
    rownames(expr) <- target_ids
  }
  if (all(target_ids %in% meta$sample_id)) {
    meta <- meta[match(target_ids, meta$sample_id), , drop = FALSE]
  } else {
    message("[Step 3] Warning: sample_id in meta does not fully match expression rownames; keeping meta order.")
  }
} else {
  message("[Step 3] Warning: meta has no 'sample_id' column; keeping meta order.")
}

# Pick a subtype column to color PCA (common in BRCA examples)
subtype_col <- NULL
candidates <- c("BRCA_Subtype_PAM50", "PAM50", "subtype", "Subtype")
for (cc in candidates) {
  if (cc %in% colnames(meta)) { subtype_col <- cc; break }
}
if (is.null(subtype_col)) {
  subtype_col <- colnames(meta)[1]
  message("[Step 3] Warning: subtype column not found; using first meta column: ", subtype_col)
}


# PCA on top variable genes
gene_sd <- colSds(expr)
top_n <- 1000
top_genes <- names(sort(gene_sd, decreasing = TRUE))[seq_len(min(top_n, length(gene_sd)))]

X <- scale(expr[, top_genes, drop = FALSE], center = TRUE, scale = TRUE)
pca <- prcomp(X, rank. = 10)

pca_df <- as.data.frame(pca$x[, 1:5, drop = FALSE])
pca_df$sample_type <- meta[[subtype_col]]

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = sample_type)) +
  geom_point(size = 2, alpha = 0.9) +
  labs(title = "PCA (top variable genes)", color = subtype_col)
ggsave(file.path(fig_dir, "pca_pc1_pc2.png"), p_pca, width = 6.5, height = 4.5)

# Clustering in PCA space (k-means + hierarchical)

#K-means clustering
set.seed(1)
k <- 3
km <- kmeans(pca$x[, 1:5, drop = FALSE], centers = k, nstart = 50)
pca_df$cluster_kmeans <- factor(km$cluster)

p_km <- ggplot(pca_df, aes(x = PC1, y = PC2, color = cluster_kmeans, shape = sample_type)) +
  geom_point(size = 2, alpha = 0.9) +
  labs(title = "k-means clustering on PCA space", color = "k-means cluster", shape = subtype_col)
ggsave(file.path(fig_dir, "kmeans_on_pca.png"), p_km, width = 6.5, height = 4.5)

#Hierarchical clustering
d <- dist(pca$x[, 1:5, drop = FALSE])
hc <- hclust(d, method = "ward.D2")
hc_clusters <- cutree(hc, k = k)
pca_df$cluster_hclust <- factor(hc_clusters)

# Save contingency tables
tab_km <- table(pca_df$cluster_kmeans, pca_df$sample_type)
tab_hc <- table(pca_df$cluster_hclust, pca_df$sample_type)
write.csv(as.data.frame.matrix(tab_km), file.path(out_dir, "contingency_kmeans_vs_subtype.csv"))
write.csv(as.data.frame.matrix(tab_hc), file.path(out_dir, "contingency_hclust_vs_subtype.csv"))

# Export sample labels for downstream DE (DE should be done on raw counts with DESeq2/edgeR/limma-voom)
clusters_out <- data.frame(
  sample = rownames(expr),
  subtype = pca_df$sample_type,
  cluster_kmeans = pca_df$cluster_kmeans,
  cluster_hclust = pca_df$cluster_hclust,
  stringsAsFactors = FALSE
)
write.csv(clusters_out, file.path(out_dir, "sample_clusters.csv"), row.names = FALSE)

sink(file.path(out_dir, "sessionInfo.txt")); print(sessionInfo()); sink()
message("Done. Results written to: ", out_dir)
