#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(tibble)
  library(readr)
  library(jsonlite)
  library(scales)
  library(cluster)
})

options(stringsAsFactors = FALSE)
set.seed(42)

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) {
  normalizePath(args[[1]], mustWork = TRUE)
} else {
  normalizePath("~/Desktop/playground/sud-pfc-sc-seq", mustWork = TRUE)
}
output_root <- if (length(args) >= 2) {
  normalizePath(args[[2]], mustWork = FALSE)
} else {
  file.path(project_root, "analysis_results", "seurat_reanalysis_2026")
}

run_id <- "mclover3_D10_12_2024-07-09"
sample_id <- "F344_SHR_M_E007_E118"
dge_dir <- file.path(
  project_root,
  "data",
  run_id,
  "raw",
  sample_id,
  "DGE_filtered"
)
scanpy_report_dir <- file.path(
  project_root,
  "data",
  run_id,
  "raw",
  sample_id,
  "report"
)

tables_dir <- file.path(output_root, "tables")
figures_dir <- file.path(output_root, "figures")
objects_dir <- file.path(output_root, "objects")
logs_dir <- file.path(output_root, "logs")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(objects_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(logs_dir, "analysis.log")
log_connection <- file(log_file, open = "wt")
sink(log_connection, type = "output", split = TRUE)
sink(log_connection, type = "message")
log_closed <- FALSE
close_log <- function(trim_whitespace = FALSE) {
  if (!log_closed) {
    sink(type = "message")
    sink(type = "output")
    close(log_connection)
    log_closed <<- TRUE
  }
  if (trim_whitespace && file.exists(log_file)) {
    log_lines <- readLines(log_file, warn = FALSE)
    writeLines(sub("[[:blank:]]+$", "", log_lines), log_file)
  }
}
on.exit(close_log(), add = TRUE)

log_step <- function(message) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), message))
}

save_plot <- function(plot, filename, width = 10, height = 7) {
  ggsave(
    filename = file.path(figures_dir, filename),
    plot = plot,
    width = width,
    height = height,
    dpi = 300,
    bg = "white"
  )
}

theme_report <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        margin = margin(b = 16)
      ),
      plot.subtitle = element_text(
        hjust = 0.5,
        color = "#475569",
        margin = margin(b = 12)
      ),
      panel.grid.minor = element_blank(),
      legend.title = element_text(face = "bold"),
      plot.margin = margin(18, 24, 18, 24)
    )
}

adjusted_rand_index <- function(labels_a, labels_b) {
  valid <- !is.na(labels_a) & !is.na(labels_b)
  labels_a <- labels_a[valid]
  labels_b <- labels_b[valid]
  contingency <- table(labels_a, labels_b)
  choose_two <- function(x) x * (x - 1) / 2
  sum_cells <- sum(choose_two(contingency))
  sum_rows <- sum(choose_two(rowSums(contingency)))
  sum_cols <- sum(choose_two(colSums(contingency)))
  total_pairs <- choose_two(sum(contingency))
  if (total_pairs == 0) {
    return(NA_real_)
  }
  expected <- sum_rows * sum_cols / total_pairs
  maximum <- (sum_rows + sum_cols) / 2
  if (maximum == expected) {
    return(1)
  }
  (sum_cells - expected) / (maximum - expected)
}

required_files <- c(
  file.path(dge_dir, "count_matrix.mtx"),
  file.path(dge_dir, "cell_metadata.csv"),
  file.path(dge_dir, "all_genes.csv"),
  file.path(scanpy_report_dir, "cluster_assignment.csv")
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required input files: ", paste(missing_files, collapse = ", "))
}

log_step("Loading Parse DGE matrix and metadata")
cell_by_gene <- readMM(file.path(dge_dir, "count_matrix.mtx"))
cell_metadata <- read_csv(
  file.path(dge_dir, "cell_metadata.csv"),
  show_col_types = FALSE
)
genes <- read_csv(
  file.path(dge_dir, "all_genes.csv"),
  show_col_types = FALSE
)

if (nrow(cell_by_gene) != nrow(cell_metadata)) {
  stop("Matrix row count does not match cell metadata")
}
if (ncol(cell_by_gene) != nrow(genes)) {
  stop("Matrix column count does not match gene metadata")
}

cell_metadata$bc_wells <- make.unique(cell_metadata$bc_wells, sep = "_dup")
genes$gene_name <- make.unique(genes$gene_name, sep = "_dup")
rownames(cell_by_gene) <- cell_metadata$bc_wells
colnames(cell_by_gene) <- genes$gene_name
counts <- as(t(cell_by_gene), "dgCMatrix")
cell_metadata_df <- as.data.frame(cell_metadata)
rownames(cell_metadata_df) <- cell_metadata_df$bc_wells

input_summary <- tibble(
  metric = c(
    "input_barcodes",
    "input_genes",
    "input_nonzero_counts"
  ),
  value = c(
    ncol(counts),
    nrow(counts),
    length(counts@x)
  )
)
write_csv(input_summary, file.path(tables_dir, "input_summary.csv"))

log_step("Creating Seurat object with min.features=300 and min.cells=5")
seurat_object <- CreateSeuratObject(
  counts = counts,
  project = "SUD_PFC_snRNA",
  min.cells = 5,
  min.features = 300,
  meta.data = cell_metadata_df
)
seurat_object$orig.ident <- "F344_SHR_M_E007_E118"
seurat_object$technical_well <- seurat_object$bc1_well

mitochondrial_pattern <- paste0(
  "^(ATP6|ATP8|COX[1-3]|CYTB|mt-Rnr[12]|ND[1-6L]?|Trn[a-y](1|2)?)$"
)
seurat_object[["percent.mt"]] <- PercentageFeatureSet(
  seurat_object,
  pattern = mitochondrial_pattern
)

imported_barcodes <- ncol(seurat_object)
imported_genes <- nrow(seurat_object)

log_step(sprintf(
  "Seurat import retained %d barcodes and %d genes",
  imported_barcodes,
  imported_genes
))

qc_before <- seurat_object@meta.data |>
  rownames_to_column("barcode") |>
  select(
    barcode,
    bc1_well,
    bc2_well,
    bc3_well,
    nCount_RNA,
    nFeature_RNA,
    percent.mt
  )
write_csv(qc_before, file.path(tables_dir, "qc_before_filtering.csv"))

qc_violin_before <- VlnPlot(
  seurat_object,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  pt.size = 0.1,
  ncol = 3
) +
  plot_annotation(
    title = "Seurat QC Before Upper-Bound Filtering",
    subtitle = sprintf(
      "%s | %s barcodes after minimum-feature import",
      sample_id,
      comma(imported_barcodes)
    ),
    theme = theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 18,
        margin = margin(b = 8)
      ),
      plot.subtitle = element_text(
        hjust = 0.5,
        color = "#475569",
        margin = margin(b = 16)
      )
    )
  )
save_plot(qc_violin_before, "01_qc_violin_before.png", 12, 6)

qc_scatter_before <- FeatureScatter(
  seurat_object,
  feature1 = "nCount_RNA",
  feature2 = "nFeature_RNA"
) +
  ggtitle("RNA Counts vs Detected Features Before Filtering") +
  theme_report()
save_plot(qc_scatter_before, "02_qc_scatter_before.png", 8, 6)

log_step("Applying legacy upper-bound QC thresholds")
seurat_object <- subset(
  seurat_object,
  subset = nFeature_RNA < 2000 &
    nCount_RNA < 3000 &
    percent.mt < 3
)

filtered_barcodes <- ncol(seurat_object)
filtered_genes <- nrow(seurat_object)
removed_barcodes <- imported_barcodes - filtered_barcodes

if (filtered_barcodes < 50) {
  stop("QC filtering retained fewer than 50 barcodes; stop before clustering")
}

log_step(sprintf(
  "QC retained %d/%d barcodes (removed %d)",
  filtered_barcodes,
  imported_barcodes,
  removed_barcodes
))

qc_after <- seurat_object@meta.data |>
  rownames_to_column("barcode") |>
  select(
    barcode,
    bc1_well,
    bc2_well,
    bc3_well,
    nCount_RNA,
    nFeature_RNA,
    percent.mt
  )
write_csv(qc_after, file.path(tables_dir, "qc_after_filtering.csv"))

qc_summary <- tibble(
  metric = c(
    "barcodes_after_minimum_feature_import",
    "genes_after_minimum_cell_import",
    "barcodes_after_upper_bound_qc",
    "barcodes_removed_by_upper_bound_qc",
    "barcode_retention_fraction",
    "median_features_after_qc",
    "median_counts_after_qc",
    "median_percent_mt_after_qc"
  ),
  value = c(
    imported_barcodes,
    imported_genes,
    filtered_barcodes,
    removed_barcodes,
    filtered_barcodes / imported_barcodes,
    median(seurat_object$nFeature_RNA),
    median(seurat_object$nCount_RNA),
    median(seurat_object$percent.mt)
  )
)
write_csv(qc_summary, file.path(tables_dir, "qc_summary.csv"))

qc_violin_after <- VlnPlot(
  seurat_object,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  pt.size = 0.1,
  ncol = 3
) +
  plot_annotation(
    title = "Seurat QC After Upper-Bound Filtering",
    subtitle = sprintf(
      "%s barcodes retained (%.1f%%)",
      comma(filtered_barcodes),
      100 * filtered_barcodes / imported_barcodes
    ),
    theme = theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 18,
        margin = margin(b = 8)
      ),
      plot.subtitle = element_text(
        hjust = 0.5,
        color = "#475569",
        margin = margin(b = 16)
      )
    )
  )
save_plot(qc_violin_after, "03_qc_violin_after.png", 12, 6)

log_step("Normalizing, selecting 2,000 variable genes, scaling, and running PCA")
seurat_object <- NormalizeData(
  seurat_object,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)
seurat_object <- FindVariableFeatures(
  seurat_object,
  selection.method = "vst",
  nfeatures = 2000,
  verbose = FALSE
)
seurat_object <- ScaleData(
  seurat_object,
  features = VariableFeatures(seurat_object),
  verbose = FALSE
)
seurat_object <- RunPCA(
  seurat_object,
  features = VariableFeatures(seurat_object),
  npcs = 50,
  seed.use = 42,
  verbose = FALSE
)

pca_stdev <- Stdev(seurat_object, reduction = "pca")
pca_variance <- pca_stdev^2 / sum(pca_stdev^2)
pca_summary <- tibble(
  pc = seq_along(pca_variance),
  variance_fraction = pca_variance,
  cumulative_variance_fraction = cumsum(pca_variance)
)
write_csv(pca_summary, file.path(tables_dir, "pca_variance.csv"))

elbow_plot <- ElbowPlot(seurat_object, ndims = 50) +
  ggtitle("PCA Elbow Plot") +
  theme_report()
save_plot(elbow_plot, "04_pca_elbow.png", 8, 6)

pcs_to_use <- 1:30
seurat_object <- FindNeighbors(
  seurat_object,
  dims = pcs_to_use,
  verbose = FALSE
)

resolutions <- c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.5)
resolution_columns <- character(length(resolutions))
resolution_cluster_counts <- integer(length(resolutions))

log_step("Running clustering resolution sweep")
for (index in seq_along(resolutions)) {
  resolution <- resolutions[[index]]
  seurat_object <- FindClusters(
    seurat_object,
    resolution = resolution,
    random.seed = 42,
    algorithm = 1,
    verbose = FALSE
  )
  column_name <- paste0(
    "seurat_clusters_res_",
    gsub("\\.", "_", format(resolution, nsmall = 1))
  )
  seurat_object[[column_name]] <- as.character(Idents(seurat_object))
  resolution_columns[[index]] <- column_name
  resolution_cluster_counts[[index]] <- nlevels(Idents(seurat_object))
  log_step(sprintf(
    "Resolution %.1f produced %d clusters",
    resolution,
    resolution_cluster_counts[[index]]
  ))
}

scanpy_assignments <- read_csv(
  file.path(scanpy_report_dir, "cluster_assignment.csv"),
  show_col_types = FALSE
) |>
  transmute(
    barcode = bc_wells,
    scanpy_cluster = as.character(cluster)
  )
scanpy_map <- setNames(
  scanpy_assignments$scanpy_cluster,
  scanpy_assignments$barcode
)
scanpy_labels <- unname(scanpy_map[colnames(seurat_object)])
pca_distance <- dist(Embeddings(seurat_object, "pca")[, pcs_to_use])

minimum_cluster_sizes <- integer(length(resolutions))
median_cluster_sizes <- numeric(length(resolutions))
clusters_below_ten <- integer(length(resolutions))
scanpy_ari_by_resolution <- numeric(length(resolutions))
mean_silhouette_by_resolution <- numeric(length(resolutions))
technical_well_cramers_v <- numeric(length(resolutions))

for (index in seq_along(resolutions)) {
  cluster_labels_at_resolution <- seurat_object@meta.data[[
    resolution_columns[[index]]
  ]]
  cluster_sizes_at_resolution <- table(cluster_labels_at_resolution)
  minimum_cluster_sizes[[index]] <- min(cluster_sizes_at_resolution)
  median_cluster_sizes[[index]] <- median(cluster_sizes_at_resolution)
  clusters_below_ten[[index]] <- sum(cluster_sizes_at_resolution < 10)
  scanpy_ari_by_resolution[[index]] <- adjusted_rand_index(
    cluster_labels_at_resolution,
    scanpy_labels
  )
  silhouette_result <- silhouette(
    as.integer(factor(cluster_labels_at_resolution)),
    pca_distance
  )
  mean_silhouette_by_resolution[[index]] <- mean(silhouette_result[, 3])
  technical_table <- table(
    cluster_labels_at_resolution,
    seurat_object$technical_well
  )
  technical_test <- suppressWarnings(chisq.test(technical_table))
  technical_well_cramers_v[[index]] <- sqrt(
    as.numeric(technical_test$statistic) /
      (
        sum(technical_table) *
          min(nrow(technical_table) - 1, ncol(technical_table) - 1)
      )
  )
}

resolution_summary <- tibble(
  resolution = resolutions,
  metadata_column = resolution_columns,
  n_clusters = resolution_cluster_counts,
  minimum_cluster_size = minimum_cluster_sizes,
  median_cluster_size = median_cluster_sizes,
  clusters_below_ten = clusters_below_ten,
  ari_vs_archived_scanpy = scanpy_ari_by_resolution,
  mean_pca_silhouette = mean_silhouette_by_resolution,
  technical_well_cramers_v = technical_well_cramers_v
)

if (length(resolutions) > 1) {
  adjacent_ari <- c(NA_real_)
  for (index in 2:length(resolutions)) {
    adjacent_ari[[index]] <- adjusted_rand_index(
      seurat_object@meta.data[[resolution_columns[[index - 1]]]],
      seurat_object@meta.data[[resolution_columns[[index]]]]
    )
  }
  resolution_summary$ari_vs_previous_resolution <- adjacent_ari
}
write_csv(resolution_summary, file.path(tables_dir, "resolution_summary.csv"))

resolution_plot <- ggplot(
  resolution_summary,
  aes(x = resolution, y = n_clusters)
) +
  geom_line(color = "#0891b2", linewidth = 1.2) +
  geom_point(color = "#0e7490", size = 3) +
  geom_text(
    aes(label = n_clusters),
    vjust = -0.8,
    fontface = "bold",
    color = "#164e63"
  ) +
  scale_x_continuous(breaks = resolutions) +
  scale_y_continuous(breaks = pretty_breaks()) +
  labs(
    title = "Cluster Count Across Seurat Resolutions",
    subtitle = "Resolution 1.5 is the legacy target; no cluster count is forced",
    x = "Resolution",
    y = "Number of clusters"
  ) +
  theme_report()
save_plot(resolution_plot, "05_resolution_sweep.png", 9, 6)

resolution_ari_plot <- ggplot(
  resolution_summary,
  aes(x = resolution, y = ari_vs_archived_scanpy)
) +
  geom_line(color = "#0e7490", linewidth = 1.2) +
  geom_point(color = "#0891b2", size = 3) +
  geom_vline(
    xintercept = c(0.6, 1.5),
    linetype = "dashed",
    color = c("#0891b2", "#f97316")
  ) +
  scale_x_continuous(breaks = resolutions) +
  scale_y_continuous(limits = c(0, 0.65), labels = number_format(accuracy = 0.01)) +
  labs(
    title = "Agreement With Archived Scanpy Clusters",
    subtitle = "ARI measures partition agreement, not biological truth",
    x = "Seurat resolution",
    y = "Adjusted Rand Index"
  ) +
  theme_report()

resolution_size_plot <- ggplot(
  resolution_summary,
  aes(x = resolution, y = minimum_cluster_size)
) +
  geom_line(color = "#c2410c", linewidth = 1.2) +
  geom_point(color = "#f97316", size = 3) +
  geom_text(
    aes(label = minimum_cluster_size),
    vjust = -0.8,
    fontface = "bold",
    color = "#9a3412"
  ) +
  scale_x_continuous(breaks = resolutions) +
  scale_y_continuous(expand = expansion(mult = c(0.04, 0.16))) +
  labs(
    title = "Smallest Cluster at Each Resolution",
    subtitle = "Resolution 1.5 creates a two-barcode outlier cluster",
    x = "Seurat resolution",
    y = "Minimum barcodes in a cluster"
  ) +
  theme_report()

resolution_diagnostics_plot <- resolution_ari_plot |
  resolution_size_plot
save_plot(
  resolution_diagnostics_plot,
  "05b_resolution_diagnostics.png",
  15,
  6
)

primary_resolution <- 0.6
legacy_resolution <- 1.5
primary_column <- resolution_columns[
  which(resolutions == primary_resolution)
]
legacy_column <- resolution_columns[
  which(resolutions == legacy_resolution)
]
seurat_object$seurat_cluster_legacy_res_1_5 <- as.character(
  seurat_object@meta.data[[legacy_column]]
)
legacy_cluster_levels <- as.character(sort(
  as.integer(unique(seurat_object$seurat_cluster_legacy_res_1_5))
))
seurat_object$seurat_cluster_legacy_res_1_5 <- factor(
  seurat_object$seurat_cluster_legacy_res_1_5,
  levels = legacy_cluster_levels
)
legacy_cluster_sizes <- table(
  seurat_object$seurat_cluster_legacy_res_1_5
)
legacy_cluster_count <- length(legacy_cluster_sizes)
legacy_minimum_cluster_size <- min(legacy_cluster_sizes)
legacy_cluster_counts <- as.data.frame(legacy_cluster_sizes) |>
  as_tibble() |>
  rename(cluster = Var1, barcodes = Freq) |>
  mutate(
    cluster = as.character(cluster),
    fraction = barcodes / sum(barcodes)
  )
write_csv(
  legacy_cluster_counts,
  file.path(tables_dir, "legacy_resolution_1_5_cluster_counts.csv")
)

Idents(seurat_object) <- seurat_object@meta.data[[primary_column]]
seurat_object$seurat_cluster_raw <- as.character(Idents(seurat_object))
primary_raw_cluster_count <- nlevels(Idents(seurat_object))

log_step(sprintf(
  paste(
    "Selected conservative resolution %.1f with %d raw clusters;",
    "legacy resolution %.1f produced %d clusters (minimum size %d)"
  ),
  primary_resolution,
  primary_raw_cluster_count,
  legacy_resolution,
  legacy_cluster_count,
  legacy_minimum_cluster_size
))

log_step("Running UMAP")
seurat_object <- RunUMAP(
  seurat_object,
  dims = pcs_to_use,
  seed.use = 42,
  n.neighbors = 30,
  min.dist = 0.3,
  verbose = FALSE
)

legacy_palette <- hue_pal(
  h = c(190, 20),
  c = 95,
  l = 55
)(legacy_cluster_count)
legacy_umap <- DimPlot(
  seurat_object,
  reduction = "umap",
  group.by = "seurat_cluster_legacy_res_1_5",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.7,
  cols = legacy_palette
) +
  labs(
    title = sprintf(
      "Legacy Seurat Target at Resolution %.1f: %d Clusters",
      legacy_resolution,
      legacy_cluster_count
    ),
    subtitle = sprintf(
      "The smallest cluster has %d barcodes; 14 clusters were not reproduced",
      legacy_minimum_cluster_size
    ),
    color = "Legacy cluster"
  ) +
  theme_report()
save_plot(legacy_umap, "06_legacy_resolution_1_5_umap.png", 10, 8)

Idents(seurat_object) <- seurat_object$seurat_cluster_raw
log_step("Building primary cluster tree and reproducible tree-based labels")
seurat_object <- BuildClusterTree(
  seurat_object,
  dims = pcs_to_use,
  reorder = TRUE,
  reorder.numeric = TRUE,
  verbose = FALSE
)
if (!"tree.ident" %in% colnames(seurat_object@meta.data)) {
  stop("BuildClusterTree did not create tree.ident")
}
Idents(seurat_object) <- seurat_object$tree.ident
seurat_object$seurat_tree_cluster <- as.character(Idents(seurat_object))
primary_tree_cluster_count <- nlevels(Idents(seurat_object))

cluster_counts <- as.data.frame(table(Idents(seurat_object))) |>
  as_tibble() |>
  rename(cluster = Var1, barcodes = Freq) |>
  mutate(
    cluster = as.character(cluster),
    fraction = barcodes / sum(barcodes)
  )
write_csv(cluster_counts, file.path(tables_dir, "cluster_counts.csv"))

cluster_palette <- hue_pal(
  h = c(190, 20),
  c = 95,
  l = 55
)(primary_tree_cluster_count)

umap_clusters <- DimPlot(
  seurat_object,
  reduction = "umap",
  group.by = "seurat_tree_cluster",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.7,
  cols = cluster_palette
) +
  labs(
    title = sprintf(
      "Conservative Seurat UMAP at Resolution %.1f: %d Clusters",
      primary_resolution,
      primary_tree_cluster_count
    ),
    subtitle = sprintf(
      "%s post-QC barcodes | selected for stable cluster size and Scanpy agreement",
      comma(filtered_barcodes)
    ),
    color = "Seurat cluster"
  ) +
  theme_report()
save_plot(umap_clusters, "06b_primary_resolution_0_6_umap.png", 10, 8)

umap_wells <- DimPlot(
  seurat_object,
  reduction = "umap",
  group.by = "technical_well",
  pt.size = 0.7,
  cols = c(D10 = "#0891b2", D11 = "#f97316", D12 = "#14b8a6")
) +
  labs(
    title = "Seurat UMAP Colored by Technical Round-1 Well",
    subtitle = "Well mixing is a technical diagnostic, not biological replication",
    color = "Well"
  ) +
  theme_report()
save_plot(umap_wells, "07_umap_technical_wells.png", 10, 8)

cluster_abundance_plot <- ggplot(
  cluster_counts,
  aes(
    x = reorder(cluster, as.numeric(cluster)),
    y = barcodes
  )
) +
  geom_col(fill = "#0e7490", width = 0.72) +
  geom_text(
    aes(label = barcodes),
    vjust = -0.4,
    size = 3.5,
    fontface = "bold"
  ) +
  labs(
    title = "Seurat Cluster Abundance",
    subtitle = sprintf(
      "%s called barcodes across %d tree-ordered clusters",
      comma(filtered_barcodes),
      primary_tree_cluster_count
    ),
    x = "Seurat tree cluster",
    y = "Barcodes"
  ) +
  theme_report()
save_plot(cluster_abundance_plot, "08_cluster_abundance.png", 10, 6)

log_step("Finding positive cluster markers")
markers <- FindAllMarkers(
  seurat_object,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  test.use = "wilcox",
  verbose = FALSE
) |>
  as_tibble()

if (nrow(markers) == 0) {
  warning("No markers passed the configured thresholds")
}

write_csv(markers, file.path(tables_dir, "markers_all.csv"))

top_markers <- markers |>
  group_by(cluster) |>
  arrange(p_val_adj, desc(avg_log2FC), .by_group = TRUE) |>
  slice_head(n = 20) |>
  ungroup()
write_csv(top_markers, file.path(tables_dir, "markers_top20.csv"))

marker_panels <- list(
  Excitatory_neuron = c(
    "Slc17a7",
    "Camk2a",
    "Satb2",
    "Tle4",
    "Slc30a3",
    "Neurod6",
    "Rorb",
    "Bcl11b"
  ),
  Inhibitory_neuron = c(
    "Gad1",
    "Gad2",
    "Slc6a1",
    "Slc32a1",
    "Erbb4",
    "Dlx1",
    "Dlx2",
    "Sst",
    "Vip"
  ),
  Astrocyte = c("Aqp4", "Gfap", "Aldoc", "Slc1a3"),
  Oligodendrocyte = c("Mbp", "Plp1", "Mobp", "Mag"),
  OPC = c("Pdgfra", "Cspg4", "Olig1", "Olig2"),
  Microglia = c("Aif1", "Cx3cr1", "P2ry12", "Csf1r"),
  Endothelial = c("Cldn5", "Flt1", "Pecam1", "Vwf"),
  Pericyte_VSMC = c("Pdgfrb", "Rgs5", "Acta2", "Kcnj8", "Abcc9"),
  Ependymal = c("Foxj1", "Pifo", "Tppp3"),
  Choroid_plexus = c("Ttr", "Klotho", "Aqp1")
)
marker_panels <- lapply(
  marker_panels,
  function(panel) panel[panel %in% rownames(seurat_object)]
)
marker_panels <- marker_panels[lengths(marker_panels) >= 2]

marker_dotplot <- DotPlot(
  seurat_object,
  features = marker_panels,
  group.by = "seurat_tree_cluster",
  dot.scale = 6
) +
  RotatedAxis() +
  labs(
    title = "Canonical Marker Evidence Across Seurat Clusters",
    subtitle = "Marker panels support provisional review; they are not automatic ground-truth labels",
    x = "Marker genes",
    y = "Seurat tree cluster",
    color = "Average expression",
    size = "Percent expressed"
  ) +
  theme_report(11) +
  theme(
    axis.text.x = element_text(angle = 55, hjust = 1),
    strip.text = element_text(face = "bold")
  )
save_plot(marker_dotplot, "09_marker_dotplot.png", 18, 9)

log_step("Calculating marker-panel scores by cluster")
normalized_data <- GetAssayData(seurat_object, assay = "RNA", layer = "data")
cluster_labels <- as.character(seurat_object$seurat_tree_cluster)
cluster_levels <- levels(Idents(seurat_object))

panel_score_rows <- list()
row_index <- 1
for (panel_name in names(marker_panels)) {
  panel_genes <- marker_panels[[panel_name]]
  for (cluster_label in cluster_levels) {
    cells_in_cluster <- which(cluster_labels == cluster_label)
    panel_matrix <- normalized_data[panel_genes, cells_in_cluster, drop = FALSE]
    mean_expression <- mean(Matrix::rowMeans(panel_matrix))
    mean_detection <- mean(Matrix::rowMeans(panel_matrix > 0))
    panel_score_rows[[row_index]] <- tibble(
      cluster = cluster_label,
      panel = panel_name,
      genes_available = length(panel_genes),
      genes = paste(panel_genes, collapse = ";"),
      mean_log_normalized_expression = mean_expression,
      mean_detection_fraction = mean_detection
    )
    row_index <- row_index + 1
  }
}

panel_scores <- bind_rows(panel_score_rows) |>
  group_by(panel) |>
  mutate(
    expression_z = if (sd(mean_log_normalized_expression) > 0) {
      as.numeric(scale(mean_log_normalized_expression))
    } else {
      0
    }
  ) |>
  ungroup()
write_csv(panel_scores, file.path(tables_dir, "marker_panel_scores.csv"))

provisional_annotations <- panel_scores |>
  group_by(cluster) |>
  arrange(desc(expression_z), desc(mean_detection_fraction), .by_group = TRUE) |>
  mutate(rank = row_number()) |>
  filter(rank <= 2) |>
  summarise(
    top_panel = panel[rank == 1],
    top_expression_z = expression_z[rank == 1],
    top_detection_fraction = mean_detection_fraction[rank == 1],
    runner_up_panel = panel[rank == 2],
    runner_up_expression_z = expression_z[rank == 2],
    score_margin = expression_z[rank == 1] - expression_z[rank == 2],
    .groups = "drop"
  ) |>
  mutate(
    candidate_label = if_else(
      top_expression_z >= 0.5 &
        score_margin >= 0.25 &
        top_detection_fraction >= 0.02,
      top_panel,
      "Unresolved"
    ),
    evidence_status = "Computational candidate; requires marker and reference review"
  ) |>
  left_join(
    select(cluster_counts, cluster, barcodes),
    by = "cluster"
  ) |>
  mutate(
    candidate_label = if_else(
      barcodes < 10,
      "Rare_outlier_unresolved",
      candidate_label
    )
  )
write_csv(
  provisional_annotations,
  file.path(tables_dir, "provisional_annotations.csv")
)

annotation_map <- setNames(
  provisional_annotations$candidate_label,
  provisional_annotations$cluster
)
seurat_object$provisional_annotation <- unname(
  annotation_map[as.character(seurat_object$seurat_tree_cluster)]
)

annotation_umap <- DimPlot(
  seurat_object,
  reduction = "umap",
  group.by = "provisional_annotation",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.7
) +
  labs(
    title = "Provisional Marker-Panel Annotation",
    subtitle = "Algorithmic candidates only; unresolved clusters remain explicit",
    color = "Candidate label"
  ) +
  theme_report()
save_plot(annotation_umap, "10_umap_provisional_annotations.png", 11, 8)

log_step("Building barcode-level Scanpy crosswalk")
cluster_assignments <- seurat_object@meta.data |>
  rownames_to_column("barcode") |>
  transmute(
    barcode,
    technical_well,
    seurat_cluster_raw,
    seurat_tree_cluster,
    seurat_cluster_legacy_res_1_5,
    provisional_annotation,
    nCount_RNA,
    nFeature_RNA,
    percent.mt
  ) |>
  left_join(scanpy_assignments, by = "barcode")
write_csv(cluster_assignments, file.path(tables_dir, "cluster_assignments.csv"))

matched_assignments <- cluster_assignments |>
  filter(!is.na(scanpy_cluster))
scanpy_ari <- adjusted_rand_index(
  matched_assignments$seurat_tree_cluster,
  matched_assignments$scanpy_cluster
)

crosswalk_counts <- as.data.frame(table(
  matched_assignments$seurat_tree_cluster,
  matched_assignments$scanpy_cluster
)) |>
  as_tibble() |>
  rename(
    seurat_cluster = Var1,
    scanpy_cluster = Var2,
    count = Freq
  ) |>
  group_by(seurat_cluster) |>
  mutate(seurat_cluster_fraction = count / sum(count)) |>
  ungroup()
write_csv(crosswalk_counts, file.path(tables_dir, "scanpy_crosswalk.csv"))

crosswalk_plot <- ggplot(
  crosswalk_counts,
  aes(
    x = factor(scanpy_cluster),
    y = factor(seurat_cluster),
    fill = seurat_cluster_fraction
  )
) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(
    aes(label = if_else(count > 0, as.character(count), "")),
    size = 3
  ) +
  scale_fill_gradient(
    low = "#ecfeff",
    high = "#0e7490",
    labels = percent_format()
  ) +
  labs(
    title = "Barcode Crosswalk: Conservative Seurat vs Archived Scanpy",
    subtitle = sprintf(
      "%s matched barcodes | Adjusted Rand Index = %.3f",
      comma(nrow(matched_assignments)),
      scanpy_ari
    ),
    x = "Archived Scanpy Leiden cluster",
    y = "Seurat tree cluster",
    fill = "Within-Seurat\nfraction"
  ) +
  coord_fixed() +
  theme_report()
save_plot(crosswalk_plot, "11_scanpy_seurat_crosswalk.png", 10, 9)

log_step("Saving final object and reproducibility metadata")
saveRDS(
  seurat_object,
  file.path(objects_dir, "seurat_reanalysis.rds"),
  compress = "xz"
)

session_lines <- sub(
  "[[:blank:]]+$",
  "",
  capture.output(sessionInfo())
)
writeLines(session_lines, file.path(logs_dir, "session_info.txt"))

primary_resolution_diagnostics <- resolution_summary |>
  filter(resolution == primary_resolution)
legacy_resolution_diagnostics <- resolution_summary |>
  filter(resolution == legacy_resolution)

analysis_summary <- list(
  analysis_name = "SUD-PFC Seurat reanalysis",
  created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  random_seed = 42,
  project_root = project_root,
  run_id = run_id,
  sample_id = sample_id,
  input = list(
    barcodes = ncol(counts),
    genes = nrow(counts),
    nonzero_counts = length(counts@x)
  ),
  qc = list(
    imported_barcodes = imported_barcodes,
    filtered_barcodes = filtered_barcodes,
    removed_barcodes = removed_barcodes,
    retention_fraction = filtered_barcodes / imported_barcodes,
    thresholds = list(
      minimum_features_at_import = 300,
      minimum_cells_per_gene_at_import = 5,
      maximum_features = 2000,
      maximum_counts = 3000,
      maximum_percent_mt = 3
    )
  ),
  dimension_reduction = list(
    variable_genes = length(VariableFeatures(seurat_object)),
    pcs_used = as.integer(pcs_to_use),
    umap_neighbors = 30,
    umap_min_dist = 0.3
  ),
  clustering = list(
    algorithm = "Louvain algorithm 1",
    primary = list(
      resolution = primary_resolution,
      raw_cluster_count = primary_raw_cluster_count,
      tree_cluster_count = primary_tree_cluster_count,
      minimum_cluster_size = primary_resolution_diagnostics$minimum_cluster_size,
      ari_vs_archived_scanpy = primary_resolution_diagnostics$ari_vs_archived_scanpy,
      mean_pca_silhouette = primary_resolution_diagnostics$mean_pca_silhouette,
      technical_well_cramers_v = primary_resolution_diagnostics$technical_well_cramers_v
    ),
    legacy_target = list(
      resolution = legacy_resolution,
      observed_cluster_count = legacy_cluster_count,
      minimum_cluster_size = legacy_minimum_cluster_size,
      ari_vs_archived_scanpy = legacy_resolution_diagnostics$ari_vs_archived_scanpy,
      mean_pca_silhouette = legacy_resolution_diagnostics$mean_pca_silhouette,
      technical_well_cramers_v = legacy_resolution_diagnostics$technical_well_cramers_v
    ),
    selection_rationale = paste(
      "Resolution 0.6 retained at least 90 barcodes per cluster and had",
      "the highest agreement with the archived Scanpy partition;",
      "resolution 1.5 was retained as a legacy sensitivity result."
    ),
    resolution_sweep = split(
      resolution_summary$n_clusters,
      resolution_summary$resolution
    )
  ),
  markers = list(
    all_marker_rows = nrow(markers),
    top_marker_rows = nrow(top_markers),
    provisional_annotation_rule = paste(
      "top panel z >= 0.5, margin >= 0.25,",
      "and mean detection >= 0.02"
    )
  ),
  scanpy_crosswalk = list(
    matched_barcodes = nrow(matched_assignments),
    adjusted_rand_index = scanpy_ari
  ),
  limitations = c(
    "One combined biological sample label; technical wells are not biological replicates.",
    "Doublet filtering was not applied because a validated platform-specific expected doublet rate was not archived.",
    "Provisional marker-panel labels require manual and reference-based review.",
    "Resolution 0.6 is a conservative computational choice, not proof of eight biological cell types.",
    "Resolution 1.5 produced 12 rather than 14 clusters and included a two-barcode outlier cluster."
  )
)
write_json(
  analysis_summary,
  file.path(output_root, "analysis_summary.json"),
  pretty = TRUE,
  auto_unbox = TRUE,
  digits = 6
)

annotation_counts <- provisional_annotations |>
  count(candidate_label, name = "clusters") |>
  arrange(desc(clusters), candidate_label)

report_lines <- c(
  "# SUD-PFC Seurat Reanalysis",
  "",
  sprintf("- Run: `%s`", run_id),
  sprintf("- Sample label: `%s`", sample_id),
  "- Random seed: `42`",
  sprintf("- Created: `%s`", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## QC",
  "",
  sprintf(
    "- Minimum-feature import retained **%s** barcodes.",
    comma(imported_barcodes)
  ),
  sprintf(
    "- Legacy upper-bound QC retained **%s** barcodes (**%.1f%%**) and removed **%s**.",
    comma(filtered_barcodes),
    100 * filtered_barcodes / imported_barcodes,
    comma(removed_barcodes)
  ),
  "- Upper-bound thresholds: nFeature_RNA < 2,000; nCount_RNA < 3,000; percent.mt < 3.",
  "",
  "## Clustering",
  "",
  sprintf(
    "- Conservative resolution **%.1f** produced **%d raw clusters** and **%d tree-ordered clusters**.",
    primary_resolution,
    primary_raw_cluster_count,
    primary_tree_cluster_count
  ),
  sprintf(
    "- This solution retained at least **%d barcodes per cluster** and had ARI **%.3f** against the archived Scanpy partition.",
    primary_resolution_diagnostics$minimum_cluster_size,
    primary_resolution_diagnostics$ari_vs_archived_scanpy
  ),
  sprintf(
    "- Legacy resolution **%.1f** produced **%d clusters**, not 14, and its smallest cluster contained **%d barcodes**.",
    legacy_resolution,
    legacy_cluster_count,
    legacy_minimum_cluster_size
  ),
  "- No cluster count was forced. Resolution 1.5 is retained as a sensitivity result rather than the primary interpretation.",
  sprintf(
    "- **%s** positive marker rows passed min.pct=0.25 and logFC=0.25.",
    comma(nrow(markers))
  ),
  "",
  "## Scanpy comparison",
  "",
  sprintf(
    "- **%s** post-QC barcodes matched the archived Scanpy assignments.",
    comma(nrow(matched_assignments))
  ),
  sprintf(
    "- Adjusted Rand Index between the conservative Seurat clusters and Scanpy Leiden clusters: **%.3f**.",
    scanpy_ari
  ),
  "",
  "## Annotation boundary",
  "",
  "- Marker-panel labels are computational candidates, not validated biological identities.",
  "- The dataset has one combined biological sample label and cannot support strain or SUD contrasts.",
  "- Doublet filtering was not applied because no validated platform-specific expected doublet rate was archived.",
  "",
  "## Outputs",
  "",
  "- `tables/`: QC, resolutions, clusters, markers, annotations, and Scanpy crosswalk.",
  "- `figures/`: QC, PCA, resolution sweep, UMAPs, marker evidence, and crosswalk.",
  "- `objects/seurat_reanalysis.rds`: completed Seurat object.",
  "- `logs/session_info.txt`: package and R versions."
)
writeLines(report_lines, file.path(output_root, "README.md"))

log_step("Analysis complete")
log_step(sprintf("Output directory: %s", output_root))
close_log(trim_whitespace = TRUE)
