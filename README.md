# SUD-PFC Single-Nucleus RNA-seq Analysis

This repository contains reproducible Seurat analyses for the SUD/PFC
single-nucleus RNA-seq dataset, selected Parse Biosciences pipeline provenance,
and planning notes for follow-up machine-learning analyses.

## Repository layout

```text
analysis_results/
  seurat_reanalysis_2026/  # Versioned figures, tables, logs, and Seurat object
provenance/
  mclover3_D10_12_2024-07-09/  # Split-pipe summary, log, and input checksums
scripts/
  seurat_reanalysis.R      # Primary reproducible Seurat workflow
scikit-learn/
  workflow.md              # Follow-up machine-learning analysis plan
seurat_mclover3_D10_12.R   # Parameterized legacy exploratory workflow
data -> /Volumes/external_1000GB_all/playground/sud-pfc-sc-seq/data
legacy_results -> /Volumes/external_1000GB_all/playground/sud-pfc-sc-seq/legacy_results
mclover3_D10_12_2024-07-09 -> /Volumes/external_1000GB_all/playground/sud-pfc-sc-seq/mclover3_D10_12_2024-07-09
```

The three symbolic links are intentionally excluded from Git. The external
drive preserves the 27 GB source-data tree, the complete original Split-pipe
report bundle, and generated outputs from the legacy workflow.

## Required inputs

The primary reanalysis expects these files below the local `data` link:

```text
data/mclover3_D10_12_2024-07-09/raw/F344_SHR_M_E007_E118/
  DGE_filtered/count_matrix.mtx
  DGE_filtered/cell_metadata.csv
  DGE_filtered/all_genes.csv
  report/cluster_assignment.csv
```

These inputs are not committed because the complete source tree is 27 GB.
The filtered matrix and its metadata remain on the external drive, while the
small Scanpy assignments are represented in the versioned crosswalk output.

## Restore the external links

Mount the external drive at `/Volumes/external_1000GB_all`, then run:

```bash
ln -s /Volumes/external_1000GB_all/playground/sud-pfc-sc-seq/data data
ln -s /Volumes/external_1000GB_all/playground/sud-pfc-sc-seq/legacy_results legacy_results
ln -s /Volumes/external_1000GB_all/playground/sud-pfc-sc-seq/mclover3_D10_12_2024-07-09 mclover3_D10_12_2024-07-09
```

For a different storage location, create equivalent links with the same local
names. The analysis code reads through the repository-relative paths rather
than embedding the external drive path.

## Run the primary analysis

From the repository root:

```bash
Rscript scripts/seurat_reanalysis.R \
  "$(pwd)" \
  "$(pwd)/analysis_results/seurat_reanalysis_2026"
```

The completed run retained 1,790 nuclei after QC and selected resolution 0.6
as the conservative eight-cluster result. Resolution 1.5 produced 12 clusters
and is retained as a sensitivity result rather than being presented as a
reproduced 14-cluster solution.

## Run the legacy workflow

The legacy script accepts the repository root and an optional output root:

```bash
Rscript seurat_mclover3_D10_12.R \
  "$(pwd)" \
  "$(pwd)/legacy_results/mclover3_D10_12_2024-07-09"
```

Its figures, tables, and RDS objects are written under `legacy_results`, which
is stored on the external drive through the symbolic link.

## Provenance policy

Git contains the small, reviewable Split-pipe summary CSV and execution log.
It also records `input_checksums.sha256` for the four files consumed by the
primary reanalysis, so the external inputs can be verified without committing
the data itself.
The duplicate ZIP bundle, large generated HTML report, and historical QC image
remain in the external archive at:

```text
/Volumes/external_1000GB_all/playground/sud-pfc-sc-seq/mclover3_D10_12_2024-07-09
```

If the source data must be shared or versioned independently, use an
artifact store, DVC, or Git LFS rather than adding the external `data` tree to
regular Git history.
