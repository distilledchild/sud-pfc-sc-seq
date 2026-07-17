# SUD-PFC Seurat Reanalysis

- Run: `mclover3_D10_12_2024-07-09`
- Sample label: `F344_SHR_M_E007_E118`
- Random seed: `42`
- Created: `2026-07-17 01:07:46 CDT`

## QC

- Minimum-feature import retained **1,838** barcodes.
- Legacy upper-bound QC retained **1,790** barcodes (**97.4%**) and removed **48**.
- Upper-bound thresholds: nFeature_RNA < 2,000; nCount_RNA < 3,000; percent.mt < 3.

## Clustering

- Conservative resolution **0.6** produced **8 raw clusters** and **8 tree-ordered clusters**.
- This solution retained at least **90 barcodes per cluster** and had ARI **0.589** against the archived Scanpy partition.
- Legacy resolution **1.5** produced **12 clusters**, not 14, and its smallest cluster contained **2 barcodes**.
- No cluster count was forced. Resolution 1.5 is retained as a sensitivity result rather than the primary interpretation.
- **1,133** positive marker rows passed min.pct=0.25 and logFC=0.25.

## Scanpy comparison

- **1,790** post-QC barcodes matched the archived Scanpy assignments.
- Adjusted Rand Index between the conservative Seurat clusters and Scanpy Leiden clusters: **0.589**.

## Annotation boundary

- Marker-panel labels are computational candidates, not validated biological identities.
- The dataset has one combined biological sample label and cannot support strain or SUD contrasts.
- Doublet filtering was not applied because no validated platform-specific expected doublet rate was archived.

## Outputs

- `tables/`: QC, resolutions, clusters, markers, annotations, and Scanpy crosswalk.
- `figures/`: QC, PCA, resolution sweep, UMAPs, marker evidence, and crosswalk.
- `objects/seurat_reanalysis.rds`: completed Seurat object.
- `logs/session_info.txt`: package and R versions.
