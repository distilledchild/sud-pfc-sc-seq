# scikit-learn Workflow: snRNA-seq Exploratory Machine Learning

## Purpose

This folder records exploratory scikit-learn analyses based on the SUD/PFC single-nucleus RNA-seq project. The goal is to gain hands-on machine-learning experience with transcriptomic data while keeping the analysis biologically interpretable.

This workflow is separate from the Hi-C regulatory-genomics workflow because the unit of analysis, preprocessing, labels, and biological interpretation are different.

## Available Data Sources

Potential inputs may include:

- Processed snRNA-seq expression matrix
- Cell or nucleus metadata
- Cell-type labels
- Strain labels
- Cluster labels
- QC metrics
- Differentially expressed genes or marker genes, if available
- Pseudobulk expression summaries, if generated

## Project 5: snRNA-seq Cell-Type Or Strain Classification

### Biological Question

Can cell-type-specific or strain-specific transcriptomic signatures be separated using simple scikit-learn models?

### Machine-Learning Framing

Classification or exploratory dimensionality reduction.

### Possible Labels

Depending on available metadata:

- Cell type
- Cluster label
- Strain
- Condition or phenotype group, if available

### Possible Input Matrices

Option 1: Cell-level matrix

- Rows: cells or nuclei
- Columns: selected genes or principal components

Option 2: Pseudobulk matrix

- Rows: sample x cell-type combinations
- Columns: gene expression summaries

Pseudobulk is often more statistically stable when sample-level biological interpretation is important.

## Candidate Features

- Highly variable genes
- Marker genes
- PCA components from normalized expression matrix
- Pseudobulk expression of selected genes
- Pathway/module scores, if available

## Candidate Methods

### Exploratory Analysis

- StandardScaler
- PCA
- UMAP outside scikit-learn, if needed
- KMeans
- AgglomerativeClustering

### Classification

- LogisticRegression
- RandomForestClassifier
- LinearSVC
- GradientBoostingClassifier

## Evaluation

- Train/test split
- Cross-validation
- Accuracy, balanced accuracy, F1 score
- Confusion matrix
- Feature importance or model coefficients

## Important Cautions

### Avoid Cell-Level Leakage

If cells from the same biological sample appear in both train and test sets, the model can look artificially accurate. When possible, split by sample or strain, not only by individual cells.

### Avoid Over-Interpreting Cell-Type Prediction

If the label is an existing cell-type annotation, a classifier may simply learn marker genes already used for annotation. That is useful for hands-on ML practice but not necessarily a new biological discovery.

### Keep Biological Interpretation Modest

This analysis should initially be described as exploratory machine-learning practice and transcriptomic signature evaluation.

## Recommended Starting Order

1. Confirm available processed snRNA-seq object and metadata.
2. Export a compact expression table or PCA embedding with labels.
3. Start with PCA visualization and simple logistic regression.
4. Compare random split vs sample-aware split, if metadata allows.
5. Add random forest or gradient boosting only after the simple baseline works.

## Practical Implementation Plan

1. Use existing Seurat/Scanpy/R objects to export a clean CSV or `.h5ad`-derived table.
2. Keep exported data in a `data/` subfolder.
3. Use scikit-learn scripts or notebooks for modeling.
4. Save metrics, confusion matrices, and feature-importance plots in a `results/` subfolder.
5. Record label definitions, feature selection, and split strategy in this workflow.

## Resume-Friendly Outcome

Possible future resume bullet after completion:

Applied scikit-learn models to exploratory snRNA-seq data to evaluate cell-type- or strain-associated transcriptomic signatures, using PCA, logistic regression, tree-based classifiers, and sample-aware validation strategies.
