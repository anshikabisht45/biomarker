# Network-Based Therapy Optimization - Biomarker Analysis

## Overview
This repository contains DESeq2 differential expression analysis for identifying biomarkers in cancer therapy optimization studies.

## Analysis Details
- **Cell Line**: BT549/BT540
- **Comparison**: Knockout (KO) vs Control (CT)
- **Dataset**: GSE300385
- **Method**: DESeq2 differential expression analysis
- **Thresholds**: 
  - Adjusted p-value < 0.05
  - |log2FoldChange| > 1

## Results Summary
- **Upregulated genes**: 131
- **Downregulated genes**: 336
- **Total differentially expressed genes**: 467

## Files

### Analysis Script
- `deseq2_analysis.R` - R script for DESeq2 differential expression analysis

### Output Files
- `upregulated_genes.txt` - List of 131 upregulated genes with statistics
- `downregulated_genes.txt` - List of 336 downregulated genes with statistics
- `deseq2_full_results.txt` - Complete DESeq2 results for all genes

### Visualizations
- `volcano_plot.pdf` - Volcano plot showing differentially expressed genes
- `ma_plot.pdf` - MA plot of fold changes vs mean expression

## Requirements
- R (version 4.0 or higher)
- DESeq2 (Bioconductor package)
- tidyverse

## Usage
```R
# Run the analysis
Rscript deseq2_analysis.R
```

## Data
The raw count data files are not included in this repository due to size constraints. 
Please download the original data from GEO (GSE300385).

## Author
Anshika Bisht (@anshikabisht45)

## Date
February 2026
