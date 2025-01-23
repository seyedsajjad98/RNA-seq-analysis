# Differential Gene Expression Analysis in Lung Adenocarcinoma

This repository contains the code and results of a differential gene expression analysis of lung adenocarcinoma RNA-seq data. By comparing normal and tumor tissues, this study identified key molecular changes contributing to lung cancer progression.

## Overview
- **Goal**: To identify differentially expressed genes between normal and tumor lung tissues.
- **Data Source**: Genomic Data Commons (GDC) Portal.
- **Samples**: 20 normal and 20 tumor tissues.
- **Methodology**:
  - **Data Processing**: DESeq2 for variance-stabilized normalization.
  - **Visualization**: Volcano plots, heatmaps, PCA.

## Key Findings
- **Significant Genes**:
  - Upregulated: 5,374
  - Downregulated: 114
- **Implications**:
  - Potential biomarkers: TSPAN6, SCYL3.
  - Therapeutic targets for lung adenocarcinoma.

## Repository Structure
- `README.md`: Project overview.
- `RNA-seq-analysis.R`: R script for data analysis.
- `deseq2_results.csv`: Results of differential expression analysis.
- `docs/`: Documentation folder containing the detailed report PDF.

## Usage
1. Download the repository from GitHub.
2. Install required R packages (`DESeq2`, `ggplot2`, etc.).
3. Run `RNA-seq-analysis.R` in RStudio or an R environment to reproduce the analysis.

## Results
- **Volcano Plot**: Highlights significantly upregulated and downregulated genes.
- **Heatmap**: Displays the top 20 differentially expressed genes across samples.
- **PCA**: Shows distinct clustering of normal and tumor samples.

## License
This project is licensed under the MIT License.

## Detailed Report
For a comprehensive description of the methodology, results, and implications, refer to the detailed report in the `docs` folder:
[Differential Gene Expression Analysis in Lung Adenocarcinoma](docs/Differential_Gene_Expression_Analysis_Lung_Adenocarcinoma.pdf)
