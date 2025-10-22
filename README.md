# Differential Gene Expression (DGE) Analysis Pipeline

## Project Summary
This repository contains a full bioinformatics pipeline for identifying **Differential Gene Expression (DGE)** between two simulated conditions (Treated vs. Control).


| Area | Tool / Concept |
| :--- | :--- |
| **Data Preparation** | **Python** (Pandas, NumPy) for raw data generation and formatting. |
| **Statistical Analysis** | **R** and **DESeq2** package. |
| **Visualization** | **R's ggplot2** for generating the specialized **Volcano Plot**. |
| **Data Filtering** | Applying statistical thresholds (Adjusted P-value < 0.05, Fold Change > 1.5). |

---

## Scientific Goal
The primary objective was to process a raw gene expression count matrix and identify statistically significant biomarkers that are either up- or downregulated in the Treated condition. The analysis successfully identified a single significant gene based on the applied stringency thresholds.

---

## Repository Contents

| File Name | Description |
| :--- | :--- |
| `dge_data_simulator.py` | **Python Script:** Generates the input files (`.csv`) for the R analysis. |
| `dge_analysis.r` | **R Script:** Performs the core DESeq2 analysis, statistical filtering, and visualization. |
| `mock_expression_counts.csv` | Input: Simulated raw count data (10,000 genes x 12 samples). |
| `sample_metadata.csv` | Input: Labels assigning each sample to the 'Treated' or 'Control' group. |
| `dge_report.md` | **Final Report:** Summary of findings and the single significant gene identified. |
| `volcano_plot.png` | **Visual Output:** The final DGE visualization plot. |

---

## How to Run the Pipeline

This project is designed to be run in two sequential steps.

### Step 1: Data Generation (Python)
1.  Ensure you have **Python 3** installed with `pandas` and `numpy`.
2.  Run the generation script from your terminal:
    ```bash
    python dge_data_simulator.py
    ```
    *(This creates the two necessary CSV files.)*

### Step 2: Analysis and Plotting (R)
1.  Ensure you have **R** and **RStudio** installed.
2.  Open RStudio, set the working directory to this folder.
3.  Open `dge_analysis.r`. Run the installation commands for `DESeq2`, `ggplot2`, and `dplyr` (if necessary).
4.  Run the entire `dge_analysis.r` script.
    *(This generates the `volcano_plot.png` and `significant_genes_report.csv`.)*

---

## Final Result

The volcano plot is the primary output, visually communicating the distribution of gene expression changes against statistical significance.
