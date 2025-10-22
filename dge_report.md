Differential Gene Expression Analysis Report

Project Summary

This project demonstrates proficiency in the standard RNA-seq data analysis pipeline. Differential Gene Expression (DGE) analysis was performed to identify genes significantly affected by a mock 'Treated' condition versus a 'Control' condition. The workflow successfully integrated Python (Pandas/NumPy) for generating expression count data and R (DESeq2/ggplot2) for core statistical modeling, filtering, and visualization.

Key Findings

A total of 9,996 genes were analyzed across 12 samples (after filtering for low counts). After filtering for an adjusted p-value (FDR) < 0.05 and a Log2(Fold Change) > 0.58 (1.5-fold change), a total of 1 gene was identified as differentially expressed.

Downregulated Gene (Treated < Control)

The only statistically significant gene identified was downregulated in the Treated condition, demonstrating high statistical significance (very low p-value) and a strong fold-change.

Gene ID

Log2(Fold Change)

Adjusted P-value (FDR)

Expression Status

Gene_4252

-1.94936398

2.78E-27

Downregulated

Visualization Output

The statistical results are visualized in a Volcano Plot (saved as volcano_plot.png). This plot shows that the one statistically significant gene (colored blue) falls outside the designated cutoffs, confirming the quality and stringency of the analysis.