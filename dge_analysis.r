# Differential Gene Expression Analysis using DESeq2

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#BiocManager::install("DESeq2") 
#install.packages("ggplot2")
#install.packages("dplyr")
#install.packages("reshape2") 

library(DESeq2)
library(ggplot2)
library(dplyr)
library(reshape2)

# --- 1. Load Data ---
print("1. Loading count and metadata files...")



# Load the count data (Gene IDs as row names)
countData <- read.csv("mock_expression_counts.csv", row.names = 1, stringsAsFactors = FALSE)

countData <- round(countData) 
countData <- as.matrix(countData)


colData <- read.csv("sample_metadata.csv", row.names = 1, stringsAsFactors = FALSE)
colData$Condition <- factor(colData$Condition) # Convert condition column to factor

# Ensure the sample names and order are identical
if(!all(rownames(colData) == colnames(countData))) {
  stop("Sample names in count data and metadata do not match or are not in the same order!")
}

# --- 2. Create DESeq2 Object and Run Analysis ---
print("2. Creating DESeq2 dataset and running analysis...")
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData,
                              design = ~ Condition)

# Filter out genes with very low counts
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Run the core DGE analysis
dds <- DESeq(dds)

# --- 3. Extract Results and Filtering ---

res <- results(dds, contrast=c("Condition", "Treated", "Control"))


res_df <- as.data.frame(res)
res_df$GeneID <- rownames(res_df)


P_THRESHOLD <- 0.05
FC_THRESHOLD <- 1.5 

res_df <- res_df %>%
  mutate(is_DGE = case_when(
    padj < P_THRESHOLD & log2FoldChange > log2(FC_THRESHOLD) ~ "Upregulated",
    padj < P_THRESHOLD & log2FoldChange < -log2(FC_THRESHOLD) ~ "Downregulated",
    TRUE ~ "Not Significant"
  ))

# --- 4. Generate Volcano Plot---
print("4. Generating Volcano Plot...")


res_df$negLog10Pval <- -log10(res_df$padj)

max_finite_pval <- max(res_df$negLog10Pval[is.finite(res_df$negLog10Pval)], na.rm = TRUE)
res_df[!is.na(res_df$negLog10Pval) & res_df$negLog10Pval == Inf, "negLog10Pval"] <- max_finite_pval

volcano_plot <- ggplot(res_df, aes(x = log2FoldChange, y = negLog10Pval, color = is_DGE)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Upregulated" = "#B31B21", "Downregulated" = "#1465AC", "Not Significant" = "gray50"),
                     name = "Expression Status") +
 
  geom_hline(yintercept = -log10(P_THRESHOLD), linetype = "dashed", color = "black", alpha = 0.5) +
  geom_vline(xintercept = c(log2(FC_THRESHOLD), -log2(FC_THRESHOLD)), linetype = "dashed", color = "black", alpha = 0.5) +
  labs(
    title = "Volcano Plot: Treated vs. Control DGE Analysis",
    x = expression(paste(Log[2], " (Fold Change)")),
    y = expression(paste("-Log[10] (Adjusted P-value)")),
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5, face = "bold"))


ggsave("volcano_plot.png", plot = volcano_plot, width = 8, height = 7, dpi = 300)
print("Volcano Plot saved as volcano_plot.png")

# --- 5. Save Filtered Results (For Final Report) ---
significant_genes <- res_df %>%
  filter(is_DGE != "Not Significant") %>%
  arrange(padj) %>%
  select(GeneID, log2FoldChange, padj, is_DGE)

write.csv(significant_genes, "significant_genes_report.csv", row.names = FALSE)
print("Significant genes report saved as significant_genes_report.csv")
