# DESeq2 Differential Expression Analysis
# BT549/BT540 Control vs Knockout

# Install required packages if not already installed
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

if (!require("DESeq2", quietly = TRUE))
    BiocManager::install("DESeq2")

if (!require("tidyverse", quietly = TRUE))
    install.packages("tidyverse")

# Load libraries
library(DESeq2)
library(tidyverse)

# Set working directory
setwd("c:/Users/Hp/Downloads/network_based_therapy_optmisation")

cat("Loading raw counts data...\n")

# Read the raw counts data
count_data <- read.table("GSE300385_raw_counts.txt", 
                         header = TRUE, 
                         row.names = 1, 
                         sep = "\t",
                         check.names = FALSE)

# Extract only BT549/BT540 samples (first 6 columns)
bt549_samples <- c("BT549_CT2", "BT540_CT3", "BT549_CT4", 
                   "BT540_KO3", "BT549_KO8", "BT549_KO11")
count_data_bt549 <- count_data[, bt549_samples]

cat("Samples included:\n")
print(colnames(count_data_bt549))
cat("\nNumber of genes:", nrow(count_data_bt549), "\n\n")

# Create sample metadata
sample_info <- data.frame(
  sample = bt549_samples,
  condition = c("Control", "Control", "Control", 
                "Knockout", "Knockout", "Knockout"),
  row.names = bt549_samples
)

# Ensure condition is a factor with Control as reference
sample_info$condition <- factor(sample_info$condition, 
                                levels = c("Control", "Knockout"))

cat("Sample metadata:\n")
print(sample_info)
cat("\n")

# Filter low count genes (at least 10 reads across all samples)
cat("Filtering low-count genes...\n")
keep <- rowSums(count_data_bt549) >= 10
count_data_filtered <- count_data_bt549[keep, ]
cat("Genes after filtering:", nrow(count_data_filtered), "\n")
cat("Genes removed:", nrow(count_data_bt549) - nrow(count_data_filtered), "\n\n")

# Create DESeq2 dataset
cat("Creating DESeq2 dataset...\n")
dds <- DESeqDataSetFromMatrix(countData = count_data_filtered,
                              colData = sample_info,
                              design = ~ condition)

# Run DESeq2 analysis
cat("Running DESeq2 differential expression analysis...\n")
cat("This may take a few minutes...\n")
dds <- DESeq(dds)

# Extract results (Knockout vs Control)
cat("\nExtracting results...\n")
res <- results(dds, contrast = c("condition", "Knockout", "Control"))

# Summary of results
cat("\nDESeq2 Results Summary:\n")
cat("========================\n")
summary(res)

# Convert to data frame and add gene IDs
res_df <- as.data.frame(res)
res_df$GeneID <- rownames(res_df)
res_df <- res_df[, c("GeneID", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]

# Remove rows with NA values in padj
res_df_clean <- res_df[!is.na(res_df$padj), ]

# Define thresholds
padj_threshold <- 0.05
log2fc_threshold <- 1

# Filter upregulated genes
upregulated <- res_df_clean %>%
  filter(log2FoldChange > log2fc_threshold & padj < padj_threshold) %>%
  arrange(padj)

# Filter downregulated genes
downregulated <- res_df_clean %>%
  filter(log2FoldChange < -log2fc_threshold & padj < padj_threshold) %>%
  arrange(padj)

# Print summary statistics
cat("\n\nDifferentially Expressed Genes:\n")
cat("================================\n")
cat("Upregulated genes (log2FC > 1, padj < 0.05):", nrow(upregulated), "\n")
cat("Downregulated genes (log2FC < -1, padj < 0.05):", nrow(downregulated), "\n")
cat("Total differentially expressed genes:", nrow(upregulated) + nrow(downregulated), "\n\n")

# Save results to files
cat("Saving results to files...\n")

# Save upregulated genes
write.table(upregulated, 
            file = "upregulated_genes.txt", 
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)
cat("Saved:", nrow(upregulated), "upregulated genes to upregulated_genes.txt\n")

# Save downregulated genes
write.table(downregulated, 
            file = "downregulated_genes.txt", 
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)
cat("Saved:", nrow(downregulated), "downregulated genes to downregulated_genes.txt\n")

# Save full results (sorted by padj)
res_df_sorted <- res_df_clean %>% arrange(padj)
write.table(res_df_sorted, 
            file = "deseq2_full_results.txt", 
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)
cat("Saved: Full results for all", nrow(res_df_sorted), "genes to deseq2_full_results.txt\n")

# Display top upregulated genes
if (nrow(upregulated) > 0) {
  cat("\n\nTop 10 Upregulated Genes:\n")
  cat("=========================\n")
  print(head(upregulated[, c("GeneID", "log2FoldChange", "padj")], 10))
}

# Display top downregulated genes
if (nrow(downregulated) > 0) {
  cat("\n\nTop 10 Downregulated Genes:\n")
  cat("===========================\n")
  print(head(downregulated[, c("GeneID", "log2FoldChange", "padj")], 10))
}

# Optional: Create visualization plots
cat("\n\nCreating visualization plots...\n")

# Volcano plot
pdf("volcano_plot.pdf", width = 10, height = 8)
plot_data <- res_df_clean
plot_data$Significant <- "Not Significant"
plot_data$Significant[plot_data$log2FoldChange > log2fc_threshold & 
                      plot_data$padj < padj_threshold] <- "Upregulated"
plot_data$Significant[plot_data$log2FoldChange < -log2fc_threshold & 
                      plot_data$padj < padj_threshold] <- "Downregulated"

plot(plot_data$log2FoldChange, 
     -log10(plot_data$padj),
     col = ifelse(plot_data$Significant == "Upregulated", "red",
                  ifelse(plot_data$Significant == "Downregulated", "blue", "gray")),
     pch = 20,
     xlab = "log2 Fold Change",
     ylab = "-log10(adjusted p-value)",
     main = "Volcano Plot: BT549/BT540 KO vs Control")
abline(v = c(-log2fc_threshold, log2fc_threshold), lty = 2, col = "black")
abline(h = -log10(padj_threshold), lty = 2, col = "black")
legend("topright", 
       legend = c(paste("Upregulated (", nrow(upregulated), ")", sep=""),
                  paste("Downregulated (", nrow(downregulated), ")", sep=""),
                  "Not Significant"),
       col = c("red", "blue", "gray"),
       pch = 20)
dev.off()
cat("Saved volcano plot to volcano_plot.pdf\n")

# MA plot
pdf("ma_plot.pdf", width = 10, height = 8)
plotMA(res, ylim = c(-5, 5), 
       main = "MA Plot: BT549/BT540 KO vs Control")
dev.off()
cat("Saved MA plot to ma_plot.pdf\n")

cat("\n\nAnalysis complete!\n")
cat("==================\n")
cat("Output files created:\n")
cat("  - upregulated_genes.txt\n")
cat("  - downregulated_genes.txt\n")
cat("  - deseq2_full_results.txt\n")
cat("  - volcano_plot.pdf\n")
cat("  - ma_plot.pdf\n")
