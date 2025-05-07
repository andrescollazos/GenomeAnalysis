setwd('~/Documents/Studies/Genome Analysis/GenomeAnalysis/data/counts/')

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace("DESeq2", quietly = TRUE)) {
  BiocManager::install("DESeq2")
}

if (!requireNamespace("apeglm", quietly = TRUE)) {
  BiocManager::install("apeglm")
}

if (!requireNamespace("pheatmap", quietly = TRUE)) {
  BiocManager::install("pheatmap")
}

library("DESeq2")
library("pheatmap")
library("ggplot2")

# Load data
counts_data <- read.delim("counts.txt", comment.char = "#", stringsAsFactors = FALSE)

# Columns: Geneid, Chr, Start, End, Strand, Length, Control_1, Control_2, Control_3, Treatment_1, ..., Treatment_3

# Remove columns: Chr, Start, End, Strand, Length
cts <- counts_data[, c(1, 7:ncol(counts_data))]

# Set gene IDs as row names
rownames(cts) <- cts$Geneid
cts <- cts[ , -1]

# Create metadata object
coldata <- data.frame(
  row.names = colnames(cts),
  condition = factor(c("control", "control", "control",
                       "treatment", "treatment", "treatment"))
)

# Construct DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
  countData = cts,
  colData = coldata,
  design = ~ condition
)

dds

# Pre filtering
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Set Factor Level
dds$condition <- relevel(dds$condition, ref = "control")

# Run DESeq
dds <- DESeq(dds)
res <- results(dds)
res
resultsNames(dds)

# Shrinkage of effect size (LFC estimates)
resLFC <- lfcShrink(dds, coef="condition_treatment_vs_control", type="apeglm")

# Summary Results
summary(res)
summary(resLFC)

# MA plot
plotMA(res)
plotMA(resLFC, ylim=c(-6,6))

  # Plot counts
plotCounts(dds, gene=which.min(res$padj), intgroup="condition")

# Convert to data frame
res_df <- as.data.frame(resLFC)

# Categorize fold change magnitude
res_df$diffexp <- cut(
  res_df$log2FoldChange,
  breaks = c(-Inf, -2, -1, 1, 2, Inf),
  labels = c(
    "Highly downregulated",
    "Downregulated",
    "Low fold change",
    "Upregulated",
    "Highly upregulated"
  )
)

ggplot(res_df, aes(x = baseMean, y = log2FoldChange, color = diffexp)) +
  geom_point(size = 1, alpha = 0.6) +
  scale_x_log10() +
  scale_color_manual(
    values = c(
      "Highly downregulated" = "blue",
      "Downregulated" = "skyblue",
      "Low fold change" = "gray40",
      "Upregulated" = "pink",
      "Highly upregulated" = "red"
    ),
    name = "Differential expression"
  ) +
  geom_hline(yintercept = 0, linetype = "solid") +
  labs(
    x = "Mean of normalized counts",
    y = "Log2 fold change (shrunken)"
  ) +
  coord_cartesian(ylim = c(-6, 6)) +
  scale_y_continuous(breaks = seq(-6, 6, by = 2)) +
  theme_minimal()


# Volcano Plot
# Add -log10(padj) column
res_df$neglog10padj <- -log10(res_df$padj)

# Define a significance label (optional)
res_df$sig <- ifelse(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, "yes", "no")
res_df$volcano_color <- ifelse(
  res_df$padj < 0.05 & res_df$log2FoldChange > 1, "up",
  ifelse(res_df$padj < 0.05 & res_df$log2FoldChange < -1, "down", "ns")
)

# Plot
ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = volcano_color)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(
    values = c(
      "up" = "red",
      "down" = "blue",
      "ns" = "gray70"
    ),
    labels = c(
      "up" = "Upregulated (FDR < 0.05)",
      "down" = "Downregulated (FDR < 0.05)",
      "ns" = "Not significant"
    ),
    name = "Significance"
  ) +
  labs(
    title = "Volcano plot",
    x = "log2FoldChange",
    y = "-log10 adjusted p-value"
  ) +
  theme_minimal()

# Heatmap
vsd <- vst(dds, blind=FALSE)

# Filter to only significant genes (padj < 0.05)
res_sig <- resLFC[!is.na(resLFC$padj) & resLFC$padj < 0.05, ]

# Order by smallest adjusted p-value
top_genes <- head(rownames(res_sig[order(res_sig$padj), ]), 50)

# Extract normalized counts (transformed)
mat <- assay(vsd)[top_genes, ]
mat_scaled <- t(scale(t(mat)))

# Annotation for columns
annotation_col <- as.data.frame(colData(dds)[, "condition", drop = FALSE])

pheatmap(
  mat_scaled,
  annotation_col = annotation_col,
  show_rownames = FALSE,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  scale = "none",
  color = colorRampPalette(c("blue", "white", "red"))(100)
)


# PCA
pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  theme_minimal()



