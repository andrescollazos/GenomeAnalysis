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













# Heat maps
ntd <- normTransform(dds)
vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)

# Select the top 30 most variable genes
select <- order(rowVars(assay(vsd)), decreasing=TRUE)

# Plot
pheatmap(assay(vsd)[select, ],
         cluster_rows=TRUE,
         show_rownames=FALSE,
         cluster_cols=TRUE,
         annotation_col=coldata,
         color = colorRampPalette(c("green", "black", "red"))(100))

