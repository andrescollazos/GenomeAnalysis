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
plotMA(resLFC, ylim=c(-2,2))

# Detect the row number of individual genes by clicking on the plot
idx <- identify(res$baseMean, res$log2FoldChange)
rownames(res)[idx]

# Plot counts
plotCounts(dds, gene=which.min(res$padj), intgroup="condition")

# Heat maps
ntd <- normTransform(dds)
vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)

# Select top 20 expressed genes
select <- order(rowMeans(counts(dds, normalized=TRUE)), decreasing=TRUE)[1:20]

# Plot heatmaps
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=coldata)

pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=coldata)

pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=coldata)

