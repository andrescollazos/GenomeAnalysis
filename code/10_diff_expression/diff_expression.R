if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace("DESeq2", quietly = TRUE)) {
  BiocManager::install("DESeq2")
}

library("DESeq2")

setwd('~/Documents/Studies/Genome Analysis/GenomeAnalysis/data/counts/')

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

dds

dds <- DESeq(dds)
res <- results(dds)
res
