setwd('~/Documents/Studies/Genome Analysis/GenomeAnalysis/data/counts/')

# Load libraries
library(ggplot2)
library(reshape2)

# Load the data
counts <- read.table("counts.txt", header = TRUE, sep = "\t", comment.char = "")

# Extract the count columns
# 1       2   3     4   5       6       7                     8
# Geneid	Chr	Start	End	Strand	Length	Control_1.sorted.bam	Control_2.sorted.bam	...
count_data <- counts[, 7:ncol(counts)]

# Compute total count per gene
counts$TotalCount <- rowSums(count_data)

# Remove zeros to avoid log10(0)
counts_nonzero <- counts[counts$TotalCount > 0, ]

# Apply log10 to total counts
counts_nonzero$Log10Count <- log10(counts_nonzero$TotalCount)

# Plot
hist(counts_nonzero$Log10Count,
     breaks = 50,
     col = "steelblue",
     border = "black",
     main = "Histogram of Gene Counts (log10)",
     xlab = "log10 Total Counts per Gene")


