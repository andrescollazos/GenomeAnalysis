if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace("rrvgo", quietly = TRUE)) {
  BiocManager::install("rrvgo")
}

library(dplyr)
library(stringr)
library(tidyr)

setwd('~/Documents/Studies/Genome Analysis/GenomeAnalysis/data/gene_enrichment/')


# ------------------------------------------------------------------------------
# Extract GO Terms from eggNOG GFF file
gff <- read.delim("n_japonicum.gff", header = FALSE, comment.char = "#", sep = "\t")
colnames(gff) <- c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")

# Filter for mRNA and extract gene ID and GO terms
go_table <- gff %>%
  filter(type == "mRNA", source == "AUGUSTUS") %>%
  mutate(
    gene_id = str_extract(attributes, "Parent=[^;]+") %>% str_remove("Parent="),
    go_terms = str_extract(attributes, "GO:[^;]+") %>% str_extract_all("GO:\\d+")
  ) %>%
  unnest(go_terms) %>%
  distinct(gene_id, go_terms)

colnames(go_table) <- c("gene", "go_id")
go_table <- go_table %>% filter(!is.na(go_id))

# Load Differential Expression Results
res_df <- read.csv("deseq2_results.csv", row.names = 1)
sig_genes <- rownames(res_df %>% filter(padj < 0.05 & abs(log2FoldChange) > 1))

# Match with GO-annotated genes
sig_go_table <- go_table %>% filter(gene %in% sig_genes)









