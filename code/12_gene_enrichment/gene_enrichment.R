if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace("rrvgo", quietly = TRUE)) {
  BiocManager::install("rrvgo")
}

if (!requireNamespace("org.At.tair.db", quietly = TRUE)) {
  BiocManager::install("org.At.tair.db")
}

if (!requireNamespace("GOstats", quietly = TRUE)) {
  BiocManager::install("GOstats")
}

if (!requireNamespace("GSEABase", quietly = TRUE)) {
  BiocManager::install("GSEABase")
}



library(dplyr)
library(stringr)
library(tidyr)
library(tibble)
library("rrvgo")
library(GOstats)
library(GSEABase)


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

# ------------------------------------------------------------------------------
# Load Differential Expression Results
res_df <- read.csv("deseq2_results.csv", row.names = 1)
sig_genes <- rownames(res_df %>% filter(padj < 0.05 & abs(log2FoldChange) > 1))

gene_universe <- rownames(res_df)
gene_ids <- sig_genes



# Filter and create clean base R data.frame
go_df <- data.frame(
  go_id = as.character(go_table$go_id),
  Evidence = "IEA",
  gene_id = as.character(go_table$gene),
  stringsAsFactors = FALSE
)
goframe <- GOFrame(go_df, organism = "Niphotrichum japonicum")
goallframe <- GOAllFrame(goframe)

# Make geneSetCollection object
gsc <- GeneSetCollection(goallframe, setType = GOCollection())

params <- GSEAGOHyperGParams(
  name = "GOstats enrichment",
  geneSetCollection = gsc,
  geneIds = gene_ids,
  universeGeneIds = gene_universe,
  ontology = "BP",
  pvalueCutoff = 0.05,
  conditional = FALSE,
  testDirection = "over"
)

go_result <- hyperGTest(params)
go_table_stats <- summary(go_result)
write.csv(go_table_stats, "gostats_enrichment.csv", row.names = FALSE)

simMatrix <- calculateSimMatrix(
  go_table_stats$GOBPID,
  orgdb = "org.At.tair.db",
  ont = "BP",
  method = "Rel"
)
scores <- setNames(-log10(go_table_stats$Pvalue), go_table_stats$GOBPID)
reducedTerms <- reduceSimMatrix(
  simMatrix,
  scores,
  orgdb = "org.At.tair.db",
  threshold = 0.7
)
scatterPlot(simMatrix, reducedTerms, addLabel = TRUE, onlyParents = TRUE)
scatterPlot(simMatrix, reducedTerms, addLabel = TRUE)

heatmapPlot(
  simMatrix,
  reducedTerms,
  annotateParent=TRUE,
  annotationLabel="parentTerm",
  fontsize=6
)

treemapPlot(reducedTerms)
wordcloudPlot(reducedTerms, min.freq=1, colors="black")












