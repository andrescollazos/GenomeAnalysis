import matplotlib.pyplot as plt

# ------------------------------------------------------------
# Gene length distribution

# Load gene lengths (See script 03_extract_genes_length.sh)
with open("gene_lengths.txt") as f:
    gene_lengths = [int(line.strip()) for line in f]

# Convert to kilobases
gene_lengths_kb = [l / 1000 for l in gene_lengths]

# Plot
plt.figure(figsize=(8, 5))
plt.hist(gene_lengths_kb, bins=30, color='#4c72b0', edgecolor='black')
plt.xlabel('Gene length (kb)')
plt.ylabel('Number of genes')
plt.title('Distribution of Gene Lengths')
plt.tight_layout()
plt.show()

# ------------------------------------------------------------
# Exon distribution

from collections import Counter

# Load exon counts
with open("exon_counts_only.txt") as f:
    exon_counts = [int(line.strip()) for line in f]

# Frequency distribution
count_distribution = Counter(exon_counts)
x = sorted(count_distribution.keys())
y = [count_distribution[i] for i in x]

# Plot
plt.figure(figsize=(8, 5))
plt.bar(x, y, color='#55a868', edgecolor='black')
plt.xlabel('Number of exons per transcript')
plt.ylabel('Number of transcripts')
plt.title('Exon Count Distribution per Transcript')
plt.tight_layout()
plt.show()

# ------------------------------------------------------------
# Chromosome space usage
# Values in bp

cds = 3651888
introns = 5449492
intergenic = 7880430

# Labels and sizes
labels = ['CDS', 'Introns', 'Intergenic']
sizes = [cds, introns, intergenic]
colors = ['#4c72b0', '#55a868', '#c44e52']

# Plot
plt.figure(figsize=(6, 6))
plt.pie(sizes, labels=labels, autopct='%1.1f%%', startangle=90, colors=colors)
plt.title('Genome Space Usage in Chromosome 3')
plt.tight_layout()
plt.show()