import matplotlib.pyplot as plt

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