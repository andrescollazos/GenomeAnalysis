import matplotlib.pyplot as plt

cog_labels = {
    'J': 'Translation, ribosomal structure and biogenesis',
    'A': 'RNA processing and modification',
    'K': 'Transcription',
    'L': 'Replication, recombination and repair',
    'B': 'Chromatin structure and dynamics',
    'D': 'Cell cycle control, cell division, chromosome partitioning',
    'Y': 'Nuclear structure',
    'V': 'Defense mechanisms',
    'T': 'Signal transduction mechanisms',
    'M': 'Cell wall/membrane/envelope biogenesis',
    'N': 'Cell motility',
    'Z': 'Cytoskeleton',
    'W': 'Extracellular structures',
    'U': 'Intracellular trafficking, secretion, vesicular transport',
    'O': 'Post-translational modification, protein turnover, chaperones',
    'C': 'Energy production and conversion',
    'G': 'Carbohydrate transport and metabolism',
    'E': 'Amino acid transport and metabolism',
    'F': 'Nucleotide transport and metabolism',
    'H': 'Coenzyme transport and metabolism',
    'I': 'Lipid transport and metabolism',
    'P': 'Inorganic ion transport and metabolism',
    'Q': 'Secondary metabolites biosynthesis, transport and catabolism',
    'R': 'General function prediction only',
    'S': 'Function unknown'
}

cog_counts = {
    'J': 72,
    'A': 48,
    'K': 115,
    'L': 67,
    'B': 16,
    'D': 8,
    'Y': 0,
    'V': 9,
    'T': 109,
    'M': 21,
    'N': 0,
    'Z': 24,
    'W': 0,
    'U': 62,
    'O': 0,
    'C': 63,
    'G': 75,
    'E': 42,
    'F': 14,
    'H': 16,
    'I': 45,
    'P': 26,
    'Q': 53,
    'R': 0,
    'S': 484,
}

# Define figure and axis
fig, ax = plt.subplots(figsize=(12, 6))

# Sort by COG category order (pre-sorted cog_counts)
codes = list(cog_counts.keys())
values = [cog_counts[c] for c in codes]
labels = [f"{c}: {cog_labels.get(c, 'Unknown')}" for c in codes]

# Bar plot
bars = ax.bar(codes, values, color=plt.cm.tab20.colors[:len(codes)], edgecolor='black')

# Set axis labels and title
ax.set_xlabel('COG Functional Categories')
ax.set_ylabel('Number of Proteins')
ax.set_title('COG Functional Classification')

# Add legend to the right
ax.legend(bars, labels, loc='center left', bbox_to_anchor=(1.01, 0.5), title='Function Class')

plt.tight_layout()
plt.show()


